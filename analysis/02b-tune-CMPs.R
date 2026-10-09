# Tune the candidate management procedures on the base case OM: maximise mean
# catch subject to the constraints of each tuning level
#
# TuneLevels:
# - T1: PGK;
# - T2: PGK and Safety.
#
# Tuned with nSimTune simulations. Run 02a-select-configs.R
# first to select the MP settings not set by tuning

source('analysis/02-define-CMPs.R')

Overwrite <- FALSE     # re-tune CMPs with a saved result for the same MP settings
DryRun    <- FALSE     # only estimate the number of projections and run time
nWorkers  <- 40        # workers for MSEtool::Project(parallel = TRUE)

Control <- MSEtool::TuneControl(Interval = TuneInterval, MaxExpand = 3, DryRun = DryRun)

# MP arguments other than tunepar; a saved result is re-used only if these are unchanged
MPSettings <- function(MP) {
  Args <- formals(MP)
  Args[names(Args) != 'tunepar']
}

# Objective and constraints; a saved result is also re-used only if these are unchanged
MetricSpec <- function(Metrics) {
  lapply(Metrics, \(m) list(Name = m$Name, Type = m$Type, Min = m$Min, Max = m$Max,
                            Args = m$Args, PM = deparse(body(m$PM))))
}
TuneSpec <- function(i) MetricSpec(c(list(Objective), CMPs$Constraints[[i]]))

# TuneMP records a failed projection as infeasible without the error
AllFailed <- function(Tuned) all(Tuned@History$FailRate >= 1)


IsCurrent <- function(Tuned, i) {
  identical(MPSettings(Tuned@MP), MPSettings(CMPs$Fn[[i]])) &&
    identical(MetricSpec(Tuned@Settings$Metrics), TuneSpec(i)) &&
    identical(Tuned@Settings$Sims, TuneSims) &&
    !AllFailed(Tuned)
}

TuneFile <- function(i) file.path(TuneDir, paste0(CMPs$Name[i], '.tune'))

ReuseConfigTuning <- function(i) {
  File <- file.path(TuneDir, paste0('Configs_', CMPs$MP[i], '.tune'))
  if (!file.exists(File)) return(FALSE)
  Tuned <- readRDS(File)
  Tuned@Settings$Sims <- TuneSims
  if (!IsCurrent(Tuned, i)) return(FALSE)
  Tuned@MPName <- CMPs$Name[i]
  saveRDS(Tuned, TuneFile(i))
  TRUE
}

TuneCMP <- function(i) {
  if (file.exists(TuneFile(i)) && !Overwrite && IsCurrent(readRDS(TuneFile(i)), i))
    return('exists')
  if (!Overwrite && ReuseConfigTuning(i))
    return('re-used from 02a-select-configs.R')
  tryCatch({
    Tuned <- MSEtool::TuneMP(HistTune, CMPs$Fn[[i]],
                             Objective   = Objective,
                             Constraints = CMPs$Constraints[[i]],
                             Control     = Control,
                             parallel    = TRUE,
                             silent      = FALSE)
    Tuned@MPName <- CMPs$Name[i]
    Tuned@Settings$Sims <- TuneSims
    saveRDS(Tuned, TuneFile(i))
    if (AllFailed(Tuned)) 'Error: every projection failed' else Tuned@Status
  }, error = function(e) paste('Error:', conditionMessage(e)))
}

# ---- Tune ----
if (DryRun) {
  # one CMP of each type
  for (i in match(unique(CMPs$MP), CMPs$MP)) {
    cli::cli_h2(CMPs$Name[i])
    MSEtool::TuneMP(HistTune, CMPs$Fn[[i]], Objective = Objective, Constraints = CMPs$Constraints[[i]],
                    Control = Control)
  }
} else {
  MSEtool::SetupParallel(workers = nWorkers)
  TuneLog <- character(nrow(CMPs))
  for (i in seq_len(nrow(CMPs))) {
    cli::cli_h1(CMPs$Name[i])
    TuneLog[i] <- TuneCMP(i)
    cli::cli_inform('{CMPs$Name[i]}: {TuneLog[i]}')
  }
  MSEtool::DisableParallel()
  print(data.frame(CMP = CMPs$Name, Result = TuneLog), row.names = FALSE)
}

# ---- Summarise ----
Tuned  <- lapply(seq_len(nrow(CMPs)), \(i) if (file.exists(TuneFile(i))) readRDS(TuneFile(i)))
IsDone <- vapply(seq_along(Tuned), \(i) !is.null(Tuned[[i]]) && IsCurrent(Tuned[[i]], i), logical(1))

if (!DryRun && any(IsDone)) {
  Summary <- CMPs[IsDone, c('Name', 'MP', 'Tuning', 'TACChange', 'TACChangeDown')]
  Summary$TuningLabel <- TuneLevelLabel[Summary$Tuning]
  Summary$Config     <- vapply(Summary$MP, \(mp) {
    Config <- TypeConfig(mp)
    if (length(Config)) ConfigLabel(Config) else 'Earlier MP settings'
  }, character(1))
  Summary$tunepar    <- vapply(Tuned[IsDone], \(x) x@Args$tunepar, numeric(1))
  Summary$TuneStatus <- vapply(Tuned[IsDone], \(x) x@Status, character(1))
  Summary[[Objective$Name]] <- vapply(Tuned[IsDone], \(x) x@Objective, numeric(1))
  # constraint values at the tuned tunepar (NA if not a constraint of the tuning level)
  for (nm in unique(unlist(lapply(TuneLevels, \(l) vapply(l, \(x) x$Name, character(1))))))
    Summary[[nm]] <- vapply(Tuned[IsDone], \(x) {
      Tab <- MSEtool::TuneTable(x, 'Constraints')
      if (nm %in% Tab$Name) Tab$Value[Tab$Name == nm][1] else NA_real_
    }, numeric(1))

  print(Summary, row.names = FALSE)
  if (any(Summary$TuneStatus == 'at_bound'))
    cli::cli_alert_warning("Catch still increasing at the end of the tuning range for {.val {Summary$Name[Summary$TuneStatus == 'at_bound']}}; widen {.arg Interval} in {.fn TuneControl}.")
  if (any(Summary$TuneStatus == 'infeasible'))
    cli::cli_alert_warning("No tuning met the constraints for {.val {Summary$Name[Summary$TuneStatus == 'infeasible']}}.")
  if (!all(IsDone))
    cli::cli_alert_warning("Not tuned with the current MP settings and tuning specification: {.val {CMPs$Name[!IsDone]}}.")

  utils::write.csv(Summary, file.path(TuneDir, 'Summary.csv'), row.names = FALSE)

  ## Tuned CMPs (feasible tunings) for the projections
  Keep <- Summary$TuneStatus != 'infeasible'
  TunedCMPs <- lapply(Tuned[IsDone][Keep], \(x) x@MP)
  names(TunedCMPs) <- Summary$Name[Keep]
  saveRDS(TunedCMPs, file.path(TuneDir, 'TunedCMPs.rds'))
}
