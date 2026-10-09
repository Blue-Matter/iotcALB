# Select the MP settings not set by tuning (ConfigCandidates in 02-define-CMPs.R)
#
# For each MP: the configuration with the highest mean catch when tuned to the
# constraints of tuning level ConfigTuning (T1: PGK). The selected settings are
# used for the CMPs of that MP at every tuning level in 02b-tune-CMPs.R.
# Selected on the base case OM only, with nSimTune simulations
source('analysis/02-define-CMPs.R')

ConfigMPs       <- names(ConfigCandidates)  # MPs to select settings for
ConfigTACChange <- 0.15                     # TAC change limit of the CMP used for the selection
nWorkers        <- 40                       # workers for MSEtool::Project(parallel = TRUE)

Control <- MSEtool::TuneControl(Interval    = TuneInterval,
                                nGridConfig = 5,
                                nSimConfig  = 40,  # simulations used to compare configurations
                                nRefine     = 1)   # best configuration tuned on all nSimTune simulations

# Configurations of an MP: all combinations of its candidate values
Candidates <- function(MP) {
  Cand <- ConfigCandidates[[MP]]
  Grid <- expand.grid(lapply(Cand, seq_along), KEEP.OUT.ATTRS = FALSE)
  lapply(seq_len(nrow(Grid)), \(i) purrr::imap(Cand, \(x, nm) x[[Grid[i, nm]]]))
}

MSEtool::SetupParallel(workers = nWorkers)

for (mp in ConfigMPs) {
  cli::cli_h1('Configurations: {mp}')
  i <- which(CMPs$MP == mp & CMPs$TACChange == ConfigTACChange & CMPs$Tuning == ConfigTuning)
  Cand <- Candidates(mp)
  Args <- lapply(Cand, \(x) ConfigArgs(mp, x))

  Tuned <- MSEtool::TuneMP(HistTune, CMPs$Fn[[i]],
                           Objective   = Objective,
                           Constraints = CMPs$Constraints[[i]],
                           Configs     = Args,
                           Control     = Control,
                           parallel    = TRUE)
  Tuned@MPName <- mp
  saveRDS(Tuned, file.path(TuneDir, paste0('Configs_', mp, '.tune')))

  Tab <- MSEtool::TuneTable(Tuned, 'Configs')
  Tab$Label <- vapply(Cand, ConfigLabel, character(1))[Tab$Config]
  print(Tab, row.names = FALSE)

  if (Tuned@Status == 'infeasible') {
    cli::cli_alert_warning("No configuration of {mp} met the constraints; settings not updated.")
  } else {
    # the candidate whose MP arguments were selected
    Selected <- which(vapply(Args, \(a) isTRUE(all.equal(Tuned@Args[names(a)], a, check.attributes = FALSE)),
                             logical(1)))
    SelectedConfigs[[mp]] <- Cand[[Selected[1]]]
    saveRDS(SelectedConfigs, ConfigsFile)
    cli::cli_alert_success("{mp}: {ConfigLabel(SelectedConfigs[[mp]])}")
  }
}

MSEtool::DisableParallel()
SelectedConfigs
