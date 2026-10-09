# Run the MSE projections of the tuned CMPs (02b-tune-CMPs.R) for each OM,
# with all simulations
library(iotcALB)

source('analysis/00-specifications.R')

TuneDir <- 'objects/Tune'
MSEDir  <- 'objects/MSE'
dir.create(MSEDir, showWarnings = FALSE)

TunedCMPs <- readRDS(file.path(TuneDir, 'TunedCMPs.rds'))
if (!length(TunedCMPs))
  cli::cli_abort("No tuned CMPs in {.file {file.path(TuneDir, 'TunedCMPs.rds')}}; run 02b-tune-CMPs.R.")

# Reference MPs
RefMPs <- list(NoFishing = 'NoFishing')

# OMs to project
ProjectOMs <- 'OM5b'

OM_DF <- data.frame(OM = OMNames) |>
  dplyr::mutate(HistFile = file.path('objects/Hist', paste0(OM, '.hist')),
                MSEFile  = file.path(MSEDir, paste0(OM, '.mse')),
                RunMSE   = OM %in% ProjectOMs & file.exists(HistFile))

nWorkers <- 40  # workers for MSEtool::Project(parallel = TRUE)

# ---- Project ----
if (any(OM_DF$RunMSE)) {
  MSEtool::SetupParallel(workers = nWorkers)
  for (i in seq_len(nrow(OM_DF))) {
    if (!OM_DF$RunMSE[i]) next
    cli::cli_h1('OM {OM_DF$OM[i]}')
    Hist <- readRDS(OM_DF$HistFile[i])
    MSE  <- MSEtool::Project(Hist, MPs = c(TunedCMPs, RefMPs), parallel = TRUE)
    MSEtool::Save(MSE, OM_DF$MSEFile[i], overwrite = TRUE)
    rm(Hist, MSE); gc()
  }
  MSEtool::DisableParallel()
}

# ---- Performance ----
PMTable <- function(MSE) {
  PMs <- list(PGK_All     = MSEtool::PM_Status(MSE),
              PGK_Tune    = MSEtool::PM_Status(MSE,
                                               Years = MSEtool::PMYears(MSE, 'last', n = 15)),
              Safety      = MSEtool::PM_Safety(MSE,
                                               Lim = 0.4,
                                               Years = MSEtool::PMYears(MSE, Skip = 5)),
              PLRP        = MSEtool::PM_SBSBlim(MSE, Lim = 0.4),
              PRed        = MSEtool::PM_Red(MSE),
              Catch_All   = MSEtool::PM_Removals(MSE),
              Catch_Short = MSEtool::PM_Removals(MSE,
                                                 Years = MSEtool::PMYears(MSE, 'first', n = 10)),
              AAVTAC      = MSEtool::PM_AAVY(MSE,
                                             Type = 'TAC',
                                             IncludeFirst = TRUE))
  out <- lapply(PMs, \(pm) colMeans(pm@Mean))
  data.frame(MP = names(out[[1]]), do.call(cbind, out), row.names = NULL)
}

Done <- OM_DF[OM_DF$OM %in% ProjectOMs & file.exists(OM_DF$MSEFile), ]
PerformanceDF <- do.call(rbind, lapply(seq_len(nrow(Done)), \(i) {
  cbind(OM = Done$OM[i], PMTable(readRDS(Done$MSEFile[i])))
}))
print(PerformanceDF, digits = 3, row.names = FALSE)
utils::write.csv(PerformanceDF, file.path(MSEDir, 'Performance.csv'),
                 row.names = FALSE)
