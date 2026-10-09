# Performance of the tuned CMPs on the base case OM with the nSimTune simulations
# used for tuning (quick check before the full projections in 03-run-projections.R).
#
# Safety: probability SB > LRP in every year of SafetyYears (tuning constraint of T2)
# PLRP: probability SB > LRP (proportion of simulations and years the MP is active)
# P_Collapse: probability SB < 0.1 SBMSY in any year the MP is active
source('analysis/02-define-CMPs.R')

TunedCMPs <- readRDS(file.path(TuneDir, 'TunedCMPs.rds'))
MSEtool::SetupParallel(workers = 24)
MSE <- MSEtool::Project(HistTune,
                        MPs = c(TunedCMPs, list(NoFishing = 'NoFishing')),
                        parallel = TRUE, silent = TRUE)
MSEtool::DisableParallel()
saveRDS(MSE, file.path(TuneDir, 'Tuned_80sims.mse'))

MinSB <- MSEtool::PM_MinStatus(MSE, Years = CatchYears)
MPs   <- colnames(MSEtool::PM_Status(MSE)@Mean)
Performance <- data.frame(
  MP          = MPs,
  PGK_Tune    = colMeans(MSEtool::PM_Status(MSE, Years = TuneYears)@Mean),
  PGK_All     = colMeans(MSEtool::PM_Status(MSE)@Mean),
  Safety      = colMeans(MSEtool::PM_Safety(MSE, Lim = LRP, Years = SafetyYears)@Mean),
  PLRP        = colMeans(MSEtool::PM_SBSBlim(MSE, Lim = LRP)@Mean),
  P_Collapse  = colMeans(MinSB@Stat[, 1, MPs] < 0.1),
  Catch_All   = colMeans(MSEtool::PM_Removals(MSE, Years = CatchYears)@Mean),
  Catch_Short = colMeans(MSEtool::PM_Removals(MSE, Years = MSEtool::PMYears(MSE, 'first', n = 10))@Mean),
  PRed        = colMeans(MSEtool::PM_Red(MSE)@Mean),
  AAVTAC      = colMeans(MSEtool::PM_AAVY(MSE, Type = 'TAC', IncludeFirst = TRUE)@Mean),
  row.names = NULL)
print(Performance, digits = 3, row.names = FALSE)
utils::write.csv(Performance, file.path(TuneDir, 'Tuned_PMs_80sims.csv'), row.names = FALSE)
