# Simulations where the stock collapses before the MPs start (2029)
#
# The interim catches (2021-2023 reported, 2024-2028 assumed) are applied to
# every simulation of the base case OM, with no fishing from 2029 (MPs cannot
# affect these years). Simulations with SB < 0.1 SBMSY in 2029 are compared
# with the other simulations on the ABC conditioning quantities, the fit to the
# CPUE used for conditioning, and the catches relative to MSY.
#
# Output (read by the TSD):
#   figures/diagnostics/interim_collapse/Collapse_Sims.csv   per-simulation diagnostics
#   figures/diagnostics/interim_collapse/SB_SBMSY.png        SB/SBMSY trajectories
#   figures/diagnostics/interim_collapse/CPUE_fit.png        ABC fit to the LL1 CPUE

library(iotcALB)

BaseOM    <- 'OM5b'
Threshold <- 0.1      # collapse: SB < Threshold SBMSY in MPStartYear
nWorkers  <- 24
OutDir    <- 'figures/diagnostics/interim_collapse'
dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)

Hist <- readRDS(file.path('objects/Hist', paste0(BaseOM, '.hist')))
MPStartYear <- Hist@OM@MPStartYear

MSEtool::SetupParallel(workers = nWorkers)
MSE <- MSEtool::Project(Hist, MPs = list(NoFishing = 'NoFishing'), parallel = TRUE, silent = TRUE)
MSEtool::DisableParallel()

# annual female SB/SBMSY and SB/SB0 (mean over the values in each calendar year)
Annual <- function(df) {
  df <- as.data.frame(df[df$Stock == 'Female', ])
  df$Year <- floor(df$Year)
  df <- stats::aggregate(Value ~ Sim + Year, df, mean)
  tapply(df$Value, list(df$Sim, df$Year), mean)
}

SBMSY <- Annual(MSEtool::SB_SBMSY(MSE, df = TRUE, Reduce = FALSE))
SB0   <- Annual(MSEtool::SB_SB0(MSE, df = TRUE, Reduce = FALSE))
Collapse <- SBMSY[, as.character(MPStartYear)] < Threshold

# ---- ABC conditioning quantities ----
Cond <- get(paste0('CondData_', BaseOM))
g    <- function(f) vapply(Cond, f, numeric(1))

ObsCPUE <- matrix(albMSE_Data@CPUE@Value[, 'LL1'], ncol = 4, byrow = TRUE)  # year x season
CPUERes <- lapply(Cond, \(x) { r <- log(ObsCPUE / x$Ihat); r - mean(r, na.rm = TRUE) })

Catch <- rowSums(ProcessNewData()$Landings)

# MSY of the OM (total removals at MSY). The MSY reported by the ABC conditioning
# model (Cmsy) is not used: it is calculated with fixed values of R0, h, and M
# rather than those of each simulation (see the TSD, Operating Models)
OMMSY <- apply(MSEtool::MSYLandings(Hist), 1, sum) + apply(MSEtool::MSYDiscards(Hist), 1, sum)
Catch <- tapply(Catch, floor(as.numeric(names(Catch))), sum)

DF <- data.frame(
  Sim             = seq_along(Cond),
  Collapse        = Collapse,
  B0              = g(\(x) x$B0),
  MSY             = OMMSY,
  M               = g(\(x) x$M),
  h               = g(\(x) x$h),
  ABC_SB_SB0      = g(\(x) utils::tail(x$dep, 1)),
  ABC_SB_SBMSY    = g(\(x) utils::tail(x$dbmsy, 1)),
  ABC_HR_HRMSY    = g(\(x) mean(utils::tail(x$hmsyrat, 5))),
  MaxH            = g(\(x) max(x$H)),
  CPUE_RMSE       = vapply(CPUERes, \(r) sqrt(mean(r^2, na.rm = TRUE)), numeric(1)),
  CPUE_Res_Recent = vapply(CPUERes, \(r) mean(r[17:21, ], na.rm = TRUE), numeric(1)),
  OM_SB_SBMSY     = SBMSY[, '2020'],
  OM_SB_SBMSY_MP  = SBMSY[, as.character(MPStartYear)],
  Catch_MSY       = mean(Catch[as.character(2021:2023)]) / OMMSY)
utils::write.csv(DF, file.path(OutDir, 'Collapse_Sims.csv'), row.names = FALSE)

cli::cli_inform('{sum(Collapse)} of {length(Collapse)} simulations ({round(100 * mean(Collapse), 1)}%) have SB < {Threshold} SBMSY in {MPStartYear}.')

# ---- Figures ----
Traj <- data.frame(Sim = rep(as.integer(rownames(SBMSY)), ncol(SBMSY)),
                   Year = rep(as.numeric(colnames(SBMSY)), each = nrow(SBMSY)),
                   Value = as.vector(SBMSY))
Traj <- Traj[Traj$Year <= MPStartYear, ]
Traj$Group <- ifelse(Collapse[as.character(Traj$Sim)], 'Collapse before MPs start', 'Other simulations')

p1 <- ggplot2::ggplot(Traj, ggplot2::aes(.data$Year, .data$Value, group = .data$Sim, colour = .data$Group)) +
  ggplot2::geom_line(data = Traj[Traj$Group == 'Other simulations', ], alpha = 0.1) +
  ggplot2::geom_line(data = Traj[Traj$Group != 'Other simulations', ], linewidth = 0.6) +
  ggplot2::geom_vline(xintercept = Hist@OM@CurrentYear, linetype = 2) +
  ggplot2::geom_hline(yintercept = c(0.4, 1), linetype = 3) +
  ggplot2::scale_colour_manual(values = c('Collapse before MPs start' = 'firebrick', 'Other simulations' = 'grey40')) +
  ggplot2::coord_cartesian(ylim = c(0, 5)) +
  ggplot2::labs(x = NULL, y = expression(SB/SB[MSY]), colour = NULL) +
  ggplot2::theme_bw() + ggplot2::theme(legend.position = 'bottom')
ggplot2::ggsave(file.path(OutDir, 'SB_SBMSY.png'), p1, width = 8, height = 4.5)

Years <- 2000:2020
CPUEFit <- do.call(rbind, lapply(seq_along(Cond), \(i) {
  Pred <- rowMeans(Cond[[i]]$Ihat)
  data.frame(Sim = i, Year = Years, Value = Pred / mean(Pred),
             Group = ifelse(Collapse[i], 'Collapse before MPs start', 'Other simulations'))
}))
Obs <- rowMeans(ObsCPUE, na.rm = TRUE)
ObsDF <- data.frame(Year = Years, Value = Obs / mean(Obs))
Ribbon <- stats::aggregate(Value ~ Year + Group, CPUEFit, \(x) stats::quantile(x, c(0.05, 0.5, 0.95)))
Ribbon <- cbind(Ribbon[, c('Year', 'Group')], as.data.frame(Ribbon$Value))
names(Ribbon)[3:5] <- c('Lower', 'Median', 'Upper')

p2 <- ggplot2::ggplot(Ribbon, ggplot2::aes(.data$Year)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$Lower, ymax = .data$Upper, fill = .data$Group), alpha = 0.2) +
  ggplot2::geom_line(ggplot2::aes(y = .data$Median, colour = .data$Group)) +
  ggplot2::geom_point(data = ObsDF, ggplot2::aes(y = .data$Value)) +
  ggplot2::scale_colour_manual(values = c('Collapse before MPs start' = 'firebrick', 'Other simulations' = 'grey40')) +
  ggplot2::scale_fill_manual(values = c('Collapse before MPs start' = 'firebrick', 'Other simulations' = 'grey40')) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::labs(x = NULL, y = 'LL1 CPUE (scaled to mean 1)', colour = NULL, fill = NULL) +
  ggplot2::theme_bw() + ggplot2::theme(legend.position = 'bottom')
ggplot2::ggsave(file.path(OutDir, 'CPUE_fit.png'), p2, width = 8, height = 4.5)
