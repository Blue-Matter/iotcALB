# Diagnostics of the surplus production model used by the SP CMPs
#
# 1. Fit to the historical data (2000-2023) with the informative priors and
#    fixed Shape/Depletion (SPSettings()), and with the MSEtool default (vague)
#    FMSY prior, for comparison
#
# 2. Convergence of the model fits in every management cycle of the tuned SP
#    CMPs (base case OM, nSimTune simulations), and the estimated vs OM stock
#    status. A management cycle where the fit does not converge holds the
#    previous TAC (OnFail = 'hold'), so a low convergence rate makes the SP
#    results uninformative about the model-based MP

source('analysis/02-define-CMPs.R')

MinConvergence <- 0.95   # minimum acceptable fraction of converged fits
nWorkers       <- 40
FigDir         <- 'figures/CMPs'
dir.create(FigDir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Historical data ----
Data       <- MSEtool::Data(Hist)[[1]][[1]]
Data@Years <- as.numeric(rownames(Data@Landings@Value))  # include the 2021-2023 data

FitArgs <- list(Data        = Data,
                Indices     = IndexName,
                IndexSource = 'Survey',
                Shape       = SPShape,
                Depletion   = SPSet$Depletion,
                Uncertainty = TRUE)
Fits <- list(Informative = do.call(MSEtool::FitSP,
                                   c(FitArgs, list(Priors = SPSet$Priors))),
             Vague       = do.call(MSEtool::FitSP,
                                   c(FitArgs, list(Priors = list(FMSY = c(0.2, 2))))
                                   )
             )

HistFitTable <- do.call(rbind, lapply(names(Fits), \(nm) {
  f <- Fits[[nm]]
  data.frame(Priors = nm, Converged = f$Converged, FMSY = f$FMSY, MSY = f$MSY, BMSY = f$BMSY,
             B_BMSY = f$Terminal[['B_BMSY']], F_FMSY = f$Terminal[['F_FMSY']],
             SE_logB_BMSY = f$SE[['logB_BMSY']], SE_logFMSY = f$SE[['logFMSY']])
}))
print(HistFitTable, digits = 3, row.names = FALSE)

# ---- 2. Tuned SP CMPs ----
TunedCMPs <- readRDS(file.path(TuneDir, 'TunedCMPs.rds'))
SPCMPs    <- TunedCMPs[CMPs$MP[match(names(TunedCMPs), CMPs$Name)] %in% c('SP', 'SPB')]
if (!length(SPCMPs))
  cli::cli_abort('No tuned SP CMPs in {.file TunedCMPs.rds}; run 02b-tune-CMPs.R.')

MSEtool::SetupParallel(workers = nWorkers)
MSE <- MSEtool::Project(HistTune, MPs = SPCMPs, parallel = TRUE)
MSEtool::DisableParallel()

Est <- MSEtool::SPEstimates(MSE)

## OM total biomass relative to total BMSY at the start of the year after the
# last data year (the estimated B_BMSY), for comparison with the SP estimate
BMSYTotal <- apply(MSEtool::BMSY(Hist), 1, sum)
OMB <- MSEtool::Biomass(MSE, df = TRUE, Reduce = FALSE)
OMB <- stats::aggregate(Value ~ Sim + Year + MP, data = as.data.frame(OMB), FUN = sum)
OMB <- OMB[OMB$Year == floor(OMB$Year), ]  # first season of each year
OMB$OM_B_BMSY <- OMB$Value / BMSYTotal[as.integer(as.character(OMB$Sim))]
Est$OM_B_BMSY <- OMB$OM_B_BMSY[match(paste(Est$MP, Est$Sim, Est$LastDataYear + 1),
                                     paste(OMB$MP, OMB$Sim, OMB$Year))]

Convergence <- do.call(rbind, lapply(split(Est, Est$MP), \(d) {
  data.frame(MP = d$MP[1], nFits = nrow(d), Converged = mean(d$Converged),
             SimsAllConverged = mean(tapply(d$Converged, d$Sim, all)),
             MedianRelErr_B_BMSY = stats::median((d$B_BMSY - d$OM_B_BMSY) / d$OM_B_BMSY, na.rm = TRUE),
             Cor_B_BMSY = suppressWarnings(stats::cor(log(d$B_BMSY), log(d$OM_B_BMSY), use = 'complete')))
}))
print(Convergence, digits = 3, row.names = FALSE)
print(table(Est$MP, Est$Message))

utils::write.csv(Convergence, file.path(TuneDir, 'SP_Convergence.csv'), row.names = FALSE)
saveRDS(list(HistFits = HistFitTable, Estimates = Est, Convergence = Convergence),
        file.path(TuneDir, 'SP_Diagnostics.rds'))

Low <- Convergence$MP[Convergence$Converged < MinConvergence]
if (length(Low)) {
  cli::cli_alert_warning("Convergence below {MinConvergence * 100}% for {.val {Low}}.")
} else {
  cli::cli_alert_success("Convergence of every SP CMP is at least {MinConvergence * 100}%.")
}

## Estimated vs OM B/BMSY
p <- ggplot2::ggplot(Est[Est$Converged, ], ggplot2::aes(x = .data$OM_B_BMSY, y = .data$B_BMSY)) +
  ggplot2::geom_point(alpha = 0.3) +
  ggplot2::geom_abline(linetype = 'dashed', colour = 'firebrick') +
  ggplot2::facet_wrap(~MP) +
  ggplot2::labs(x = 'OM B/BMSY (total biomass)', y = 'Estimated B/BMSY') +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(FigDir, 'SP_B_BMSY_est_vs_OM.png'), p, width = 8, height = 4)
