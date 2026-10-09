# Define the candidate management procedures (CMPs) and the tuning specifications
#
# CMPs are tuned on the base case OM (OM5b) to maximise mean catch over the
# years the MP is active, at two tuning levels:
#   T1: P(Kobe green) (PGK) >= 0.6 over the last 15 years of the projections
#   T2: T1, and P(SB > 0.4 SBMSY in every year 2034-2058) >= 0.9 (safety; the
#       first 5 MP years are excluded so the stock can recover from the interim
#       catches: with no fishing, P = 0.94 over 2029-2058 and 0.98 over 2034-2058),
#       with TAC decreases of up to 30% (increases up to 15%)
#
# The MP settings (ConfigCandidates) are selected for each MP
# in 02a-select-configs.R; the CMPs are tuned in 02b-tune-CMPs.R.
library(iotcALB)

source('analysis/00-specifications.R')

BaseOM   <- 'OM5b'
HistFile <- file.path('objects/Hist', paste0(BaseOM, '.hist'))
Hist     <- readRDS(HistFile)

# ---- Reference Levels ----
# IOTC Resolution 15/10: target SBMSY and FMSY, limit 0.4 SBMSY for albacore
LRP       <- 0.4                 # limit reference point: SB / SBMSY
PGKMin    <- 0.6                 # minimum P(Kobe green) when tuning
PSafety   <- 0.9                 # minimum P(SB > LRP in every safety year) when tuning (T2)
TACChange <- 0.15                # maximum TAC change (up or down) between management cycles

# Maximum TAC decrease of each tuning level.
TACChangeDown <- c(T1 = 0.15, T2 = 0.30)

TuneYears   <- MSEtool::PMYears(Hist@OM, 'last', n = 15)  # 2044-2058
CatchYears  <- MSEtool::PMYears(Hist@OM)                   # 2029-2058
SafetyYears <- MSEtool::PMYears(Hist@OM, Skip = 5)         # 2034-2058

# ---- Indices ----
# LL1 (NW) CPUE: the index the base case OM is conditioned on, used in the
# earlier MPs.
IndexName <- 'LL1'

## Index level at SBMSY: mean annual index / median mean SB/SBMSY over
# IndexCalYears (base case OM), so the HCR control points of the index MPs are
# in units of SB/SBMSY
IndexCalYears <- 2016:2020
IndexSBMSY    <- IndexAtSBMSY(Hist, Index = IndexName, Years = IndexCalYears)

# ---- Surplus Production Model ----
# Schaefer model (shape n = 2, BMSY/K = 0.5), with informative priors on FMSY
# and MSY and the initial depletion (2000) fixed at the base case OM value
SPSet   <- SPSettings(Hist)
SPShape <- 2
BMSY_K  <- function(Shape) if (abs(Shape - 1) < 1e-6) exp(-1) else Shape^(1 / (1 - Shape))

# Buffer HCR of the earlier work (BUF): limit 0.25, buffers 0.35 and 0.60, and
# slope ratio 0.10 in units of the LL1 index, scaled so that the limit is the LRP
BUFScale <- LRP / 0.25
BUFArgs  <- list(Limit = 0.25 * BUFScale, LowerBuffer = 0.35 * BUFScale,
                 UpperBuffer = 0.60 * BUFScale, SlopeRatio = 0.10, nYears = 4, Metric = 'wmean')

cli::cli_inform(c(
  'i' = 'Index at SBMSY ({IndexName}, {min(IndexCalYears)}-{max(IndexCalYears)}): {signif(IndexSBMSY, 3)}',
  'i' = 'SP: Schaefer (BMSY/K = {BMSY_K(SPShape)}), Depletion = {signif(SPSet$Depletion, 3)}',
  'i' = 'SP priors: FMSY = {signif(SPSet$Priors$FMSY[1], 3)} (CV {signif(SPSet$Priors$FMSY[2], 2)}), MSY = {round(SPSet$Priors$MSY[1])} (CV {signif(SPSet$Priors$MSY[2], 2)})'
))

# ---- Tuning Specifications ----
TuneDir <- file.path(normalizePath('objects'), 'Tune')
dir.create(TuneDir, showWarnings = FALSE)

nSimTune <- 80   # simulations used for tuning (500 for the final projections)

# Simulations used for tuning: a stratified random sample, so that the tuning
# set has the same distribution of stock productivity and status as all
# simulations.
TuneSeed   <- 1
nMSYGroups <- 8
TuneSims <- local({
  MSY <- apply(MSEtool::MSYLandings(Hist), 1, sum) + apply(MSEtool::MSYDiscards(Hist), 1, sum)
  Dep <- MSEtool::SB_SB0(Hist, df = TRUE, Reduce = FALSE)
  Dep <- Dep[Dep$Stock == 'Female' & floor(Dep$Year) == OMSpecs$CurrentYear, ]
  Dep <- tapply(Dep$Value, as.integer(as.character(Dep$Sim)), mean)
  Sims <- as.integer(names(Dep))
  MSYGroup <- cut(rank(MSY[Sims], ties.method = 'first'), nMSYGroups, labels = FALSE)
  DepGroup <- stats::ave(Dep, MSYGroup,
                         FUN = \(x) cut(rank(x, ties.method = 'first'), nSimTune / nMSYGroups, labels = FALSE))
  Cells <- split(Sims, paste(MSYGroup, DepGroup))
  withr::with_seed(TuneSeed, sort(unname(vapply(Cells, \(s) s[sample.int(length(s), 1)], integer(1)))))
})
HistTune <- MSEtool::Subset(Hist, Sims = TuneSims)

# PM_Removals: mean annual catch (t; removals summed over the seasons of each
# calendar year) over the years the MP is active
Objective <- MSEtool::TuneObjective(MSEtool::PM_Removals,
                                    Years = CatchYears,
                                    Name  = 'MeanCatch')

# Constraints of each tuning level:
# - PGK, PM_Status: P(SB > SBMSY and F < FMSY) over the last 15 projection
#   years (proportion of simulations and years)
# - Safety, PM_Safety: proportion of simulations in which SB > LRP in every
#   year of SafetyYears (a single year below the LRP is a breach)
PGKConstraint <- MSEtool::TuneConstraint(MSEtool::PM_Status,
                                         Definition = 'SBiomass',
                                         Years      = TuneYears,
                                         Min        = PGKMin,
                                         Name       = 'PGK')

SafetyConstraint <- MSEtool::TuneConstraint(MSEtool::PM_Safety,
                                            Lim        = LRP,
                                            Definition = 'SBiomass',
                                            Years      = SafetyYears,
                                            Min        = PSafety,
                                            Name       = 'Safety')

TuneLevels <- list(T1 = list(PGKConstraint),
                   T2 = list(PGKConstraint, SafetyConstraint))
TuneLevelLabel <- c(T1 = sprintf('P(Kobe green, %d-%d) >= %.2g',
                                 min(TuneYears), max(TuneYears), PGKMin),
                    T2 = sprintf('P(Kobe green, %d-%d) >= %.2g and P(SB > %.1f SBMSY in every year %d-%d) >= %.2g, with TAC decreases of up to %d%%',
                                 min(TuneYears), max(TuneYears), PGKMin, LRP, min(SafetyYears),
                                 max(SafetyYears), PSafety, round(100 * TACChangeDown[['T2']])))

# tuning level used to select the MP settings in 02a-select-configs.R
ConfigTuning <- 'T1'

# initial range of tunepar; extended by a factor of 3 up to MaxExpand times in
# 02b-tune-CMPs.R if needed
TuneInterval <- c(0.1, 3)

# ---- MP Settings Selected by 02a-select-configs.R ----
# Candidate values of the settings not set by tuning, for the single-index
# CMPs.
#
# - RecentYears: years averaged for the current index (the earlier MP used 4)
# - Smooth: loess-smooth the index before averaging
# - Buffer: lower and upper buffers of the buffer HCR (SB/SBMSY)
# - HCRThreshold: B/BMSY below which F is reduced (SP)

ConfigCandidates <- list(IR  = list(RecentYears  = c(1, 3),
                                    Smooth       = c(TRUE, FALSE)),
                         ITB = list(RecentYears  = c(1, 3),
                                    Buffer       = list(c(0.8, 1.2), c(1.0, 1.5))),
                         SP  = list(HCRThreshold = c(0.8, 1)))

DefaultConfig <- list(IR  = list(RecentYears     = 3,
                                Smooth           = FALSE),
                      ITB = list(RecentYears     = 3,
                                 Buffer          = c(0.8, 1.2)),
                      SP  = list(HCRThreshold    = 1))

# CMP types: the settings each uses (NA: none; BUF replicates the earlier MP)
MPTypes <- data.frame(MP           = c('IR', 'ITB', 'BUF', 'SP'),
                      ConfigFamily = c('IR', 'ITB', NA, 'SP'))

ConfigsFile     <- file.path(TuneDir, 'Configs.rds')
SelectedConfigs <- if (file.exists(ConfigsFile)) readRDS(ConfigsFile) else list()
SelectedConfigs <- lapply(stats::setNames(nm = names(ConfigCandidates)),
                          \(mp) utils::modifyList(DefaultConfig[[mp]], SelectedConfigs[[mp]] %||% list()))

ConfigLabel <- function(Config) {
  paste(names(Config), vapply(Config, \(x) paste(x, collapse = '-'), character(1)), sep = ' = ', collapse = ', ')
}
for (mp in names(ConfigCandidates))
  cli::cli_inform("{mp} settings: {ConfigLabel(SelectedConfigs[[mp]])}")

# Settings of a CMP type, from its config family
TypeConfig <- function(MP) {
  Family <- MPTypes$ConfigFamily[match(MP, MPTypes$MP)]
  if (is.na(Family)) list() else SelectedConfigs[[Family]]
}

# ---- Candidate MPs ----
#
# - IR: IndexRate. TAC = catch/index rate (calibrated over the last 2 historical
#   years, x tunepar) x current index. The rate is reduced linearly from the
#   full rate at SBMSY to zero at the LRP (index / IndexSBMSY)
#
# - ITB: IndexTarget with the buffer HCR of the earlier albacore MSE
#   (BufferHCR()). TAC = reference catch (mean of the last 5 historical years,
#   x tunepar) x HCR(index / IndexSBMSY): the reference catch between the
#   buffers, decreasing below the lower buffer (more rapidly below the LRP),
#   and increasing slowly above the upper buffer
#
# - BUF: the CPUE + buffer HCR MP of the earlier albacore MSE (BufferMP()): the
#   HCR multiplies the previous TAC; tunepar moves the upper buffer
#
# - SP: SurplusProduction. F = tunepar x FMSY (estimated), reduced linearly
#   from the threshold (B/BMSY) to zero at the LRP; TAC is the mean catch at
#   that F over the 3-year management interval

# MP arguments of a configuration
ConfigArgs <- function(MP, Config) {
  Family <- MPTypes$ConfigFamily[match(MP, MPTypes$MP)]
  if (is.na(Family)) return(list())
  Config <- utils::modifyList(DefaultConfig[[Family]], Config)
  switch(MP,
         IR  = list(RecentYears           = Config$RecentYears,
                    Smooth                = Config$Smooth,
                    HCRControlPointsIndex = c(LRP, 1),
                    HCRControlPointsRate  = c(0, 1)),
         ITB = c(list(RecentYears = Config$RecentYears, Smooth = FALSE),
                 BufferHCR(Config$Buffer, Limit = LRP)),
         SP  = list(Shape                   = SPShape,
                    HCRControlPointsBiomass = c(LRP, Config$HCRThreshold),
                    HCRControlPointsRate    = c(0, 1)))
}

MakeCMP <- function(MP, TACChange, TACChangeDown = TACChange, Config = TypeConfig(MP)) {
  Delta     <- c(0.01, TACChange)
  DeltaDown <- c(0.01, TACChangeDown)
  SPBase <- list(MSEtool::SurplusProduction, IndexSource = 'Survey',
                 Depletion = SPSet$Depletion, Priors = SPSet$Priors, OnFail = 'hold')
  Base <- switch(MP,
                 IR  = list(MSEtool::IndexRate,
                            Indices     = IndexName,
                            IndexSource = 'Survey',
                            IndexTarget = unname(IndexSBMSY[IndexName])),
                 ITB = list(MSEtool::IndexTarget,
                            Indices     = IndexName,
                            IndexSource = 'Survey',
                            IndexTarget = unname(IndexSBMSY[IndexName])),
                 BUF = c(list(BufferMP,
                              Indices     = IndexName,
                              IndexSource = 'Survey',
                              IndexRef    = unname(IndexSBMSY[IndexName])), BUFArgs),
                 SP  = c(SPBase, list(Indices = IndexName)))

  if (MP == 'BUF') Delta[1] <- DeltaDown[1] <- 0
  do.call(MSEtool::SetMPArgs, c(Base, list(DeltaDown = DeltaDown, DeltaUp = Delta), ConfigArgs(MP, Config)))
}

# Each CMP type with its selected settings, at each tuning level (e.g. IR_T1).
# TACChange: maximum increase; TACChangeDown: maximum decrease
CMPs               <- expand.grid(MP     = MPTypes$MP,
                                  Tuning = names(TuneLevels),
                                  stringsAsFactors = FALSE)
CMPs$TACChange     <- TACChange
CMPs$TACChangeDown <- unname(TACChangeDown[CMPs$Tuning])
CMPs$Name          <- paste(CMPs$MP, CMPs$Tuning, sep = '_')
CMPs               <- CMPs[order(CMPs$Tuning, match(CMPs$MP, MPTypes$MP)), ]
CMPs$ConfigFamily  <- MPTypes$ConfigFamily[match(CMPs$MP, MPTypes$MP)]
rownames(CMPs)     <- NULL

CMPs$Fn            <- purrr::pmap(list(CMPs$MP, CMPs$TACChange, CMPs$TACChangeDown),
                                  MakeCMP)
CMPs$Constraints   <- unname(TuneLevels[CMPs$Tuning])
