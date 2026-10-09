# Build the Slick object for results visualisation
#
# Reads the MSE results of 03-run-projections.R and builds the Slick object with
# MSEtool::MSE2Slick(), from the MSEtool performance metric (PM) functions and
# annual time series. Also saves figures and tables for reports (FigDir).
library(iotcALB)

MSEDir  <- 'objects/MSE'
TuneDir <- 'objects/Tune'
FigDir  <- 'figures/Slick'

# Reference levels (IOTC Resolution 15/10 for albacore: target SBMSY and FMSY,
# limits 0.4 SBMSY and 1.4 FMSY) and the tuning objective
SBLim  <- 0.4
FLim   <- 1.4
PGKMin <- 0.6        # tuning (T1, T2): P(Kobe green) over the last 15 projection years
PSafety <- 0.9       # tuning (T2): P(SB > SBLim SBMSY in every year 2034-2058)

# CMPs in the Slick object (NULL for all tuned CMPs). The reference MPs
# projected by 03-run-projections.R (NoFishing) are not included
SelectedCMPs <- NULL

# OMs in the Slick object: only the Base Case OM is projected (see
# 03-run-projections.R); an older .mse of another OM is not included
SlickOMs <- 'OM5b'

TunedCMPs   <- readRDS(file.path(TuneDir, 'TunedCMPs.rds'))
TuneSummary <- utils::read.csv(file.path(TuneDir, 'Summary.csv'))

# ---- OMs with MSE results ----
OM_Design <- data.frame(OM          = c('OM5b', 'OM5a', 'OM6b'),
                        Type        = c('Base Case', 'Robustness', 'Robustness'),
                        Description = c('Base Case: LL1 (NW) CPUE as the index of abundance',
                                        'Robustness: SW CPUE as the index of abundance',
                                        'Robustness: 1% annual increase in longline catchability')) |>
  dplyr::mutate(MSEFile = file.path(MSEDir, paste0(OM, '.mse'))) |>
  dplyr::filter(OM %in% SlickOMs, file.exists(MSEFile))

if (!nrow(OM_Design))
  cli::cli_abort("No MSE results in {.file {MSEDir}}; run 03-run-projections.R.")
OM_Design$OM_Name <- paste0(OM_Design$OM, ': ', OM_Design$Type)
nOM <- nrow(OM_Design)

# ---- MSE results (selected CMPs only) ----
MP_all <- names(TunedCMPs)
if (!is.null(SelectedCMPs)) {
  if (length(setdiff(SelectedCMPs, MP_all)))
    cli::cli_abort("Not tuned CMPs: {.val {setdiff(SelectedCMPs, MP_all)}}.")
  MP_all <- SelectedCMPs
}

MSEList <- lapply(seq_len(nOM), \(i) {
  MSE <- readRDS(OM_Design$MSEFile[i])
  Missing <- setdiff(MP_all, names(MSE@MPs))
  if (length(Missing))
    cli::cli_abort("OM {OM_Design$OM[i]} has no results for {.val {Missing}}; re-run 03-run-projections.R.")
  MSEtool::Subset(MSE, MPs = MP_all)
})
names(MSEList) <- OM_Design$OM_Name
MSE1 <- MSEList[[1]]

# ---- MP metadata ----
MP_all <- names(MSE1@MPs)   # MP order of the MSE objects
nMP    <- length(MP_all)
MPInfo <- TuneSummary[match(MP_all, TuneSummary$Name), ]

MP_Type <- c(IR  = 'Index Rate',
             ITB = 'Index Target Buffer',
             BUF = 'Buffer MP',
             SP  = 'Surplus Production')

MP_Class <- c(IR = 'Model-free', ITB = 'Model-free', BUF = 'Model-free', SP = 'Model-based')

TuningName <- c(T1 = 'T1: P(Kobe green) only',
                T2 = 'T2: P(Kobe green) and Safety')

MP_TypeDescription <- c(
  IR  = 'The TAC is a catch/index rate multiplied by the current LL1 (NW) CPUE. The rate is reduced linearly from SBMSY to zero at the limit reference point (0.4 SBMSY), with stock status estimated as the index relative to its level at SBMSY.',
  ITB = 'The TAC is a reference catch scaled by a buffer HCR on the LL1 (NW) CPUE: the reference catch while the index is between the buffers, decreasing below the lower buffer (more rapidly below the limit reference point), and increasing slowly above the upper buffer.',
  BUF = 'A replica of the CPUE + buffer HCR MP of the earlier albacore MSE work. The TAC is the previous TAC multiplied by a buffer HCR on the 4-year weighted mean of the LL1 (NW) CPUE: unchanged while the index is between the buffers, decreased below the lower buffer (more rapidly below the limit), and increased slowly above the upper buffer. Tuned by the position of the upper buffer.',
  SP  = 'A Schaefer surplus production model is fitted to the catch and the LL1 (NW) CPUE, with informative priors on FMSY and MSY. F is a multiple of the estimated FMSY, reduced linearly from the HCR threshold to zero at 0.4 B/BMSY.'
)

MP_Code  <- MP_all
MP_Label <- sprintf('%s (%s)', MP_Type[MPInfo$MP], MPInfo$Tuning)
MP_Description <- sprintf(
  '%s TAC increases are limited to %d%% and decreases to %d%% between management cycles.',
  MP_TypeDescription[MPInfo$MP], round(100 * MPInfo$TACChange), round(100 * MPInfo$TACChangeDown))


MP_Description[5] <- 'Same as IR_T1 but TAC increases are limited to 15% and decreases to 30% between management cycles.'
MP_Description[6] <- 'Same as ITB_T1 but TAC increases are limited to 15% and decreases to 30% between management cycles.'
MP_Description[7] <- 'Same as SP_T1 but TAC increases are limited to 15% and decreases to 30% between management cycles.'

# , MPInfo$Config,
  # MPInfo$Tuning, MPInfo$TuningLabel, MPInfo$tunepar)

# ---- Performance indicators ----
# IOTC MSE categories: Status, Safety, Yield, Stability. The probability that
# the TAC changes by no more than 15% is not included: every CMP limits TAC increases to 15%
PM_Codes <- c('PGK_Tune', 'PGK_All', 'PGK_Short', 'PGK_Medium', 'PGK_Long',
              'PRed_All', 'SB_SBMSY_All', 'F_FMSY_All', 'SB_SBMSY_Min',
              'Safety_SB', 'PLim_SB', 'PLim_F',
              'Yield_All', 'Yield_Short', 'Yield_Medium', 'Yield_Long',
              'AAV_TAC', 'P_TACLimited')
nPI_pm <- length(PM_Codes)

PM_Category <- c(PGK_Tune = 'Status', PGK_All = 'Status', PGK_Short = 'Status', PGK_Medium = 'Status',
                 PGK_Long = 'Status', PRed_All = 'Status', SB_SBMSY_All = 'Status', F_FMSY_All = 'Status',
                 SB_SBMSY_Min = 'Status', Safety_SB = 'Safety', PLim_SB = 'Safety', PLim_F = 'Safety',
                 Yield_All = 'Yield', Yield_Short = 'Yield', Yield_Medium = 'Yield', Yield_Long = 'Yield',
                 AAV_TAC = 'Stability', P_TACLimited = 'Stability')

PM_Label <- c(
  PGK_Tune     = 'P(Kobe green) 2044-2058',
  PGK_All      = 'P(Kobe green)',
  PGK_Short    = 'P(Kobe green) (years 1-10)',
  PGK_Medium   = 'P(Kobe green) (years 11-20)',
  PGK_Long     = 'P(Kobe green) (years 21-30)',
  PRed_All     = 'P(Kobe red)',
  SB_SBMSY_All = 'Mean SB/SBMSY',
  F_FMSY_All   = 'Mean F/FMSY',
  SB_SBMSY_Min = 'Minimum SB/SBMSY',
  Safety_SB    = paste0('P(SB > ', SBLim, ' SBMSY every year) 2034-2058'),
  PLim_SB      = paste0('P(SB > ', SBLim, ' SBMSY)'),
  PLim_F       = paste0('P(F < ', FLim, ' FMSY)'),
  Yield_All    = 'Mean catch',
  Yield_Short  = 'Mean catch (years 1-10)',
  Yield_Medium = 'Mean catch (years 11-20)',
  Yield_Long   = 'Mean catch (years 21-30)',
  AAV_TAC      = 'Mean TAC change',
  P_TACLimited = 'P(TAC change limited)',
  Rel_Yield    = 'Relative mean catch'
)

PM_Description <- c(
  PGK_Tune     = 'Probability of being in the green quadrant of the Kobe plot (SB > SBMSY and F < FMSY) over 2044-2058, the last 15 projection years (the tuning objective).',
  PGK_All      = 'Probability of being in the green quadrant of the Kobe plot (SB > SBMSY and F < FMSY) over all years the MP is active (2029-2058).',
  PGK_Short    = 'Probability of being in the green quadrant of the Kobe plot over the first 10 years the MP is active (2029-2038).',
  PGK_Medium   = 'Probability of being in the green quadrant of the Kobe plot over the middle 10 years the MP is active (2039-2048).',
  PGK_Long     = 'Probability of being in the green quadrant of the Kobe plot over the last 10 years the MP is active (2049-2058).',
  PRed_All     = 'Probability of being in the red quadrant of the Kobe plot (SB < SBMSY and F > FMSY) over all years the MP is active.',
  SB_SBMSY_All = 'Mean annual female spawning biomass relative to SBMSY over all years the MP is active.',
  F_FMSY_All   = 'Mean annual fishing mortality relative to FMSY over all years the MP is active.',
  SB_SBMSY_Min = 'Lowest annual SB/SBMSY reached in the years the MP is active.',
  Safety_SB    = paste0('Probability that spawning biomass is above the limit reference point (', SBLim, ' SBMSY) in every year 2034-2058, i.e. the proportion of simulations that never fall below the limit (the safety tuning constraint of T2). The first 5 years the MP is active are excluded so the stock can recover from the interim catches.'),
  PLim_SB      = paste0('Probability that spawning biomass is above the limit reference point (', SBLim, ' SBMSY) over all years the MP is active (proportion of simulations and years).'),
  PLim_F       = paste0('Probability that fishing mortality is below the limit reference point (', FLim, ' FMSY) over all years the MP is active.'),
  Yield_All    = 'Mean annual catch (t) over all years the MP is active (2029-2058).',
  Yield_Short  = 'Mean annual catch (t) over the first 10 years the MP is active (2029-2038).',
  Yield_Medium = 'Mean annual catch (t) over the middle 10 years the MP is active (2039-2048).',
  Yield_Long   = 'Mean annual catch (t) over the last 10 years the MP is active (2049-2058).',
  AAV_TAC      = 'Mean absolute proportional change in the TAC between management cycles.',
  P_TACLimited = 'Probability that the change in the TAC between management cycles is at the maximum allowed by the MP (the TAC change limit is binding).',
  Rel_Yield    = 'Mean annual catch (all years) relative to the highest-yielding MP in the same OM.'
)

Prob01_Codes   <- c('PGK_Tune', 'PGK_All', 'PGK_Short', 'PGK_Medium', 'PGK_Long', 'PRed_All',
                    'Safety_SB', 'PLim_SB', 'PLim_F', 'P_TACLimited')
Headline_Codes <- c('PGK_Tune', 'PGK_All', 'PRed_All', 'Safety_SB', 'PLim_SB', 'PLim_F', 'Yield_All',
                    'Yield_Short', 'AAV_TAC', 'P_TACLimited')
Spider_Codes   <- c('PGK_Tune', 'PGK_All', 'Safety_SB', 'PLim_F', 'Rel_Yield')
Tradeoff_Codes <- c('Yield_All', 'PGK_All', 'PGK_Tune', 'Safety_SB', 'PLim_SB', 'SB_SBMSY_All', 'AAV_TAC')


PMSpec <- function(Code, PM, Args = list())
  list(PM = PM, Args = Args, Code = Code, Label = PM_Label[[Code]], Description = PM_Description[[Code]])

PMSpecs <- list(
  PGK_Tune     = PMSpec('PGK_Tune', MSEtool::PM_Status, list(Years = MSEtool::PMYears(MSE1, 'last', n = 15))),
  PGK_All      = PMSpec('PGK_All', MSEtool::PM_Status),
  PGK_Short    = PMSpec('PGK_Short', MSEtool::PM_Status, list(Years = MSEtool::PMYears(MSE1, 'first', n = 10))),
  PGK_Medium   = PMSpec('PGK_Medium', MSEtool::PM_Status, list(Years = MSEtool::PMYears(MSE1, 'middle', n = 10))),
  PGK_Long     = PMSpec('PGK_Long', MSEtool::PM_Status, list(Years = MSEtool::PMYears(MSE1, 'last', n = 10))),
  PRed_All     = PMSpec('PRed_All', MSEtool::PM_Red),
  SB_SBMSY_All = PMSpec('SB_SBMSY_All', MSEtool::PM_SBSBMSY, list(Ref = NULL)),
  F_FMSY_All   = PMSpec('F_FMSY_All', MSEtool::PM_FFMSY, list(Ref = NULL)),
  SB_SBMSY_Min = PMSpec('SB_SBMSY_Min', MSEtool::PM_MinStatus),
  Safety_SB    = PMSpec('Safety_SB', MSEtool::PM_Safety, list(Lim = SBLim, Years = MSEtool::PMYears(MSE1, Skip = 5))),
  PLim_SB      = PMSpec('PLim_SB', MSEtool::PM_SBSBlim, list(Lim = SBLim)),
  PLim_F       = PMSpec('PLim_F', MSEtool::PM_FFMSY, list(Ref = FLim)),
  Yield_All    = PMSpec('Yield_All', MSEtool::PM_Removals),
  Yield_Short  = PMSpec('Yield_Short', MSEtool::PM_Removals, list(Years = MSEtool::PMYears(MSE1, 'first', n = 10))),
  Yield_Medium = PMSpec('Yield_Medium', MSEtool::PM_Removals, list(Years = MSEtool::PMYears(MSE1, 'middle', n = 10))),
  Yield_Long   = PMSpec('Yield_Long', MSEtool::PM_Removals, list(Years = MSEtool::PMYears(MSE1, 'last', n = 10))),
  AAV_TAC      = PMSpec('AAV_TAC', MSEtool::PM_AAVY, list(Type = 'TAC', IncludeFirst = TRUE)),
  P_TACLimited = PMSpec('P_TACLimited', MSEtool::PM_TACLimited),
  Rel_Yield    = PMSpec('Rel_Yield', MSEtool::PM_Removals)
)
Quilt_Codes <- setdiff(PM_Codes, c('SB_SBMSY_Min', 'PRed_All'))

# Time series: catch is the removals summed over the seasons of each year, and
# SB/SB0 the mean over the seasons of each year
TS_Codes <- c('SB_SBMSY', 'F_FMSY', 'SB_SB0', 'Removals')
TS_Label <- c(SB_SBMSY = 'SB/SBMSY', F_FMSY = 'F/FMSY', SB_SB0 = 'SB/SB0', Removals = 'Catch (t)')
TS_Description <- c(
  SB_SBMSY = 'Female spawning biomass relative to SBMSY.',
  F_FMSY   = 'Annual apical fishing mortality relative to FMSY.',
  SB_SB0   = 'Female spawning biomass relative to the equilibrium unfished spawning biomass.',
  Removals = 'Total annual catch (t) of all fleets.'
)
TS_Target <- c(SB_SBMSY = 1,     F_FMSY = 1,    SB_SB0 = NA, Removals = NA)
TS_Limit  <- c(SB_SBMSY = SBLim, F_FMSY = FLim, SB_SB0 = NA, Removals = NA)
TS_File   <- c(SB_SBMSY = 'SB_SBMSY', F_FMSY = 'F_FMSY', SB_SB0 = 'SB_SB0', Removals = 'Yield')  # figure file names

CategoryPreset <- function(Codes) {
  c(list(Summary = match(intersect(Headline_Codes, Codes), Codes),
         All     = seq_along(Codes)),
    lapply(split(seq_along(Codes), PM_Category[Codes]), identity),
    list(ShortTerm = match(grep('_Short$',  Codes, value = TRUE), Codes),
         MediumTerm = match(grep('_Medium$', Codes, value = TRUE), Codes),
         LongTerm = match(grep('_Long$',   Codes, value = TRUE), Codes)))
}

# ---- Assemble and save ----
Introduction <- "
Preliminary results of the management strategy evaluation (MSE) for Indian Ocean albacore (*Thunnus alalunga*).

The results are for the Base Case operating model (OM), conditioned with the Approximate Bayesian Computation (ABC) approach in the earlier phase of the albacore MSE work on data up to 2020, with the reported catches for 2021-2023 and the CMPs using the catch and CPUE data of the 2025 stock assessment. The OMs may be reconditioned to the updated data, and these results should be considered preliminary until then.

The performance indicators follow the IOTC categories of Status (Kobe quadrants, SB/SBMSY, F/FMSY), Safety (probability of being above the interim limit reference points, 0.4 SBMSY and 1.4 FMSY), Yield (mean catch), and Stability (TAC variability).

**These results are preliminary and do not necessarily reflect the views of the IOTC or its members.**
"

slick <- MSEtool::MSE2Slick(
  MSEList,
  Title          = 'Indian Ocean Albacore',
  Author         = 'Adrian Hordyk',
  Email          = 'adrian@bluematterscience.com',
  Introduction   = Introduction,
  PMs            = unname(PMSpecs[PM_Codes]),
  QuiltPMs       = unname(PMSpecs[Quilt_Codes]),
  SpiderPMs      = unname(PMSpecs[Spider_Codes]),
  TradeoffPMs    = unname(PMSpecs[Tradeoff_Codes]),
  TimeseriesCode        = TS_Codes,
  TimeseriesLabel       = unname(TS_Label[TS_Codes]),
  TimeseriesDescription = unname(TS_Description[TS_Codes]),
  TimeseriesTarget      = TS_Target[!is.na(TS_Target)],
  TimeseriesLimit       = TS_Limit[!is.na(TS_Limit)],
  SeasonalBasis  = 'Mean',
  KobeLimit      = c(SBLim, FLim),
  MPCode         = MP_all,
  MPLabel        = MP_Label,
  MPDescription  = MP_Description,
  Design         = data.frame(OM = OM_Design$OM),
  AutoPreset     = FALSE,
  OMsPreset      = c(if ('OM5b' %in% OM_Design$OM) list(`Base Case` = list(match('OM5b', OM_Design$OM))),
                     if (any(OM_Design$Type == 'Robustness')) list(Robustness = list(which(OM_Design$Type == 'Robustness'))),
                     list(All = list(seq_len(nOM)))),
  MPsPreset      = c(list(All = seq_len(nMP)),
                     split(seq_len(nMP), TuningName[MPInfo$Tuning]),
                     split(seq_len(nMP), MP_Class[MPInfo$MP])),
  BoxplotPreset  = CategoryPreset(PM_Codes),
  QuiltPreset    = CategoryPreset(Quilt_Codes),
  TradeoffPreset = list(All            = seq_along(Tradeoff_Codes),
                        CatchStatus    = match(c('Yield_All', 'PGK_All'), Tradeoff_Codes),
                        CatchSafety    = match(c('Yield_All', 'Safety_SB'), Tradeoff_Codes),
                        CatchStability = match(c('Yield_All', 'AAV_TAC'), Tradeoff_Codes))
)
MPStartYear <- MSE1@OM@MPStartYear
rm(MSEList, MSE1)
gc()

# OM descriptions
OMs_obj <- Slick::OMs(slick)
Slick::Factors(OMs_obj) <- data.frame(Factor = 'OM', Level = OM_Design$OM, Description = OM_Design$Description)
rownames(Slick::Design(OMs_obj)) <- OM_Design$OM_Name
Slick::OMs(slick) <- OMs_obj

# Dimension names of the per-simulation values
Time <- Slick::Time(Slick::Timeseries(slick))
PM_Value <- Slick::Value(Slick::Boxplot(slick))
dimnames(PM_Value) <- list(Sim = seq_len(dim(PM_Value)[1]), OM = OM_Design$OM_Name, MP = MP_all, PI = PM_Codes)
TS_Value <- Slick::Value(Slick::Timeseries(slick))
dimnames(TS_Value) <- list(Sim = seq_len(dim(TS_Value)[1]), OM = OM_Design$OM_Name, MP = MP_all, PI = TS_Codes,
                           Time = Time)
Boxplot_obj <- Slick::Boxplot(slick);       Slick::Value(Boxplot_obj) <- PM_Value;   Slick::Boxplot(slick) <- Boxplot_obj
Timeseries_obj <- Slick::Timeseries(slick); Slick::Value(Timeseries_obj) <- TS_Value; Slick::Timeseries(slick) <- Timeseries_obj

print(Slick::Check(slick))

MSEtool::Save(slick, 'objects/slick.rds', overwrite = TRUE)
# MSEtool::Save(slick, '../SlickLibrary/Slick_Objects/Indian_Ocean_Albacore.slick', overwrite = TRUE)

# ---- Run App ----
if (interactive())
  Slick::App(slick = slick)

# ---- Make manual plots ----
# The Base Case OM, the robustness OMs, and all OMs
OM_Sets <- c(if ('OM5b' %in% OM_Design$OM) list(BaseCase = match('OM5b', OM_Design$OM)),
             if (any(OM_Design$Type == 'Robustness')) list(Robustness = which(OM_Design$Type == 'Robustness')),
             list(All = seq_len(nOM)))

## ---- Timeseries ----
for (om_set in names(OM_Sets)) {
  om_idx <- OM_Sets[[om_set]]
  for (code in TS_Codes) {
    p <- Slick::plotTimeseries(slick, PI = match(code, TS_Codes), OMs = om_idx,
                               byMP = TRUE, includeHist = FALSE)
    ggplot2::ggsave(file.path(FigDir, 'Timeseries', om_set, paste0(TS_File[[code]], '.png')),
                    p, width = 10, height = 6, dpi = 300, create.dir = TRUE)
  }
}

## ---- Kobe and Kobe-time ----
for (om_set in names(OM_Sets)) {
  om_idx <- OM_Sets[[om_set]]

  p_kobe <- Slick::plotKobe(slick, xPI = 1, yPI = 2, Time = FALSE, OMs = om_idx, percentile = NULL)
  ggplot2::ggsave(file.path(FigDir, 'Kobe', paste0(om_set, '.png')),
                  p_kobe, width = 8, height = 7, dpi = 300, create.dir = TRUE)

  p_kobe_time <- Slick::plotKobe(slick, xPI = 1, yPI = 2, Time = TRUE, OMs = om_idx)
  ggplot2::ggsave(file.path(FigDir, 'Kobe', paste0(om_set, '_Time.png')),
                  p_kobe_time, width = 10, height = 5, dpi = 300)
}

## ---- Quilt tables ----
Quilt_Primary_Codes  <- intersect(Headline_Codes, Quilt_Codes)
Quilt_Extended_Codes <- c(Quilt_Primary_Codes,
                          'PGK_Short', 'PGK_Medium', 'PGK_Long',
                          'Yield_Medium', 'Yield_Long',
                          'SB_SBMSY_All', 'F_FMSY_All')

Save_Quilt_Table <- function(pi_codes, om_idx, file_stub) {
  slick_sub <- Slick::FilterSlick(slick, OMs = om_idx, PIs = match(pi_codes, Quilt_Codes), plot = 'Quilt')
  tbl       <- Slick::plotQuilt(slick_sub, kable = TRUE)

  prob_cols <- intersect(Prob01_Codes, pi_codes)
  if (length(prob_cols))
    tbl <- flextable::colformat_double(tbl, j = prob_cols, digits = 2)

  dir <- file.path(FigDir, 'Quilt')
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  htmltools::save_html(flextable::htmltools_value(tbl), file = file.path(dir, paste0(file_stub, '.html')))
  flextable::save_as_docx(tbl, path = file.path(dir, paste0(file_stub, '.docx')))
}

for (om_set in names(OM_Sets)) {
  Save_Quilt_Table(Quilt_Primary_Codes,  OM_Sets[[om_set]], paste0(om_set, '_Primary'))
  Save_Quilt_Table(Quilt_Extended_Codes, OM_Sets[[om_set]], paste0(om_set, '_Extended'))
}

## ---- Comparison across OMs ----
OM_Colours <- scales::hue_pal()(nOM)
OM_Colours[OM_Design$OM == 'OM5b'] <- 'black'
names(OM_Colours) <- OM_Design$OM_Name

# Median time series in each OM, one panel per MP
TS_Median <- apply(TS_Value, c(2, 3, 4, 5), stats::median, na.rm = TRUE) |>
  as.data.frame.table(responseName = 'Value', stringsAsFactors = FALSE)
TS_Median$Time <- as.numeric(TS_Median$Time)
TS_Median$MP   <- factor(TS_Median$MP, levels = MP_all)
TS_Median$OM   <- factor(TS_Median$OM, levels = OM_Design$OM_Name)

for (code in TS_Codes) {
  p <- ggplot2::ggplot(TS_Median[TS_Median$PI == code, ], ggplot2::aes(Time, Value, colour = OM)) +
    ggplot2::geom_vline(xintercept = MPStartYear, linetype = 3, colour = 'grey40') +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::facet_wrap(~MP, ncol = 3, scales = 'free_y') +
    ggplot2::scale_colour_manual(values = OM_Colours) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = NULL, y = TS_Label[[code]], colour = NULL,
                  caption = 'Median over simulations. Dotted line: first year of MP advice.') +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = 'bottom')
  if (!is.na(TS_Target[[code]]))
    p <- p + ggplot2::geom_hline(yintercept = TS_Target[[code]], linetype = 2, colour = 'darkgreen')
  if (!is.na(TS_Limit[[code]]))
    p <- p + ggplot2::geom_hline(yintercept = TS_Limit[[code]], linetype = 2, colour = 'red')
  ggplot2::ggsave(file.path(FigDir, 'Timeseries', 'ByOM', paste0(TS_File[[code]], '.png')),
                  p, width = 12, height = 8, dpi = 300, create.dir = TRUE)
}

# Headline performance indicators in each OM
PI_OM <- apply(PM_Value, c(2, 3, 4), mean, na.rm = TRUE) |>
  as.data.frame.table(responseName = 'Value', stringsAsFactors = FALSE)

d <- PI_OM[PI_OM$PI %in% Headline_Codes, ]
d$PI <- factor(PM_Label[d$PI], levels = PM_Label[Headline_Codes])
d$MP <- factor(d$MP, levels = MP_all)
d$OM <- factor(d$OM, levels = OM_Design$OM_Name)
Min  <- data.frame(PI = factor(PM_Label[c('PGK_Tune', 'Safety_SB')], levels = levels(d$PI)),
                   Value = c(PGKMin, PSafety))

p <- ggplot2::ggplot(d, ggplot2::aes(MP, Value, colour = OM)) +
  ggplot2::geom_hline(data = Min, ggplot2::aes(yintercept = Value), linetype = 2, colour = 'red') +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.7), size = 2) +
  ggplot2::facet_wrap(~PI, ncol = 3, scales = 'free_y', labeller = ggplot2::label_wrap_gen(35)) +
  ggplot2::scale_colour_manual(values = OM_Colours) +
  ggplot2::labs(x = NULL, y = NULL, colour = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = 'bottom', axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(file.path(FigDir, 'Summary', 'PI_by_OM.png'), p, width = 11, height = 10, dpi = 300,
                create.dir = TRUE)

# Table: one row per PI and CMP, one column per OM, and the mean and minimum over OMs
Tab <- stats::reshape(PI_OM, idvar = c('PI', 'MP'), timevar = 'OM', direction = 'wide')
names(Tab) <- sub('^Value[.]', '', names(Tab))
Tab <- Tab[order(match(Tab$PI, PM_Codes), match(Tab$MP, MP_all)), ]
Tab$Mean <- rowMeans(Tab[OM_Design$OM_Name])
Tab$Min  <- apply(Tab[OM_Design$OM_Name], 1, min)
utils::write.csv(Tab, file.path(FigDir, 'Summary', 'PI_by_OM.csv'), row.names = FALSE)
