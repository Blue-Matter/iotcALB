#' Index Level Corresponding to SBMSY
#'
#' Calibrates relative abundance indices to stock status, so that the harvest
#' control rules of the index-based MPs ([MSEtool::IndexRate()],
#' [MSEtool::IndexTarget()], [BufferMP()]) can be expressed in units of
#' `SB/SBMSY`. The index level at `SBMSY` is the mean annual index over `Years`
#' divided by the median (over simulations) of the annual `SB/SBMSY`
#' averaged over the same years (see [AnnualStatus()]).
#'
#' @param Hist A `hist` object (the base case OM).
#' @param Index Character vector. Names of the indices (in the `Survey` slot
#'   after [MSEtool::CombineFleets()], else `CPUE`). Default `'LL1'`.
#' @param Years Numeric. Calendar years used for the calibration. Must be
#'   historical years of `Hist`. Default `2016:2020`.
#'
#' @return Named numeric vector (one value per index): the index level at
#'   `SBMSY`, with attributes `Index` (the mean index over `Years`) and
#'   `SB_SBMSY` (the median status over `Years`).
#' @export
IndexAtSBMSY <- function(Hist, Index = 'LL1', Years = 2016:2020) {
  Data <- MSEtool::AnnualData(Hist@Data[[1]][[1]])
  Slot <- .IndexSlot(Data, Index)
  Val  <- methods::slot(Data, Slot)@Value
  MeanIndex <- colMeans(Val[as.character(Years), Index, drop = FALSE], na.rm = TRUE)

  Status <- AnnualStatus(Hist)
  Status <- Status[Status$Year %in% Years, ]
  SimStatus <- tapply(Status$SB_SBMSY, Status$Sim, mean)
  MedStatus <- stats::median(SimStatus)

  structure(MeanIndex / MedStatus, names = Index, Index = MeanIndex, SB_SBMSY = MedStatus)
}

#' Inverse-Variance Weights of Indices from the Conditioned Observation Error
#'
#' Weights for combining several indices in the index-based MPs, proportional
#' to `1 / SD^2`, where `SD` is the median (over simulations) standard deviation
#' of the log-scale observation error of each index conditioned in the OM.
#'
#' @param Hist A `hist` object (the base case OM).
#' @param Index Character vector. Names of the indices.
#'
#' @return Named numeric vector of weights, summing to 1, with attribute `SD`.
#' @export
IndexWeights <- function(Hist, Index) {
  Obs <- Hist@OM@Obs[[1]]
  SD  <- vapply(Index, \(ix) stats::median(Obs[[ix]]@Survey@Stats$SD), numeric(1))
  W   <- 1 / SD^2
  structure(W / sum(W), names = Index, SD = SD)
}

#' Surplus Production Model Settings from the Operating Model
#'
#' Derives the fixed parameters and informative priors of the surplus
#' production model used by [MSEtool::SurplusProduction()] from the
#' operating model, so that the assessment model in the MP matches the
#' productivity and initial depletion of the base case OM.
#'
#' The surplus production model has a single biomass pool, so the OM
#' quantities are calculated for the total biomass of both sexes:
#'
#' - `MSY`: equilibrium removals at MSY (`MSYLandings + MSYDiscards`).
#' - `FMSY`: the harvest rate at MSY, `MSY / BMSY`, where `BMSY` is the total
#'   equilibrium biomass at MSY. This is the `FMSY` of the surplus production
#'   model, not the apical fishing mortality of the OM.
#' - `Shape`: the Pella-Tomlinson shape parameter `n` with
#'   `n^(1/(1-n)) = BMSY/B0`, the median over simulations.
#' - `Depletion`: total biomass relative to unfished at the start of the
#'   first historical year (`B/B0`), the median over simulations.
#'
#' The priors are lognormal, with median equal to the median over simulations
#' and CV equal to the CV of the OM values over simulations or `MinCV`,
#' whichever is larger.
#'
#' @param Hist A `hist` object (the base case OM).
#' @param MinCV Named numeric vector. Minimum CV of the priors on `FMSY` and
#'   `MSY`. Default `c(FMSY = 0.3, MSY = 0.4)`.
#'
#' @return A list with `Priors` (a list of `c(median, CV)` for `FMSY` and
#'   `MSY`, see [MSEtool::FitSP()]), `Shape`, `Depletion`, and `OM` (a
#'   `data.frame` of the per-simulation OM values).
#' @export
SPSettings <- function(Hist, MinCV = c(FMSY = 0.3, MSY = 0.4)) {
  MSY  <- apply(MSEtool::MSYLandings(Hist), 1, sum) +
    rep_len(apply(MSEtool::MSYDiscards(Hist), 1, sum), MSEtool::nSim(Hist))
  BMSY <- apply(MSEtool::BMSY(Hist), 1, sum)

  B   <- MSEtool::Biomass(Hist, df = TRUE, Reduce = FALSE)
  B   <- B[B$Year == min(B$Year), ]
  BB0 <- MSEtool::B_B0(Hist, df = TRUE, Reduce = FALSE)
  BB0 <- BB0[BB0$Year == min(BB0$Year), ]
  B0  <- tapply(B$Value / BB0$Value[match(paste(B$Sim, B$Stock), paste(BB0$Sim, BB0$Stock))],
                B$Sim, sum)
  B1  <- tapply(B$Value, B$Sim, sum)

  OM <- data.frame(Sim = seq_along(MSY), MSY = unname(MSY), BMSY = unname(BMSY),
                   B0 = as.numeric(B0), FMSY = unname(MSY / BMSY),
                   BMSY_B0 = unname(BMSY) / as.numeric(B0),
                   Depletion = as.numeric(B1 / B0))

  LogNormalCV <- function(x) sqrt(exp(stats::sd(log(x))^2) - 1)
  Prior <- function(x, nm) c(stats::median(x), max(LogNormalCV(x), MinCV[[nm]]))

  BMSY_K <- stats::median(OM$BMSY_B0)
  Shape  <- stats::uniroot(\(n) n^(1 / (1 - n)) - BMSY_K, c(1 + 1e-6, 50))$root

  list(Priors    = list(FMSY = Prior(OM$FMSY, 'FMSY'), MSY = Prior(OM$MSY, 'MSY')),
       Shape     = Shape,
       Depletion = stats::median(OM$Depletion),
       OM        = OM)
}

#' Harvest Control Rule Control Points of a Buffer HCR
#'
#' Converts the buffer harvest control rule used in the earlier albacore MSE
#' work (Mosqueira and Hillary 2025, IOTC-2025-WPM16-11) to the control
#' points of [MSEtool::IndexTarget()]. The TAC is the reference catch while
#' the index is between the lower and upper buffers; it decreases below the
#' lower buffer, more rapidly below the limit, and increases above the upper
#' buffer at `SlopeRatio` times the rate of decrease between the limit and
#' the lower buffer, up to `MaxIndex`.
#'
#' @param Buffer Numeric, length 2. Lower and upper buffers, in units of
#'   `status / IndexTarget`.
#' @param Limit Numeric. Limit, in units of `status / IndexTarget`. Default
#'   `0.4`.
#' @param LimitRate Numeric. TAC multiplier at the limit. Default `0.5`.
#' @param Zero Numeric. Index level at which the TAC is zero. Default `0.2`.
#' @param SlopeRatio Numeric. Slope of the increase above the upper buffer
#'   relative to the slope of the decrease below the lower buffer. Default
#'   `0.25`.
#' @param MaxIndex Numeric. Index level above which the TAC no longer
#'   increases. Default `2`.
#'
#' @return A list with `HCRControlPointsIndex` and `HCRControlPointsRate`.
#' @export
BufferHCR <- function(Buffer, Limit = 0.4, LimitRate = 0.5, Zero = 0.2,
                      SlopeRatio = 0.25, MaxIndex = 2) {
  SlopeDown <- (1 - LimitRate) / (Buffer[1] - Limit)
  list(HCRControlPointsIndex = c(Zero, Limit, Buffer[1], Buffer[2], MaxIndex),
       HCRControlPointsRate  = c(0, LimitRate, 1, 1,
                                 1 + SlopeRatio * SlopeDown * (MaxIndex - Buffer[2])))
}

#' Buffer Harvest Control Rule Management Procedure
#'
#' A model-free management procedure that replicates the CPUE + buffer HCR MP
#' of the earlier albacore MSE work (`cpue.ind` and `buffer.hcr` of the FLR
#' `mse`/`msemodules` packages; Mosqueira and Hillary 2025, IOTC-2025-WPM16-11,
#' IOTC-2026-WPM17(MSE)-03). Unlike [MSEtool::IndexTarget()] with a buffer HCR
#' ([BufferHCR()]), the HCR output is a multiplier on the previous TAC, so the
#' TAC keeps changing in each management cycle while the index is outside the
#' buffer.
#'
#' ## Procedure
#'
#' 1. The data are aggregated to calendar years with [MSEtool::AnnualData()]
#'    (each index averaged over the seasons with data).
#' 2. The index metric is the mean of each index over the last `nYears` data
#'    years: a weighted mean (`Metric = 'wmean'`), with weight 0.5 for the last
#'    year and the other 0.5 shared among the earlier years in proportion to
#'    their position (e.g. 1/12, 2/12, 3/12 for `nYears = 4`), or an
#'    unweighted mean (`Metric = 'mean'`).
#' 3. Stock status is the metric relative to `IndexRef` (e.g. the index level
#'    at `SBMSY`, [IndexAtSBMSY()]), combined over indices with `IndexWeight`.
#' 4. The TAC multiplier is, for status `x`, limit `L`, lower buffer `B1` and
#'    upper buffer `B2`:
#'    - `(x / L)^2 / 2` for `x <= L`;
#'    - `0.5 (1 + (x - L) / (B1 - L))` for `L < x <= B1`;
#'    - `1` for `B1 < x < B2`;
#'    - `1 + SlopeRatio (0.5 / (B1 - L)) (x - B2)` for `x >= B2`.
#' 5. The TAC is the previous TAC ([MSEtool::LastTAC()]) times the multiplier,
#'    with the change constrained by `DeltaDown` and `DeltaUp`
#'    ([MSEtool::ConstrainTAC()]).
#'
#' ## Tuning
#'
#' In the earlier work the MP was tuned by the position of the upper buffer.
#' Here `tunepar` sets the upper buffer to
#' `LowerBuffer + (UpperBuffer - LowerBuffer) / tunepar`, so that larger values
#' of `tunepar` move the upper buffer towards the lower buffer and increase the
#' catch. `tunepar = 1` uses `UpperBuffer`.
#'
#' @param Data A `data` object.
#' @param Indices Character vector. Names of the indices in `IndexSource`.
#'   Default `'LL1'`.
#' @param IndexSource Character. `'Survey'` (default) or `'CPUE'`.
#' @param IndexRef Numeric vector, one per index. The index level that
#'   corresponds to a status of 1 (e.g. the index at `SBMSY`).
#' @param IndexWeight `NULL` (default, equal weights) or a numeric vector of
#'   weights, one per index.
#' @param nYears Positive integer. Number of data years in the metric. Default
#'   `4`.
#' @param Metric Character. `'wmean'` (default) or `'mean'`.
#' @param Limit,LowerBuffer,UpperBuffer Numeric. The limit and the lower and
#'   upper buffers, in units of status (index / `IndexRef`). Defaults `0.4`,
#'   `0.56`, and `0.96`: the values of the earlier work (0.25, 0.35, and 0.60
#'   in units of the index) scaled so that the limit is at 0.4 `SBMSY`.
#' @param SlopeRatio Numeric. Slope of the increase above the upper buffer
#'   relative to the slope of the decrease between the limit and the lower
#'   buffer. Default `0.10`.
#' @param DeltaDown,DeltaUp Numeric vectors, length 2. Minimum and maximum
#'   proportional decrease and increase in the TAC (see
#'   [MSEtool::ConstrainTAC()]). Default `c(0, 0.15)`.
#' @param TACType Character. `'Removals'` (default) or `'Landings'`.
#' @param tunepar Positive number. Tuning parameter; see Tuning.
#'
#' @return An `advice` object.
#' @export
BufferMP <- function(Data,
                     Indices     = 'LL1',
                     IndexSource = 'Survey',
                     IndexRef,
                     IndexWeight = NULL,
                     nYears      = 4,
                     Metric      = c('wmean', 'mean'),
                     Limit       = 0.4,
                     LowerBuffer = 0.56,
                     UpperBuffer = 0.96,
                     SlopeRatio  = 0.10,
                     DeltaDown   = c(0, 0.15),
                     DeltaUp     = c(0, 0.15),
                     TACType     = c('Removals', 'Landings'),
                     tunepar     = 1) {

  Metric  <- match.arg(Metric)
  TACType <- match.arg(TACType)
  if (!is.numeric(tunepar) || length(tunepar) != 1 || !is.finite(tunepar) || tunepar <= 0)
    cli::cli_abort('{.arg tunepar} must be a single positive number.')
  if (length(IndexRef) != length(Indices))
    cli::cli_abort('{.arg IndexRef} must have one value for each of {.arg Indices}.')
  if (!(Limit < LowerBuffer && LowerBuffer < UpperBuffer))
    cli::cli_abort('{.arg Limit}, {.arg LowerBuffer}, and {.arg UpperBuffer} must be increasing.')
  if (is.null(IndexWeight)) IndexWeight <- rep(1, length(Indices))

  AD  <- MSEtool::AnnualData(Data)
  Val <- methods::slot(AD, IndexSource)@Value[, Indices, drop = FALSE]
  Val <- utils::tail(Val[as.numeric(rownames(Val)) <= max(AD@Years), , drop = FALSE], nYears)
  Val[!is.finite(Val) | Val <= 0] <- NA

  n   <- nrow(Val)
  Wts <- if (Metric == 'wmean' && n > 1) c(0.5 * seq_len(n - 1) / sum(seq_len(n - 1)), 0.5) else rep(1 / n, n)
  Met <- apply(Val, 2, \(v) if (all(is.na(v))) NA_real_ else sum(Wts * v, na.rm = TRUE) / sum(Wts[!is.na(v)]))

  Status <- stats::weighted.mean(Met / IndexRef, IndexWeight, na.rm = TRUE)

  B2   <- LowerBuffer + (UpperBuffer - LowerBuffer) / tunepar
  Mult <- if (!is.finite(Status)) {
    1
  } else if (Status <= Limit) {
    (Status / Limit)^2 / 2
  } else if (Status <= LowerBuffer) {
    0.5 * (1 + (Status - Limit) / (LowerBuffer - Limit))
  } else if (Status < B2) {
    1
  } else {
    1 + SlopeRatio * (0.5 / (LowerBuffer - Limit)) * (Status - B2)
  }

  PrevTAC <- MSEtool::LastTAC(Data, TACType)
  TAC <- MSEtool::ConstrainTAC(PrevTAC, Mult, DeltaDown, DeltaUp, c(0, Inf))

  Adv <- MSEtool::Advice(TAC = TAC, TACType = TACType, TACUnit = Data@Landings@Units)
  Adv@Misc$BufferMP <- list(Status = Status, Multiplier = Mult, UpperBuffer = B2)
  Adv
}
class(BufferMP) <- 'mp'
