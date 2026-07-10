#' Import an MSEtool OM from IOTC Albacore conditioning output
#'
#' Constructs a multi-stock, multi-fleet MSEtool [MSEtool::OM()] from the MCMC
#' posterior samples produced by the IOTC Albacore conditioning model. The OM
#' is quarterly (4 seasons) with two sexes (Female, Male; from
#' [StockObjectList]) and six fleets.
#'
#' Natural mortality and stock-recruit parameters are drawn directly from
#' `CondData_<OMName>` (e.g. [CondData_OM5b]). `object$R0` is total
#' (both-sex) unfished recruitment and is split evenly between the two
#' stocks; recruitment is restricted to the recruitment season (season 4)
#' via a mask on `R0`, and `RecDevInit`/`RecDevHist` are derived from the
#' conditioning model's own numbers-at-age and recruitment deviations.
#' `R0` is passed to [MSEtool::SRR()] in thousands (`Units = 1000`), so
#' [MSEtool::Simulate()] output (e.g. `Hist@Number`) is in thousands of fish,
#' not the raw individual counts used by `CondData_<OMName>`.
#'
#' Historical fishing mortality is recovered from the MCMC harvest rates
#' (`object$H`), and fleet selectivity from `object$sela`.
#'
#' @param OMName Character scalar. Suffix identifying which conditioning
#'   dataset to load, via `get(paste0('CondData_', OMName))` (e.g. `'OM5b'`
#'   loads [CondData_OM5b]). Also used as the `Name` of the returned OM.
#'
#' @return An [MSEtool::OM()] object with stocks and fleets fully populated
#'   from the MCMC posterior, ready to pass to [MSEtool::Simulate()].
#'
#' @export
ImportOM <- function(OMName = 'OM5b') {

  # ---- Load MCMC output ----
  object <- get(paste0('CondData_', OMName))

  # ---- Initialize OM ----
  Seasons     <- 4
  nSim        <- length(object)
  CurrentYear <- 2020
  nYear       <- length(2000:CurrentYear)

  nFleet     <- 6
  FleetNames <- c(paste0('LL', 1:4), 'PS', 'Other')

  OM <- MSEtool::OM(Name        = OMName,
                    Agency      = 'ITOC',
                    nSim        = nSim,
                    nYear       = nYear,
                    CurrentYear = CurrentYear,
                    Seasons     = Seasons)

  HistYears <- MSEtool::Years(OM, 'H')

  StockList  <- StockObjectList
  StockNames <- names(StockList)

  # ---- Update Stocks with Natural Mortality & SRR ----
  for (st in seq_along(StockNames)) {

    stock <- StockList[[st]]

    AgeClasses <- MSEtool::Classes(stock)
    nAge       <- length(AgeClasses)

    AnnualAges <- unique(as.integer(AgeClasses))
    nAnnualAge  <- length(AnnualAges)
    nTS        <- nYear * Seasons

    flat_age_idx  <- c(1L, rep(seq(2L, nAnnualAge), each = 4L), rep(nAnnualAge, 3L))
    flat_seas_idx <- c(4L, rep(1:4, nAnnualAge - 1L), 1:3)


    ## ---- Natural Mortality ----
    natural_mortality <- purrr::map_dbl(object, 'M')

    M_array <- array(natural_mortality, c(nSim, nAge),
                     dimnames = list(
                       Sim = seq_len(nSim),
                       Age = AgeClasses)
                     )

    MSEtool::NaturalMortality(stock) <- MSEtool::NaturalMortality(
      MeanAtAge = MSEtool::AddDimension(M_array, 'Year', HistYears[1])
      )


    ## ---- Stock-Recruit ----
    psi <- 0.5
    R0s <- purrr::map_dbl(object, 'R0') * psi
    R0_array <- array(R0s, dim = c(nSim, nYear * Seasons),
                      dimnames = list(
                        Sim = seq_len(nSim),
                        Year = HistYears)
    )

    steepness <- purrr::map_dbl(object, 'h')
    sigmaR    <- purrr::map_dbl(object, 'sigmar')
    ac        <- purrr::map_dbl(object, 'rho')

    spawn_season <- 3L
    rec_season   <- 4L
    rec_devs <- purrr::map(object, 'epsrx')

    # R0 is only nonzero in the recruitment season
    rec_mask <- matrix(0L, nSim, nYear * Seasons)
    rec_mask[, seq(rec_season, nYear * Seasons, by = Seasons)] <- 1L
    R0_thousands <- R0_array * rec_mask / 1000

    RecDevHist <- array(0,
                        dim = c(nSim, nTS),
                        dimnames = list(
                          Sim = seq_len(nSim),
                          Year = HistYears)
                        )

    rec_ts <- seq(from = rec_season, by = Seasons, length.out = nYear )
    RecDevHist[,rec_ts] <- exp(cbind(0,do.call('rbind', rec_devs)))

    n_primary <- 1L + (nAnnualAge - 1L) * 4L
    init_pos  <- 2:nAge
    N_obs_init <- do.call(rbind, purrr::map(object, \(it) {
      vapply(init_pos, \(pos) {
        if (pos > n_primary) return(0)
        aidx <- flat_age_idx[pos]
        seas <- flat_seas_idx[pos]
        if (seas == 1L && aidx <= nAnnualAge) it$N[1L, aidx, 1L, st] else 0
      }, numeric(1L))
    }))

    N0_eq <- R0s *  t(exp(-outer(seq_along(AgeClasses) - 1L, natural_mortality)))
    RecDevInit <- N_obs_init / N0_eq[,-1]

    MSEtool::SRR(stock) <- MSEtool::SRR(Pars       = list(h = steepness),
                                        R0         = R0_thousands,
                                        SD         = sigmaR,
                                        AC         = ac,
                                        SPFrom     = StockNames[1],
                                        RecDevInit = RecDevInit,
                                        RecDevHist = RecDevHist,
                                        Units      = 1000,
                                        SpawnLag   = rec_season - spawn_season)


    StockList[[st]] <- stock
  }
  MSEtool::Stock(OM) <- StockList

  # ---- Fleets ----
  FleetList <- MSEtool::MakeNamedList(StockNames, MSEtool::MakeNamedList(FleetNames))


  F_ref_sims <- purrr::imap(object, \(it, sim) {
    F_ref_ts <- matrix(0, nTS, nFleet)
    for (ts in seq_len(nTS)) {
      yr <- ((ts - 1L) %/% Seasons) + 1L
      s  <- ((ts - 1L) %% Seasons) + 1L
      F_ref_ts[ts, ] <- -log(1 - pmax(it$H[yr, s, ], 0))
    }
    F_ref_ts
  })

  for (st in seq_along(StockNames)) {
    for (fl in seq_along(FleetNames)) {

      fleet <- MSEtool::Fleet(Name = FleetNames[fl])

      ## ---- Selectivity ----
      sel_age <- purrr::map(object, \(it) {
        vapply(seq_len(nAge), \(a) {
          s     <- flat_seas_idx[a]
          ai    <- flat_age_idx[a]
          denom <- max(it$sela[, s, st, fl])
          if (denom > 0) it$sela[ai, s, st, fl] / denom else 0
        }, numeric(1L))
      })

      ## ---- Effort ----
      Effort_mat <- do.call(rbind, purrr::imap(object, \(it, sim) {
        vapply(seq_len(nTS), \(ts) {
          s        <- ((ts - 1L) %% Seasons) + 1L
          max_sela <- max(it$sela[, s, st, fl])
          F_ref_sims[[sim]][ts, fl] * max_sela
        }, numeric(1L))
      }))

      dimnames(Effort_mat)  <- list(Sim = seq_len(nSim), Year = HistYears)
      sel_age_mat <- do.call(rbind, sel_age)
      dimnames(sel_age_mat) <- list(Sim = seq_len(nSim), Age = AgeClasses)

      MSEtool::Effort(fleet)       <- MSEtool::Effort(Effort = Effort_mat)
      MSEtool::Catchability(fleet) <- MSEtool::Catchability(Efficiency = 1)
      MSEtool::Selectivity(fleet)  <- MSEtool::Selectivity(
        MeanAtAge = MSEtool::AddDimension(sel_age_mat, 'Year', HistYears[1]))

      FleetList[[st]][[fl]] <- fleet
    }
  }

  MSEtool::Fleet(OM) <- FleetList

  # ---- Data ----
  OM@Data <- list(Combined = iotcALB::albMSE_Data)

  Complexes(OM) <- list('Combined' = 1:2)

  OM
}



# This is code to import the data shared by Rich, saved in object `Cond_BaseCase`
# It has now been updated to load the OM conditioning data saved in
# `data/CondData_OMxx.rda` and built in `data-raw/albMSE_CondData.R`

# #' Import an MSEtool OM from IOTC Albacore conditioning output
# #'
# #' Constructs a multi-stock, multi-fleet MSEtool [MSEtool::OM()] from the MCMC
# #' posterior samples produced by the IOTC Albacore conditioning model. The OM
# #' is quarterly (4 seasons) with two sexes (Female, Male) and six fleets.
# #'
# #' Biological schedules (weight, maturity, natural mortality) and stock-recruit
# #' parameters are drawn directly from the conditioning object. Historical
# #' fishing mortality is recovered from the MCMC numbers-at-age via N-at-age
# #' transitions.
# #'
# #' @param object A named list of conditioning output, as produced by the IOTC
# #'   Albacore conditioning model. Defaults to [Cond_BaseCase]. Must contain:
# #'   * `mcmcvars`: list of MCMC iterations, each with arrays `N`, `H`, `sela`
# #'   * `steepness`, `natural_mortality`, `R0`, `sigmaR`: numeric vectors of
# #'     length `nSim`
# #'   * `maturity`, `weight`: arrays `[annual_age, season, sex]`
# #' @param Name Character scalar. Name assigned to the returned OM object.
# #'
# #' @return An [MSEtool::OM()] object with stocks and fleets fully populated
# #'   from the MCMC posterior, ready to pass to [MSEtool::Simulate()].
# #'
# #' @export
# ImportOM <- function(object = Cond_BaseCase, Name = 'Base Case') {
#
#   # ---- Initialize OM ----
#   Seasons     <- 4
#   nSim        <- length(object$mcmcvars)
#   CurrentYear <- 2020
#   nYear       <- length(2000:CurrentYear)
#
#   StockNames <- c('Female', 'Male')
#   nFleet     <- 6
#   FleetNames <- c(paste0('LL', 1:4), 'PS', 'Other')
#
#   OM <- MSEtool::OM(Name        = Name,
#                     Agency      = 'ITOC',
#                     nSim        = nSim,
#                     nYear       = nYear,
#                     CurrentYear = CurrentYear,
#                     Seasons     = Seasons)
#
#   HistYears <- MSEtool::Years(OM, 'H')
#
#   spawning_season <- 3L
#
#   # flat_age_idx[p]  : annual age (1-based) for OM position p
#   # flat_seas_idx[p] : reference season for OM position p
#   # Sequence starts at annual age 1 / S4 (recruit entry) and advances at each S4→S1 boundary.
#   nAnnAge_r <- dim(object$mcmcvars[[1]]$N)[2]
#   nTS       <- nYear * Seasons
#
#   flat_age_idx  <- c(1L, rep(seq(2L, nAnnAge_r), each = 4L), rep(nAnnAge_r, 3L))
#   flat_seas_idx <- c(4L, rep(1:4, nAnnAge_r - 1L), 1:3)
#
#   # ---- Import Stocks ----
#   StockList <- MSEtool::MakeNamedList(StockNames)
#
#   for (st in seq_along(StockNames)) {
#
#     stock <- MSEtool::Stock(Name = StockNames[st])
#
#     ## ---- Age structure ----
#     MSEtool::Ages(stock) <- MSEtool::Ages(MinAge    = 0,
#                                           MaxAge    = 59,
#                                           Units     = 'quarter',
#                                           PlusGroup = FALSE)
#
#     AgeClasses <- MSEtool::Classes(stock)
#     nAge       <- MSEtool::nAge(stock)
#     nAnnualAge <- dim(object$mcmcvars[[1]]$N)[2]
#
#     ## ---- Weight ----
#     wght <- do.call(cbind, lapply(seq_len(Seasons), \(s)
#       object$weight[flat_age_idx, s, st]
#     ))
#     MSEtool::Weight(stock) <- MSEtool::Weight(MeanAtAge = array(wght,
#                                                                 dim = c(1, dim(wght)),
#                                                                 dimnames = list(
#                                                                   Sim  = 1,
#                                                                   Age  = AgeClasses,
#                                                                   Year = HistYears[seq_len(ncol(wght))])))
#
#     ## ---- Natural Mortality ----
#     M_array <- matrix(object$natural_mortality, nSim, nAge)
#     dimnames(M_array) <- list(Sim = seq_len(nSim), Age = AgeClasses)
#     MSEtool::NaturalMortality(stock) <- MSEtool::NaturalMortality(MeanAtAge =
#                                           MSEtool::AddDimension(M_array, 'Year', HistYears[1]))
#
#     ## ---- Maturity ----
#     maturity <- do.call(cbind, lapply(seq_len(Seasons), \(s)
#       object$maturity[flat_age_idx, s, st]
#     ))
#     MSEtool::Maturity(stock) <- MSEtool::Maturity(MeanAtAge = array(maturity,
#                                                                     dim = c(1, dim(maturity)),
#                                                                     dimnames = list(
#                                                                       Sim  = 1,
#                                                                       Age  = AgeClasses,
#                                                                       Year = HistYears[seq_len(ncol(maturity))])))
#
#     ## ---- Fecundity ----
#     fecundity <- matrix(0, nrow = nAge, ncol = Seasons)
#     fecundity[, spawning_season] <- object$maturity[flat_age_idx, spawning_season, st] *
#                                     object$weight[flat_age_idx, spawning_season, st]
#     MSEtool::Fecundity(stock) <- MSEtool::Fecundity(MeanAtAge = array(fecundity,
#                                                                       dim = c(1, dim(fecundity)),
#                                                                       dimnames = list(
#                                                                         Sim  = 1,
#                                                                         Age  = AgeClasses,
#                                                                         Year = HistYears[seq_len(ncol(fecundity))])))
#
#     ## ---- Stock-Recruit ----
#     R0_array <- array(object$R0, dim = c(nSim, nYear * Seasons),
#                       dimnames = list(Sim = seq_len(nSim), Year = HistYears))
#
#     rec_mask <- matrix(0L, nSim, nYear * Seasons)
#     rec_mask[, seq(4L, nYear * Seasons, by = Seasons)] <- 1L
#     R0_array <- R0_array * rec_mask
#
#     AnnualYears <- unique(floor(as.numeric(HistYears)))
#     year_idx    <- rep(seq_len(nYear), each = Seasons)
#
#     annual_recruits <- do.call(rbind, purrr::map(object$mcmcvars, \(it)
#       it$N[, 1L, 4L, st]
#     ))
#     dimnames(annual_recruits) <- list(Sim = seq_len(nSim), Year = AnnualYears)
#
#     SSB_mat <- do.call(rbind, purrr::map(object$mcmcvars, \(it) it$SSB))
#     dimnames(SSB_mat) <- list(Sim = seq_len(nSim), Year = AnnualYears)
#
#     pos_s3 <- which(flat_seas_idx == spawning_season)
#     q_s3   <- pos_s3 - 1L
#     sw_s3  <- object$maturity[flat_age_idx[pos_s3], spawning_season, 1L] *
#                object$weight[flat_age_idx[pos_s3], spawning_season, 1L]
#     SSB0 <- vapply(seq_len(nSim), \(sim) {
#       sum(object$R0[sim] * exp(-object$natural_mortality[sim] * q_s3) * sw_s3)
#     }, numeric(1))
#
#     R_pred   <- (4 * object$steepness * object$R0 * SSB_mat) /
#                 (SSB0 * (1 - object$steepness) + SSB_mat * (5 * object$steepness - 1))
#     rec_devs <- annual_recruits / R_pred
#     dimnames(rec_devs) <- list(Sim = seq_len(nSim), Year = AnnualYears)
#
#     RecDevHist <- rec_devs[, year_idx]
#     dimnames(RecDevHist) <- list(Sim = seq_len(nSim), Year = HistYears)
#
#     n_primary <- 1L + (nAnnAge_r - 1L) * 4L
#     init_pos  <- 2:nAge
#     N_obs_init <- do.call(rbind, purrr::map(object$mcmcvars, \(it) {
#       vapply(init_pos, \(pos) {
#         if (pos > n_primary) return(0)
#         aidx <- flat_age_idx[pos]; seas <- flat_seas_idx[pos]
#         if (seas == 1L && aidx <= nAnnualAge) it$N[1L, aidx, 1L, st] else 0
#       }, numeric(1L))
#     }))
#     dimnames(N_obs_init) <- list(Sim = seq_len(nSim), Age = AgeClasses[init_pos])
#
#     N0_eq <- sweep(exp(-outer(object$natural_mortality, as.numeric(init_pos) - 1L)), 1, object$R0, "*")
#     dimnames(N0_eq) <- list(Sim = seq_len(nSim), Age = AgeClasses[init_pos])
#
#     RecDevInit <- N_obs_init / N0_eq
#
#     MSEtool::SRR(stock) <- MSEtool::SRR(Pars       = list(h = object$steepness),
#                                          R0         = R0_array,
#                                          SD         = object$sigmaR,
#                                          SPFrom     = StockNames[1],
#                                          RecDevInit = RecDevInit,
#                                          RecDevHist = RecDevHist,
#                                          SpawnLag   = 1)
#
#     StockList[[st]] <- stock
#   }
#
#   MSEtool::Stock(OM) <- StockList
#
#   # ---- Import Fleets ----
#   FleetList <- MSEtool::MakeNamedList(StockNames, MSEtool::MakeNamedList(FleetNames))
#
#   F_ref_sims <- purrr::imap(object$mcmcvars, \(it, sim) {
#     F_ref_ts <- matrix(0, nTS, nFleet)
#     for (ts in seq_len(nTS)) {
#       yr <- ((ts - 1L) %/% Seasons) + 1L
#       s  <- ((ts - 1L) %% Seasons) + 1L
#       F_ref_ts[ts, ] <- -log(1 - pmax(it$H[yr, s, ], 0))
#     }
#     F_ref_ts
#   })
#
#   for (st in seq_along(StockNames)) {
#     for (fl in seq_along(FleetNames)) {
#
#       fleet <- MSEtool::Fleet(Name = FleetNames[fl])
#
#       ## ---- Selectivity ----
#       sel_age <- purrr::map(object$mcmcvars, \(it) {
#         vapply(seq_len(nAge), \(a) {
#           s     <- flat_seas_idx[a]
#           ai    <- flat_age_idx[a]
#           denom <- max(it$sela[, s, st, fl])
#           if (denom > 0) it$sela[ai, s, st, fl] / denom else 0
#         }, numeric(1L))
#       })
#
#       ## ---- Effort ----
#       Effort_mat <- do.call(rbind, purrr::imap(object$mcmcvars, \(it, sim) {
#         vapply(seq_len(nTS), \(ts) {
#           s        <- ((ts - 1L) %% Seasons) + 1L
#           max_sela <- max(it$sela[, s, st, fl])
#           F_ref_sims[[sim]][ts, fl] * max_sela
#         }, numeric(1L))
#       }))
#
#       dimnames(Effort_mat)  <- list(Sim = seq_len(nSim), Year = HistYears)
#       sel_age_mat <- do.call(rbind, sel_age)
#       dimnames(sel_age_mat) <- list(Sim = seq_len(nSim), Age = AgeClasses)
#
#       MSEtool::Effort(fleet)       <- MSEtool::Effort(Effort = Effort_mat)
#       MSEtool::Catchability(fleet) <- MSEtool::Catchability(Efficiency = 1)
#       MSEtool::Selectivity(fleet)  <- MSEtool::Selectivity(
#         MeanAtAge = MSEtool::AddDimension(sel_age_mat, 'Year', HistYears[1]))
#
#       FleetList[[st]][[fl]] <- fleet
#     }
#   }
#
#   MSEtool::Fleet(OM) <- FleetList
#
#   Complexes(OM) <- list('Combined' = 1:2)
#
#   # ---- Data ----
#   data <- Data()
#   data@YearLH <- OM@CurrentYear
#   data@Years  <- HistYears
#
#   data@CPUE@Name <- FleetNames[1:4]
#   data@CPUE@Value <- array(NA, dim=c(length(HistYears), 4),
#                            dimnames = list(Year = HistYears,
#                                            Fleet =   data@CPUE@Name))
#   for (fl in seq_along(data@CPUE@Name)) {
#     index <- as.vector(t(object$Index[,,fl]))
#     data@CPUE@Value[,fl] <- index/mean(index, na.rm = TRUE)
#   }
#
#   OM@Data <- list(Combined = data)
#
#   # ---- Observation ----
#   # Assume no observation error on catches
#   # simulated value = data values
#   # can replace with real observed catches when available
#   # but not the discrepancy due to HR vs F approaches
#   OM@Obs <- MakeNamedList('Combined',
#                           MakeNamedList(FleetNames,
#                                         Obs(Landings = CatchObs(CV = 0),
#                                             Discards = CatchObs(CV = 0))
#                                         )
#   )
#
#
#   OM
# }
#
