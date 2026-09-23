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

