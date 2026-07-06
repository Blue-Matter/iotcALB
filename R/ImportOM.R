
#' Import an MSEtool OM from IOTC Albacore conditioning output
#'
#' Constructs a multi-stock, multi-fleet MSEtool [MSEtool::OM()] from the MCMC
#' posterior samples produced by the IOTC Albacore conditioning model. The OM
#' is quarterly (4 seasons) with two sexes (Female, Male) and six fleets.
#'
#' Biological schedules (weight, maturity, natural mortality) and stock-recruit
#' parameters are drawn directly from the conditioning object. Historical
#' fishing mortality is recovered from the MCMC numbers-at-age via N-at-age
#' transitions, which exactly reproduces the reference model's survival without
#' requiring knowledge of the internal H parameterisation.
#'
#' @param object A named list of conditioning output, as produced by the IOTC
#'   Albacore conditioning model. Defaults to [Cond_BaseCase]. Must contain:
#'   * `mcmcvars`: list of MCMC iterations, each with arrays `N`, `H`, `sela`
#'   * `steepness`, `natural_mortality`, `R0`, `sigmaR`: numeric vectors of
#'     length `nSim`
#'   * `maturity`, `weight`: arrays `[annual_age, season, sex]`
#' @param Name Character scalar. Name assigned to the returned OM object.
#'
#' @return An [MSEtool::OM()] object with stocks and fleets fully populated
#'   from the MCMC posterior, ready to pass to [MSEtool::Simulate()].
#'
#' @export
ImportOM <- function(object = Cond_BaseCase, Name = 'Base Case') {

  # ---- Initialize OM ----
  Seasons     <- 4                        # quarterly model
  nSim        <- length(object$mcmcvars)
  CurrentYear <- 2020                     # last historical year included in the model
  nYear       <- length(2000:CurrentYear) # model starts in 2000


  StockNames <- c('Female', 'Male')

  nFleet <- 6
  FleetNames <- paste('Fleet', seq_len(nFleet))

  OM <- MSEtool::OM(Name        = Name,
                    Agency      = 'ITOC',
                    nSim        = nSim,
                    nYear       = nYear,
                    CurrentYear = CurrentYear,
                    Seasons     = Seasons
  )

  HistYears <- MSEtool::Years(OM, 'H')

  # ---- Import Stocks ----

  StockList  <- MSEtool::MakeNamedList(StockNames)

  for (st in seq_along(StockNames)) {

    ## ---- Create `stock` object ----
    stock <- MSEtool::Stock(Name = StockNames[st])

    ## ---- Age structure ----
    MSEtool::Ages(stock) <- MSEtool::Ages(MinAge    = 0,
                                          MaxAge    = 59,
                                          Units     = 'quarter',
                                          PlusGroup = FALSE)

    AgeClasses <- MSEtool::Classes(stock)
    nAge       <- MSEtool::nAge(stock)
    nAnnualAge <- dim(object$mcmcvars[[1]]$N)[2] # nAge/Seasons

    ## ---- Length -----
    # only have ALK, not mean-at-length, classes, or CV
    # not strictly needed

    ## ---- Weight ----
    wght <- ExpandAgeSeason(object$weight[,,st], AgeClasses)
    MSEtool::Weight(stock) <- MSEtool::Weight(MeanAtAge = array(wght,
                                                                dim = c(1, dim(wght)),
                                                                dimnames = list(
                                                                  Sim  = 1,
                                                                  Age  = AgeClasses,
                                                                  Year = HistYears[seq_len(ncol(wght))]))
    )

    ## ---- Natural Mortality ----
    M_array <- matrix(object$natural_mortality, nSim, nAge)
    dimnames(M_array) <- list(Sim = seq_len(nSim),
                              Age = AgeClasses)
    MSEtool::NaturalMortality(stock) <- MSEtool::NaturalMortality(MeanAtAge =
                                                                    MSEtool::AddDimension(M_array, 'Year', HistYears[1]))

    ## ---- Maturity ----
    maturity <- ExpandAgeSeason(object$maturity[,,st], AgeClasses)
    MSEtool::Maturity(stock) <- MSEtool::Maturity(MeanAtAge = array(maturity,
                                                                    dim = c(1, dim(maturity)),
                                                                    dimnames = list(
                                                                      Sim  = 1,
                                                                      Age  = AgeClasses,
                                                                      Year = HistYears[seq_len(ncol(maturity))]))
    )

    ## ---- Stock-Recruit ----
    R0_array <- array(object$R0, dim = c(nSim, nYear * Seasons),
                      dimnames = list(
                        Sim  = seq_len(nSim),
                        Year = HistYears))

    # seasonal recruitment
    recruits <- purrr::map(object$mcmcvars, \(it) {
      as.vector(t(it$N[,1,,st]))
    })

    recruits_mat <- do.call(rbind, recruits)
    dimnames(recruits_mat) <- list(Sim = seq_len(nSim), Year = HistYears)

    rec_season <- recruits_mat
    rec_season[rec_season > 0] <- 1

    R0_array <- R0_array * rec_season

    # back-calculate recruitment deviations from annual SSB
    AnnualYears <- unique(floor(as.numeric(HistYears)))
    year_idx    <- rep(seq_len(nYear), each = Seasons)
    annual_recruits <- t(apply(recruits_mat, 1, \(x) tapply(x, year_idx, sum)))
    dimnames(annual_recruits) <- list(Sim = seq_len(nSim), Year = AnnualYears)

    # SSB at S4 (spawning season) — always use female (st=1) to match SPFrom='Female'
    sw      <- object$maturity[, 4, 1L] * object$weight[, 4, 1L]
    SSB_mat <- do.call(rbind, purrr::map(object$mcmcvars, \(it) drop(it$N[, , 4, 1L] %*% sw)))
    dimnames(SSB_mat) <- list(Sim = seq_len(nSim), Year = AnnualYears)

    # Unfished equilibrium SSB at S4 per sim
    SSB0 <- vapply(seq_len(nSim), \(sim) {
      N0_S4 <- object$R0[sim] * exp(-object$natural_mortality[sim] * 4 * seq(0, nAnnualAge - 1))
      sum(N0_S4 * sw)
    }, numeric(1))

    R_pred <- (4 * object$steepness * object$R0 * SSB_mat) / (SSB0 * (1 - object$steepness) + SSB_mat * (5 * object$steepness - 1))
    rec_devs <- annual_recruits / R_pred
    dimnames(rec_devs) <- list(Sim = seq_len(nSim), Year = AnnualYears)

    RecDevHist <- rec_devs[, year_idx]
    dimnames(RecDevHist) <- list(Sim = seq_len(nSim), Year = HistYears)

    # Initial recruitment deviations: ratio of observed N at Year 1 Season 1 to
    # unfished equilibrium N at each quarterly age class.
    AgeClasses_init <- AgeClasses[-1]
    seas_from_q     <- round((AgeClasses_init %% 1) * Seasons)
    seas_from_q     <- ifelse(seas_from_q == 0L, Seasons, seas_from_q)
    age_arr_idx     <- floor(AgeClasses_init) + 2L

    # observed N at initial year S1
    N_obs_init <- do.call(rbind, purrr::map(object$mcmcvars, \(it) {
      mapply(\(aidx, seas) {
        if (seas == 1L && aidx <= nAnnualAge) it$N[1, aidx, 1L, st] else 0
      }, age_arr_idx, seas_from_q)
    }))
    dimnames(N_obs_init) <- list(Sim = seq_len(nSim), Age = AgeClasses_init)

    # unfished equilibrium N: R0 * exp(-M_quarterly * q_quarters)
    N0_eq <- sweep(exp(-outer(object$natural_mortality, AgeClasses_init * 4)), 1, object$R0, "*")
    dimnames(N0_eq) <- list(Sim = seq_len(nSim), Age = AgeClasses_init)

    RecDevInit <- N_obs_init / N0_eq

    MSEtool::SRR(stock) <- MSEtool::SRR(Pars = list(h = object$steepness),
                                        R0 = R0_array,
                                        SD = object$sigmaR,
                                        SPFrom = StockNames[1],
                                        RecDevInit = RecDevInit,
                                        RecDevHist = RecDevHist
    )


    StockList[[st]] <- stock
  }

  MSEtool::Stock(OM) <- StockList


  # ---- Import Fleets ----
  FleetList  <- MSEtool::MakeNamedList(StockNames,
                                       MSEtool::MakeNamedList(FleetNames))

  nAnnAge_r <- dim(object$mcmcvars[[1]]$N)[2]
  nTS       <- nYear * Seasons

  # Index vectors mapping each OM quarterly age class to (annual_age_idx, season)
  flat_age_idx  <- c(1L,
                     rep(seq(2L, nAnnAge_r), each = 4L),
                     rep(nAnnAge_r, 3L))
  flat_seas_idx <- c(4L,
                     rep(1:4, nAnnAge_r - 1L),
                     1:3)

  # Per-season apical: last (oldest) populated OM age class for each season
  apical_by_seas <- vapply(seq_len(Seasons), \(s) {
    if (s == Seasons) max(seq(1L, by = 4L, length.out = nAnnAge_r))
    else              max(seq(s + 1L, by = 4L, length.out = nAnnAge_r - 1L))
  }, integer(1L))

  # determine stock & fleet-specific apical Fs
  nAA_within <- nAnnAge_r - 1L   # 14 ages for S1–S3
  nAA_across <- nAnnAge_r - 2L   # 13 ages for S4

  ## ---- Seasonal Stock/Fleet apical Fs ----
  F_ref_sims <- purrr::imap(object$mcmcvars, \(it, sim) {
    m  <- object$natural_mortality[sim]
    Nf <- it$N[, , , 1L]   # [yr, ann_age, seas]  female
    Nm <- it$N[, , , 2L]   # male

    F_ref_ts <- matrix(0, nTS, nFleet)

    for (ts in seq_len(nTS - 1L)) {
      yr  <- ((ts - 1L) %/% Seasons) + 1L
      s   <- ((ts - 1L) %% Seasons) + 1L
      nA  <- if (s < Seasons) nAA_within else nAA_across
      nRow <- nA * 2L

      f_vec    <- numeric(nRow)
      sela_mat <- matrix(0, nRow, nFleet)

      for (sex in 1:2) {
        Nsex <- if (sex == 1L) Nf else Nm
        row0 <- (sex - 1L) * nA
        for (ai in seq_len(nA)) {
          n0 <- Nsex[yr, ai, s]
          n1 <- if (s < Seasons) Nsex[yr, ai, s + 1L] else Nsex[yr + 1L, ai + 1L, 1L]
          f_vec[row0 + ai]      <- if (n0 > 0 && n1 > 0) max(-log(n1 / n0) - m, 0) else 0
          sela_mat[row0 + ai, ] <- it$sela[ai, s, sex, ]
        }
      }

      F_ref_ts[ts, ] <- nnls::nnls(sela_mat, f_vec)$x
    }

    F_ref_ts[nTS, ] <- F_ref_ts[nTS - 1L, ]
    F_ref_ts   # [nTS, nFleet]
  })

  for (st in seq_along(StockNames)) {
    for (fl in seq_along(FleetNames)) {

      fleet <- MSEtool::Fleet(Name = FleetNames[fl])

      ## ---- Selectivity ----
      sel_age <- purrr::map(object$mcmcvars, \(it) {
        vapply(seq_len(nAge), \(a) {
          s     <- flat_seas_idx[a]
          ai    <- flat_age_idx[a]
          denom <- max(it$sela[, s, st, fl])
          if (denom > 0) it$sela[ai, s, st, fl] / denom else 0
        }, numeric(1L))
      })

      ## ---- Effort ----
      # Effort = F_ref[fl] × max_sela[fl, s, st] per timestep.
      # F_ref is stock-independent (NNLS uses both sexes); Effort differs per
      # stock only because each sex has its own max_sela scaling.
      Effort_mat <- do.call(rbind, purrr::imap(object$mcmcvars, \(it, sim) {
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
        MeanAtAge = MSEtool::AddDimension(sel_age_mat, 'Year', HistYears[1])
      )

      FleetList[[st]][[fl]] <- fleet
    }
  }

  MSEtool::Fleet(OM) <- FleetList

  Complexes(OM) <- list('Combined' = 1:2)

  OM
}

