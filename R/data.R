#' Base-case operating model conditioning data for Indian Ocean Albacore
#'
#' A list of MCMC posterior samples, biological parameters, biological
#' schedules, and observed data from the Indian Ocean Albacore
#' base-case OM conditioning, used to condition the MSEtool operating model
#' (`OM_BaseCase`).
#'
#' @format A named list with the following elements:
#'
#' * `mcmcvars` — List of 500 MCMC iterations. Each element contains derived
#'   quantities: `N` (numbers-at-age), `H` (harvest rates), `sela`
#'   (selectivity-at-age), `SSB`, `dep`, and others.
#' * `steepness` — Numeric vector (500). Beverton-Holt steepness *h*.
#' * `natural_mortality` — Numeric vector (500). Quarterly natural mortality *M*.
#' * `R0` — Numeric vector (500). Unfished recruitment.
#' * `sigmaR` — Numeric vector (500). Recruitment standard deviation.
#' * `maturity` — Array `[annual_age, season, sex]`. Maturity-at-age by
#'   season and sex.
#' * `weight` — Array `[annual_age, season, sex]`. Mean weight-at-age (kg) by
#'   season and sex.
#' * `ALK` — Array `[length, annual_age, season, sex]`. Age-length key
#'   (proportion-at-length given age).
#' * `Index` — Array `[year, season, fleet]`. Longline CPUE indices.
#' * `LengthComp` — Length composition data by year, class, season and fleet.
#'
"Cond_BaseCase"
