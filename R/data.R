#' Base-case operating model conditioning data for Indian Ocean Albacore
#'
#' A list of MCMC posterior samples, biological parameters, biological
#' schedules, and observed data from the Indian Ocean Albacore
#' base-case OM conditioning, used to condition the MSEtool operating model
#' (`OM_BaseCase`).
#'
#' This was built using data sent from Rich Hillary June 2026
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


#' Raw biological schedules from the albMSE conditioning data
#'
#' Inputs used to build [StockObjectList]'s length, weight, and maturity
#' schedules.
#'
#' @format A named list:
#' * `ages` — Numeric vector of annual age classes.
#' * `mula`, `sdla` — Arrays `[annual_age, season, sex]`. Mean and SD of length-at-age.
#' * `mulbins` — Numeric vector of length bin midpoints.
#' * `wta` — Array `[annual_age, season, sex]`. Mean weight-at-age (kg).
#' * `mata` — Array `[annual_age, season, sex]`. Maturity-at-age.
#'
#' @source `albMSE::boot/data/alb_abcdata.rda`
"albMSE_Biology"

#' Albacore operating model stock schedules
#'
#' Fixed biological schedules for the IOTC albacore operating model, built
#' from `albMSE` conditioning data
#' (`https://github.com/iagomosqueira/albMSE/tree/main/boot/data/alb_abcdata.rda`).
#'
#' Contains one [MSEtool::Stock()] object per sex, with seasonal age structure and
#' length, weight, maturity, and fecundity schedules. Natural mortality and the
#' stock-recruitrelationship are not set here; they are populated from MCMC
#' posterior samples.
#'
#' @format A named list of length 2 (`"Female"`, `"Male"`), each element a
#'   [MSEtool::stock-class] object with:
#'   - `Ages`: seasonal age classes spanning `min(ages)` to `max(ages) * 4`, no plus group
#'   - `Length`: mean length-at-age and CV-at-age by season
#'   - `Weight`: mean weight-at-age by season
#'   - `Maturity`: mean maturity-at-age by season
#'   - `Fecundity`: maturity-at-age times weight-at-age, zero outside the spawning season (season 3)
#'
#' @source `albMSE::boot/data/alb_abcdata.rda`
"StockObjectList"

#' MCMC-derived variables for OM5a (SW CPUE reference OM)
#'
#' Simulated population and observation quantities for each retained MCMC/ABC
#' posterior draw of the OM5a operating model
#' OM5a conditions on the south-west (SW) CPUE index
#' (`fcpue = 3`) rather than the base case north-west (NW) index used in OM5b.
#'
#' @format A list of length equal to the number of retained iterations. Each
#'   element is itself a list with components including `N` (numbers-at-age),
#'   `Rtot`, `SSB`, `dep`, `dbmsy`, `Bmsy`, `Cmsy`, `hmsy`, `hmsyrat`, `H`,
#'   `Ihat`, `LFhat`, `B0`, `R0`, `M`, `h`, and `sela`.
#'
#' @source `albMSE::data/om5a/mcvars_abc5a.rda`
"CondData_OM5a"

#' MCMC-derived variables for OM5b (base case reference OM)
#'
#' Simulated population and observation quantities for each retained MCMC/ABC
#' posterior draw of the OM5b operating model.
#' OM5b is the reference operating model: conditioned on
#' the north-west (NW) LL CPUE index, with SSB priors, recruitment
#' variability, and the overfishing penalty.
#'
#' @format A list of length equal to the number of retained iterations. Each
#'   element is itself a list with components including `N` (numbers-at-age),
#'   `Rtot`, `SSB`, `dep`, `dbmsy`, `Bmsy`, `Cmsy`, `hmsy`, `hmsyrat`, `H`,
#'   `Ihat`, `LFhat`, `B0`, `R0`, `M`, `h`, and `sela`.
#'
#' @source `albMSE::data/om5b/mcvars_abc5b.rda`
"CondData_OM5b"

#' MCMC-derived variables for OM6b (1% effort creep robustness OM)
#'
#' Simulated population and observation quantities for each retained MCMC/ABC
#' posterior draw of the OM6b operating model.
#' OM6b is a robustness variant of the reference OM
#' (OM5b) with an added 1% annual catchability trend (`qtrend = TRUE`).
#'
#' @format A list of length equal to the number of retained iterations. Each
#'   element is itself a list with components including `N` (numbers-at-age),
#'   `Rtot`, `SSB`, `dep`, `dbmsy`, `Bmsy`, `Cmsy`, `hmsy`, `hmsyrat`, `H`,
#'   `Ihat`, `LFhat`, `B0`, `R0`, `M`, `h`, and `sela`.
#'
#' @source `albMSE::data/om6b/mcvars_abc6b.rda`
"CondData_OM6b"


#' MSEtool Data object for IOTC albacore
#'
#' Observed catch and CPUE data for Indian Ocean albacore tuna (*Thunnus
#' alalunga*), reformatted from the `albMSE` ABC conditioning data
#' (`boot/data/alb_abcdata.rda`) into an [MSEtool::data-class] object for use
#' with candidate management procedures (CMPs). Covers 2000-2020 at
#' quarterly (4-season) resolution across 6 fleets.
#'
#' @format A [MSEtool::data-class] object with:
#'   - `Landings`: a [MSEtool::catchdata-class] object, quarterly catch (t)
#'     by fleet (`LL1`-`LL4`, `PS`, `Other`), `nTS = 84` timesteps
#'   - `CPUE`: a [MSEtool::indicesdata-class] object, quarterly CPUE by
#'     fleet, standardized to a mean of 1 per fleet
#'
#'
#' @source `albMSE::boot/data/alb_abcdata.rda` (`C`, `I` objects)
"albMSE_Data"


