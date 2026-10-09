#' Annual Stock Status Relative to MSY Reference Points
#'
#' Annual spawning biomass and fishing mortality relative to the MSY reference
#' points, from [MSEtool::SB_SBMSY()] and [MSEtool::F_FMSY()].
#'
#' The OMs are seasonal (4 seasons). The MSY reference points are annual, and
#' MSEtool reports the ratios on the same basis, one value per calendar year:
#' `F_FMSY` is the apical fishing mortality summed over the seasons of the
#' year, and `SB_SBMSY` the spawning biomass in the reference season (see
#' [MSEtool::RefSeason()]). This requires a version of MSEtool that reports
#' annual ratios for seasonal OMs; an error is given if the ratios are per
#' time step.
#'
#' @param object An `mse` or `hist` object.
#' @param Stock Character. Name of the spawning stock. Default `'Female'`.
#'
#' @return A `data.frame` with columns `Sim`, `Year` (calendar year),
#'   `Period`, `MP`, `SB_SBMSY`, and `F_FMSY`. Historical rows have
#'   `MP = 'Historical'`.
#' @export
AnnualStatus <- function(object, Stock = 'Female') {
  Prep <- function(df) {
    df <- as.data.frame(df)
    if (!'MP' %in% names(df)) df$MP <- 'Historical'
    if (any(abs(df$Year - round(df$Year)) > 1e-6))
      cli::cli_abort(c('{.fn MSEtool::{unique(df$Variable)}} returned values per time step.',
                       'i' = 'Update MSEtool to a version that reports annual MSY ratios for seasonal OMs.'))
    df$Sim  <- as.integer(as.character(df$Sim))
    df$Year <- round(df$Year)
    df[, c('Sim', 'Year', 'Period', 'MP', 'Value')]
  }

  SB <- MSEtool::SB_SBMSY(object, df = TRUE, Reduce = FALSE)
  SB <- Prep(SB[SB$Stock %in% Stock, ])
  FF <- Prep(MSEtool::F_FMSY(object, df = TRUE, Reduce = FALSE))

  out <- merge(SB, FF, by = c('Sim', 'Year', 'Period', 'MP'), suffixes = c('_SB', '_F'))
  names(out)[names(out) == 'Value_SB'] <- 'SB_SBMSY'
  names(out)[names(out) == 'Value_F']  <- 'F_FMSY'
  out[order(out$MP, out$Sim, out$Year), ]
}
