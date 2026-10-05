#' Process new catch and CPUE data from an SS3 data file
#'
#' Reads an SS3 data file with [r4ss::SS_readdat()] and aggregates the
#' quarterly catch and CPUE observations from 2000 onwards into the fleet
#' structure used by the operating models.
#'
#' Catches are summed across the SS3 fleets within each OM fleet:
#' `LL1` (fleets 1-4), `LL2` (5-8), `LL3` (9-12), `LL4` (13-16), `PS` (19),
#' and `Other` (17-18, 20-23).
#'
#' CPUE indices are taken from SS3 index fleets 24-27 (`LL1`), 28-31 (`LL2`),
#' 32-35 (`LL3`), and 36-39 (`LL4`). Seasons with no observation are `NA`.
#'
#' Time steps are labelled with decimal years from [MSEtool::CalcYears()]
#' (4 seasons per year, current year 2023).
#'
#' @param file Path to the SS3 data file. Defaults to
#'   `'data-raw/new_data.dat'`, relative to the package source directory.
#'
#' @return A named list with two matrices, each with dimensions
#'   `c(Year, Fleet)`:
#' \describe{
#'   \item{Landings}{Total catch by season and fleet (`LL1`, `LL2`, `LL3`,
#'   `LL4`, `PS`, `Other`).}
#'   \item{CPUE}{CPUE index by season and longline fleet (`LL1`-`LL4`).}
#' }
#'
#' @export
ProcessNewData <- function(file = 'data-raw/new_data.dat') {
  newdata <- r4ss::SS_readdat(file, verbose = FALSE)

  # Catch
  FleetGroups <- list(LL1 = 1:4,
                      LL2 = 5:8,
                      LL3 = 9:12,
                      LL4 = 13:16,
                      PS  = 19,
                      Other = c(17:18, 20:23)
  )

  yrs <- newdata$catch$year
  years <- unique(yrs[yrs>=2000])
  Years <- MSEtool::CalcYears(nYear = length(years),
                              pYear = 0,
                              CurrentYear = 2023,
                              Seasons = 4)


  catch <- purrr::imap(FleetGroups, \(fleets, i) {

    mat <- newdata$catch |> dplyr::filter(.data$fleet %in% fleets,
                                          .data$year %in% years) |>
      dplyr::arrange(.data$year, .data$seas) |>
      dplyr::group_by(.data$year, .data$seas) |>
      dplyr::summarise(C = sum(.data$catch), .groups = 'drop') |>
      dplyr::select('C') |>
      as.matrix()
    mat
  }) |> abind::abind(along = 2, use.dnns = TRUE)

  dimnames(catch) <- list(
    Year = Years,
    Fleet = names(FleetGroups)
  )

  # CPUE
  FleetGroups <- list(LL1 = 24:27,
                      LL2 = 28:31,
                      LL3 = 32:35,
                      LL4 = 36:39)

  index <- matrix(NA, nrow = length(Years), ncol = 4,
                  dimnames = list(
                    Year = Years,
                    Fleet = names(FleetGroups)
                  )
  )

  months <- unique(newdata$CPUE$month)

  for (i in seq_along(FleetGroups)) {
    fleets <- FleetGroups[[i]]
    for (j in seq_along(years)) {
      yr <- years[j]
      yr_ind <- match(yr, Years)
      yr_ind <- yr_ind:(yr_ind+3)
      mat <- newdata$CPUE |> dplyr::filter(.data$index %in% fleets,
                                           .data$year %in% yr) |>
        dplyr::arrange(.data$year, .data$month)

      val <- matrix(NA, 4,1, dimnames = list(Year = Years[yr_ind],
                                             Fleet = names(FleetGroups)[i])
      )
      val[match(mat$month, months),1] <- mat$obs
      abind::afill(index) <- val
    }
  }

  list(Landings = catch, CPUE = index)
}
