#' Update a simulated Hist object with new catch and CPUE data
#'
#' Replaces the historical landings and CPUE index data in a `Hist` object
#' (from [MSEtool::Simulate()]) with the updated observations returned by
#' [ProcessNewData()].
#'
#' This must be called *after* the observation model (`Obs`) has been
#' populated in [MSEtool::Simulate()] using the catch and index data that
#' were used in conditioning. Otherwise the differences between the
#' conditioned and updated data would be interpreted as observation error
#' in the observed catches.
#'
#' The data are replaced in both the OM data (`Hist@OM@Data`) and the
#' simulated data (`Hist@Data`) for the combined stock.
#'
#' If the fleets of the OM have been combined with [MSEtool::CombineFleets()],
#' the new landings are summed over fleets, and the CPUE indices, which
#' [MSEtool::CombineFleets()] moves to the `Survey` slot, are replaced there.
#'
#' @param Hist A `Hist` object returned by [MSEtool::Simulate()].
#' @param plot Logical. If `TRUE`, plots comparing the conditioned and updated
#'   data are saved to `figures/data/`:
#' \describe{
#'   \item{`catch_seasonal_fleet.png`}{Seasonal catch by fleet.}
#'   \item{`catch_annual_fleet.png`}{Annual catch by fleet.}
#'   \item{`catch_annual_overall.png`}{Annual total catch.}
#'   \item{`index_seasonal.png`}{Seasonal CPUE index by fleet, standardized
#'   to the mean over the conditioned years.}
#'   \item{`index_annual_mean.png`}{Annual mean of the standardized CPUE
#'   index by fleet.}
#' }
#'
#' @return The `Hist` object with updated landings and CPUE data.
#'
#' @export
UpdateData <- function(Hist, plot = FALSE) {

  Histcopy <- Hist
  NewData <- ProcessNewData()

  Landings   <- .MatchLandingsFleets(NewData$Landings, Hist@OM@Data$Combined)
  IndexSlot  <- .IndexSlot(Hist@OM@Data$Combined, colnames(NewData$CPUE))

  Hist@OM@Data$Combined@Landings@Value <- Landings
  methods::slot(Hist@OM@Data$Combined, IndexSlot)@Value <- NewData$CPUE

  Hist@Data$`1`$Combined@Landings@Value <- Landings
  methods::slot(Hist@Data$`1`$Combined, IndexSlot)@Value <- NewData$CPUE

  if (plot)
    .PlotDataUpdate(Histcopy@Data$`1`$Combined, Hist@Data$`1`$Combined, IndexSlot)

  Hist
}

# Plots comparing the conditioned and updated catch and CPUE data (two
# MSEtool `data` objects with the same fleets), saved to figures/data/. The
# updated index is scaled to its mean over the years of the conditioned index
.PlotDataUpdate <- function(Conditioned, Updated, IndexSlot) {

  # ---- Catch Data ----
  catch_df <- dplyr::bind_rows(
    MSEtool::Array2DF(Conditioned@Landings@Value) |> dplyr::mutate(Data = 'Conditioned'),
    MSEtool::Array2DF(Updated@Landings@Value) |> dplyr::mutate(Data = 'Updated')
  )

  # Seasonal by Fleet
  c1 <- ggplot2::ggplot(catch_df, ggplot2::aes(x = .data$Year, y = .data$Value, color = .data$Data)) +
    ggplot2::facet_wrap(~Fleet, scales = 'free_y') +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::labs(y = 'Catch (t)') +
    ggplot2::expand_limits(y = 0)

  ggplot2::ggsave('figures/data/catch_seasonal_fleet.png', c1, width = 14, height = 7, create.dir = TRUE)

  # Annual by Fleet
  catch_df2 <- catch_df |>
    dplyr::mutate(Year = floor(.data$Year)) |>
    dplyr::group_by(.data$Year, .data$Fleet, .data$Data) |>
    dplyr::summarize(Value = sum(.data$Value), .groups = 'drop_last')

  c2 <- ggplot2::ggplot(catch_df2, ggplot2::aes(x = .data$Year, y = .data$Value, color = .data$Data)) +
    ggplot2::facet_wrap(~Fleet, scales = 'free_y') +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::labs(y = 'Catch (t)') +
    ggplot2::expand_limits(y = 0)

  ggplot2::ggsave('figures/data/catch_annual_fleet.png', c2, width = 14, height = 7, create.dir = TRUE)

  # Annual Overall
  catch_df3 <- catch_df |>
    dplyr::mutate(Year = floor(.data$Year)) |>
    dplyr::group_by(.data$Year, .data$Data) |>
    dplyr::summarize(Value = sum(.data$Value), .groups = 'drop_last')

  c3 <- ggplot2::ggplot(catch_df3, ggplot2::aes(x = .data$Year, y = .data$Value, color = .data$Data)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::labs(y = 'Catch (t)') +
    ggplot2::expand_limits(y = 0)

  ggplot2::ggsave('figures/data/catch_annual_overall.png', c3, width = 9, height = 5, create.dir = TRUE)

  # ---- Index Data ----
  conditioned_index <- MSEtool::Array2DF(methods::slot(Conditioned, IndexSlot)@Value) |>
    dplyr::mutate(Data = 'Conditioned')
  updated_index <- MSEtool::Array2DF(methods::slot(Updated, IndexSlot)@Value) |>
    dplyr::mutate(Data = 'Updated')

  # standardize to mean over the conditioned years
  years <- unique(conditioned_index$Year)
  std_conditioned <- conditioned_index |>
    dplyr::filter(.data$Year %in% years) |>
    dplyr::group_by(.data$Fleet) |>
    dplyr::mutate(mean = mean(.data$Value, na.rm = TRUE)) |>
    dplyr::mutate(StValue = .data$Value / .data$mean)

  # all years of the updated index, scaled to its mean over the conditioned years
  std_updated <- updated_index |>
    dplyr::group_by(.data$Fleet) |>
    dplyr::mutate(mean = mean(.data$Value[.data$Year %in% years], na.rm = TRUE)) |>
    dplyr::mutate(StValue = .data$Value / .data$mean)

  index_df <- dplyr::bind_rows(std_conditioned, std_updated)

  # Seasonal by Fleet
  i1 <- ggplot2::ggplot(index_df, ggplot2::aes(x = .data$Year, y = .data$StValue, color = .data$Data)) +
    ggplot2::facet_wrap(~Fleet, scales = 'free_y') +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::labs(y = 'Index') +
    ggplot2::expand_limits(y = 0)

  ggplot2::ggsave('figures/data/index_seasonal.png', i1, width = 14, height = 7, create.dir = TRUE)

  # Annual by Fleet
  index_df2 <- index_df |>
    dplyr::mutate(Year = floor(.data$Year)) |>
    dplyr::group_by(.data$Year, .data$Fleet, .data$Data) |>
    dplyr::summarize(StValue = mean(.data$StValue, na.rm = TRUE), .groups = 'drop_last')

  i2 <- ggplot2::ggplot(index_df2, ggplot2::aes(x = .data$Year, y = .data$StValue, color = .data$Data)) +
    ggplot2::facet_wrap(~Fleet, scales = 'free_y') +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::labs(y = 'Index') +
    ggplot2::expand_limits(y = 0)

  ggplot2::ggsave('figures/data/index_annual_mean.png', i2, width = 14, height = 7, create.dir = TRUE)

  invisible(NULL)
}

# Sum the new landings over fleets when the OM fleets have been combined
.MatchLandingsFleets <- function(Landings, Data) {
  OMFleets <- colnames(Data@Landings@Value)
  if (identical(OMFleets, colnames(Landings)))
    return(Landings)
  if (length(OMFleets) != 1)
    cli::cli_abort('New landings have fleets {.val {colnames(Landings)}} but the OM data has {.val {OMFleets}}.')
  array(rowSums(Landings, na.rm = TRUE), dim = c(nrow(Landings), 1),
        dimnames = list(Year = rownames(Landings), Fleet = OMFleets))
}

# Slot holding the CPUE indices: `CPUE`, or `Survey` after [MSEtool::CombineFleets()]
.IndexSlot <- function(Data, IndexNames) {
  SurveyNames <- c(Data@Survey@Name, colnames(Data@Survey@Value))
  if (all(IndexNames %in% SurveyNames)) 'Survey' else 'CPUE'
}
