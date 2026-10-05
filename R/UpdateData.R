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

  Hist@OM@Data$Combined@Landings@Value <- NewData$Landings
  Hist@OM@Data$Combined@CPUE@Value <- NewData$CPUE

  Hist@Data$`1`$Combined@Landings@Value <- NewData$Landings
  Hist@Data$`1`$Combined@CPUE@Value <- NewData$CPUE

  if (plot) {

    # ---- Catch Data ----
    catch_df <- dplyr::bind_rows(
      MSEtool::Array2DF(Histcopy@Data$`1`$Combined@Landings@Value) |> dplyr::mutate(Data = 'Conditioned'),
      MSEtool::Array2DF(Hist@Data$`1`$Combined@Landings@Value) |> dplyr::mutate(Data = 'Updated')
    )

    # Seasonal by Fleet
    c1 <- ggplot2::ggplot(catch_df, ggplot2::aes(x = .data$Year, y = .data$Value, color = .data$Data)) +
      ggplot2::facet_wrap(~Fleet, scales = 'free_y') +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::labs(y = 'Catch (t)') +
      ggplot2::expand_limits(y = 0)

    ggplot2::ggsave('figures/data/catch_seasonal_fleet.png', c1, create.dir = TRUE)

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

    ggplot2::ggsave('figures/data/catch_annual_fleet.png', c2, create.dir = TRUE)

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

    ggplot2::ggsave('figures/data/catch_annual_overall.png', c3, create.dir = TRUE)

    # ---- Index Data ----
    conditioned_index <- MSEtool::Array2DF(Histcopy@Data$`1`$Combined@CPUE@Value) |>
      dplyr::mutate(Data = 'Conditioned')
    updated_index <- MSEtool::Array2DF(Hist@Data$`1`$Combined@CPUE@Value) |>
      dplyr::mutate(Data = 'Updated')

    # standardize to mean over the conditioned years
    years <- unique(conditioned_index$Year)
    std_conditioned <- conditioned_index |>
      dplyr::filter(.data$Year %in% years) |>
      dplyr::group_by(.data$Fleet) |>
      dplyr::mutate(mean = mean(.data$Value, na.rm = TRUE)) |>
      dplyr::mutate(StValue = .data$Value / .data$mean)

    std_updated <- updated_index |>
      dplyr::filter(.data$Year %in% years) |>
      dplyr::group_by(.data$Fleet) |>
      dplyr::mutate(mean = mean(.data$Value, na.rm = TRUE)) |>
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

    ggplot2::ggsave('figures/data/index_seasonal.png', i1, create.dir = TRUE)

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

    ggplot2::ggsave('figures/data/index_annual_mean.png', i2, create.dir = TRUE)

  }

  Hist
}
