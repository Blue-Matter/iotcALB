
#' Compare OM reference quantities against the conditioning model
#'
#' Compares an [MSEtool::Hist()] simulation against the reference quantities
#' recovered from the MCMC posterior of the conditioning model: total
#' numbers, spawning biomass, catch-by-fleet, and numbers-at-age.
#'
#' For total numbers, spawning biomass, and catch, a ribbon plot (median with a
#' 10th-90th percentile band, OM vs. reference) and a per-simulation ratio
#' plot (OM / reference) are produced.
#'
#' Numbers-at-age is compared by mapping
#' each MCMC iteration's N-at-age into the OM's quarterly age x timestep
#' layout; the distribution of OM/reference ratios is summarised as a ribbon
#' by quarterly age class and season, and cell-level diagnostics are printed.
#'
#' All plots are optionally saved to disk.
#'
#' @param Hist A [MSEtool::Hist()] object produced by [MSEtool::Simulate()].
#'   The conditioning data used to build `Hist@OM` is loaded via
#'   `get(paste0('CondData_', Hist@OM@Name))` (same object used by
#'   [ImportOM()]).
#' @param tol Numeric. Ratio threshold outside `[1-tol, 1+tol]` that triggers
#'   a printed warning for the N-at-age comparison. Defaults to `0.02` (2%).
#' @param min_pct Numeric. Minimum fraction of N-at-age cells (0-1) that must
#'   exceed `tol` before the worst-offenders table is printed. Defaults to
#'   `0.05` (5%).
#' @param alpha Numeric. Transparency used for ribbon fills and
#'   per-simulation ratio lines. Defaults to `0.2`.
#' @param probs Numeric vector of length 2 giving the lower/upper quantiles
#'   used for ribbon plots. Defaults to `c(0.1, 0.9)`.
#' @param save_plots Logical. If `TRUE` (default), diagnostic plots are
#'   written as PNG files to `file.path(outdir, Hist@OM@Name)`.
#' @param outdir Character. Base directory for saved plots. Defaults to
#'   `"figures/diagnostics/OM"`.
#' @param width,height Numeric. Width/height (inches) passed to
#'   [ggplot2::ggsave()]. If `NULL` (default), each plot is sized
#'   automatically from its number of facet panels: single-panel plots
#'   (Number, SSB) use 6 x 4; the Catch-by-fleet and N-at-age plots, which
#'   facet by fleet/season, scale up accordingly. Set either to a number to
#'   use that fixed size for every saved plot instead.
#' @param verbose Logical. If `TRUE` (default), print N-at-age match
#'   diagnostics (success/warning messages and the worst-offenders table).
#'
#' @return Invisibly returns a named list:
#'   * `n`, `ssb`, `catch`: each a list with `ref` and `om` data frames
#'   * `natage`: long data frame of OM vs reference N-at-age with columns
#'     `Stock`, `Sim`, `Age`, `Timestep`, `Season`, `Year`, `OM`, `Ref`,
#'     `Ratio`
#'   * `plots`: a list of the diagnostic `ggplot` objects — `number_ribbon`,
#'     `number_ratio`, `ssb_ribbon`, `ssb_ratio`, `catch_ribbon`,
#'     `catch_ratio`, `natage_ratio`
#'
#' @export
ValidateOM <- function(Hist,
                       tol = 0.05, min_pct = 0.02,
                       alpha = 0.2, probs = c(0.1, 0.9), save_plots = TRUE,
                       outdir = 'figures/diagnostics/OM', width = NULL, height = NULL,
                       verbose = TRUE) {

  if (!inherits(Hist, 'hist'))
    cli::cli_abort('`Hist` must be an {.help MSEtool::Hist} object')

  OMName <- Hist@OM@Name
  object <- get(paste0('CondData_', OMName))

  req_mcmc_fields <- c('N', 'H', 'sela', 'SSB')
  missing_fields  <- setdiff(req_mcmc_fields, names(object[[1]]))
  if (length(missing_fields))
    cli::cli_abort('Each element of {.arg object} must contain: {.field {req_mcmc_fields}}')

  Seasons    <- 4L
  Years      <- MSEtool::Years(Hist, 'H')
  CalYears   <- floor(Years) |> unique()
  nYear      <- MSEtool::nYear(Hist)
  AgeClasses <- MSEtool::Classes(Hist@OM@Stock$Female)
  nAge       <- length(AgeClasses)
  nAnnAge    <- floor(AgeClasses) |> unique() |> length()
  nTS        <- length(Years)
  nSim       <- min(MSEtool::nSim(Hist), length(object))
  StockNames <- MSEtool::StockNames(Hist)
  FleetNames <- MSEtool::FleetNames(Hist)
  nSex       <- length(StockNames)
  nFleet     <- length(FleetNames)

  figdir <- file.path(outdir, OMName)

  # `ImportOM()` feeds `SRR(R0 = ..., Units = 1000)`, so `Hist@Number` (and
  # anything derived from it: SProduction, Removals) is in thousands of
  # fish, with biomass (SSB, Catch) coming out directly in tonnes because
  # weight-at-age is in kg (thousands of fish x kg = tonnes). The reference
  # `object` is in raw individual fish counts, so its N (and anything
  # derived from N, i.e. catch) is rescaled to match; SSB needs no rescaling
  # since it's already reported in tonnes.
  n_scale <- 1 / 1000

  # ---- OM quantities ----
  om_n <- MSEtool::Number(Hist) |>
    dplyr::mutate(Sim = as.integer(Sim),
                  Year = floor(as.numeric(as.character(Year)))) |>
    dplyr::group_by(Sim, Year) |>
    dplyr::summarise(N = sum(Value), .groups = 'drop')

  om_ssb <- MSEtool::SProduction(Hist) |>
    dplyr::filter(Stock == 'Female') |>
    dplyr::mutate(Year = floor(Year)) |>
    dplyr::group_by(Sim, Year) |>
    dplyr::summarise(SSB = sum(Value), .groups = 'drop')

  om_catch <- MSEtool::Removals(Hist) |>
    dplyr::mutate(Year = floor(Year)) |>
    dplyr::group_by(Sim, Year, Fleet) |>
    dplyr::summarise(Catch = sum(Value), .groups = 'drop')

  # ---- Reference quantities across all MCMC sims ----
  ref_list <- purrr::imap(object, \(it, sim) {

    it$N <- it$N * n_scale

    n   <- data.frame(Sim = sim, Year = CalYears, N = apply(it$N, 1, sum))
    ssb <- data.frame(Sim = sim, Year = CalYears, SSB = it$SSB)

    catch_fl <- matrix(0, nYear, nFleet)
    for (s in seq_len(Seasons))
      for (sex in seq_len(nSex)) {
        sw       <- it$sela[, s, sex, ] * albMSE_Biology$wta[, s, sex]
        catch_fl <- catch_fl + it$N[, , s, sex] %*% sw * it$H[, s, ]
      }

    catch <- data.frame(Sim   = sim,
                        Year  = rep(CalYears, nFleet),
                        Fleet = rep(FleetNames, each = nYear),
                        Catch = as.vector(catch_fl))

    list(n = n, ssb = ssb, catch = catch)
  })

  ref_n     <- do.call(rbind, lapply(ref_list, `[[`, 'n'))
  ref_ssb   <- do.call(rbind, lapply(ref_list, `[[`, 'ssb'))
  ref_catch <- do.call(rbind, lapply(ref_list, `[[`, 'catch'))

  # ---- Ribbon / ratio plots: Number, SSB, Catch ----
  n_diag <- .validateOM_diag_plots(ref_n, om_n, value_col = 'N', by = 'Year',
                                   ylab = "Numbers ('000s)", label = 'Total Numbers',
                                   probs = probs, alpha = alpha)

  ssb_diag <- .validateOM_diag_plots(ref_ssb, om_ssb, value_col = 'SSB', by = 'Year',
                                     ylab = 'SSB (t)', label = 'Spawning Biomass',
                                     probs = probs, alpha = alpha)

  catch_diag <- .validateOM_diag_plots(ref_catch, om_catch, value_col = 'Catch',
                                       by = c('Year', 'Fleet'), ylab = 'Catch (t)',
                                       label = 'Catch by Fleet', probs = probs,
                                       alpha = alpha, facet = 'Fleet',
                                       facet_levels = FleetNames)

  if (save_plots) {
    .validateOM_save(n_diag$ribbon,     figdir, 'Number_ribbon.png', width, height, n_diag$dims)
    .validateOM_save(n_diag$ratio,      figdir, 'Number_ratio.png',  width, height, n_diag$dims)
    .validateOM_save(ssb_diag$ribbon,   figdir, 'SSB_ribbon.png',    width, height, ssb_diag$dims)
    .validateOM_save(ssb_diag$ratio,    figdir, 'SSB_ratio.png',     width, height, ssb_diag$dims)
    .validateOM_save(catch_diag$ribbon, figdir, 'Catch_ribbon.png',  width, height, catch_diag$dims)
    .validateOM_save(catch_diag$ratio,  figdir, 'Catch_ratio.png',   width, height, catch_diag$dims)
  }

  # ---- Numbers-at-age ----

  # Map reference N[yr, annual_age, seas, st] into OM layout [age_pos, ts]
  ref_mat <- function(it, st) {
    N_ref <- it$N[, , , st] * n_scale
    mat   <- matrix(0, nAge, nTS)
    for (yr in seq_len(nYear)) {
      for (s in seq_len(Seasons)) {
        ts <- (yr - 1L) * Seasons + s
        if (s == Seasons) {
          sind          <- seq(1L, by = 4L, length.out = nAnnAge)
          mat[sind, ts] <- N_ref[yr, , s]
        } else {
          sind          <- seq(s + 1L, by = 4L, length.out = nAnnAge - 1L)
          mat[sind, ts] <- N_ref[yr, 2:nAnnAge, s]
        }
      }
    }
    mat
  }

  cli::cli_progress_bar("Comparing N-at-age", total = nSim * length(StockNames))

  # Build long data frame across all sims and both stocks
  df_list <- vector('list', nSim * length(StockNames))
  idx <- 1L
  for (sim in seq_len(nSim)) {
    it <- object[[sim]]
    for (sti in seq_along(StockNames)) {
      om_arr <- Hist@Number[[StockNames[sti]]][sim, , , 1L]  # [age, ts]
      ref    <- ref_mat(it, sti)
      keep   <- which(ref > 0, arr.ind = TRUE)
      df_list[[idx]] <- data.frame(
        Stock    = StockNames[sti],
        Sim      = sim,
        Age      = AgeClasses[keep[, 1L]],
        Timestep = as.integer(keep[, 2L]),
        OM       = om_arr[keep],
        Ref      = ref[keep]
      )
      idx <- idx + 1L
      cli::cli_progress_update()
    }
  }
  cli::cli_progress_done()

  df          <- do.call(rbind, df_list)
  df$Season   <- ((df$Timestep - 1L) %% Seasons) + 1L
  df$Year     <- 1999L + ((df$Timestep - 1L) %/% Seasons) + 1L
  df$Ratio    <- df$OM / df$Ref
  df          <- df[is.finite(df$Ratio), ]

  #  Diagnostics
  bad <- df[abs(df$Ratio - 1) > tol, ]

  frac_bad <- nrow(bad) / nrow(df)
  pct_bad  <- round(frac_bad * 100, 2)

  if (verbose) {
    if (frac_bad < min_pct) {
      cli::cli_alert_success(
        "N-at-age match: {100 - pct_bad}% of cells within {tol * 100}% of reference \\
        across all {nSim} simulation{?s}."
      )
    } else {
      cli::cli_alert_warning(
        "{pct_bad}% of age×timestep×sim cells deviate >{tol * 100}% from reference \\
        ({nrow(bad)} of {nrow(df)})."
      )

      bad$AbsDev <- abs(bad$Ratio - 1)
      agg_bad <- aggregate(
        cbind(MeanRatio = Ratio, MeanAbsDev = AbsDev) ~ Stock + Season + Age,
        data = bad,
        FUN  = mean
      )
      agg_bad <- agg_bad[order(agg_bad$MeanAbsDev, decreasing = TRUE), ]
      agg_bad$MeanRatio  <- round(agg_bad$MeanRatio,  4)
      agg_bad$MeanAbsDev <- round(agg_bad$MeanAbsDev, 4)

      cli::cli_h2("Worst offenders (top 20 by mean |ratio - 1|)")
      print(utils::head(agg_bad, 20L), row.names = FALSE)
    }
  }

  # plot
  agg <- aggregate(Ratio ~ Stock + Age + Season + Year + Sim, data = df, FUN = mean)

  ribbon <- do.call(rbind, Filter(
    Negate(is.null),
    lapply(
      split(agg, list(agg$Stock, agg$Age, agg$Season), drop = TRUE),
      function(x) {
        if (nrow(x) == 0L) return(NULL)
        data.frame(
          Stock  = x$Stock[1L],
          Age    = x$Age[1L],
          Season = paste0('S', x$Season[1L]),
          lo     = stats::quantile(x$Ratio, probs[1]),
          mid    = stats::quantile(x$Ratio, 0.50),
          hi     = stats::quantile(x$Ratio, probs[2])
        )
      }
    )
  ))
  ribbon$Stock <- factor(ribbon$Stock, levels = StockNames)

  max_dev <- max(abs(c(ribbon$lo, ribbon$hi) - 1), tol, na.rm = TRUE)
  y_pad   <- max_dev * 0.15
  y_lim   <- 1 + c(-1, 1) * (max_dev + y_pad)

  p <- ggplot2::ggplot(ribbon, ggplot2::aes(x = Age)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = alpha, fill = 'steelblue', colour = NA) +
    ggplot2::geom_line(ggplot2::aes(y = mid), colour = 'steelblue') +
    ggplot2::geom_hline(yintercept = 1, linetype = 'dashed', colour = 'grey30') +
    ggplot2::geom_hline(yintercept = c(1 - tol, 1 + tol),
                        linetype = 'dotted', colour = 'firebrick', alpha = 0.6) +
    ggplot2::facet_grid(Stock ~ Season) +
    ggplot2::coord_cartesian(ylim = y_lim) +
    ggplot2::labs(
      title    = paste0('OM / ABC N-at-age  (', nSim, ' simulation', if (nSim > 1) 's' else '', ')'),
      subtitle = paste0('Line = median  |  Ribbon = ', probs[1] * 100, '-', probs[2] * 100,
                        'th percentile across sims  |  Dotted = ±', tol * 100, '% tolerance'),
      x        = 'Quarterly age class',
      y        = 'OM / ABC ratio'
    ) +
    ggplot2::theme_bw()

  if (save_plots)
    .validateOM_save(p, figdir, 'NatAge_ratio.png', width, height,
                     list(nrow = length(StockNames), ncol = Seasons))

  invisible(list(
    n      = list(ref = ref_n,     om = om_n),
    ssb    = list(ref = ref_ssb,   om = om_ssb),
    catch  = list(ref = ref_catch, om = om_catch),
    natage = df,
    plots  = list(
      number_ribbon = n_diag$ribbon,
      number_ratio  = n_diag$ratio,
      ssb_ribbon    = ssb_diag$ribbon,
      ssb_ratio     = ssb_diag$ratio,
      catch_ribbon  = catch_diag$ribbon,
      catch_ratio   = catch_diag$ratio,
      natage_ratio  = p
    )
  ))
}


.validateOM_diag_plots <- function(ref, om, value_col, by, ylab, label, probs, alpha,
                                   facet = NULL, facet_levels = NULL) {

  join_cols <- c('Sim', by)

  summarise_q <- function(d, source) {
    d |>
      dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
      dplyr::summarise(
        lo  = stats::quantile(.data[[value_col]], probs[1]),
        mid = stats::median(.data[[value_col]]),
        hi  = stats::quantile(.data[[value_col]], probs[2]),
        .groups = 'drop'
      ) |>
      dplyr::mutate(Source = source)
  }

  ribbon_df <- dplyr::bind_rows(summarise_q(ref, 'ABC'), summarise_q(om, 'OM'))

  ratio_df <- dplyr::inner_join(ref, om, by = join_cols, suffix = c('_ref', '_om'))
  ratio_df$Ratio <- ratio_df[[paste0(value_col, '_om')]] / ratio_df[[paste0(value_col, '_ref')]]
  ratio_df <- ratio_df[is.finite(ratio_df$Ratio), ]

  if (!is.null(facet)) {
    ribbon_df[[facet]] <- factor(ribbon_df[[facet]], levels = facet_levels)
    ratio_df[[facet]]  <- factor(ratio_df[[facet]],  levels = facet_levels)
  }

  ribbon_p <- ggplot2::ggplot(ribbon_df, ggplot2::aes(x = Year, colour = Source, fill = Source)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = alpha, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(y = mid)) +
    ggplot2::scale_colour_manual(values = c(ABC = 'firebrick', OM = 'steelblue')) +
    ggplot2::scale_fill_manual(values   = c(ABC = 'firebrick', OM = 'steelblue')) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(title = paste0(label, ': ABC vs OM'),
                  subtitle = paste0('Line = median  |  Ribbon = ', probs[1] * 100, '-',
                                    probs[2] * 100, 'th percentile'),
                  y = ylab, x = NULL) +
    ggplot2::theme_bw()

  half_width <- max(stats::quantile(abs(ratio_df$Ratio - 1), 0.99, na.rm = TRUE), 1e-3)
  ylim       <- 1 + c(-1, 1) * half_width * 1.05

  ratio_p <- ggplot2::ggplot(ratio_df, ggplot2::aes(x = Year, y = Ratio, group = Sim)) +
    ggplot2::geom_line(alpha = alpha, linewidth = 0.3, colour = 'steelblue') +
    ggplot2::geom_hline(yintercept = 1, colour = 'firebrick', linetype = 'dashed') +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(title = paste0(label, ': OM / ABC by simulation'),
                  y = 'Ratio (OM / ABC)', x = NULL) +
    ggplot2::theme_bw()

  dims <- if (is.null(facet)) list(nrow = 1L, ncol = 1L) else .validateOM_wrap_dims(length(facet_levels))

  if (!is.null(facet)) {
    ribbon_p <- ribbon_p + ggplot2::facet_wrap(stats::reformulate(facet), scales = 'free_y',
                                               nrow = dims$nrow, ncol = dims$ncol)
    ratio_p  <- ratio_p  + ggplot2::facet_wrap(stats::reformulate(facet),
                                               nrow = dims$nrow, ncol = dims$ncol)
  }

  list(ribbon = ribbon_p, ratio = ratio_p, dims = dims)
}

# Near-square (nrow, ncol) grid for `n` facet panels, matching ggplot2's own
# facet_wrap layout heuristic.
.validateOM_wrap_dims <- function(n) {
  ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n / ncol)
  list(nrow = nrow, ncol = ncol)
}

# Plot size (inches) for a `dims$nrow` x `dims$ncol` panel grid: 6 x 4 for a
# single panel, scaling up per additional row/column so faceted plots (e.g.
# Catch by fleet, N-at-age by stock/season) don't come out squashed.
.validateOM_auto_dims <- function(dims) {
  list(width  = 3.5 + dims$ncol * 2.5,
       height = 2   + dims$nrow * 2)
}

.validateOM_save <- function(plot, dir, filename, width, height, dims = list(nrow = 1L, ncol = 1L)) {
  if (is.null(width) || is.null(height)) {
    auto <- .validateOM_auto_dims(dims)
    if (is.null(width))  width  <- auto$width
    if (is.null(height)) height <- auto$height
  }
  ggplot2::ggsave(file.path(dir, filename), plot, width = width, height = height, create.dir = TRUE)
  invisible(NULL)
}
