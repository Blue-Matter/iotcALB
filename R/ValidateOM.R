
#' Compare OM N-at-age against conditioning model reference
#'
#' Maps each MCMC iteration's N-at-age into the OM's quarterly age × timestep
#' layout and compares it against the simulated N-at-age from a
#' [MSEtool::Hist()] object. Diagnostics are computed across all simulations;
#' the plot summarises the distribution of OM/Reference ratios as a ribbon.
#'
#' @param hist A [MSEtool::Hist()] object produced by [MSEtool::Simulate()].
#' @param object A named list of conditioning output (same object passed to
#'   [ImportOM()]). Defaults to [Cond_BaseCase].
#' @param tol Numeric. Ratio threshold outside `[1-tol, 1+tol]` that triggers
#'   a printed warning. Defaults to `0.02` (2%).
#' @param min_pct Numeric. Minimum fraction of cells (0–1) that must exceed
#'   `tol` before the worst-offenders table is printed. Defaults to `0.02`
#'   (2%).
#'
#' @return Invisibly returns a data frame of OM vs reference N-at-age with
#'   columns `Stock`, `Sim`, `Age`, `Timestep`, `Season`, `Year`, `OM`, `Ref`,
#'   `Ratio`.
#'
#' @export
ValidateOM <- function(hist, object = Cond_BaseCase, tol = 0.02, min_pct = 0.02) {

  Seasons <- 4L
  nYear   <- dim(object$mcmcvars[[1]]$N)[1]
  nAnnAge <- dim(object$mcmcvars[[1]]$N)[2]
  nTS     <- nYear * Seasons
  nSim    <- min(MSEtool::nSim(hist), length(object$mcmcvars))

  AgeClasses <- as.numeric(dimnames(hist@Number$Female)[[2]])
  nAge       <- length(AgeClasses)

  StockNames <- c('Female', 'Male')

  # Map reference N[yr, annual_age, seas, st] into OM layout [age_pos, ts]
  ref_mat <- function(it, st) {
    N_ref <- it$N[, , , st]
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
    it <- object$mcmcvars[[sim]]
    for (sti in seq_along(StockNames)) {
      om_arr <- hist@Number[[StockNames[sti]]][sim, , , 1L]  # [age, ts]
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
  df$Season   <- factor(paste0('S', ((df$Timestep - 1L) %% Seasons) + 1L))
  df$Year     <- ((df$Timestep - 1L) %/% Seasons) + 1L
  df$Ratio    <- df$OM / df$Ref
  df          <- df[is.finite(df$Ratio), ]

  #  Diagnostics
  bad <- df[abs(df$Ratio - 1) > tol, ]

  frac_bad <- nrow(bad) / nrow(df)
  pct_bad  <- round(frac_bad * 100, 2)

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

  # plot
  agg <- aggregate(Ratio ~ Stock + Age + Season + Sim, data = df, FUN = mean)

  ribbon <- do.call(rbind, Filter(
    Negate(is.null),
    lapply(
      split(agg, list(agg$Stock, agg$Age, agg$Season), drop = TRUE),
      function(x) {
        if (nrow(x) == 0L) return(NULL)
        data.frame(
          Stock  = x$Stock[1L],
          Age    = x$Age[1L],
          Season = x$Season[1L],
          lo     = stats::quantile(x$Ratio, 0.10),
          mid    = stats::quantile(x$Ratio, 0.50),
          hi     = stats::quantile(x$Ratio, 0.90)
        )
      }
    )
  ))
  ribbon$Stock <- factor(ribbon$Stock, levels = StockNames)

  max_dev <- max(abs(c(ribbon$lo, ribbon$hi) - 1), na.rm = TRUE)
  y_pad   <- max(max_dev * 0.15, tol * 0.5)
  y_lim   <- 1 + c(-1, 1) * (max_dev + y_pad)

  p <- ggplot2::ggplot(ribbon, ggplot2::aes(x = Age)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2, fill = 'steelblue', colour = NA) +
    ggplot2::geom_line(ggplot2::aes(y = mid), colour = 'steelblue') +
    ggplot2::geom_hline(yintercept = 1, linetype = 'dashed', colour = 'grey30') +
    ggplot2::geom_hline(yintercept = c(1 - tol, 1 + tol),
                        linetype = 'dotted', colour = 'firebrick', alpha = 0.6) +
    ggplot2::facet_grid(Stock ~ Season) +
    ggplot2::coord_cartesian(ylim = y_lim) +
    ggplot2::labs(
      title    = paste0('OM / Reference N-at-age  (', nSim, ' simulation', if (nSim > 1) 's' else '', ')'),
      subtitle = paste0('Line = median  |  Ribbon = 10–90th percentile across sims  |  ',
                        'Dotted = ±', tol * 100, '% tolerance'),
      x        = 'Quarterly age class',
      y        = 'OM / Reference ratio'
    ) +
    ggplot2::theme_bw()

  print(p)
  invisible(df)
}
