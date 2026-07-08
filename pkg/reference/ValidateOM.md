# Compare OM reference quantities against the conditioning model

Compares an
[`MSEtool::Hist()`](https://msetool.openmse.com/reference/Hist.html)
simulation against the reference quantities recovered from the MCMC
posterior of the conditioning model: total numbers, spawning biomass,
catch-by-fleet, and numbers-at-age.

## Usage

``` r
ValidateOM(
  Hist,
  object = Cond_BaseCase,
  tol = 0.05,
  min_pct = 0.02,
  alpha = 0.2,
  probs = c(0.1, 0.9),
  save_plots = TRUE,
  outdir = "figures/diagnostics/OM",
  width = NULL,
  height = NULL,
  verbose = TRUE
)
```

## Arguments

- Hist:

  A [`MSEtool::Hist()`](https://msetool.openmse.com/reference/Hist.html)
  object produced by
  [`MSEtool::Simulate()`](https://msetool.openmse.com/reference/Simulate.html)

- object:

  A named list of conditioning output (same object passed to
  [`ImportOM()`](https://iotcalb.bluematterscience.com/pkg/reference/ImportOM.md)).
  Defaults to
  [Cond_BaseCase](https://iotcalb.bluematterscience.com/pkg/reference/Cond_BaseCase.md).
  Must contain `mcmcvars` (a list of MCMC iterations, each with arrays
  `N`, `H`, `sela`, `SSB`) and `weight`.

- tol:

  Numeric. Ratio threshold outside `[1-tol, 1+tol]` that triggers a
  printed warning for the N-at-age comparison. Defaults to `0.02` (2%).

- min_pct:

  Numeric. Minimum fraction of N-at-age cells (0-1) that must exceed
  `tol` before the worst-offenders table is printed. Defaults to `0.05`
  (5%).

- alpha:

  Numeric. Transparency used for ribbon fills and per-simulation ratio
  lines. Defaults to `0.2`.

- probs:

  Numeric vector of length 2 giving the lower/upper quantiles used for
  ribbon plots. Defaults to `c(0.1, 0.9)`.

- save_plots:

  Logical. If `TRUE` (default), diagnostic plots are written as PNG
  files to `file.path(outdir, Hist@OM@Name)`.

- outdir:

  Character. Base directory for saved plots. Defaults to
  `"figures/diagnostics/OM"`.

- width, height:

  Numeric. Width/height (inches) passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).
  If `NULL` (default), each plot is sized automatically from its number
  of facet panels: single-panel plots (Number, SSB) use 6 x 4; the
  Catch-by-fleet and N-at-age plots, which facet by fleet/season, scale
  up accordingly. Set either to a number to use that fixed size for
  every saved plot instead.

- verbose:

  Logical. If `TRUE` (default), print N-at-age match diagnostics
  (success/warning messages and the worst-offenders table).

## Value

Invisibly returns a named list:

- `n`, `ssb`, `catch`: each a list with `ref` and `om` data frames

- `natage`: long data frame of OM vs reference N-at-age with columns
  `Stock`, `Sim`, `Age`, `Timestep`, `Season`, `Year`, `OM`, `Ref`,
  `Ratio`

## Details

For total numbers, spawning biomass, and catch, a ribbon plot (median
with a 10th-90th percentile band, OM vs. reference) and a per-simulation
ratio plot (OM / reference) are produced.

Numbers-at-age is compared by mapping each MCMC iteration's N-at-age
into the OM's quarterly age x timestep layout; the distribution of
OM/reference ratios is summarised as a ribbon by quarterly age class and
season, and cell-level diagnostics are printed.

All plots are optionally saved to disk.
