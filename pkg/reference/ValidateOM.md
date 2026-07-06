# Compare OM N-at-age against conditioning model reference

Maps each MCMC iteration's N-at-age into the OM's quarterly age ×
timestep layout and compares it against the simulated N-at-age from a
[`MSEtool::Hist()`](https://msetool.openmse.com/reference/Hist.html)
object. Diagnostics are computed across all simulations; the plot
summarises the distribution of OM/Reference ratios as a ribbon.

## Usage

``` r
ValidateOM(hist, object = Cond_BaseCase, tol = 0.02, min_pct = 0.02)
```

## Arguments

- hist:

  A [`MSEtool::Hist()`](https://msetool.openmse.com/reference/Hist.html)
  object produced by
  [`MSEtool::Simulate()`](https://msetool.openmse.com/reference/Simulate.html).

- object:

  A named list of conditioning output (same object passed to
  [`ImportOM()`](https://iotcalb.bluematterscience.com/pkg/reference/ImportOM.md)).
  Defaults to
  [Cond_BaseCase](https://iotcalb.bluematterscience.com/pkg/reference/Cond_BaseCase.md).

- tol:

  Numeric. Ratio threshold outside `[1-tol, 1+tol]` that triggers a
  printed warning. Defaults to `0.02` (2%).

- min_pct:

  Numeric. Minimum fraction of cells (0–1) that must exceed `tol` before
  the worst-offenders table is printed. Defaults to `0.02` (2%).

## Value

Invisibly returns a data frame of OM vs reference N-at-age with columns
`Stock`, `Sim`, `Age`, `Timestep`, `Season`, `Year`, `OM`, `Ref`,
`Ratio`.
