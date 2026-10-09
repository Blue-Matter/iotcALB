# Surplus Production Model Settings from the Operating Model

Derives the fixed parameters and informative priors of the surplus
production model used by
[`MSEtool::SurplusProduction()`](https://msetool.openmse.com/reference/SurplusProduction.html)
from the operating model, so that the assessment model in the MP matches
the productivity and initial depletion of the base case OM.

## Usage

``` r
SPSettings(Hist, MinCV = c(FMSY = 0.3, MSY = 0.4))
```

## Arguments

- Hist:

  A `hist` object (the base case OM).

- MinCV:

  Named numeric vector. Minimum CV of the priors on `FMSY` and `MSY`.
  Default `c(FMSY = 0.3, MSY = 0.4)`.

## Value

A list with `Priors` (a list of `c(median, CV)` for `FMSY` and `MSY`,
see
[`MSEtool::FitSP()`](https://msetool.openmse.com/reference/FitSP.html)),
`Shape`, `Depletion`, and `OM` (a `data.frame` of the per-simulation OM
values).

## Details

The surplus production model has a single biomass pool, so the OM
quantities are calculated for the total biomass of both sexes:

- `MSY`: equilibrium removals at MSY (`MSYLandings + MSYDiscards`).

- `FMSY`: the harvest rate at MSY, `MSY / BMSY`, where `BMSY` is the
  total equilibrium biomass at MSY. This is the `FMSY` of the surplus
  production model, not the apical fishing mortality of the OM.

- `Shape`: the Pella-Tomlinson shape parameter `n` with
  `n^(1/(1-n)) = BMSY/B0`, the median over simulations.

- `Depletion`: total biomass relative to unfished at the start of the
  first historical year (`B/B0`), the median over simulations.

The priors are lognormal, with median equal to the median over
simulations and CV equal to the CV of the OM values over simulations or
`MinCV`, whichever is larger.
