# Import an MSEtool OM from IOTC Albacore conditioning output

Constructs a multi-stock, multi-fleet MSEtool
[`MSEtool::OM()`](https://msetool.openmse.com/reference/OM.html) from
the MCMC posterior samples produced by the IOTC Albacore conditioning
model. The OM is quarterly (4 seasons) with two sexes (Female, Male) and
six fleets.

## Usage

``` r
ImportOM(object = Cond_BaseCase, Name = "Base Case")
```

## Arguments

- object:

  A named list of conditioning output, as produced by the IOTC Albacore
  conditioning model. Defaults to
  [Cond_BaseCase](https://iotcalb.bluematterscience.com/pkg/reference/Cond_BaseCase.md).
  Must contain:

  - `mcmcvars`: list of MCMC iterations, each with arrays `N`, `H`,
    `sela`

  - `steepness`, `natural_mortality`, `R0`, `sigmaR`: numeric vectors of
    length `nSim`

  - `maturity`, `weight`: arrays `[annual_age, season, sex]`

- Name:

  Character scalar. Name assigned to the returned OM object.

## Value

An [`MSEtool::OM()`](https://msetool.openmse.com/reference/OM.html)
object with stocks and fleets fully populated from the MCMC posterior,
ready to pass to
[`MSEtool::Simulate()`](https://msetool.openmse.com/reference/Simulate.html).

## Details

Biological schedules (weight, maturity, natural mortality) and
stock-recruit parameters are drawn directly from the conditioning
object. Historical fishing mortality is recovered from the MCMC
numbers-at-age via N-at-age transitions, which exactly reproduces the
reference model's survival without requiring knowledge of the internal H
parameterisation.
