# Update a simulated Hist object with new catch and CPUE data

Replaces the historical landings and CPUE index data in a `Hist` object
(from
[`MSEtool::Simulate()`](https://msetool.openmse.com/reference/Simulate.html))
with the updated observations returned by
[`ProcessNewData()`](https://iotcalb.bluematterscience.com/pkg/reference/ProcessNewData.md).

## Usage

``` r
UpdateData(Hist, plot = FALSE)
```

## Arguments

- Hist:

  A `Hist` object returned by
  [`MSEtool::Simulate()`](https://msetool.openmse.com/reference/Simulate.html).

- plot:

  Logical. If `TRUE`, plots comparing the conditioned and updated data
  are saved to `figures/data/`:

  `catch_seasonal_fleet.png`

  :   Seasonal catch by fleet.

  `catch_annual_fleet.png`

  :   Annual catch by fleet.

  `catch_annual_overall.png`

  :   Annual total catch.

  `index_seasonal.png`

  :   Seasonal CPUE index by fleet, standardized to the mean over the
      conditioned years.

  `index_annual_mean.png`

  :   Annual mean of the standardized CPUE index by fleet.

## Value

The `Hist` object with updated landings and CPUE data.

## Details

This must be called *after* the observation model (`Obs`) has been
populated in
[`MSEtool::Simulate()`](https://msetool.openmse.com/reference/Simulate.html)
using the catch and index data that were used in conditioning. Otherwise
the differences between the conditioned and updated data would be
interpreted as observation error in the observed catches.

The data are replaced in both the OM data (`Hist@OM@Data`) and the
simulated data (`Hist@Data`) for the combined stock.

If the fleets of the OM have been combined with
[`MSEtool::CombineFleets()`](https://msetool.openmse.com/reference/CombineFleets.html),
the new landings are summed over fleets, and the CPUE indices, which
[`MSEtool::CombineFleets()`](https://msetool.openmse.com/reference/CombineFleets.html)
moves to the `Survey` slot, are replaced there.
