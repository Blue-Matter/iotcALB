# Annual Stock Status Relative to MSY Reference Points

Annual spawning biomass and fishing mortality relative to the MSY
reference points, from
[`MSEtool::SB_SBMSY()`](https://msetool.openmse.com/reference/relative_ref.html)
and
[`MSEtool::F_FMSY()`](https://msetool.openmse.com/reference/relative_ref.html).

## Usage

``` r
AnnualStatus(object, Stock = "Female")
```

## Arguments

- object:

  An `mse` or `hist` object.

- Stock:

  Character. Name of the spawning stock. Default `'Female'`.

## Value

A `data.frame` with columns `Sim`, `Year` (calendar year), `Period`,
`MP`, `SB_SBMSY`, and `F_FMSY`. Historical rows have
`MP = 'Historical'`.

## Details

The OMs are seasonal (4 seasons). The MSY reference points are annual,
and MSEtool reports the ratios on the same basis, one value per calendar
year: `F_FMSY` is the apical fishing mortality summed over the seasons
of the year, and `SB_SBMSY` the spawning biomass in the reference season
(see
[`MSEtool::RefSeason()`](https://msetool.openmse.com/reference/RefSeason.html)).
This requires a version of MSEtool that reports annual ratios for
seasonal OMs; an error is given if the ratios are per time step.
