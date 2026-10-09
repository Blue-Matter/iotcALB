# Index Level Corresponding to SBMSY

Calibrates relative abundance indices to stock status, so that the
harvest control rules of the index-based MPs
([`MSEtool::IndexRate()`](https://msetool.openmse.com/reference/IndexMPs.html),
[`MSEtool::IndexTarget()`](https://msetool.openmse.com/reference/IndexMPs.html),
[`BufferMP()`](https://iotcalb.bluematterscience.com/pkg/reference/BufferMP.md))
can be expressed in units of `SB/SBMSY`. The index level at `SBMSY` is
the mean annual index over `Years` divided by the median (over
simulations) of the annual `SB/SBMSY` averaged over the same years (see
[`AnnualStatus()`](https://iotcalb.bluematterscience.com/pkg/reference/AnnualStatus.md)).

## Usage

``` r
IndexAtSBMSY(Hist, Index = "LL1", Years = 2016:2020)
```

## Arguments

- Hist:

  A `hist` object (the base case OM).

- Index:

  Character vector. Names of the indices (in the `Survey` slot after
  [`MSEtool::CombineFleets()`](https://msetool.openmse.com/reference/CombineFleets.html),
  else `CPUE`). Default `'LL1'`.

- Years:

  Numeric. Calendar years used for the calibration. Must be historical
  years of `Hist`. Default `2016:2020`.

## Value

Named numeric vector (one value per index): the index level at `SBMSY`,
with attributes `Index` (the mean index over `Years`) and `SB_SBMSY`
(the median status over `Years`).
