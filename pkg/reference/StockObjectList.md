# Albacore operating model stock schedules

Fixed biological schedules for the IOTC albacore operating model, built
from `albMSE` conditioning data
(`https://github.com/iagomosqueira/albMSE/tree/main/boot/data/alb_abcdata.rda`).

## Usage

``` r
StockObjectList
```

## Format

A named list of length 2 (`"Female"`, `"Male"`), each element a
[MSEtool::stock](https://msetool.openmse.com/reference/Stock-class.html)
object with:

- `Ages`: seasonal age classes spanning `min(ages)` to `max(ages) * 4`,
  no plus group

- `Length`: mean length-at-age and CV-at-age by season

- `Weight`: mean weight-at-age by season

- `Maturity`: mean maturity-at-age by season

- `Fecundity`: maturity-at-age times weight-at-age, zero outside the
  spawning season (season 3)

## Source

`albMSE::boot/data/alb_abcdata.rda`

## Details

Contains one
[`MSEtool::Stock()`](https://msetool.openmse.com/reference/Stock.html)
object per sex, with seasonal age structure and length, weight,
maturity, and fecundity schedules. Natural mortality and the
stock-recruitrelationship are not set here; they are populated from MCMC
posterior samples.
