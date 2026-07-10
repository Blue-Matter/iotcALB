# Raw biological schedules from the albMSE conditioning data

Inputs used to build
[StockObjectList](https://iotcalb.bluematterscience.com/pkg/reference/StockObjectList.md)'s
length, weight, and maturity schedules.

## Usage

``` r
albMSE_Biology
```

## Format

A named list:

- `ages` — Numeric vector of annual age classes.

- `mula`, `sdla` — Arrays `[annual_age, season, sex]`. Mean and SD of
  length-at-age.

- `mulbins` — Numeric vector of length bin midpoints.

- `wta` — Array `[annual_age, season, sex]`. Mean weight-at-age (kg).

- `mata` — Array `[annual_age, season, sex]`. Maturity-at-age.

## Source

`albMSE::boot/data/alb_abcdata.rda`
