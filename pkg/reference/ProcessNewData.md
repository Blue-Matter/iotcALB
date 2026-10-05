# Process new catch and CPUE data from an SS3 data file

Reads an SS3 data file with
[`r4ss::SS_readdat()`](https://rdrr.io/pkg/r4ss/man/SS_readdat.html) and
aggregates the quarterly catch and CPUE observations from 2000 onwards
into the fleet structure used by the operating models.

## Usage

``` r
ProcessNewData(file = "data-raw/new_data.dat")
```

## Arguments

- file:

  Path to the SS3 data file. Defaults to `'data-raw/new_data.dat'`,
  relative to the package source directory.

## Value

A named list with two matrices, each with dimensions `c(Year, Fleet)`:

- Landings:

  Total catch by season and fleet (`LL1`, `LL2`, `LL3`, `LL4`, `PS`,
  `Other`).

- CPUE:

  CPUE index by season and longline fleet (`LL1`-`LL4`).

## Details

Catches are summed across the SS3 fleets within each OM fleet: `LL1`
(fleets 1-4), `LL2` (5-8), `LL3` (9-12), `LL4` (13-16), `PS` (19), and
`Other` (17-18, 20-23).

CPUE indices are taken from SS3 index fleets 24-27 (`LL1`), 28-31
(`LL2`), 32-35 (`LL3`), and 36-39 (`LL4`). Seasons with no observation
are `NA`.

Time steps are labelled with decimal years from
[`MSEtool::CalcYears()`](https://msetool.openmse.com/reference/Years.html)
(4 seasons per year, current year 2023).
