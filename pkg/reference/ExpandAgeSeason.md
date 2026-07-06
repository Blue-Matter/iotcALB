# Expand an annual-age-by-season schedule onto quarterly age classes

The conditioning model tabulates life-history schedules (weight,
maturity, etc.) by annual age (rows, starting at annual age 0) and
season (columns). An MSEtool `Stock` with a quarterly age structure
instead indexes age by quarterly age classes (`0, 0.25, 0.5, ...`). This
expands the annual-age rows onto the quarterly age classes, holding the
value constant across the sub-year age classes within each annual age,
while keeping the season dimension intact.

## Usage

``` r
ExpandAgeSeason(x, AgeClasses)
```

## Arguments

- x:

  A matrix with dimensions `[annual_age, season]`, or an array with
  additional trailing dimensions (e.g. `[annual_age, season, stock]`),
  indexed from annual age 0 in row 1.

- AgeClasses:

  Numeric vector of quarterly age classes to expand onto, e.g.
  `Classes(stock)`.

## Value

An array with the first dimension replaced by `length(AgeClasses)`,
keeping the season dimension and any trailing dimensions of `x`.
