# Buffer Harvest Control Rule Management Procedure

A model-free management procedure that replicates the CPUE + buffer HCR
MP of the earlier albacore MSE work (`cpue.ind` and `buffer.hcr` of the
FLR `mse`/`msemodules` packages; Mosqueira and Hillary 2025,
IOTC-2025-WPM16-11, IOTC-2026-WPM17(MSE)-03). Unlike
[`MSEtool::IndexTarget()`](https://msetool.openmse.com/reference/IndexMPs.html)
with a buffer HCR
([`BufferHCR()`](https://iotcalb.bluematterscience.com/pkg/reference/BufferHCR.md)),
the HCR output is a multiplier on the previous TAC, so the TAC keeps
changing in each management cycle while the index is outside the buffer.

## Usage

``` r
BufferMP(
  Data,
  Indices = "LL1",
  IndexSource = "Survey",
  IndexRef,
  IndexWeight = NULL,
  nYears = 4,
  Metric = c("wmean", "mean"),
  Limit = 0.4,
  LowerBuffer = 0.56,
  UpperBuffer = 0.96,
  SlopeRatio = 0.1,
  DeltaDown = c(0, 0.15),
  DeltaUp = c(0, 0.15),
  TACType = c("Removals", "Landings"),
  tunepar = 1
)
```

## Arguments

- Data:

  A `data` object.

- Indices:

  Character vector. Names of the indices in `IndexSource`. Default
  `'LL1'`.

- IndexSource:

  Character. `'Survey'` (default) or `'CPUE'`.

- IndexRef:

  Numeric vector, one per index. The index level that corresponds to a
  status of 1 (e.g. the index at `SBMSY`).

- IndexWeight:

  `NULL` (default, equal weights) or a numeric vector of weights, one
  per index.

- nYears:

  Positive integer. Number of data years in the metric. Default `4`.

- Metric:

  Character. `'wmean'` (default) or `'mean'`.

- Limit, LowerBuffer, UpperBuffer:

  Numeric. The limit and the lower and upper buffers, in units of status
  (index / `IndexRef`). Defaults `0.4`, `0.56`, and `0.96`: the values
  of the earlier work (0.25, 0.35, and 0.60 in units of the index)
  scaled so that the limit is at 0.4 `SBMSY`.

- SlopeRatio:

  Numeric. Slope of the increase above the upper buffer relative to the
  slope of the decrease between the limit and the lower buffer. Default
  `0.10`.

- DeltaDown, DeltaUp:

  Numeric vectors, length 2. Minimum and maximum proportional decrease
  and increase in the TAC (see
  [`MSEtool::ConstrainTAC()`](https://msetool.openmse.com/reference/ConstrainTAC.html)).
  Default `c(0, 0.15)`.

- TACType:

  Character. `'Removals'` (default) or `'Landings'`.

- tunepar:

  Positive number. Tuning parameter; see Tuning.

## Value

An `advice` object.

## Details

### Procedure

1.  The data are aggregated to calendar years with
    [`MSEtool::AnnualData()`](https://msetool.openmse.com/reference/AnnualData.html)
    (each index averaged over the seasons with data).

2.  The index metric is the mean of each index over the last `nYears`
    data years: a weighted mean (`Metric = 'wmean'`), with weight 0.5
    for the last year and the other 0.5 shared among the earlier years
    in proportion to their position (e.g. 1/12, 2/12, 3/12 for
    `nYears = 4`), or an unweighted mean (`Metric = 'mean'`).

3.  Stock status is the metric relative to `IndexRef` (e.g. the index
    level at `SBMSY`,
    [`IndexAtSBMSY()`](https://iotcalb.bluematterscience.com/pkg/reference/IndexAtSBMSY.md)),
    combined over indices with `IndexWeight`.

4.  The TAC multiplier is, for status `x`, limit `L`, lower buffer `B1`
    and upper buffer `B2`:

    - `(x / L)^2 / 2` for `x <= L`;

    - `0.5 (1 + (x - L) / (B1 - L))` for `L < x <= B1`;

    - `1` for `B1 < x < B2`;

    - `1 + SlopeRatio (0.5 / (B1 - L)) (x - B2)` for `x >= B2`.

5.  The TAC is the previous TAC
    ([`MSEtool::LastTAC()`](https://msetool.openmse.com/reference/DataHelpers.html))
    times the multiplier, with the change constrained by `DeltaDown` and
    `DeltaUp`
    ([`MSEtool::ConstrainTAC()`](https://msetool.openmse.com/reference/ConstrainTAC.html)).

### Tuning

In the earlier work the MP was tuned by the position of the upper
buffer. Here `tunepar` sets the upper buffer to
`LowerBuffer + (UpperBuffer - LowerBuffer) / tunepar`, so that larger
values of `tunepar` move the upper buffer towards the lower buffer and
increase the catch. `tunepar = 1` uses `UpperBuffer`.
