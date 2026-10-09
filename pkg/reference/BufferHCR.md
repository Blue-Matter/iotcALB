# Harvest Control Rule Control Points of a Buffer HCR

Converts the buffer harvest control rule used in the earlier albacore
MSE work (Mosqueira and Hillary 2025, IOTC-2025-WPM16-11) to the control
points of
[`MSEtool::IndexTarget()`](https://msetool.openmse.com/reference/IndexMPs.html).
The TAC is the reference catch while the index is between the lower and
upper buffers; it decreases below the lower buffer, more rapidly below
the limit, and increases above the upper buffer at `SlopeRatio` times
the rate of decrease between the limit and the lower buffer, up to
`MaxIndex`.

## Usage

``` r
BufferHCR(
  Buffer,
  Limit = 0.4,
  LimitRate = 0.5,
  Zero = 0.2,
  SlopeRatio = 0.25,
  MaxIndex = 2
)
```

## Arguments

- Buffer:

  Numeric, length 2. Lower and upper buffers, in units of
  `status / IndexTarget`.

- Limit:

  Numeric. Limit, in units of `status / IndexTarget`. Default `0.4`.

- LimitRate:

  Numeric. TAC multiplier at the limit. Default `0.5`.

- Zero:

  Numeric. Index level at which the TAC is zero. Default `0.2`.

- SlopeRatio:

  Numeric. Slope of the increase above the upper buffer relative to the
  slope of the decrease below the lower buffer. Default `0.25`.

- MaxIndex:

  Numeric. Index level above which the TAC no longer increases. Default
  `2`.

## Value

A list with `HCRControlPointsIndex` and `HCRControlPointsRate`.
