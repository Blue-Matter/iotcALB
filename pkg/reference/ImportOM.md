# Import an MSEtool OM from IOTC Albacore conditioning output

Constructs a multi-stock, multi-fleet MSEtool
[`MSEtool::OM()`](https://msetool.openmse.com/reference/OM.html) from
the MCMC posterior samples produced by the IOTC Albacore conditioning
model. The OM is quarterly (4 seasons) with two sexes (Female, Male;
from
[StockObjectList](https://iotcalb.bluematterscience.com/pkg/reference/StockObjectList.md))
and six fleets.

## Usage

``` r
ImportOM(OMName = "OM5b")
```

## Arguments

- OMName:

  Character scalar. Suffix identifying which conditioning dataset to
  load, via `get(paste0('CondData_', OMName))` (e.g. `'OM5b'` loads
  [CondData_OM5b](https://iotcalb.bluematterscience.com/pkg/reference/CondData_OM5b.md)).
  Also used as the `Name` of the returned OM.

## Value

An [`MSEtool::OM()`](https://msetool.openmse.com/reference/OM.html)
object with stocks and fleets fully populated from the MCMC posterior,
ready to pass to
[`MSEtool::Simulate()`](https://msetool.openmse.com/reference/Simulate.html).

## Details

Natural mortality and stock-recruit parameters are drawn directly from
`CondData_<OMName>` (e.g.
[CondData_OM5b](https://iotcalb.bluematterscience.com/pkg/reference/CondData_OM5b.md)).
`object$R0` is total (both-sex) unfished recruitment and is split evenly
between the two stocks; recruitment is restricted to the recruitment
season (season 4) via a mask on `R0`, and `RecDevInit`/`RecDevHist` are
derived from the conditioning model's own numbers-at-age and recruitment
deviations. `R0` is passed to
[`MSEtool::SRR()`](https://msetool.openmse.com/reference/SRR.html) in
thousands (`Units = 1000`), so
[`MSEtool::Simulate()`](https://msetool.openmse.com/reference/Simulate.html)
output (e.g. `Hist@Number`) is in thousands of fish, not the raw
individual counts used by `CondData_<OMName>`.

Historical fishing mortality is recovered from the MCMC harvest rates
(`object$H`), and fleet selectivity from `object$sela`.
