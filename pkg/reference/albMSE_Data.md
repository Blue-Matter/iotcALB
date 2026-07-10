# MSEtool Data object for IOTC albacore

Observed catch and CPUE data for Indian Ocean albacore tuna (*Thunnus
alalunga*), reformatted from the `albMSE` ABC conditioning data
(`boot/data/alb_abcdata.rda`) into an
[MSEtool::data](https://msetool.openmse.com/reference/data-class.html)
object for use with candidate management procedures (CMPs). Covers
2000-2020 at quarterly (4-season) resolution across 6 fleets.

## Usage

``` r
albMSE_Data
```

## Format

A [MSEtool::data](https://msetool.openmse.com/reference/data-class.html)
object with:

- `Landings`: a
  [MSEtool::catchdata](https://msetool.openmse.com/reference/catchdata-class.html)
  object, quarterly catch (t) by fleet (`LL1`-`LL4`, `PS`, `Other`),
  `nTS = 84` timesteps

- `CPUE`: a
  [MSEtool::indicesdata](https://msetool.openmse.com/reference/indicesdata-class.html)
  object, quarterly CPUE by fleet, standardized to a mean of 1 per fleet

## Source

`albMSE::boot/data/alb_abcdata.rda` (`C`, `I` objects)
