# MCMC-derived variables for OM6b (1% effort creep robustness OM)

Simulated population and observation quantities for each retained
MCMC/ABC posterior draw of the OM6b operating model. OM6b is a
robustness variant of the reference OM (OM5b) with an added 1% annual
catchability trend (`qtrend = TRUE`).

## Usage

``` r
CondData_OM6b
```

## Format

A list of length equal to the number of retained iterations. Each
element is itself a list with components including `N` (numbers-at-age),
`Rtot`, `SSB`, `dep`, `dbmsy`, `Bmsy`, `Cmsy`, `hmsy`, `hmsyrat`, `H`,
`Ihat`, `LFhat`, `B0`, `R0`, `M`, `h`, and `sela`.

## Source

`albMSE::data/om6b/mcvars_abc6b.rda`
