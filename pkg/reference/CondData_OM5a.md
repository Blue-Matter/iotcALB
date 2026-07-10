# MCMC-derived variables for OM5a (SW CPUE reference OM)

Simulated population and observation quantities for each retained
MCMC/ABC posterior draw of the OM5a operating model OM5a conditions on
the south-west (SW) CPUE index (`fcpue = 3`) rather than the base case
north-west (NW) index used in OM5b.

## Usage

``` r
CondData_OM5a
```

## Format

A list of length equal to the number of retained iterations. Each
element is itself a list with components including `N` (numbers-at-age),
`Rtot`, `SSB`, `dep`, `dbmsy`, `Bmsy`, `Cmsy`, `hmsy`, `hmsyrat`, `H`,
`Ihat`, `LFhat`, `B0`, `R0`, `M`, `h`, and `sela`.

## Source

`albMSE::data/om5a/mcvars_abc5a.rda`
