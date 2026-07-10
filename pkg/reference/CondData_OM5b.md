# MCMC-derived variables for OM5b (base case reference OM)

Simulated population and observation quantities for each retained
MCMC/ABC posterior draw of the OM5b operating model. OM5b is the
reference operating model: conditioned on the north-west (NW) LL CPUE
index, with SSB priors, recruitment variability, and the overfishing
penalty.

## Usage

``` r
CondData_OM5b
```

## Format

A list of length equal to the number of retained iterations. Each
element is itself a list with components including `N` (numbers-at-age),
`Rtot`, `SSB`, `dep`, `dbmsy`, `Bmsy`, `Cmsy`, `hmsy`, `hmsyrat`, `H`,
`Ihat`, `LFhat`, `B0`, `R0`, `M`, `h`, and `sela`.

## Source

`albMSE::data/om5b/mcvars_abc5b.rda`
