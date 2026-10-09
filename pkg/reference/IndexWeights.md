# Inverse-Variance Weights of Indices from the Conditioned Observation Error

Weights for combining several indices in the index-based MPs,
proportional to `1 / SD^2`, where `SD` is the median (over simulations)
standard deviation of the log-scale observation error of each index
conditioned in the OM.

## Usage

``` r
IndexWeights(Hist, Index)
```

## Arguments

- Hist:

  A `hist` object (the base case OM).

- Index:

  Character vector. Names of the indices.

## Value

Named numeric vector of weights, summing to 1, with attribute `SD`.
