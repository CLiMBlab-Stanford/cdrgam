# Information criteria for block and sparse fits

[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) of block and sparse fits
correct the effective degrees of freedom for smoothing-parameter
uncertainty, as `mgcv` does for `gam` fits through `logLik.gam()` (Wood,
Pya, and Saefken 2016). Uncorrected, the criterion favours overly
complex models, especially ones with random effects. The correction
needs dense matrices of the size of the model, one per smoothing
parameter, so a fit with more coefficients than
`getOption("cdrgam.edf2_max_coefficients", 3000)` gives an error rather
than an uncorrected value. It needs a REML or ML fit.

## Usage

``` r
# S3 method for class 'cdrgam_block'
AIC(object, ..., k = 2)

# S3 method for class 'cdrgam_sparse'
AIC(object, ..., k = 2)

# S3 method for class 'cdrgam_block'
BIC(object, ...)

# S3 method for class 'cdrgam_sparse'
BIC(object, ...)
```

## Arguments

- object:

  A block or sparse `cdrgam` fit.

- ...:

  Further fits to compare.

- k:

  The penalty per effective degree of freedom; 2 gives AIC.

## Value

A number for one fit; for several, a data frame with `df` and the
criterion, one row per fit.

## Details

Criteria are comparable only between models fitted to the same
responses, with the same response variable on the same scale.

## References

Wood, S. N., Pya, N., and Saefken, B. (2016). Smoothing parameter and
model selection for general smooth models. *Journal of the American
Statistical Association*, 111(516), 1548–1563.
