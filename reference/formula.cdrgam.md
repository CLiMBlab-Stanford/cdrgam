# Extract a CDR-GAM formula

Returns the raw user formula, the normalized CDR formula with structural
defaults, the effective formula after automatic simplification, or the
native mgcv export formula.

## Usage

``` r
# S3 method for class 'cdrgam_design'
formula(x, type = c("raw", "normalized", "effective", "mgcv"), ...)

# S3 method for class 'cdrgam'
formula(x, type = c("raw", "normalized", "effective", "mgcv"), ...)
```

## Arguments

- x:

  A fitted `cdrgam` object or prepared `cdrgam_design`.

- type:

  Formula representation: the raw formula supplied by the user, the
  normalized CDR formula including structural defaults, or the effective
  CDR formula after automatic simplification, or `"mgcv"`, the native
  formula exported by
  [`as_mgcv_spec()`](https://climblab.org/cdrgam/reference/as_mgcv_spec.md).
  Its named matrix terms receive their penalties through the rest of
  that specification.

- ...:

  Unused.

## Value

A formula, or a parameter-named formula list for a distributional model.
