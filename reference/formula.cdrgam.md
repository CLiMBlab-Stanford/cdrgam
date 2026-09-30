# Extract a CDR-GAM formula

Returns the raw user formula, the normalized CDR formula with structural
defaults, the effective formula after automatic simplification, or the
translated mgcv formula.

## Usage

``` r
# S3 method for class 'cdrgam_design'
formula(
  x,
  type = c("raw", "normalized", "effective"),
  ...
)

# S3 method for class 'cdrgam'
formula(
  x,
  type = c("raw", "normalized", "effective", "mgcv"),
  ...
)
```

## Arguments

- x:

  A fitted `cdrgam` object or prepared `cdrgam_design`.

- type:

  Formula representation to return. The `"mgcv"` translation is
  available for fitted models.

- ...:

  Unused.

## Value

The selected formula, or a parameter-named formula list for a
distributional model.
