# Construct a supported CDR-GAM response family

`cdrgam_family()` is the validated family registry used by the
command-line harness. Family objects may still be passed directly to
[`cdrgam()`](https://climblab.org/cdrgam/reference/cdrgam.md) and
[`cdrgam.fit()`](https://climblab.org/cdrgam/reference/cdrgam.md). The
native `mgcv` backend supports every combination returned here. `gaulss`
selects mgcv's two-predictor Gaussian location–scale family. The dense
block backend supports a joint location–scale reference solver plus
`gaussian(identity)`, `binomial(logit)`, `poisson(log)`, and
estimated-dispersion `Gamma(log)`. The sparse backend supports those
single-predictor combinations with streamed PIRLS and an exact Laplace
score, plus joint location–scale LAML with an exact outer score.

## Usage

``` r
cdrgam_family(
  family = c("gaussian", "binomial", "poisson", "Gamma", "gaulss"),
  link = NULL
)
```

## Arguments

- family:

  A supported family name or an existing family object. The default is
  `"gaussian"`. A family object is returned unchanged.

- link:

  Optional link name. It may be supplied only when `family` is a name
  other than `"gaulss"`, whose links are fixed (identity for the
  location and log for the scale).

## Value

A standard R family object.
