# Summarize a fitted CDR-GAM

The default output follows
[`mgcv::summary.gam()`](https://rdrr.io/pkg/mgcv/man/summary.gam.html):
it reports parametric coefficients and one row per smooth term instead
of expanded smooth and random-effect coefficients. Fully penalized
random-effect and grouped IRF terms report effective degrees of freedom
and an explicit "not computed" test status. The printed output reports
automatic formula simplifications directly; normalized and effective
formulas remain available through
[`formula()`](https://rdrr.io/r/stats/formula.html). This status does
not indicate fit failure. Distributional Gaussian summaries report
deviance, null deviance, and deviance explained without an adjusted
R-squared statistic. Set `all.coefficients = TRUE` to request expanded
coefficients.

## Usage

``` r
# S3 method for class 'cdrgam'
summary(
  object,
  dispersion = NULL,
  freq = FALSE,
  re.test = TRUE,
  all.coefficients = FALSE,
  ...
)

# S3 method for class 'cdrgam_block'
summary(
  object,
  dispersion = NULL,
  freq = FALSE,
  re.test = TRUE,
  all.coefficients = FALSE,
  ...
)

# S3 method for class 'cdrgam_sparse'
summary(
  object,
  dispersion = NULL,
  freq = FALSE,
  re.test = TRUE,
  all.coefficients = FALSE,
  ...
)

# S3 method for class 'summary.cdrgam'
print(x, ...)
```

## Arguments

- object:

  A fitted `cdrgam` model.

- dispersion:

  Optional known dispersion.

- freq, re.test:

  Controls passed to
  [`mgcv::summary.gam()`](https://rdrr.io/pkg/mgcv/man/summary.gam.html)
  for a native mgcv fit. Custom backends retain these arguments for
  interface compatibility.

- all.coefficients:

  Include the expanded smooth and random-effect coefficient table. The
  default is false.

- ...:

  Additional arguments passed to native
  [`mgcv::summary.gam()`](https://rdrr.io/pkg/mgcv/man/summary.gam.html)
  or the print method.

- x:

  A CDR-GAM summary object.

## Value

A `summary.gam`-style object containing parametric and smooth-term
tables, the raw and translated mgcv formulas, and recorded automatic
simplifications. Native fits retain the complete `summary.gam` result.
Custom backends provide its commonly used coefficient, smooth-term,
fit-statistic, and formula fields. The `s.test` field states whether
each smooth-term test was computed.
