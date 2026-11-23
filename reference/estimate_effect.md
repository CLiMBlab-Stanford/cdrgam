# Evaluate fitted effects on declarative axis grids

Computes selected impulse-response effects and their pointwise
uncertainty as joint linear functionals of the fitted coefficients. Axes
may vary over numeric values, quantiles, or the fitted domain, or be
fixed at a numeric value or a training-data mean or median with an
optional standard deviation offset.

## Usage

``` r
estimate_effect(
  object,
  terms = NULL,
  axes = list(),
  impulses = NULL,
  responses = NULL,
  composition = c("term", "total"),
  grouping = c("population", "deviation", "conditional"),
  groups = NULL,
  reference = list(),
  se = TRUE,
  unconditional = FALSE,
  level = 0.95,
  n = 200L,
  boundary = c("estimate", "bound"),
  points = NULL
)
```

## Arguments

- object:

  A fitted `cdrgam` model.

- terms:

  Term labels, or a selector list with optional fields `names` (term
  labels), `predictors` (predictor names, or `"*"` for any term with
  predictors), `match` (`"contains"`, the default, or `"exact"`), and
  `grouped` (logical). `NULL` selects every impulse-response term.

- axes:

  Named axis requests. Predictor requests are nested under `predictors`;
  `"*"` supplies a fallback for every predictor axis.

- impulses, responses:

  Optional training streams used to recover axis summaries from models
  fitted before summaries were stored.

- composition:

  Return the selected term alone or add its hierarchy of fitted marginal
  terms.

- grouping:

  Return population, grouped-deviation, or conditional effects.

- groups:

  Optional grouping levels.

- reference:

  Optional named list of reference values for evaluation axes, using the
  same request forms as `axes`. Each request must resolve to one value.
  A predictor may instead be `NA` or `"default"` to leave it at the
  model default. Names are axis names such as `lag` or a predictor.

- se:

  Include pointwise standard errors.

- unconditional:

  Include smoothing-parameter uncertainty when supported. Smoothing
  parameters classified as on the boundary of their range contribute no
  uncertainty, because the criterion has no local curvature along them.
  Sparse fits stop with an error when the outer Hessian is not positive
  definite after those parameters are excluded.

- level:

  Confidence level for `lower` and `upper`.

- n:

  Default number of points for fitted-domain grids.

- boundary:

  `"estimate"` (the default) uses the fitted smoothing parameters.
  `"bound"` repeats the calculation with each boundary smoothing
  parameter set to its profile bound, the largest variance component the
  data still support, while the other smoothing parameters stay at their
  estimates. Compare the two to check whether conclusions depend on
  treating those variance components as zero. `"bound"` requires a
  sparse Gaussian REML fit with
  `sparse_control$boundary_sensitivity = TRUE` and
  `unconditional = FALSE`.

- points:

  Optional data frame of evaluation points, one per row, with columns
  named by axis, such as `lag` or a predictor. Its columns are paired
  rather than crossed, so the result has one row per point. Each axis
  that is not a column takes its value from `axes`, and that request
  must resolve to one value. Use it to evaluate along a path through the
  axes, such as each predictor value at its own lag.

## Value

A `cdrgam_effect_grid` data frame with evaluation axes, estimates,
standard errors, and pointwise confidence limits on the linear-predictor
scale. Each reference axis adds a `<axis>_reference` column, which is
`NA` for a model default.

## Details

With `composition = "total"`, the selected term is summed with every
ungrouped term in its hierarchy: terms on the same `by` covariate whose
predictors are a subset of the selected term's predictors, with the same
linear or smooth treatment of each shared predictor. For
`irf(x, z, k_p = list(NULL, 3))` this adds `irf(x)`, `irf(z, k_p = 3)`,
and `irf(1)` when they are present. A time-invariant term is included
below a time-varying term with the same predictors. Under
`grouping = "conditional"`, the grouped deviations of those terms at the
same level are also added.

`reference` subtracts the same functional evaluated at reference values
of some axes. Axes not named in `reference` keep their evaluated values.
The standard errors account for the covariance of the two evaluations. A
predictor whose reference is `NA` or `"default"` is left at the model
default: every term that involves it contributes zero to the reference.
For a centered linear predictor this equals its center. For a smooth
predictor it is the training-average contribution implied by the term's
identifiability constraint. With every predictor at its default, the
reference under `composition = "total"` reduces to `irf(1)`, the
response to an average event.

For example, with `composition = "total"`, `x` at its mean plus one
standard deviation, and
`reference = list(x = list(at = list(summary = "mean")))`, the result is
the effect of a one-SD increase in `x` at each requested value of `z`. A
reference value of `z` gives differences between the interaction curves.

## See also

[`test_effect()`](https://climblab.org/cdrgam/reference/test_effect.md)
for joint Wald tests of the same functionals.
