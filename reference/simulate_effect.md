# Simulate fitted effects from the coefficient posterior

Draws coefficient vectors from the approximate posterior \\N(\hat\beta,
V)\\ and evaluates the functionals that
[`estimate_effect()`](https://climblab.org/cdrgam/reference/estimate_effect.md)
would return for each draw. Use the draws for uncertainty in nonlinear
summaries of an effect, such as the lag at which a response peaks, that
pointwise standard errors do not determine. The arguments select the
same terms, axes, compositions, groupings, and references.

## Usage

``` r
simulate_effect(
  object,
  terms = NULL,
  axes = list(),
  impulses = NULL,
  responses = NULL,
  composition = c("term", "total"),
  grouping = c("population", "deviation", "conditional"),
  groups = NULL,
  reference = list(),
  unconditional = FALSE,
  draws = 1000L,
  n = 200L,
  seed = NULL,
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

- unconditional:

  Include smoothing-parameter uncertainty when supported. Smoothing
  parameters classified as on the boundary of their range contribute no
  uncertainty, because the criterion has no local curvature along them.
  Sparse fits stop with an error when the outer Hessian is not positive
  definite after those parameters are excluded.

- draws:

  Number of posterior draws.

- n:

  Default number of points for fitted-domain grids.

- seed:

  Optional integer seed. When given, the global random-number state is
  restored on exit.

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

A list with `grid`, a data frame like the one returned by
[`estimate_effect()`](https://climblab.org/cdrgam/reference/estimate_effect.md)
with `estimate` but no interval columns, and `draws`, a matrix with one
row per row of `grid` and one column per draw.

## Details

Draws for each term and grouping level are sampled separately, so draws
are joint across the rows of one functional but not across functionals.
`V` is the conditional Bayesian covariance, or the covariance corrected
for smoothing-parameter uncertainty when `unconditional = TRUE`. Draws
are on the linear-predictor scale. A rank-deficient `V` is handled
through its eigendecomposition, and directions with negative eigenvalues
contribute nothing.

## See also

[`estimate_effect()`](https://climblab.org/cdrgam/reference/estimate_effect.md)
for pointwise intervals of the same functionals.
