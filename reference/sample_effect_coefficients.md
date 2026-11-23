# Sample coefficients shared by several effect simulations

Draws once from the approximate coefficient posterior \\N(\hat\beta,
V)\\, restricted to the coefficients that the selected terms use under
`composition` and `grouping`. Passing the result to
[`simulate_effect()`](https://climblab.org/cdrgam/reference/simulate_effect.md)
as `coefficients` evaluates every call on the same draws, so summaries
computed from different calls are joint: for example, the lag at which
one contrast peaks in each draw and another effect at that lag in the
same draw.

## Usage

``` r
sample_effect_coefficients(
  object,
  terms = NULL,
  impulses = NULL,
  responses = NULL,
  composition = c("term", "total"),
  grouping = c("population", "deviation", "conditional"),
  groups = NULL,
  unconditional = FALSE,
  draws = 1000L,
  seed = NULL,
  boundary = c("estimate", "bound")
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

- unconditional:

  Include smoothing-parameter uncertainty when supported. Smoothing
  parameters classified as on the boundary of their range contribute no
  uncertainty, because the criterion has no local curvature along them.
  Sparse fits stop with an error when the outer Hessian is not positive
  definite after those parameters are excluded.

- draws:

  Number of posterior draws.

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

## Value

An object of class `cdrgam_coefficient_draws`: a list with `index`, the
coefficient positions; `draws`, a matrix with one row per coefficient in
`index` and one column per draw; `unconditional`; and `boundary`.

## See also

[`simulate_effect()`](https://climblab.org/cdrgam/reference/simulate_effect.md).
