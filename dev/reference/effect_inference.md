# Inference for summary statistics of fitted effects

Estimates, interval-estimates, and tests summary statistics of one or
more fitted effects, such as contrasts between groups, areas under an
impulse-response function, averages over observed covariates, or peak
lags. Every effect is evaluated as in
[`estimate_effect()`](https://climblab.org/cdrgam/dev/reference/estimate_effect.md),
and all of them are evaluated on one joint distribution of the
coefficients they use, so summaries that combine effects with shared
coefficients have the right variance.

## Usage

``` r
effect_inference(
  object,
  effects,
  estimand = NULL,
  null = 0,
  method = c("auto", "wald", "simulation", "delta"),
  level = 0.95,
  draws = 2000L,
  seed = NULL,
  unconditional = FALSE,
  boundary = c("estimate", "bound"),
  impulses = NULL,
  responses = NULL,
  keep_draws = FALSE
)
```

## Arguments

- object:

  A fitted `cdrgam` model.

- effects:

  A named list of effect specifications. Each is a list with any of
  `terms`, `axes`, `composition`, `grouping`, `groups`, `reference`,
  `n`, and `points`, as in
  [`estimate_effect()`](https://climblab.org/cdrgam/dev/reference/estimate_effect.md).

- estimand:

  `NULL`, a matrix with one row per estimand and one column per effect
  row (effects stacked in order), a named list whose elements are lists
  of weights named by effect (a weight of length one is recycled over
  that effect's rows), or a function. A function receives a named list
  of numeric effect values and a named list of effect grids, and returns
  a numeric vector, named to label the estimands.

- null:

  The value of each estimand under the null hypothesis, of length one or
  one per estimand.

- method:

  `"auto"` uses Wald inference for linear estimands and simulation for
  function estimands. `"delta"` applies only to function estimands.

- level:

  Confidence level for `lower` and `upper`.

- draws:

  Number of coefficient draws for simulation, and of Gaussian draws for
  simultaneous critical values.

- seed:

  Optional seed. The previous random state is restored.

- unconditional:

  Include smoothing-parameter uncertainty when supported. Smoothing
  parameters classified as on the boundary of their range contribute no
  uncertainty, because the criterion has no local curvature along them.
  Sparse fits stop with an error when the outer Hessian is not positive
  definite after those parameters are excluded.

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

- impulses, responses:

  Optional training streams used to recover axis summaries from models
  fitted before summaries were stored.

- keep_draws:

  If `TRUE`, store the estimand draws of a simulation as attribute
  `draws`, with one row per estimand.

## Value

A `cdrgam_effect_inference` data frame with one row per estimand:
`estimand`, `estimate`, `std_error`, `lower`, `upper`, `null`,
`statistic` (the standardized difference from `null`), `p_value`,
`simultaneous_lower`, `simultaneous_upper`, `adjusted_p_value`, and
`method`. Attributes hold `covariance`, the estimand covariance;
`joint_test`, a Wald test that every estimand equals its null; `grids`,
the effect grids with their estimates; `level`; and `draws` when kept.

## Details

The effects are a linear map `D beta` of the coefficients. A linear
estimand `L D beta` has covariance `L D V D' L'`, and Wald inference is
exact given the Gaussian approximation to the coefficient posterior.
Give it as a matrix or a list of weights, or leave `estimand` `NULL` to
make every effect row an estimand.

A function estimand is nonlinear in general. `method = "simulation"`
draws `beta ~ N(beta-hat, V)` and applies the function to each draw. The
estimate is the function at `beta-hat`, the standard error is the
standard deviation of the draws, the interval uses their quantiles, and
the p-value is the two-sided posterior tail area beyond `null`, whose
resolution is `1 / draws`. `method = "delta"` instead linearizes the
function at the estimate by central differences. It is invalid for
functions that are not differentiable there, such as maxima and
argmaxima.

Simultaneous intervals and adjusted p-values cover all estimands of one
call. They use the maximum over estimands of the standardized absolute
deviation, from Gaussian draws with the estimand correlation for Wald
and delta inference, or from the coefficient draws for simulation. The
adjusted p-values are single-step max-modulus p-values, estimated from
those draws and never below the unadjusted p-value.

## See also

[`estimate_effect()`](https://climblab.org/cdrgam/dev/reference/estimate_effect.md),
[`test_effect()`](https://climblab.org/cdrgam/dev/reference/test_effect.md),
[`simulate_effect()`](https://climblab.org/cdrgam/dev/reference/simulate_effect.md)
