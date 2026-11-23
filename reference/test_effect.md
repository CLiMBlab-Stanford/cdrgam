# Test fitted effects jointly over declarative axis grids

Computes a Wald test of the null hypothesis that a functional returned
by
[`estimate_effect()`](https://climblab.org/cdrgam/reference/estimate_effect.md)
is zero at every grid point. The arguments select the same terms, axes,
compositions, groupings, and references. One test is returned per
selected term and grouping level.

## Usage

``` r
test_effect(
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
  n = 50L,
  tolerance = 1e-08,
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

- n:

  Default number of points for fitted-domain grids.

- tolerance:

  Relative eigenvalue threshold for the rank of `L V L'`.

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

A `cdrgam_effect_test` data frame with `term_id`, `term`, `group`,
`composition`, `grouping`, `statistic`, `df`, `p_value`, and `points`.
`df` is zero and `statistic` and `p_value` are missing when the
functional is identically zero.

## Details

The statistic is \\\hat\beta^\top L^\top (L V L^\top)^{-} L \hat\beta\\,
where `L` stacks the grid rows of the functional and `V` is the
coefficient covariance. The pseudo-inverse keeps eigenvalues above
`tolerance` times the largest, and the reference distribution is
chi-squared with that rank as its degrees of freedom. Use a grid dense
enough to span the functional; the rank stops growing once it does.

Penalized coefficients make this test approximate. Wood (2013) reports
that such tests can reject too rarely when the null holds on the
boundary of the smoothing-parameter space.

## References

Wood, S. N. (2013). On p-values for smooth components of an extended
generalized additive model. *Biometrika*, 100(1), 221–228.

## See also

[`estimate_effect()`](https://climblab.org/cdrgam/reference/estimate_effect.md)
