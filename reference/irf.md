# Specify an impulse-response term

`irf()` is used inside a
[`cdrgam()`](https://climblab.org/cdrgam/reference/cdrgam.md) formula.
Its predictor arguments name numeric columns in the impulse stream, or
the sole predictor is literal `1` for the deconvolutional intercept.
Ordinary formula terms retain their usual `mgcv` meaning and are
evaluated against the response stream. Use `irf()` only as a top-level
additive term.

## Usage

``` r
irf(
  predictor,
  ...,
  window = NULL,
  knots_l = NULL,
  k_l = 10,
  k_t = NULL,
  k_p = NULL,
  bs_l = "cr",
  bs_t = "cr",
  bs_p = "cr",
  by = NULL,
  group = NULL,
  group_penalty = c("re", "t2", "fs"),
  center = NULL,
  k = NULL,
  bs = NULL,
  nonlinear = NULL,
  varying = NULL
)
```

## Arguments

- predictor:

  Unquoted numeric impulse-stream column, or literal `1` for the
  deconvolutional intercept.

- ...:

  Additional unquoted numeric impulse-stream predictors. Multiple
  predictors define one interaction IRF.

- window:

  Inclusive lag window as `c(minimum, maximum)`, where lag is response
  time minus impulse time. Negative lags include impulses after the
  response time. `NULL` inherits the model-level `window`; without
  either, the default is `c(0, Inf)`. The minimum must be finite; the
  maximum may be infinite.

- knots_l:

  Optional strictly increasing lag-basis construction points in the
  original lag units. Supply exactly `k_l` values spanning the linked
  training lags. This control is unavailable for `bs_l="ps"`, whose knot
  contract differs from the other supported numeric marginals.

- k_l:

  Lag-basis dimension. This axis is always smooth and cannot be `NULL`.

- k_t:

  Optional response-time basis dimension. `NULL` makes the IRF
  stationary; a number models nonstationarity against `response_time`.

- k_p:

  Predictor-axis dimensions. Supply one value to recycle it or one entry
  per predictor. `NULL` means that predictor enters linearly, while a
  number gives a smooth marginal. Use a list, such as `list(NULL, 4)`,
  when mixing linear and smooth predictors because base R reduces
  `c(NULL, 4)` to `4`.

- bs_l, bs_t:

  Lag and response-time marginal basis names accepted by `mgcv`.

- bs_p:

  Predictor marginal basis names. Supply one value to recycle it or one
  entry per predictor.

- by:

  Optional unquoted response-aligned numeric covariate multiplying the
  complete IRF contribution.

- group:

  Optional unquoted response-stream factor defining grouped IRF
  deviations. Add a corresponding ungrouped `irf()` term for a
  population-plus-group-specific model.

- group_penalty:

  Penalty construction for a grouped IRF. `"re"` keeps the parsimonious
  tensor marginal penalties plus one full-rank group shrinkage penalty.
  `"t2"` uses the grouped-null-space construction from
  `mgcv::t2(full=FALSE)`. `"fs"` separates every null-space component as
  in `mgcv::t2(full=TRUE)` and a univariate `bs="fs"` factor smooth.
  `"t2"` and `"fs"` introduce progressively more smoothing parameters as
  the number of smooth IRF axes increases.

- center:

  Centering of linear predictors. `NULL` inherits the model-level
  `center_predictors` setting. `TRUE` or `FALSE` centers all or none of
  this term's linear predictors. For interactions, a logical vector with
  one entry per predictor, such as `c(FALSE, TRUE)`, or a character
  vector of predictor names selects the predictors to center. A centered
  predictor enters as its value minus its mean over the training history
  links, so the term is zero for an event at that mean. Predictors with
  a single linked value are not centered. Smooth predictor axes are
  always centered by their own constraint and are unaffected. `irf(1)`
  has no predictors to center.

- k, bs, nonlinear, varying:

  Compatibility controls for the original single-predictor API. They
  cannot be mixed with explicitly supplied axis-specific controls.

## Value

An internal `cdrgam_irf_spec` when evaluated by the formula compiler.
Before constraints, a tensor term's coefficient count is the product of
its axis dimensions.

## Details

By default each predictor enters linearly as a convolution weight.
Numeric `k_p` entries add smooth predictor marginals, and `k_t` adds a
response-time marginal. Smooth non-lag marginals are centered, so a
tensor term represents the interaction beyond its lower-order terms.
Linear predictors are centered at their training means by default; see
`center` and the `center_predictors` argument of
[`cdrgam()`](https://climblab.org/cdrgam/reference/cdrgam.md).

Before constraints, a tensor term's coefficient count is the product of
its axis dimensions. Each added axis therefore multiplies preparation,
fitting, and inference costs.

With `group`, the term is expanded into factor-specific deviations, and
`group_penalty` selects their shared penalty decomposition.

Unless removed with `- irf(1)`, the formula compiler inserts an implicit
`irf(1)`. An implicit rate term whose design is zero after the required
centering constraint is removed with a warning. The same condition is an
error when `irf(1)` was requested explicitly.

Requested axis basis dimensions are upper bounds. When the linked data
contain fewer distinct values, preparation reduces the dimension,
subject to the marginal basis minimum, and records the change in the
prepared design.
