# Specify an impulse-response term

Use `irf()` as a top-level additive term inside a
[`cdrgam()`](https://climblab-stanford.github.io/cdrgam/dev/reference/cdrgam.md)
formula. Ordinary terms retain their usual R and mgcv meanings. By
default each predictor enters linearly as a convolution weight. Numeric
`k_p` entries add smooth predictor marginals and `k_t` adds a
response-time marginal. Smooth non-lag marginals are centered, so a
tensor term represents the interaction beyond its lower-order terms.
Before constraints, a tensor term's coefficient count grows as the
product of its axis dimensions. Larger bases are supported, but each
added axis therefore increases preparation, fitting, and inference costs
multiplicatively. With `group`, the term is expanded into
factor-specific deviations with shared roughness penalties and an
additional full-rank shrinkage penalty. Unless explicitly removed with
`- irf(1)`, the formula compiler inserts an implicit `irf(1)`. An
implicit rate term that has zero design after the required centering
constraint is removed with a warning; the same situation is an error
when `irf(1)` was explicitly requested. Requested axis basis dimensions
are upper bounds. The stream compiler reduces a dimension when the
linked data contain fewer distinct values, subject to the marginal basis
minimum, and records the change in the prepared design.

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
  either, the default is `c(0, Inf)`.

- knots_l:

  Optional strictly increasing lag-basis construction points in the
  original lag units. Supply exactly `k_l` values spanning the linked
  training lags. This control is unavailable for `bs_l = "ps"`.

- k_l:

  Lag-basis dimension. This axis is always smooth and cannot be `NULL`.

- k_t:

  Optional response-time basis dimension. `NULL` makes the IRF
  stationary; a number models nonstationarity against the configured
  response time.

- k_p:

  Predictor-axis dimensions. Supply one value to recycle it or one entry
  per predictor. `NULL` means linear and a number gives a smooth
  marginal. Use a list such as `list(NULL, 4)` for mixed axes because R
  reduces `c(NULL, 4)` to `4`.

- bs_l, bs_t:

  Lag and response-time marginal basis names accepted by mgcv.

- bs_p:

  Predictor marginal basis names; one value is recycled or supply one
  per predictor.

- k, bs, nonlinear, varying:

  Compatibility controls for the original single-predictor API. They
  cannot be mixed with explicitly supplied axis-specific controls.

- by:

  Optional unquoted response-aligned numeric covariate multiplying the
  complete IRF contribution.

- group:

  Optional unquoted response-stream factor defining grouped IRF
  deviations. Add a corresponding ungrouped `irf()` term for a
  population-plus-group-specific model.

## Value

An internal `cdrgam_irf_spec` used by the formula compiler.
