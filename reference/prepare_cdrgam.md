# Compile untiled impulse and response streams for CDR-GAM fitting

`prepare_cdrgam()` validates and aligns untiled streams, selects or
applies a history layout, and constructs response-level spline designs
in bounded-memory chunks. The returned design can be fitted with
[`cdrgam.fit()`](https://climblab.org/cdrgam/reference/cdrgam.md).

## Usage

``` r
prepare_cdrgam(
  formula,
  impulses,
  responses,
  window = NULL,
  knots_l = NULL,
  k_l = NULL,
  k_t = NULL,
  k_p = NULL,
  bs_l = NULL,
  bs_t = NULL,
  bs_p = NULL,
  series = character(),
  impulse_time = "time",
  response_time = "time",
  history = c("auto", "dense", "ragged"),
  chunk_size = 10000,
  rescale_predictors = TRUE,
  center_predictors = TRUE,
  quiet = FALSE,
  drop.unused.levels = TRUE
)
```

## Arguments

- formula:

  Extended model formula containing ordinary `mgcv` terms and optional
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) terms. An
  implicit `irf(1)` is added unless suppressed. For `family="gaulss"`,
  this may instead be a named list with `location` and `scale` formulas.
  The location formula supplies the response; the scale formula may be
  one-sided.

- impulses:

  Data frame with one row per impulse.

- responses:

  Data frame with one row per response.

- window:

  Default inclusive lag window inherited by every
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) term that does
  not define its own `window`. `NULL` retains the term-level default
  `c(0, Inf)`. Every impulse in the applicable series and window is
  included in the response-level design.

- knots_l:

  Optional model-level default for
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) lag-basis
  construction points. Term-level `knots_l`, including `NULL`, takes
  precedence.

- k_l, k_t, k_p:

  Optional model-level defaults for the corresponding
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) axis
  dimensions. An argument supplied by an individual term, including an
  explicit `NULL`, overrides its model-level default.

- bs_l, bs_t, bs_p:

  Optional model-level defaults for the corresponding
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) marginal basis
  names. Term-level arguments take precedence.

- series:

  Character vector of columns identifying independent series.

- impulse_time:

  Name of the impulse-time column.

- response_time:

  Name of the response-time column.

- history:

  One of `"auto"`, `"dense"`, or `"ragged"`.

- chunk_size:

  Maximum history links transformed at once in ragged mode, or responses
  transformed at once in dense mode.

- rescale_predictors:

  Divide continuous numeric predictors and internal lag coordinates by
  training-data standard deviations. Scaling is enabled by default; set
  `FALSE` to fit in native units. Binary, constant, categorical,
  grouping, and series variables are not changed. Lag and explicit time
  coordinates share one divisor, the standard deviation of the untiled
  impulse-time column, applied after history links are built in native
  time units.

- center_predictors:

  Center linear IRF predictors, including binary ones, at their means
  over the training history links. Every term is then defined relative
  to an average event: removing `irf(1)` removes the response to an
  average event, and omitting a lower-order term of an interaction
  removes that effect at the other predictors' means. The prepared
  design's `centering` table records each predictor's center. With
  centering, `irf(1)` is the response to an event whose linear
  predictors are at those means, and each linear term is the deviation
  from that response. Lower- order terms of a linear interaction are
  effects with the other predictors at their means. Effects and
  predictions still use source predictor values. `FALSE` restores
  uncentered predictors, where `irf(1)` is the response at zero.
  Individual terms can override this with `irf(center=)`.
  Response-aligned `by` covariates are never centered.

- quiet:

  Suppress the preparation plan message.

- drop.unused.levels:

  Drop factor levels without training observations. The default matches
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html). Retained
  levels remain part of the fitted factor vocabulary. For
  [`cdrgam.fit()`](https://climblab.org/cdrgam/reference/cdrgam.md), use
  `NULL` to inherit the choice already compiled into `design`; a
  supplied value must match it.

## Value

A reusable `cdrgam_design` object. Its `simplifications` data frame
records automatic reductions of axis basis dimensions and
lower-dimensional representations of discrete predictor axes.

## Details

Requested axis basis dimensions are reduced when the linked data support
fewer distinct values and the marginal basis minimum remains. Nonlinear
IRFs with one or two linked predictor values are represented as linear
IRFs using their original numeric coding.

Optional scaling is an upstream data transformation and does not modify
the raw, normalized, or effective model formula. Its training-data
divisors and exclusions are stored in the returned object's `scaling`
component.
