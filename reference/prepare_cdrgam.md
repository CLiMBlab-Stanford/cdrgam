# Compile untiled streams for CDR-GAM fitting

Validates and aligns untiled streams, selects or applies a history
layout, and constructs response-level spline designs in bounded-memory
chunks. Requested axis basis dimensions are reduced when the linked data
support fewer distinct values and the marginal basis minimum remains.
Nonlinear IRFs with one or two linked predictor values are represented
as linear IRFs using their original numeric coding. Optional scaling is
an upstream data transformation and does not modify the user,
normalized, or effective model formula. Its training-data divisors and
exclusions are stored in the returned object's `scaling` component.

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
  rescale_predictors = FALSE,
  quiet = FALSE,
  drop.unused.levels = TRUE
)
```

## Arguments

- formula:

  Extended formula containing ordinary terms and optional
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) terms. An
  implicit `irf(1)` is added unless suppressed.

- impulses:

  Data frame with one row per impulse.

- responses:

  Data frame with one row per response.

- window:

  Default inclusive lag window inherited by each
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) term that does
  not supply its own window. `NULL` retains the term-level default
  `c(0, Inf)`. Every impulse in the applicable series and window is
  included in the response-level design.

- knots_l:

  Optional model-level default for
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) lag-basis
  construction points. A term-level value, including `NULL`, takes
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

  Character vector naming independent-series columns.

- impulse_time:

  Name of the impulse-time column.

- response_time:

  Name of the response-time column.

- history:

  One of `"auto"`, `"dense"`, or `"ragged"`.

- chunk_size:

  Maximum responses or links transformed in a chunk.

- rescale_predictors:

  If `TRUE`, divide continuous numeric predictors by their training-data
  standard deviations without centering. Binary, constant, categorical,
  grouping, and series variables are unchanged. Lag and explicit time
  predictors share the standard deviation of the untiled impulse-time
  column, after history links have been constructed in native time
  units.

- drop.unused.levels:

  Drop factor levels without training observations. The default matches
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html). Retained
  levels remain part of the fitted factor vocabulary.

- quiet:

  Suppress preparation-plan messages.

## Value

A reusable `cdrgam_design` object. Its `simplifications` data frame
records automatic axis-basis reductions and lower-dimensional
representations of discrete predictor axes.
