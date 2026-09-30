# Evaluate fitted impulse-response function terms

Evaluates each fitted IRF term using its stored marginal basis, fitted
coefficients, and covariance matrix.

## Usage

``` r
estimate_irf(
  object,
  term = NULL,
  lag = NULL,
  n = 200,
  predictor = NULL,
  n_predictor = 25,
  at = list(),
  group = NULL,
  se = TRUE,
  unconditional = FALSE
)
```

## Arguments

- object:

  A fitted `cdrgam` model.

- term:

  Term name or one-based term number; omit for every IRF.

- lag:

  Optional finite numeric lag grid.

- n:

  Number of default grid points.

- predictor:

  Optional predictor-value grid for nonlinear IRFs.

- n_predictor:

  Number of default predictor-value grid points.

- at:

  Named list of evaluation values for smooth response-time and predictor
  axes. `predictor` remains an alias for the first non-lag axis;
  additional axes default to their fitted-grid medians.

- group:

  Optional grouping levels for grouped IRF deviations. By default all
  fitted levels are returned.

- se:

  Include pointwise standard errors when covariance is available.

- unconditional:

  Include the first-order smoothing-parameter uncertainty correction
  when supported by the fitted backend.

## Value

A data frame containing term, grouping level, lag, predictor value,
estimate, and standard error.
