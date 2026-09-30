# Construct response-history scenarios

Expands response–impulse incidences from a
[`marginaleffects_view()`](https://climblab.org/cdrgam/dev/reference/marginaleffects_view.md)
over a Cartesian grid of lag and impulse-predictor values. Each output
row retains the identity of its reference incidence, so prediction
replaces one event's contribution while holding the rest of that
response's history fixed.

## Usage

``` r
response_scenarios(object, at = list(), rows = NULL)
```

## Arguments

- object:

  A `cdrgam_marginaleffects` view.

- at:

  Named list of values for `lag` or impulse predictors. `lag` is stored
  in the adapter column `.cdrgam_lag`.

- rows:

  Optional integer or logical subset of reference incidences.

## Value

A data frame suitable for
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`estimate_response_effect()`](https://climblab.org/cdrgam/dev/reference/estimate_response_effect.md),
or direct use as `newdata` in marginaleffects.

## Details

Changing `lag` re-evaluates the selected impulse-response function
windows. A term stops contributing when the scenario lag lies outside
its window.
