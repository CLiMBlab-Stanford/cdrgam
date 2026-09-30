# Generate response-effect plot data

Evaluates predictions, slopes, or comparisons on response-history
scenarios through the marginaleffects extension API. This preserves each
response's remaining impulse history and returns point estimates with
pointwise uncertainty in a tidy table.

## Usage

``` r
estimate_response_effect(
  object,
  newdata = NULL,
  estimand = c("prediction", "slope", "comparison"),
  variables = NULL,
  by = NULL,
  type = c("response", "link"),
  vcov = TRUE,
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A `cdrgam_marginaleffects` view.

- newdata:

  Scenario data, usually returned by
  [`response_scenarios()`](https://climblab-stanford.github.io/cdrgam/dev/reference/response_scenarios.md).
  The default uses the unmodified incidences in the view.

- estimand:

  Quantity to estimate: fitted response values, numeric slopes, or
  finite comparisons.

- variables:

  Predictor specification passed to marginaleffects. It is required for
  slopes and comparisons.

- by:

  Optional columns defining averages. `NULL` retains one result per
  scenario row.

- type:

  Evaluate on the response or linear-predictor scale.

- vcov:

  Include fitted-coefficient uncertainty, or supply a covariance matrix
  or function accepted by
  [`marginaleffects_view()`](https://climblab-stanford.github.io/cdrgam/dev/reference/marginaleffects_view.md).

- level:

  Confidence level for pointwise intervals.

- ...:

  Additional arguments passed to the selected marginaleffects function.

## Value

A `cdrgam_response_effect_grid` data frame. Adapter identity columns and
scenario axes are retained when the requested aggregation permits it.

## Details

Grouping columns are passed to marginaleffects. Use `.cdrgam_scenario`
to average each scenario over its reference incidences, and add
`.cdrgam_response` or `.cdrgam_impulse` when a response- or
event-weighted result is required.
