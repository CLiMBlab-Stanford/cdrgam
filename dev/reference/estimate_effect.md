# Evaluate fitted effects on declarative axis grids

Evaluates selected impulse-response terms as joint linear functionals of
the fitted coefficients. Axes may vary over numeric values, quantiles,
or the fitted domain, or be fixed at a training-data summary with an
optional standard deviation offset.

## Usage

``` r
estimate_effect(
  object,
  terms = NULL,
  axes = list(),
  impulses = NULL,
  responses = NULL,
  composition = c("term", "deviation", "total"),
  grouping = c("population", "deviation", "conditional"),
  groups = NULL,
  se = TRUE,
  unconditional = FALSE,
  level = 0.95,
  n = 200L
)
```

## Arguments

- object:

  A fitted `cdrgam` model.

- terms:

  Term labels or a selector with `names`, `predictors`, `match`, and
  `grouped` fields.

- axes:

  Named axis requests. Predictor requests are nested below `predictors`;
  `"*"` supplies a fallback.

- impulses, responses:

  Optional training streams used to recover axis summaries from models
  fitted before summaries were stored.

- composition:

  Whether to return the selected term or add its fitted marginal
  hierarchy.

- grouping:

  Whether grouped terms represent population effects, deviations, or
  their conditional sum.

- groups:

  Optional grouping levels to evaluate.

- se:

  Whether to calculate pointwise standard errors.

- unconditional:

  Whether to include smoothing-parameter uncertainty when the fitted
  backend supports it.

- level:

  Confidence level for `lower` and `upper`.

- n:

  Default number of points for fitted-domain grids.

## Value

A `cdrgam_effect_grid` data frame containing evaluation axes, estimates,
standard errors, and pointwise confidence limits on the linear-predictor
scale.
