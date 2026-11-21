# List the estimable impulse-response effects in a fitted model

Returns the term and axis metadata used to validate and construct effect
evaluation grids. Axis ranges and summaries are reported in the source
units used by callers.

## Usage

``` r
effect_catalog(object, impulses = NULL, responses = NULL)
```

## Arguments

- object:

  A fitted `cdrgam` model.

- impulses, responses:

  Optional training streams used to recover axis summaries from models
  fitted before summaries were stored.

## Value

A data frame with one row per impulse-response term. `predictors`,
`axes`, `axis_summaries`, and `centers` are list columns. `centers`
gives the source-unit value at which each centered linear predictor
makes the term zero; it is empty for terms without centered linear
predictors.
