# List estimable impulse-response effects

Returns the term and source-scale axis metadata needed to construct
effect evaluation grids.

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

A data frame with one row per impulse-response term. The `predictors`,
`axes`, `axis_summaries`, and `centers` columns are list columns.
`centers` gives the source-unit value at which each centered linear
predictor makes the term zero; it is empty for terms without centered
linear predictors.
