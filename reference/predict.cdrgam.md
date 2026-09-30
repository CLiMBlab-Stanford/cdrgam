# Predict from a fitted CDR-GAM

Rebuilds the response-level design from new untiled impulse and response
streams. Stored training-data scaling is applied before prediction.
Unseen random-effect and grouped-IRF levels produce a warning and have
their deviations set to zero.

## Usage

``` r
# S3 method for class 'cdrgam'
predict(object, newdata = NULL, ...)

# S3 method for class 'cdrgam_block'
predict(object, newdata = NULL, ...)

# S3 method for class 'cdrgam_sparse'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted `cdrgam` model.

- newdata:

  A named list containing `impulses` and `responses` data frames. Native
  mgcv fits also accept a response-side data frame.

- ...:

  Prediction controls passed to the CDR stream predictor, including
  `type`, `se.fit`, `chunk_size`, and `unconditional`.

## Value

A prediction vector, linear-predictor matrix, or a list containing `fit`
and `se.fit`.
