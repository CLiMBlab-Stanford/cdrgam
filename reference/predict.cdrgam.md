# Predict from a fitted CDR-GAM

Supply new untiled streams as
`newdata=list(impulses=impulses, responses=responses)`. Stored
training-data scaling is applied before the response-level design is
rebuilt. Native `mgcv` fits also accept a response-side data frame.
Factor labels are mapped to the fitted training vocabulary. Unseen
random-effect and grouped-IRF levels produce a warning and contribute
zero deviation.

## Usage

``` r
# S3 method for class 'cdrgam_block'
predict(object, newdata = NULL, ...)

# S3 method for class 'cdrgam'
predict(object, newdata = NULL, ...)

# S3 method for class 'cdrgam_distributional_block'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "lpmatrix"),
  se.fit = FALSE,
  chunk_size = 10000,
  unconditional = FALSE,
  ...
)

# S3 method for class 'cdrgam_distributional_sparse'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "lpmatrix"),
  se.fit = FALSE,
  chunk_size = 10000,
  unconditional = FALSE,
  ...
)

# S3 method for class 'cdrgam_sparse'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted `cdrgam` model.

- newdata:

  A named list containing `impulses` and `responses` data frames, or a
  response-side data frame for a native `mgcv` fit.

- ...:

  Prediction controls, including `type`, `se.fit`, `chunk_size`, and
  `unconditional`.

## Value

A prediction vector, linear-predictor matrix, or list containing `fit`
and `se.fit`.
