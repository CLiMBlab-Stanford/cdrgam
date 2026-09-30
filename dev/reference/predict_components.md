# Predict Response and Link Components

`predict_components()` returns a stable representation of response- and
link-scale predictions for artifact writers and other downstream tools.
Ordinary models return vectors. Distributional models additionally
return matrices with one column per response-distribution parameter.

## Usage

``` r
predict_components(object, newdata)
```

## Arguments

- object:

  A fitted `cdrgam` model.

- newdata:

  Prediction data accepted by
  [`predict.cdrgam()`](https://climblab-stanford.github.io/cdrgam/dev/reference/predict.cdrgam.md).

## Value

A named list containing primary response- and link-scale estimates and
standard errors. Distributional fits also include parameter names and
the complete parameter matrices on both scales.
