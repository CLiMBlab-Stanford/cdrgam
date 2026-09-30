# Extract smoothing variance components

Reports the variance components implied by the smoothing parameters of
the expanded mgcv model. Native fits are passed directly to
[`mgcv::gam.vcomp()`](https://rdrr.io/pkg/mgcv/man/gam.vcomp.html). For
block-backend fits, the same calculation is applied to the stored REML
scale, penalty rescaling metadata, and outer Hessian. The function does
not estimate covariance parameters beyond those present in the expanded
model.

## Usage

``` r
variance_components(object, rescale = TRUE, conf.lev = 0.95)
```

## Arguments

- object:

  A fitted `cdrgam` model.

- rescale:

  Apply the penalty rescaling used in the original smooth specification,
  as in
  [`mgcv::gam.vcomp()`](https://rdrr.io/pkg/mgcv/man/gam.vcomp.html).

- conf.lev:

  Confidence level for intervals when an outer REML Hessian is
  available.

## Value

The value returned by
[`mgcv::gam.vcomp()`](https://rdrr.io/pkg/mgcv/man/gam.vcomp.html):
normally a matrix of standard deviations and confidence intervals, or a
named vector/list when intervals are unavailable.
