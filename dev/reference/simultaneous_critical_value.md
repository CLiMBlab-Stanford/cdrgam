# Simultaneous critical values

Returns the multiplier `c` such that the interval estimate plus or minus
`c` standard errors covers every coordinate of a group jointly with
probability `level`: the `level` quantile of the largest absolute
standardized deviation within the group. This is the max-modulus value
that
[`effect_inference()`](https://climblab.org/cdrgam/dev/reference/effect_inference.md)
uses for its simultaneous intervals.

## Usage

``` r
simultaneous_critical_value(
  covariance = NULL,
  draws = NULL,
  estimate = NULL,
  groups = NULL,
  level = 0.95,
  simulations = 10000L,
  seed = NULL
)
```

## Arguments

- covariance:

  A covariance matrix of the coordinates.

- draws:

  A matrix with one row per coordinate and one column per posterior
  draw.

- estimate:

  The point estimate of each row of `draws`.

- groups:

  Optional group of each coordinate. Each group gets its own critical
  value; a covariance is then split into its diagonal blocks.

- level:

  Joint coverage probability.

- simulations:

  Number of Gaussian draws for a covariance.

- seed:

  Optional seed for those draws. The previous random state is restored,
  so a seed makes the value reproducible without disturbing the caller's
  random stream.

## Value

A numeric vector with the critical value of each coordinate's group, one
element per coordinate.

## Details

Give either `covariance`, for a Gaussian vector with that covariance, or
`draws` and `estimate`, for posterior draws. A covariance gives a value
estimated from `simulations` Gaussian draws and is never below the
pointwise normal quantile. Draws are standardized by each row's draw
standard deviation. Coordinates with no variance are ignored.

## See also

[`effect_inference()`](https://climblab.org/cdrgam/dev/reference/effect_inference.md)
