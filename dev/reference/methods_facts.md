# Describe estimation procedures for methods reporting

`methods_facts()` resolves the smoothing-parameter criterion that
[`cdrgam.fit()`](https://climblab.org/cdrgam/dev/reference/cdrgam.md)
uses for a configuration and returns short descriptions of the
estimation procedure and the significance tests, each with the keys of
the references that support it. Report writers can use these facts
instead of restating the defaults and backend behavior of this package.

## Usage

``` r
methods_facts(
  backend = c("mgcv", "block", "sparse"),
  engine = c("bam", "gam"),
  method = NULL,
  family = stats::gaussian()
)
```

## Arguments

- backend, engine, method:

  As in
  [`cdrgam.fit()`](https://climblab.org/cdrgam/dev/reference/cdrgam.md).
  `method = NULL` selects the default criterion.

- family:

  A family object or family name, as in
  [`cdrgam.fit()`](https://climblab.org/cdrgam/dev/reference/cdrgam.md).

## Value

A list with elements `backend`, `engine` (`NULL` unless the backend is
`"mgcv"`), `method` (the resolved criterion), `criterion` and
`estimator` (sentence fragments), `citations` (a list of reference keys
with elements `criterion` and `estimator`), and `term_tests`, a data
frame with the columns `method` (the test labels used by
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`test_effect()`](https://climblab.org/cdrgam/dev/reference/test_effect.md)),
`description`, and `citation`. The reference keys name entries in the
bibliography that ships with `cdrgam.cli`.

## Details

The criterion follows the fitting code: with `method = NULL`, the native
`mgcv` backend uses the `mgcv` default (`"fREML"` for `bam()`, `"REML"`
for `gam()` with an extended or general family, and `"GCV.Cp"`
otherwise), and the custom backends use REML.
