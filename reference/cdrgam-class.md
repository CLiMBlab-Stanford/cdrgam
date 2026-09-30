# CDR-GAM model class helpers

Fits from the native backend subclass their mgcv `gam` or `bam` result.
Custom-backend fits implement the corresponding standard model methods
without inheriting from `gam`. These helpers test the common CDR-GAM
class or remove its subclass from a native-backend fit.

## Usage

``` r
is_cdrgam(x)

as_gam(x)
```

## Arguments

- x:

  An R object; for `as_gam()`, a fitted `cdrgam` object from the native
  mgcv backend.

## Value

`is_cdrgam()` returns a logical scalar. `as_gam()` returns the same
fitted model with the `cdrgam` subclass removed.
