# CDR-GAM model class helpers

Fits from the native backend subclass their `mgcv` `gam` or `bam`
result. Custom-backend fits implement the corresponding standard model
methods without inheriting from `gam`. `is_cdrgam()` tests for the
common CDR-GAM class. `as_gam()` removes that subclass from a
native-backend fit.

## Usage

``` r
is_cdrgam(x)

as_gam(x)
```

## Arguments

- x:

  An R object; for `as_gam()`, a fitted `cdrgam` object from the native
  `mgcv` backend.

## Value

`is_cdrgam()` returns a logical scalar. `as_gam()` returns the same
fitted model with the `cdrgam` subclass removed, leaving only its native
`mgcv` classes.

## Details

Calling `as_gam()` is rarely necessary: standard `mgcv` methods already
dispatch through the inherited `gam` or `bam` class. It can be useful
for code that checks exact class names. `as_gam()` signals an error for
objects that are not CDR-GAM fits and for custom-backend fits.
