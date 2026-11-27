# Export an exact native mgcv fitting specification

Convert a prepared CDR-GAM design into arguments accepted by
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html). Impulse-response
terms become meaningfully named matrix terms, and native penalty
interfaces carry their exact penalties. Scalar-response specifications
use `paraPen`. Distributional specifications use `mgcv`'s prefit `G`
interface because [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html)
does not accept `paraPen` with a formula list. The export does not use
CDR-GAM's private compressed smooth constructor.

## Usage

``` r
as_mgcv_spec(x, family = stats::gaussian(), method = NULL, ...)
```

## Arguments

- x:

  A `cdrgam_design` returned by
  [`prepare_cdrgam()`](https://climblab.org/cdrgam/dev/reference/prepare_cdrgam.md).
  Boundary-reduced fitted models are not supported.

- family:

  A family accepted by
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html). Distributional
  designs require their corresponding general family.

- method:

  Optional smoothing-parameter estimation method.

- ...:

  Additional arguments to include in the
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) call.

## Value

A named argument list suitable for `do.call(mgcv::gam, spec)`.

## Details

The returned data materialize grouped impulse-response designs. This is
an interoperability and validation interface, not the memory-efficient
fitting path used by the sparse backend.
