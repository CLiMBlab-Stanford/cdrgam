# Export an exact native mgcv fitting specification

Converts a prepared CDR-GAM design into native mgcv arguments.
Impulse-response functions become meaningfully named matrix terms.
Scalar models carry their exact penalties through `paraPen`.
Distributional models use mgcv's prefit `G` interface because
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) does not accept
`paraPen` with a formula list. The export does not use CDR-GAM's private
compressed smooth constructor.

## Usage

``` r
as_mgcv_spec(x, family = stats::gaussian(), method = NULL, ...)
```

## Arguments

- x:

  A `cdrgam_design` returned by
  [`prepare_cdrgam()`](https://climblab.org/cdrgam/reference/prepare_cdrgam.md).

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

The returned data materialize grouped impulse-response designs. This
interface supports interoperability and numerical validation; it is not
the memory-efficient fitting path used by the sparse backend.
Boundary-reduced fitted models are not supported.
