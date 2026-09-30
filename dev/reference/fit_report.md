# Extract Stable Fit Metadata and Summary Tables

`fit_report()` returns the fitted-model information needed by artifact
writers and other downstream tools without requiring them to inspect the
internal representation of a `cdrgam` object.

## Usage

``` r
fit_metadata(object)

fit_report(object)
```

## Arguments

- object:

  A fitted `cdrgam` model.

## Value

`fit_metadata()` returns fitting metadata, formulas, distributional
parameter names, and default plotting metadata. `fit_report()` adds
printable summary text, parametric and smooth summary tables, and
convergence diagnostics.
