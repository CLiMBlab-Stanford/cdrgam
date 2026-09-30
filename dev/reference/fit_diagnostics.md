# Report standardized fit convergence diagnostics

Report standardized fit convergence diagnostics

## Usage

``` r
fit_diagnostics(object)
```

## Arguments

- object:

  A fitted `cdrgam` model.

## Value

A named list containing at least `converged`, `code`, and `message`.
Custom backends also report available objective counts, gradient norms,
boundary information, and Hessian diagnostics.
