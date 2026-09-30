# Suggest model simplifications from a fitted CDR-GAM

The report ranks terms whose fitted effective degrees of freedom have
collapsed toward their penalized limit or whose unresolved outer
derivative favors stronger smoothing. For a nonconverged distributional
fit, it also identifies a parameter submodel when the
simplification-directed score is concentrated there. The function does
not change or refit the model.

## Usage

``` r
suggest_simplifications(object, max_candidates = 10L, conservatism = 1)
```

## Arguments

- object:

  A fitted `cdrgam` model.

- max_candidates:

  Maximum number of candidates to return. Use `Inf` to retain every
  candidate.

- conservatism:

  Strength of the ranking penalty for model degrees of freedom removed
  by a candidate. Zero ranks only by diagnostic evidence.

## Value

A `cdrgam_simplification_report`. Use
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to
extract its ranked candidate table. The `patches` member contains exact,
preconditioned formula replacements for candidates that are safe to
automate mechanically; selecting and refitting them remains the caller's
responsibility.
