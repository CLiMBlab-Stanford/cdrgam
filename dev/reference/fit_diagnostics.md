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

A named list containing at least `outcome`, `reasons`, `converged`,
`code`, and `message`. `outcome` is `"success"`, `"warning"`, or
`"nonconverged"`. A warning means that optimization converged but a
recorded numerical or inferential diagnostic requires review. Custom
backends also report available objective counts, gradient norms,
boundary information, conditioning, and Hessian diagnostics. `notes`
lists findings that do not change `outcome`. A sparse Gaussian REML fit
adds a `boundary_variance_components` note when smoothing parameters are
classified as on the boundary of their range; their profile bounds are
in `hessian_boundary`. The `optimizer_stopped_short` reason means the
reported optimum is not one: either a negative-curvature direction of
the interior Hessian lies along a single smoothing parameter (listed in
`hessian_negative_directions`), or, when the analytic Hessian was
computed, a Newton step would still lower the -2 REML criterion by more
than 0.1 (`predicted_decrease`). Negative curvature that mixes several
smoothing parameters is reported as `indefinite_outer_hessian`. Such a
fit also reports `newton_step_max`, the largest coordinate of that
Newton step on the log smoothing-parameter scale, and
`predicted_decrease_floor`, the smallest criterion change resolvable at
the fitted criterion (`1e-10` times its magnitude, at least `1e-8`). A
dense generalized block fit whose optimizer stalls on a flat criterion
plateau counts as converged only when every gradient component of the -2
LAML criterion is below 0.01 and at least one log smoothing parameter
exceeds 8 in absolute value; it then reports outcome `"warning"`, the
parameters in `plateau`, and the reason `smoothing_parameter_plateau`.
`gradient_norm` is the largest gradient component over all smoothing
parameters and `exact_gradient_maximum` the largest over the free ones.
Along a well-curved direction either can exceed the gradient tolerance
while the Newton step and `predicted_decrease` stay negligible;
convergence is judged from `predicted_decrease` against
`predicted_decrease_floor`.
