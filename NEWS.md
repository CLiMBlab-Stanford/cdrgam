# cdrgam (development version)

- Memory-limited parallel plans now charge each forked worker half of the
  parent process's resident memory in addition to its own estimated
  allocations, because a child's garbage collector causes inherited pages to
  be copied. For a 126,058-coefficient model the previous estimate allowed
  two finite-difference workers that together added about 100 GB and
  exceeded a 160 GB allocation. The sparse fit log reports the selected
  gradient workers and the memory charged to each.

- Chunked cross products now multiply columns that are nonzero in most rows
  of a chunk, such as population IRF terms, as a dense block with BLAS, and
  only the remaining active columns with the sparse product. For a 200-subject
  version of a 126,058-coefficient reading model, one 10,000-row chunk took
  0.48 s instead of 9.7 s, with results equal to within 4e-15 relative. This
  is the dominant cost of forming the penalized system for such models.

- The dense location-scale inner fit now counts a stalled line search as
  converged when even its smallest step could not change the objective by
  more than its numerical floor. Its absolute gradient tolerance could lie
  below that floor after predictor rescaling, so small rescaled fits ended
  with "Dense gaulss LAML optimization ended at an invalid fit".

- Sparse backends now add chunked cross products pairwise when forming the
  penalized system and the PIRLS and location-scale information matrices.
  Adding each chunk to a running total rewrote the whole accumulated matrix
  once per chunk; for a 155,698-coefficient model with 205 million system
  nonzeros and 355 chunks, cross-product accumulation took about 2.7 hours.
  On a smaller model with the same structure, pairwise accumulation was 7.4
  times faster.

- A sparse location-scale fit keeps its optimizer endpoint unless an earlier
  evaluation is lower by more than `1e-9` times the criterion. Previously any
  lower evaluation, including one that differed only by inner-solve noise,
  replaced an endpoint that had met the gradient tolerance and marked the fit
  as nonconverged.

- `rescale_predictors` now defaults to `TRUE` in `prepare_cdrgam()` and
  `cdrgam()`, so continuous predictors and lag coordinates are divided by
  their training standard deviations unless `rescale_predictors = FALSE`.
  Effect estimates and plots are still reported in native units.
  `center_predictors` already defaulted to `TRUE`.

- The exact phase of sparse location-scale fits now solves inner fits to
  `sparse_control$exact_inner_relative_tolerance` (default the tighter of
  `1e-11` and `inner_relative_tolerance`). At the previous `1e-8` relative
  tolerance, inner-fit error of about 0.002 on a criterion near 10^6 exceeded
  the remaining trust-region improvements, so fits stalled and were reported
  as nonconverged.

- The safeguarded trust-region optimizer used by the sparse backends starts
  assessing convergence as soon as three successive rejected steps, each no
  longer than `1e-5` times the largest absolute parameter (at least
  `1e-5`), lose more criterion than the model predicted. When the caller
  supplies a curvature assessment, it certifies the result as it did after
  the radius collapsed; otherwise the optimizer reports "converged to
  criterion precision". It previously kept shrinking the radius until it
  reported "trust radius became too small".

- Sparse Gaussian location-scale REML fits warm up with extended
  Fellner-Schall smoothing-parameter updates (Wood and Fasiolo, 2017) instead
  of stochastic scores whenever the hybrid gradient is selected. Each update
  needs one inner fit and penalty traces but no outer score or line search,
  and smoothing parameters still moving toward a limit start the exact phase
  on the search bound, where the projected trust region treats them as
  inactive. `sparse_control$warmup` chooses `"efs"`, `"stochastic"`, or
  `"none"`, and `efs_maxit` limits the updates.

- Generalized sparse REML fits now compute an outer Hessian by central
  differences of their exact score after optimization (`sparse_control$hessian
  = "gradient"`, the default; `"none"` skips it). It enables the boundary
  classification and profile bounds, the `optimizer_stopped_short` and
  `indefinite_outer_hessian` diagnostics, and unconditional inference for
  generalized fits, which previously had no outer Hessian. The cost is two
  score evaluations per outer parameter.

- The exact generalized GCV gradient now streams batches of inverse columns
  through the same routine as the Gaussian analytic Hessian. It previously
  stored a dense selected inverse over the union of penalty supports, whose
  memory grew with the square of that union.

- `fit_diagnostics()` distinguishes an optimizer that stopped short from a
  genuinely indefinite outer Hessian. Negative curvature along a single
  smoothing parameter, or a Newton step that would still lower the -2 REML
  criterion by more than 0.1, is reported as `optimizer_stopped_short`, with
  the parameters in `hessian_negative_directions` and the estimate in
  `predicted_decrease`; the fit warning suggests `final_newton = TRUE`.
  Negative curvature that mixes several smoothing parameters is still
  reported as `indefinite_outer_hessian`.

- Sparse Gaussian REML fits accept `sparse_control$final_newton = TRUE`,
  which refines the smoothing-parameter optimum with safeguarded Newton steps
  using the exact gradient and analytic outer Hessian after the outer
  optimizer stops. It holds boundary smoothing parameters fixed, descends
  along negative curvature, and stops at a positive-definite optimum or after
  `newton_max_steps`. It is on by default; it is skipped, and recorded as
  `skipped_memory`, when the analytic Hessian does not fit the automatic
  memory budget. Set `final_newton = FALSE` to keep the outer optimizer's
  result.

- The penalty log pseudo-determinant and its derivatives now use a stable
  reparameterization (Wood 2011) within each penalty block. Penalties are
  grouped by scaled size, each group's range is fixed from its unit-normalized
  sum, and smaller penalties are carried into the remaining null space. The
  previous single eigendecomposition dropped the eigenvalues of a penalty more
  than about 10^14 times smaller than another in the same block, which
  changed the computed rank and made the criterion discontinuous.

- Unconditional sparse inference (`unconditional = TRUE` in `vcov()`,
  `estimate_effect()`, and `test_effect()`) now stops with an error when the
  outer Hessian is not positive definite after boundary smoothing parameters
  are excluded. It previously inverted an indefinite block whenever the
  inversion succeeded numerically, which could add negative or unbounded
  smoothing-parameter variance to the covariance without warning.

- The sparse condition indicator now measures collinearity rather than
  scale. It is `min_i L_ii^2 / A_ii` over the Cholesky pivots of the
  penalized system `A`, which is one minus the largest squared multiple
  correlation of a coefficient with those eliminated before it. It previously
  squared the ratio of the smallest to the largest pivot, so a single very
  large smoothing parameter could trigger `ill_conditioned_penalized_system`
  with no collinearity present. Gaussian REML fits also leave out the
  coefficients of penalties classified as on the boundary, because a very
  large rank-deficient penalty makes its own coefficients nearly dependent;
  `fit$sparse$condition_indicator_full` keeps the value over all
  coefficients. For one 3.5-million-observation fit the reported indicator
  moved from 1.8e-16 to 3.9e-4.

- Saved sparse fits no longer contain a copy of the prepared design. The
  stored expanded `mgcv` formula took its environment from the function that
  built it, which held the design, so `saveRDS()` wrote the full design with
  every fit. For a 3.5-million-observation fit this was 6.1 of 6.6 GB. The
  formula now uses the environment of the design's ordinary formula.

- Sparse Gaussian REML fits now separate smoothing parameters on the
  boundary of their range from genuine outer-Hessian problems. A parameter
  is on the boundary when its absolute log value is at least
  `sparse_control$boundary_hessian_log_sp` (default `8`) and its whole
  outer-Hessian row is at most `boundary_curvature_tolerance` (default
  `0.05`). Such parameters are listed in
  `fit$sparse$convergence$hessian_boundary`, excluded from the
  positive-definiteness check and from unconditional covariance, and
  reported as a `boundary_variance_components` note in `fit_diagnostics()`
  that does not change `outcome`. Previously they made the Hessian
  indefinite and produced a warning. Flat directions that mix several
  smoothing parameters are still reported as an indefinite Hessian.
  `hessian_min_eigenvalue` now refers to the interior block;
  `hessian_full_min_eigenvalue` keeps the full-matrix value.

- For each smoothing parameter at its upper boundary, a variance component
  estimated as zero, the fit records a profile bound: the smallest log
  smoothing parameter whose -2 REML criterion stays within
  `qchisq(2 * boundary_profile_level - 1, 1)` of the optimum, with the other
  smoothing parameters fixed. With `sparse_control$boundary_sensitivity =
  TRUE`, the fit at the bounds is stored, and `estimate_effect()` and
  `test_effect()` accept `boundary = "bound"` to check whether conclusions
  depend on treating those components as zero.

- The sparse backend's analytic outer Hessian and exact GCV gradient use much
  less memory. They previously kept `A^{-1} S_j` for every smoothing
  parameter at once, which took memory proportional to the number of
  penalized coefficients times the summed penalty supports. They now solve
  for one batch of inverse columns at a time, accumulate every trace that
  uses the batch, and discard it. Batches follow the block structure of the
  penalties, so grouped penalties are split by group, and penalties that
  share a support share one solve. Peak memory scales with the number of
  penalized coefficients times the batch width, which is usually
  `trace_chunk_size` or the size of the largest penalty block.
  `cdrgam.max_analytic_hessian_elements` now counts the elements of one
  batch. The `hessian = "auto"` memory estimate counts every array the
  batched computation allocates, and the analytic Hessian is selected when
  that estimate is at most 75% of available memory and leaves at least 1 GiB
  free. The previous rule allowed at most half of available memory.

- `estimate_effect(composition = "total")` now sums the full hierarchy of
  ungrouped terms below the selected term: terms on the same `by` covariate
  whose predictors are a subset of its predictors, with the same linear or
  smooth treatment. Previously only a single-predictor smooth term or a
  time-varying term found a parent, so `irf(x, z, k_p = list(NULL, 3))` was
  returned without `irf(x)`, `irf(z)`, or `irf(1)`. The total effect of a
  linear single-predictor term now includes `irf(1)`, as it already did for a
  smooth one. Under `grouping = "conditional"`, the grouped deviations of the
  hierarchy terms are added at the same level.

- `estimate_effect()` gains `reference`, which subtracts the functional at
  reference axis values and carries the joint covariance into the standard
  errors. A reference value of a linear predictor gives the effect of
  changing that predictor at each value of its moderators; a reference value
  of a moderator gives differences between interaction curves. A predictor
  set to `NA` or `"default"` is left at the model default, so the terms that
  involve it drop out of the reference.

- The `deviation` composition is removed. It returned the same values as
  `term`; grouped deviations remain available through `grouping`.

- New `test_effect()` computes joint Wald tests that an `estimate_effect()`
  functional is zero over its grid, using a rank-truncated pseudo-inverse.

- The sparse backend now reports exactly collinear IRF terms as rank
  deficient regardless of rounding. Its identifiability check previously
  passed whenever CHOLMOD returned a factor, including after warning that the
  system was not positive definite or when rounding left a near-zero positive
  pivot, so some rank-deficient models were fitted without error. A CHOLMOD
  warning now counts as a failed factorization, and a pivot that is a
  negligible fraction of its diagonal entry, below `rank_tol`, is reported as
  rank deficiency. The ratio does not depend on predictor scaling.

- Linear IRF predictors, including binary ones, are now centered by default at
  their means over the training history links. `irf(1)` is then the response
  to an event with average linear predictors, linear terms are deviations from
  it, and removing `irf(1)` removes the average-event response. Previously `irf(1)` was the response at zero on the coded
  scale, which for a strictly positive predictor is an extrapolation, and the
  two terms were strongly collinear. New fits of existing formulas therefore
  estimate a different `irf(1)` and different lower-order interaction terms;
  total predicted responses are essentially unchanged, except in models that
  omit `irf(1)` or lower-order interaction terms, whose meaning changes. Use
  `center_predictors = FALSE` to restore uncentered predictors, or
  `irf(center = )` per term. Fits made before this change keep their
  uncentered meaning. `effect_catalog()` reports the stored centers.

- Solver checkpoints resume when the stored floating-point summaries of the
  model differ from the current ones only by rounding error, as happens when a
  fit resumes on a node with a different CPU or BLAS configuration. Structural
  fields must still match exactly, and the mismatch error now names the
  fields that differ.

- The sparse backend builds `ti()` terms with one random-effect factor margin,
  such as `ti(time, subject, bs=c("cr", "re"))`, as sparse per-level blocks.
  These terms were previously expanded to a dense design, which could exhaust
  memory with many factor levels. The design and penalties equal those of
  `mgcv::ti()` up to column order. Other `ti()` forms still use the dense
  ordinary design.

# cdrgam 0.2.0

- Establishes the first public development release of the current CDR-GAM
  implementation.
- Provides native, dense block, and sparse fitting backends for continuous-time
  deconvolutional generalized additive models.
- Supports Gaussian, generalized, and Gaussian location-scale models with
  stream-aware prediction, inference, diagnostics, and effect estimation.
- Adds bounded-memory history construction, checkpointed optimization, and
  resource-aware parallel execution.
