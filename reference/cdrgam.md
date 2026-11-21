# Fit a continuous-time deconvolutional GAM

`cdrgam.fit()` fits an already prepared `cdrgam_design`, such as one
returned by
[`prepare_cdrgam()`](https://climblab.org/cdrgam/reference/prepare_cdrgam.md),
without rebuilding its stream histories. It does not parse formulas.

`cdrgam()` compiles untiled impulse and response streams and fits the
resulting model. Use
[`prepare_cdrgam()`](https://climblab.org/cdrgam/reference/prepare_cdrgam.md)
followed by `cdrgam.fit()` when a compiled design will be reused across
fits.

## Usage

``` r
cdrgam.fit(
  design,
  family = stats::gaussian(),
  method = NULL,
  engine = c("bam", "gam"),
  backend = c("mgcv", "block", "sparse"),
  checkpoint = NULL,
  solver_trace = FALSE,
  sparse_control = list(),
  rank_action = c("error", "minimum_norm", "drop", "penalize"),
  rank_tol = NULL,
  rank_penalty = NULL,
  drop.unused.levels = NULL,
  ...
)

cdrgam(
  formula,
  impulses,
  responses,
  window = NULL,
  knots_l = NULL,
  k_l = NULL,
  k_t = NULL,
  k_p = NULL,
  bs_l = NULL,
  bs_t = NULL,
  bs_p = NULL,
  series = character(),
  impulse_time = "time",
  response_time = "time",
  history = c("auto", "dense", "ragged"),
  chunk_size = 10000,
  rescale_predictors = TRUE,
  center_predictors = TRUE,
  family = stats::gaussian(),
  method = NULL,
  engine = c("bam", "gam"),
  backend = c("mgcv", "block", "sparse"),
  checkpoint = NULL,
  solver_trace = FALSE,
  sparse_control = list(),
  rank_action = c("error", "minimum_norm", "drop", "penalize"),
  rank_tol = NULL,
  rank_penalty = NULL,
  drop.unused.levels = TRUE,
  ...
)
```

## Arguments

- design:

  A reusable `cdrgam_design` returned by
  [`prepare_cdrgam()`](https://climblab.org/cdrgam/reference/prepare_cdrgam.md).

- family:

  A standard family object or a family name accepted by
  [`cdrgam_family()`](https://climblab.org/cdrgam/reference/cdrgam_family.md).
  `gaulss` rejects non-unit prior weights because mgcv's family accepts
  but does not use them.

- method:

  Smoothing-parameter estimation method. The sparse Gaussian backend
  accepts REML (the default, also spelled `"fREML"`) and `"GCV.Cp"`.
  Sparse GCV uses the exact effective degrees of freedom in the GCV
  score and accepts mgcv's positive `gamma` multiplier through `...`.
  Its initial implementation uses L-BFGS-B with a memory-gated exact
  gradient and a parallel finite-difference fallback. It omits
  smoothing-parameter uncertainty and does not support automatic
  boundary reduction or rank regularization. The single-predictor
  generalized sparse solver also accepts `"GCV.Cp"`: it minimizes UBRE
  for binomial and Poisson fits and GCV for Gamma fits using streamed
  PIRLS and parallel finite differences. The sparse `gaulss` solver
  accepts REML (the default) and `"QNCV"`. QNCV uses a streamed
  quadratic leave-one-response-out approximation and an exact analytic
  smoothing-parameter gradient. It accepts mgcv's positive `gamma`
  multiplier through `...`; custom neighborhoods are not yet supported.

- engine:

  Either `"bam"` or `"gam"` for the native backend.

- backend:

  Fitting backend. `"mgcv"` uses native fitting, `"block"` selects the
  dense reference solver (Gaussian REML, generalized LAML, or joint
  Gaussian location–scale LAML), and `"sparse"` uses streamed sparse
  Gaussian, supported generalized, or joint Gaussian location–scale
  optimization. Under REML, the distributional sparse solver ends with
  exact outer scores; large systems use matrix-free stochastic scores
  for a warm-up by default, then use safeguarded trust-region BFGS for
  exact refinement. Distributional QNCV uses exact gradients throughout.
  Select the Gaussian trust-region optimizer with
  `sparse_control=list(gradient="exact", outer_optimizer="bfgs_trust")`.

- checkpoint:

  Optional checkpoint path for a custom backend. Checkpoints are written
  atomically after periodic smoothing-criterion evaluations and contain
  validated optimization/restart state. The trust optimizer checkpoints
  its dense curvature approximation, trust radius, counters, parameters,
  score, and criterion so interruption resumes the same optimization
  trajectory; older checkpoints without that state remain parameter-only
  warm starts. Reusing the path resumes an interrupted fit or skips an
  already completed outer optimization. Sparse distributional fits
  checkpoint accepted outer states and resume from the best retained
  smoothing parameters; the joint factorization is reconstructed once on
  resume.

- solver_trace:

  Custom-backend progress reporting. `FALSE` or `0` is silent; `TRUE` or
  `1` reports phases, improving solutions, and trust-optimizer
  diagnostics; `2` reports every objective evaluation; and `3` adds
  low-level chunk/factorization events. Trust diagnostics include
  projected-gradient and tolerance ratios, step norms, trust radius,
  predicted and actual improvement, acceptance statistics, and expected
  post-fit Hessian work. At the numerical trust-radius floor, the sparse
  trust optimizer uses the resolved Hessian strategy to either certify
  practical convergence or perform a bounded positive-curvature reset.
  Automatic selection uses the same process memory budget as the
  post-fit Hessian calculation. A function receives the same progress
  events as named lists.

- sparse_control:

  Named control list for the sparse backend. The single-predictor
  generalized sparse solver accepts `crossprod_chunk_size`,
  `supernodal`, `optimizer_maxit`, `optimizer_gradient_tolerance`,
  `optimizer_trust_radius`, `cores`, `gradient`, `score_workers`,
  `score_batch_size`, `finite_difference_step`, `hessian`,
  `hessian_step`, and the boundary classification options
  `boundary_curvature_tolerance`, `boundary_hessian_log_sp`,
  `boundary_profile`, and `boundary_profile_level`. Generalized REML
  uses an exact gradient and safeguarded trust-region optimization.
  After optimization, `hessian = "gradient"` (the default)
  central-differences the exact score with step `hessian_step` (default
  `1e-2`), two score evaluations per outer parameter, to obtain the
  outer Hessian used for boundary classification, profile bounds,
  optimizer diagnostics, and unconditional inference; `hessian = "none"`
  skips it. Generalized `"GCV.Cp"` uses L-BFGS-B and accepts
  `gradient="auto"` (the default), `"exact"`, or `"finite"`. Automatic
  selection compares the memory and structural cost of a
  selected-inverse analytic gradient with two conditional PIRLS fits per
  smoothing parameter. The exact calculation differentiates the
  deviance, working weights, and effective degrees of freedom through
  the converged PIRLS equations. The remaining controls below apply to
  Gaussian sparse fits. The distributional sparse solver accepts
  `crossprod_chunk_size`, `supernodal`, `optimizer_maxit`,
  `optimizer_gradient_tolerance`, `optimizer_trust_radius`,
  `inner_tolerance`, `inner_relative_tolerance`, `inner_maxit`,
  `inner_algorithm`, `inner_schur_rank`, `cores`, `gradient`,
  `gradient_probes`, `gradient_workers`, `score_batch_size`, and
  `qncv_batch_size`. Its `gradient` may be `"auto"` (the default),
  `"exact"`, or `"hybrid"`. Automatic selection uses exact scores for
  small joint systems and a deterministic matrix-free stochastic warm-up
  followed by exact-score refinement for large systems. QNCV resolves
  `"auto"` to `"exact"` and rejects `"hybrid"` because its analytic
  score is evaluated in streamed response batches.
  `inner_relative_tolerance` defaults to `1e-8` and accepts a
  conditional solve when both the actual and locally predicted objective
  reductions, divided by one plus the objective magnitude, fall below
  that threshold. This complements the absolute penalized-score
  `inner_tolerance` for large coefficient systems.
  `inner_algorithm="schur"` enables the experimental matrix-free Schur
  solver. It eliminates the larger location or scale coefficient block,
  preconditions the reduced system with the remaining sparse block
  factor, and certifies the reconstructed step against the full joint
  system. It augments that factor with a deterministic low-rank
  approximation to the Schur correction. `inner_schur_rank` controls the
  approximation rank and defaults to `6`. Hybrid warm-up stops early
  when rejected conditional fits exhaust a budget based on recent fit
  times, then refines from the best valid point. Exact refinement
  requests a score only after accepting an objective improvement. At
  exhausted trust-region recovery it compares the exact score with a
  directional finite difference to diagnose numerical stagnation.
  `gradient_probes` defaults to 64. Both generalized solvers warm-start
  each inner fit from the preceding valid outer evaluation and retry
  cold if it fails. They compute the exact score only when the outer
  optimizer requests a gradient, so objective-only line-search probes do
  not perform trace solves. `cores` is the total process and BLAS core
  budget. By default it is inferred from scheduler, container, affinity,
  and operating-system limits. Serial factorization phases use that
  budget for BLAS threads. Independent score or finite-difference tasks
  divide it between processes and BLAS threads without nested
  oversubscription. `score_workers` for the single-predictor generalized
  solver and `gradient_workers` for the distributional solver optionally
  override the number of exact-score processes or finite-difference
  prediction-error directions; exact generalized GCV gradients instead
  use one process and the full BLAS thread budget. The resolved worker
  value is capped by the core budget and available task count. By
  default, single-predictor score phases use the integer square root of
  the core budget as workers. Distributional exact scores use the full
  core budget when a memory estimate permits a shared full-design cache;
  workers divide likelihood-curvature directions and inverse columns
  while sharing read-only inputs through copy-on-write memory. Otherwise
  the likelihood derivatives retain the streamed serial path. Each
  process receives the remaining BLAS thread budget without nested
  oversubscription. Process parallelism falls back to one worker on
  Windows. Exact score directions are evaluated together so they share
  multi-right-hand-side solves. `score_batch_size` optionally caps that
  group size; by default the process uses its available-memory estimate
  to retain the largest safe shared batch. The inverse-column chunk
  width also grows with available memory, up to a fixed cap, to improve
  sparse-solve throughput. `qncv_batch_size` optionally caps the number
  of responses evaluated at once by distributional QNCV. Its default is
  calibrated from available process memory and the joint coefficient
  dimension. `gradient` entry may be `"auto"` (the default), `"finite"`,
  `"exact"`, `"stochastic"`, or `"hybrid"`; `gradient_probes` controls
  the fixed Rademacher trace probes used by the latter two methods;
  `"hybrid"` uses a warm start and then refines with exact scores. For
  Gaussian location-scale REML, `warmup` selects that warm start:
  `"efs"` runs up to `efs_maxit` (default `30`) extended Fellner-Schall
  updates (Wood and Fasiolo, 2017), each needing one inner fit and
  penalty traces but no outer score or line search, and sets smoothing
  parameters still moving toward a limit to the search bound before the
  exact phase; `"stochastic"` uses stochastic scores as before; `"none"`
  skips the warm start; and `"auto"` (the default) uses `"efs"` whenever
  `"hybrid"` would be selected and no warm start otherwise. The exact
  phase solves each inner fit to `exact_inner_relative_tolerance`
  (default the tighter of `1e-11` and `inner_relative_tolerance`), so
  that criterion noise from the inner fit stays below the trust-region
  improvements; `"auto"` compares predicted exact-score time and memory
  with the cost of finite differences after the first objective
  factorization; `gradient_workers` optionally overrides the
  automatically selected number of central finite-difference worker
  processes. The resolved value is capped by `cores` and the number of
  directions and falls back to one on Windows; `finite_difference_step`
  defaults to `1e-3`; sparse `"GCV.Cp"` uses the same cost-aware choice
  between exact and finite gradients. For Gaussian GCV, `"auto"` uses
  safeguarded Newton when the exact gradient fits the memory budget and
  there are at most six smoothing parameters, and L-BFGS-B otherwise.
  The `cdrgam.gcv_newton_max_parameters` option changes that threshold;
  `outer_optimizer` may be `"auto"` (the default), `"lbfgsb"`,
  `"bfgs_trust"`, or `"newton"`; `"bfgs_trust"` applies to REML and
  `"newton"` applies to Gaussian GCV. Explicitly requesting either
  curvature optimizer with `gradient="auto"` forces an exact gradient;
  `optimizer_maxit` defaults to an initial budget of `200`; when that
  default is exhausted, analytic curvature recovery may add up to two
  more 200-iteration blocks, while an explicitly supplied value is a
  hard limit; `hessian` selects the post-fit outer-Hessian calculation:
  `"auto"` (the default) uses the analytic Hessian when its estimated
  peak memory fits a conservative fraction of the process's available
  memory and otherwise uses `"gradient"`; `"gradient"` differences the
  exact REML score; `"analytic"` differentiates the Gaussian REML
  criterion using the fitted factor and a selected inverse, `"profiled"`
  retains the slower reference calculation that differences the profiled
  scalar criterion, `"optimhess"` directly differences the criterion,
  `"defer"` retains the sparse workspace and computes the gradient-based
  Hessian on the first unconditional-inference request, and `"none"`
  skips smoothing-parameter uncertainty; the automatic memory check
  honors the effective Linux cgroup or SLURM allocation when available.
  On other systems, automatic selection uses the gradient Hessian when
  available memory cannot be determined; `cdrgam.memory_limit_bytes` can
  provide an explicit process limit; sparse Gaussian `"GCV.Cp"` uses a
  post-fit Hessian as a local stability diagnostic. It differentiates
  the exact sparse gradient when that gradient fits memory, otherwise it
  differences the scalar criterion. GCV curvature is not reported as
  smoothing-parameter covariance; `hessian_step` controls the
  finite-difference step and defaults to `1e-2`; `supernodal` optionally
  overrides the automatically selected CHOLMOD factorization form;
  `trace_method` may be `"auto"` (the default), `"solve"`, or
  `"schur_inverse"`; the latter contracts penalties with block-local
  inverse entries and requires `schur="always"`; `trace_chunk_size`
  bounds exact-score solves; and `boundary_action` controls empirically
  inactive penalty subspaces: `"reduce"` (the default) removes certified
  zero-variance IRF terms, projects other certified IRF boundaries to
  their joint null space, and reoptimizes until no newly certified
  boundary remains; `"report"` records the same diagnostics without
  changing the fitted design; and `"error"` turns optimizer
  nonconvergence into an error. `boundary_log_sp` is the minimum log
  smoothing parameter considered for reduction and defaults to `12`.
  `boundary_edf_tolerance` is the maximum effective degrees of freedom
  that the penalized subspace may retain before reduction and defaults
  to `0.01`. This scale-invariant check prevents a numerically large
  smoothing parameter from removing active curvature. After the outer
  Hessian is computed, a Gaussian REML fit classifies a smoothing
  parameter as on the boundary when its absolute log value is at least
  `boundary_hessian_log_sp` (default `8`) and every entry of its
  outer-Hessian row, which is on the negative log restricted-likelihood
  scale, is at most `boundary_curvature_tolerance` (default `0.05`) in
  absolute value. The criterion is then flat along that parameter alone.
  A flat direction that mixes several smoothing parameters is not
  classified this way and still makes the Hessian indefinite. Boundary
  parameters are listed in `fit$sparse$convergence$hessian_boundary`,
  are excluded from the positive-definiteness check and from
  unconditional covariance, and do not by themselves trigger a warning.
  For a parameter at its upper boundary, a variance component estimated
  as zero, `boundary_profile = TRUE` (the default) records a profile
  bound: the smallest log smoothing parameter whose -2 REML criterion
  stays within `qchisq(2 * boundary_profile_level - 1, 1)` of the
  optimum, with `boundary_profile_level` defaulting to `0.95`. This
  one-sided level follows the 50:50 chi-squared mixture for a variance
  component on its boundary. The other smoothing parameters stay at
  their estimates during the search, so the bound can understate the
  largest supported variance when parameters trade off.
  `boundary_sensitivity = TRUE` also stores the fit at the bounds, which
  roughly doubles the stored factorization, so that
  [`estimate_effect()`](https://climblab.org/cdrgam/reference/estimate_effect.md)
  and
  [`test_effect()`](https://climblab.org/cdrgam/reference/test_effect.md)
  accept `boundary = "bound"`. For Gaussian fits, `final_newton = TRUE`
  (the default) refines the optimum after the outer optimizer stops. For
  REML each step uses the exact gradient and analytic outer Hessian and
  holds boundary smoothing parameters fixed. For GCV.Cp it differences
  the exact gradient, holds no parameter fixed, and scales
  `newton_tolerance` by the criterion divided by the number of
  observations. Either way the step replaces negative eigenvalues by
  their absolute values so the step still descends, follows the most
  negative eigenvector when the gradient is negligible but curvature is
  negative, and is accepted only if the criterion decreases. It stops
  when the predicted decrease falls below `newton_tolerance` (default
  `1e-3` on the -2 REML scale) with a positive-definite Hessian, or
  after `newton_max_steps` steps (default `5`). Each step costs about
  one analytic Hessian, so the phase is skipped with a warning when that
  Hessian exceeds the automatic memory budget unless
  `hessian = "analytic"`; the warning appears only when `final_newton`
  is set explicitly. The outcome is recorded in
  `fit$sparse$convergence$newton_refinement`. Reduced-model inference is
  conditional on the recorded data-selected boundary reduction. Every
  applied reduction records its source penalty, effective dimension, fit
  stage, and reduced-model objective check. `schur = "always"` enables
  the experimental response-group Schur solver instead of the default
  `"never"`. It batches block elimination, exact traces, and solves
  through BLAS; a threaded BLAS implementation improves large-core
  runtime. `crossprod_chunk_size` bounds the number of response rows
  used while accumulating sparse normal equations; `restarts` requests
  additional deterministic smoothing-parameter starts. Generalized
  sparse fits currently accept `crossprod_chunk_size`, `supernodal`,
  `optimizer_maxit`, `optimizer_gradient_tolerance`,
  `optimizer_trust_radius`, `hessian`, `hessian_step`, and the
  `boundary_curvature_tolerance`, `boundary_hessian_log_sp`,
  `boundary_profile`, and `boundary_profile_level` options described
  above; their exact score and safeguarded trust optimizer are
  mandatory. For generalized REML, `hessian = "gradient"` (the default)
  central- differences the exact score with step `hessian_step` (default
  `1e-2`) after optimization, which costs two score evaluations per
  smoothing parameter. The result supports boundary classification and
  profile bounds, the optimizer diagnostics in
  [`fit_diagnostics()`](https://climblab.org/cdrgam/reference/fit_diagnostics.md),
  and unconditional inference; `hessian = "none"` skips it. Generalized
  fits do not support `boundary = "bound"`.

- rank_action:

  Policy for scientifically meaningful nonidentifiability in a custom
  backend: `"error"` (the default), `"minimum_norm"`, `"drop"`, or
  `"penalize"`. Exact aliases among ordinary unpenalized parametric
  columns are always removed and reported. `"drop"` does not arbitrarily
  delete directions spanning smooth or IRF terms.

- rank_tol:

  Relative numerical rank tolerance. The default is
  `sqrt(.Machine$double.eps)`.

- rank_penalty:

  Relative fixed ridge used by `rank_action="penalize"`. The default is
  `1e-6`.

- drop.unused.levels:

  Drop factor levels without training observations. The default matches
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html). Retained
  levels remain part of the fitted factor vocabulary. For
  `cdrgam.fit()`, use `NULL` to inherit the choice already compiled into
  `design`; a supplied value must match it.

- ...:

  Additional arguments passed to
  [`mgcv::bam()`](https://rdrr.io/pkg/mgcv/man/bam.html) or
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html).

- formula:

  Extended model formula containing ordinary `mgcv` terms and optional
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) terms. An
  implicit `irf(1)` is added unless suppressed. For `family="gaulss"`,
  this may instead be a named list with `location` and `scale` formulas.
  The location formula supplies the response; the scale formula may be
  one-sided.

- impulses:

  Data frame with one row per impulse.

- responses:

  Data frame with one row per response.

- window:

  Default inclusive lag window inherited by every
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) term that does
  not define its own `window`. `NULL` retains the term-level default
  `c(0, Inf)`. Every impulse in the applicable series and window is
  included in the response-level design.

- knots_l:

  Optional model-level default for
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) lag-basis
  construction points. Term-level `knots_l`, including `NULL`, takes
  precedence.

- k_l, k_t, k_p:

  Optional model-level defaults for the corresponding
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) axis
  dimensions. An argument supplied by an individual term, including an
  explicit `NULL`, overrides its model-level default.

- bs_l, bs_t, bs_p:

  Optional model-level defaults for the corresponding
  [`irf()`](https://climblab.org/cdrgam/reference/irf.md) marginal basis
  names. Term-level arguments take precedence.

- series:

  Character vector of columns identifying independent series.

- impulse_time:

  Name of the impulse-time column.

- response_time:

  Name of the response-time column.

- history:

  One of `"auto"`, `"dense"`, or `"ragged"`.

- chunk_size:

  Maximum history links transformed at once in ragged mode, or responses
  transformed at once in dense mode.

- rescale_predictors:

  Divide continuous numeric predictors and internal lag coordinates by
  training-data standard deviations. Scaling is enabled by default; set
  `FALSE` to fit in native units. Binary, constant, categorical,
  grouping, and series variables are not changed. Lag and explicit time
  coordinates share one divisor, the standard deviation of the untiled
  impulse-time column, applied after history links are built in native
  time units.

- center_predictors:

  Center linear IRF predictors, including binary ones, at their means
  over the training history links. Every term is then defined relative
  to an average event: removing `irf(1)` removes the response to an
  average event, and omitting a lower-order term of an interaction
  removes that effect at the other predictors' means. The prepared
  design's `centering` table records each predictor's center. With
  centering, `irf(1)` is the response to an event whose linear
  predictors are at those means, and each linear term is the deviation
  from that response. Lower- order terms of a linear interaction are
  effects with the other predictors at their means. Effects and
  predictions still use source predictor values. `FALSE` restores
  uncentered predictors, where `irf(1)` is the response at zero.
  Individual terms can override this with `irf(center=)`.
  Response-aligned `by` covariates are never centered.

## Value

The native backend's fit also inherits from its selected mgcv fitting
class (`"bam"` or `"gam"`); custom backends provide corresponding
methods.

A fitted `cdrgam` object.
