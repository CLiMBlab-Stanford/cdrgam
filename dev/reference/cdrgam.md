# Fit a continuous-time deconvolutional GAM

`cdrgam()` compiles untiled impulse and response streams and fits the
expanded model. `cdrgam.fit()` fits an already prepared `cdrgam_design`
without rebuilding its stream histories.

## Usage

``` r
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
  rescale_predictors = FALSE,
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
```

## Arguments

- formula:

  An extended formula containing ordinary mgcv terms and optional
  [`irf()`](https://climblab.org/cdrgam/dev/reference/irf.md) terms. An
  implicit `irf(1)` is added unless suppressed with `- irf(1)`. For
  `family="gaulss"`, this may instead be a named list with `location`
  and `scale` formulas. The location formula supplies the response; the
  scale formula may be one-sided.

- impulses:

  Data frame with one row per impulse.

- responses:

  Data frame with one row per response.

- window:

  Default inclusive lag window inherited by each
  [`irf()`](https://climblab.org/cdrgam/dev/reference/irf.md) term that
  does not supply its own window. `NULL` retains the term-level default
  `c(0, Inf)`. Every impulse in the applicable series and window is
  included in the response-level design.

- knots_l:

  Optional model-level default for
  [`irf()`](https://climblab.org/cdrgam/dev/reference/irf.md) lag-basis
  construction points. A term-level value, including `NULL`, takes
  precedence.

- k_l, k_t, k_p:

  Optional model-level defaults for the corresponding
  [`irf()`](https://climblab.org/cdrgam/dev/reference/irf.md) axis
  dimensions. An argument supplied by an individual term, including an
  explicit `NULL`, overrides its model-level default.

- bs_l, bs_t, bs_p:

  Optional model-level defaults for the corresponding
  [`irf()`](https://climblab.org/cdrgam/dev/reference/irf.md) marginal
  basis names. Term-level arguments take precedence.

- series:

  Character vector naming independent-series columns.

- impulse_time:

  Name of the impulse-time column.

- response_time:

  Name of the response-time column.

- history:

  History layout: `"auto"`, `"dense"`, or `"ragged"`.

- chunk_size:

  Maximum responses or impulse-response links transformed in one chunk.

- rescale_predictors:

  If `TRUE`, divide continuous numeric predictors by their training-data
  standard deviations before model construction. No centering is
  performed. Binary, constant, categorical, grouping, and series
  variables are unchanged. Lag and explicit time predictors share the
  standard deviation of the untiled impulse-time column; history windows
  are still constructed in native time units. Stored scaling metadata
  are reused for prediction, while CDR effect estimates and plots are
  reported in native source units.

- drop.unused.levels:

  Drop factor levels without training observations. The default matches
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html). Retained
  levels remain part of the fitted factor vocabulary. For
  `cdrgam.fit()`, use `NULL` to inherit the choice already compiled into
  `design`; a supplied value must match it.

- family:

  A standard family object or a family name accepted by
  [`cdrgam_family()`](https://climblab.org/cdrgam/dev/reference/cdrgam_family.md).
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

  Either `"bam"` or `"gam"`.

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
  validated optimization and restart state. The trust optimizer
  checkpoints its dense curvature approximation, trust radius, counters,
  parameters, score, and criterion so interruption resumes the same
  optimization trajectory; older checkpoints without that state remain
  parameter-only warm starts. Reusing the path resumes an interrupted
  fit or skips an already completed outer optimization. Sparse
  distributional fits checkpoint accepted outer states and resume from
  the best retained smoothing parameters; the joint factorization is
  reconstructed once on resume.

- solver_trace:

  Custom-backend progress reporting. `FALSE` or `0` is silent; `TRUE` or
  `1` reports phases, improving solutions, and trust-optimizer
  diagnostics; `2` reports every objective evaluation; and `3` adds
  low-level chunk and factorization events. Trust diagnostics include
  projected-gradient and tolerance ratios, step norms, trust radius,
  predicted and actual improvement, acceptance statistics, and the
  expected post-fit Hessian workload. At the numerical trust-radius
  floor, the sparse trust optimizer uses the resolved Hessian strategy
  to either certify practical convergence or perform a bounded
  positive-curvature reset. Automatic selection uses the same process
  memory budget as the post-fit Hessian calculation. A function receives
  the progress events as named lists.

- sparse_control:

  Named control list for the sparse backend. The single-predictor
  generalized sparse solver accepts `crossprod_chunk_size`,
  `supernodal`, `optimizer_maxit`, `optimizer_gradient_tolerance`,
  `optimizer_trust_radius`, `cores`, `gradient`, `score_workers`,
  `score_batch_size`, and `finite_difference_step`. Generalized REML
  uses an exact gradient and safeguarded trust-region optimization.
  Generalized `"GCV.Cp"` uses L-BFGS-B and accepts `gradient="auto"`
  (the default), `"exact"`, or `"finite"`. Automatic selection compares
  the memory and structural cost of a selected-inverse analytic gradient
  with two conditional PIRLS fits per smoothing parameter. The exact
  calculation differentiates the deviance, working weights, and
  effective degrees of freedom through the converged PIRLS equations.
  The remaining controls below apply to Gaussian sparse fits. The
  distributional sparse solver accepts `crossprod_chunk_size`,
  `supernodal`, `optimizer_maxit`, `optimizer_gradient_tolerance`,
  `inner_tolerance`, `inner_maxit`, `cores`, `gradient`,
  `gradient_probes`, `gradient_workers`, `score_batch_size`, and
  `qncv_batch_size`. Its `gradient` may be `"auto"` (the default),
  `"exact"`, or `"hybrid"`. Automatic selection uses exact scores for
  small joint systems and a deterministic matrix-free stochastic warm-up
  followed by exact-score refinement for large systems. QNCV resolves
  `"auto"` to `"exact"` and rejects `"hybrid"` because its analytic
  score is evaluated in streamed response batches. Hybrid warm-up stops
  early when rejected conditional fits exhaust a budget based on recent
  fit times, then refines from the best valid point. `gradient_probes`
  defaults to 64. Both generalized solvers warm-start each inner fit
  from the preceding valid outer evaluation and retry cold if it fails.
  They compute the exact score only when the outer optimizer requests a
  gradient, so objective-only line-search probes do not perform trace
  solves. `cores` is the total process and BLAS core budget. By default
  it is inferred from scheduler, container, affinity, and
  operating-system limits. Serial factorization phases use that budget
  for BLAS threads. Independent score or finite-difference tasks divide
  it between processes and BLAS threads without nested oversubscription.
  `score_workers` for the single-predictor generalized solver and
  `gradient_workers` for the distributional solver optionally override
  the number of exact-score processes or generalized finite-difference
  prediction-error directions. Exact generalized GCV gradients instead
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
  `"exact"`, `"stochastic"`, or `"hybrid"`. `gradient_probes` controls
  the fixed Rademacher trace probes used by the latter two methods;
  `"hybrid"` uses stochastic scores for a warm start and then refines
  with exact scores. `"auto"` compares predicted exact-score time and
  memory with the cost of finite differences after the first objective
  factorization. `gradient_workers` optionally overrides the
  automatically selected number of central finite-difference worker
  processes. The resolved value is capped by `cores` and the number of
  directions and falls back to one on Windows. `finite_difference_step`
  defaults to `1e-3`. Sparse `"GCV.Cp"` uses the same cost-aware choice
  between exact and finite gradients but always uses L-BFGS-B.
  `outer_optimizer` may be `"auto"` (the default), `"lbfgsb"`, or the
  experimental safeguarded `"bfgs_trust"`. Automatic selection uses
  trust-region BFGS for exact gradients and L-BFGS-B otherwise.
  Explicitly requesting `"bfgs_trust"` with `gradient="auto"` forces an
  exact gradient. `optimizer_maxit` defaults to an initial budget of
  `200`. When that default is exhausted, analytic curvature recovery may
  add up to two more 200-iteration blocks. An explicitly supplied value
  is a hard limit. `optimizer_gradient_tolerance` and
  `optimizer_trust_radius` control its projected-gradient tolerance and
  initial trust radius. `hessian` selects the post-fit outer-Hessian
  calculation: `"auto"` (the default) uses the analytic Hessian when its
  estimated peak memory fits a conservative fraction of the process's
  available memory and otherwise uses `"gradient"`. `"gradient"`
  differences the exact REML score, `"analytic"` differentiates the
  Gaussian REML criterion using the fitted factor and a selected
  inverse, `"profiled"` retains the slower reference calculation that
  differences the profiled scalar criterion, `"optimhess"` directly
  differences the unprofiled criterion, `"defer"` retains the sparse
  workspace and computes the gradient-based Hessian on the first
  unconditional-inference request, and `"none"` skips
  smoothing-parameter uncertainty. The automatic memory check honors the
  effective Linux cgroup or SLURM allocation when available. On other
  systems, automatic selection uses the gradient Hessian when available
  memory cannot be determined; the `cdrgam.memory_limit_bytes` option
  can provide an explicit process limit. Sparse `"GCV.Cp"` resolves
  `hessian="auto"` to `"none"` because smoothing-parameter uncertainty
  is not yet implemented. `hessian_step` controls the finite-difference
  step and defaults to `1e-2`. `supernodal` optionally overrides the
  automatically selected CHOLMOD factorization form. `trace_method` may
  be `"auto"` (the default), `"solve"`, or `"schur_inverse"`; the latter
  contracts penalties with block-local inverse entries and requires
  `schur="always"`. `trace_chunk_size` bounds exact-score solves.
  `boundary_action` controls empirically inactive penalty subspaces:
  `"report"` (the default) records them, `"reduce"` projects certified
  IRF curvature boundaries to their joint null space and reoptimizes
  until no newly certified boundary remains, and `"error"` turns
  optimizer nonconvergence into an error. `boundary_log_sp` is the
  minimum log smoothing parameter considered for reduction and defaults
  to `12`. `boundary_edf_tolerance` is the maximum effective degrees of
  freedom that the penalized subspace may retain before reduction and
  defaults to `0.01`. This scale-invariant check prevents a numerically
  large smoothing parameter from removing active curvature.
  Reduced-model inference is conditional on the recorded data-selected
  boundary reduction. Full-rank boundaries that would remove an entire
  term are reported but require user confirmation. `schur = "always"`
  enables the experimental response-group Schur solver instead of the
  default `"never"`. It batches block elimination, exact traces, and
  solves through BLAS; a threaded BLAS implementation improves
  large-core runtime. `crossprod_chunk_size` bounds the number of
  response rows used while accumulating sparse normal equations, and
  `restarts` requests additional deterministic smoothing-parameter
  starts.

- rank_action:

  Policy for nonidentifiability in a custom backend: `"error"` (the
  default), `"minimum_norm"`, `"drop"`, or `"penalize"`. Exact aliases
  among ordinary unpenalized parametric columns are automatically
  removed and reported. `"drop"` will not arbitrarily remove directions
  spanning smooth or IRF terms.

- rank_tol:

  Relative numerical rank tolerance; defaults to
  `sqrt(.Machine$double.eps)`.

- rank_penalty:

  Relative fixed ridge used by `rank_action="penalize"`; defaults to
  `1e-6`.

- ...:

  Additional arguments passed to
  [`mgcv::bam()`](https://rdrr.io/pkg/mgcv/man/bam.html) or
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html).

- design:

  A reusable `cdrgam_design` returned by
  [`prepare_cdrgam()`](https://climblab.org/cdrgam/dev/reference/prepare_cdrgam.md).

## Value

A `cdrgam` object. The native backend also inherits from its selected
mgcv fitting class; custom backends provide corresponding methods.
