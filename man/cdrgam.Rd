% Generated manually from R/formula.R and R/compressed.R.
\name{cdrgam}
\alias{cdrgam}
\alias{cdrgam.fit}
\title{Fit a continuous-time deconvolutional GAM}
\usage{
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
}
\arguments{
\item{formula}{An extended formula containing ordinary \pkg{mgcv} terms and
optional \code{irf()} terms. An implicit
\code{irf(1)} is added unless suppressed with \code{- irf(1)}. For
\code{family="gaulss"}, this may instead be a named list with \code{location}
and \code{scale} formulas. The location formula supplies the response; the
scale formula may be one-sided.}
\item{impulses}{Data frame with one row per impulse.}
\item{responses}{Data frame with one row per response.}
\item{window}{Default inclusive lag window inherited by each \code{irf()}
term that does not supply its own window. \code{NULL} retains the term-level
default \code{c(0, Inf)}. Every impulse in the applicable series and window is
included in the response-level design.}
\item{knots_l}{Optional model-level default for \code{irf()} lag-basis
construction points. A term-level value, including \code{NULL}, takes
precedence.}
\item{k_l, k_t, k_p}{Optional model-level defaults for the corresponding
\code{irf()} axis dimensions. An argument supplied by an individual term,
including an explicit \code{NULL}, overrides its model-level default.}
\item{bs_l, bs_t, bs_p}{Optional model-level defaults for the corresponding
\code{irf()} marginal basis names. Term-level arguments take precedence.}
\item{series}{Character vector naming independent-series columns.}
\item{impulse_time}{Name of the impulse-time column.}
\item{response_time}{Name of the response-time column.}
\item{history}{History layout: \code{"auto"}, \code{"dense"}, or
\code{"ragged"}.}
\item{chunk_size}{Maximum responses or impulse-response links transformed in
one chunk.}
\item{rescale_predictors}{If \code{TRUE}, divide continuous numeric
predictors by their training-data standard deviations before model
construction. No centering is performed. Binary, constant, categorical,
grouping, and series variables are unchanged. Lag and explicit time
predictors share the standard deviation of the untiled impulse-time column;
history windows are still constructed in native time units. Stored scaling
metadata are reused for prediction, while CDR effect estimates and plots are
reported in native source units.}
\item{center_predictors}{If \code{TRUE}, the default, center linear IRF
predictors, including binary ones, at their means over the training history
links. Every term is then defined relative to an average event: removing
\code{irf(1)} removes the response to an average event, and omitting a
lower-order term of an interaction removes that effect at the other
predictors' means. The prepared design's \code{centering} table records each
predictor's center. \code{irf(1)} is then the response to an event whose linear predictors
are at those means, each linear term is the deviation from that response, and
lower-order terms of a linear interaction are effects with the other
predictors at their means. Effects and predictions still take source
predictor values. \code{FALSE} leaves linear predictors uncentered, so
\code{irf(1)} is the response at zero. Terms can override the setting with
\code{irf(center = )}. Response-aligned \code{by} covariates are never
centered.}
\item{drop.unused.levels}{Drop factor levels without training observations.
The default matches \code{mgcv::gam()}. Retained levels remain part of the
fitted factor vocabulary. For \code{cdrgam.fit()}, use \code{NULL} to inherit
the choice already compiled into \code{design}; a supplied value must match
it.}
\item{family}{A standard family object or a family name accepted by
\code{cdrgam_family()}. \code{gaulss} rejects non-unit prior weights because
mgcv's family accepts but does not use them.}
\item{method}{Smoothing-parameter estimation method. The sparse Gaussian
backend accepts REML (the default, also spelled \code{"fREML"}) and
\code{"GCV.Cp"}. Sparse GCV uses the exact effective degrees of freedom in the
GCV score and accepts mgcv's positive \code{gamma} multiplier through
\code{...}. Its initial implementation uses L-BFGS-B with a memory-gated exact
gradient and a parallel finite-difference fallback. It omits
smoothing-parameter uncertainty and does not support automatic boundary
reduction or rank regularization. The single-predictor generalized sparse
solver also accepts \code{"GCV.Cp"}: it minimizes UBRE for binomial and
Poisson fits and GCV for Gamma fits using streamed PIRLS and parallel finite
differences. The sparse \code{gaulss} solver accepts REML (the default) and
\code{"QNCV"}. QNCV uses a streamed quadratic leave-one-response-out
approximation and an exact analytic smoothing-parameter gradient. It accepts
mgcv's positive \code{gamma} multiplier through \code{...}; custom
neighborhoods are not yet supported.}
\item{engine}{Either \code{"bam"} or \code{"gam"}.}
\item{backend}{Fitting backend. \code{"mgcv"} uses native fitting,
\code{"block"} selects the dense reference solver (Gaussian REML,
generalized LAML, or joint Gaussian location--scale LAML), and
\code{"sparse"} uses streamed sparse Gaussian, supported generalized, or
joint Gaussian location--scale optimization. Under REML, the distributional
sparse solver ends with exact outer scores; large systems use matrix-free
stochastic scores for a warm-up by default, then use safeguarded trust-region
BFGS for exact refinement. Distributional QNCV uses exact gradients
throughout. Select the Gaussian trust-region
optimizer with
\code{sparse_control=list(gradient="exact", outer_optimizer="bfgs_trust")}.}
\item{checkpoint}{Optional checkpoint path for a custom backend. Checkpoints
are written atomically after periodic smoothing-criterion evaluations and
contain validated
optimization and restart state. The trust optimizer checkpoints its dense
curvature approximation, trust radius, counters, parameters, score, and
criterion so interruption resumes the same optimization trajectory; older
checkpoints without that state remain parameter-only warm starts. Reusing the
path resumes an interrupted fit or skips an already completed outer
optimization. Sparse distributional fits checkpoint accepted outer states and
resume from the best retained smoothing parameters; the joint factorization is
reconstructed once on resume.}
\item{solver_trace}{Custom-backend progress reporting. \code{FALSE} or
\code{0} is silent; \code{TRUE} or \code{1} reports phases, improving
solutions, and trust-optimizer diagnostics; \code{2} reports every objective
evaluation; and \code{3} adds low-level chunk and factorization events. Trust
diagnostics include projected-gradient and tolerance ratios, step norms,
trust radius, predicted and actual improvement, acceptance statistics, and
the expected post-fit Hessian workload. At the numerical trust-radius floor,
the sparse trust optimizer uses the resolved Hessian strategy to either
certify practical convergence or perform a bounded positive-curvature reset.
Automatic selection uses the same process memory budget as the post-fit
Hessian calculation.
A function receives the progress events as named lists.}
\item{sparse_control}{Named control list for the sparse backend. The
single-predictor generalized sparse solver accepts
\code{crossprod_chunk_size}, \code{supernodal}, \code{optimizer_maxit},
\code{optimizer_gradient_tolerance}, \code{optimizer_trust_radius},
\code{cores}, \code{gradient}, \code{score_workers},
\code{score_batch_size}, \code{finite_difference_step}, \code{hessian},
\code{hessian_step}, and the boundary classification options
\code{boundary_curvature_tolerance}, \code{boundary_hessian_log_sp},
\code{boundary_profile}, and \code{boundary_profile_level}. Generalized REML
uses an exact gradient and safeguarded trust-region optimization. After
optimization, \code{hessian = "gradient"} (the default) central-differences
the exact score with step \code{hessian_step} (default \code{1e-2}), two
score evaluations per outer parameter, to obtain the outer Hessian used for
boundary classification, profile bounds, optimizer diagnostics, and
unconditional inference; \code{hessian = "none"} skips it. Generalized
\code{"GCV.Cp"} uses L-BFGS-B and accepts \code{gradient="auto"} (the
default), \code{"exact"}, or \code{"finite"}. Automatic selection compares
the memory and structural cost of a selected-inverse analytic gradient with
two conditional PIRLS fits per smoothing parameter. The exact calculation
differentiates the deviance, working weights, and effective degrees of freedom
through the converged PIRLS equations. The
remaining controls below apply to Gaussian sparse fits. The distributional
sparse solver accepts \code{crossprod_chunk_size}, \code{supernodal},
\code{optimizer_maxit}, \code{optimizer_gradient_tolerance},
\code{inner_tolerance}, \code{inner_relative_tolerance},
\code{inner_maxit}, \code{inner_algorithm}, \code{inner_schur_rank},
\code{cores}, \code{gradient},
\code{gradient_probes}, \code{gradient_workers}, \code{score_batch_size}, and
\code{qncv_batch_size}. Its \code{gradient} may be \code{"auto"} (the
default), \code{"exact"}, or \code{"hybrid"}. Automatic selection uses exact
scores for small joint systems and a deterministic matrix-free stochastic
warm-up followed by exact-score refinement for large systems.
QNCV resolves \code{"auto"} to \code{"exact"} and rejects \code{"hybrid"}
because its analytic score is evaluated in streamed response batches.
\code{inner_relative_tolerance} defaults to \code{1e-8} and accepts a
conditional solve when both the actual and locally predicted objective
reductions, divided by one plus the objective magnitude, fall below that
threshold. This complements the absolute penalized-score
\code{inner_tolerance} for large coefficient systems.
\code{inner_algorithm="schur"} enables the experimental matrix-free Schur
solver. It eliminates the larger location or scale coefficient block,
preconditions the reduced system with the remaining sparse block factor, and
certifies the reconstructed step against the full joint system.
It augments that factor with a deterministic low-rank approximation to the
Schur correction. \code{inner_schur_rank} controls the approximation rank and
defaults to \code{6}.
Hybrid warm-up stops early when rejected conditional fits exhaust a budget
based on recent fit times, then refines from the best valid point.
\code{gradient_probes} defaults to 64. Both generalized solvers warm-start each inner fit from
the preceding valid outer evaluation and retry cold if it fails.
They compute the exact score only when the outer optimizer requests a gradient,
so objective-only line-search probes do not perform trace solves.
\code{cores} is the total process and BLAS core budget. By default it is
inferred from scheduler, container, affinity, and operating-system limits.
Serial factorization phases use that budget for BLAS threads. Independent
score or finite-difference tasks divide it between processes and BLAS threads
without nested oversubscription. \code{score_workers} for the single-predictor
generalized solver and \code{gradient_workers} for the distributional solver
optionally override the number of exact-score processes or generalized
finite-difference prediction-error directions. Exact generalized GCV
gradients instead use one process and the full BLAS thread budget. The
resolved worker value is capped by the core budget and available task count.
By default, single-predictor score phases use the
integer square root of the core budget as workers. Distributional exact
scores use the full core budget when a memory estimate permits a shared
full-design cache; workers divide
likelihood-curvature directions and inverse columns while sharing read-only
inputs through copy-on-write memory. Otherwise the likelihood derivatives
retain the streamed serial path. Each process receives the remaining BLAS
thread budget without nested oversubscription. Process parallelism falls back
to one worker on Windows. Exact score directions are evaluated together so
they share multi-right-hand-side solves.
\code{score_batch_size} optionally caps that group size; by default the process
uses its available-memory estimate to retain the largest safe shared batch.
The inverse-column chunk width also grows with available memory, up to a fixed
cap, to improve sparse-solve throughput.
\code{qncv_batch_size} optionally caps the number of responses evaluated at
once by distributional QNCV. Its default is calibrated from available process
memory and the joint coefficient dimension.
\code{gradient} entry may be \code{"auto"} (the default), \code{"finite"},
\code{"exact"}, \code{"stochastic"}, or \code{"hybrid"}.
\code{gradient_probes} controls the fixed Rademacher trace probes used by the
latter two methods; \code{"hybrid"} uses stochastic scores for a warm start
and then refines with exact scores. \code{"auto"} compares predicted
exact-score time and memory with the cost of finite differences after the
first objective factorization. \code{gradient_workers} optionally overrides the
automatically selected number of central finite-difference worker processes.
The resolved value is capped by \code{cores} and the number of directions and
falls back to one on Windows.
\code{finite_difference_step} defaults to \code{1e-3}.
Sparse \code{"GCV.Cp"} uses the same cost-aware choice between exact and finite
gradients. For Gaussian GCV, \code{"auto"} uses safeguarded Newton when the
exact gradient fits the memory budget and there are at most six smoothing
parameters, and L-BFGS-B otherwise. The
\code{cdrgam.gcv_newton_max_parameters} option changes that threshold.
\code{outer_optimizer} may be \code{"auto"} (the default), \code{"lbfgsb"},
the safeguarded \code{"bfgs_trust"}, or \code{"newton"}.
\code{"bfgs_trust"} applies to REML and \code{"newton"} applies to Gaussian
GCV. Explicitly requesting either curvature optimizer with
\code{gradient="auto"} forces an exact gradient.
\code{optimizer_maxit} defaults to an initial budget of \code{200}.
When that default is exhausted, analytic curvature recovery may add up to two
more 200-iteration blocks. An explicitly supplied value is a hard limit.
\code{optimizer_gradient_tolerance} and \code{optimizer_trust_radius}
control its projected-gradient tolerance and initial trust radius.
\code{hessian} selects the post-fit outer-Hessian calculation:
\code{"auto"} (the default) uses the analytic Hessian when its estimated peak
memory fits a conservative fraction of the process's available memory and
otherwise uses \code{"gradient"}. \code{"gradient"} differences the exact
REML score, \code{"analytic"} differentiates the Gaussian REML criterion
using the fitted factor and a selected inverse, \code{"profiled"} retains the slower reference
calculation that differences the profiled scalar criterion,
\code{"optimhess"} directly differences the unprofiled criterion,
\code{"defer"} retains the sparse workspace and computes
the gradient-based Hessian on the first unconditional-inference request, and
\code{"none"} skips smoothing-parameter uncertainty. The automatic memory
check honors the effective Linux cgroup or SLURM allocation when available.
On other systems, automatic selection uses the gradient Hessian when available
memory cannot be determined; the \code{cdrgam.memory_limit_bytes} option can
provide an explicit process limit. Sparse Gaussian \code{"GCV.Cp"} uses a
post-fit Hessian as a local stability diagnostic. It differentiates the exact
sparse gradient when that gradient fits memory, otherwise it differences the
scalar criterion. GCV curvature is not reported as smoothing-parameter
covariance. \code{hessian_step} controls the
finite-difference step and defaults
to \code{1e-2}.
\code{supernodal} optionally overrides the automatically selected CHOLMOD
factorization form. \code{trace_method} may be \code{"auto"} (the default),
\code{"solve"}, or \code{"schur_inverse"}; the latter contracts penalties
with block-local inverse entries and requires \code{schur="always"}.
\code{trace_chunk_size} bounds exact-score solves. \code{boundary_action}
controls empirically inactive
penalty subspaces: \code{"reduce"} (the default) removes certified
zero-variance IRF terms, projects other certified IRF boundaries to their
joint null space, and reoptimizes until no newly certified boundary remains;
\code{"report"} records the same diagnostics without changing the fitted
design; and
\code{"error"} turns optimizer nonconvergence
into an error. \code{boundary_log_sp} is the minimum log smoothing parameter
considered for reduction and defaults to \code{12}.
\code{boundary_edf_tolerance} is the maximum effective degrees of freedom that
the penalized subspace may retain before reduction and defaults to \code{0.01}.
This scale-invariant check prevents a numerically large smoothing parameter
from removing active curvature.
After the outer Hessian is computed, a Gaussian REML fit classifies a
smoothing parameter as on the boundary when its absolute log value is at
least \code{boundary_hessian_log_sp} (default \code{8}) and every entry of its
outer-Hessian row, which is on the negative log restricted-likelihood scale,
is at most \code{boundary_curvature_tolerance} (default \code{0.05}) in
absolute value. The criterion is then flat along that parameter alone. A flat
direction that mixes several smoothing parameters is not classified this way
and still makes the Hessian indefinite. Boundary parameters are listed in
\code{fit$sparse$convergence$hessian_boundary}, are excluded from the
positive-definiteness check and from unconditional covariance, and do not by
themselves trigger a warning. For a parameter at its upper boundary, a
variance component estimated as zero, \code{boundary_profile = TRUE} (the
default) records a profile bound: the smallest log smoothing parameter whose
-2 REML criterion stays within \code{qchisq(2 * boundary_profile_level - 1, 1)}
of the optimum, with \code{boundary_profile_level} defaulting to \code{0.95}.
This one-sided level follows the 50:50 chi-squared mixture for a variance
component on its boundary. The other smoothing parameters stay at their
estimates during the search, so the bound can understate the largest
supported variance when parameters trade off. \code{boundary_sensitivity =
TRUE} also stores the fit at the bounds, which roughly doubles the stored
factorization, so that \code{\link{estimate_effect}} and
\code{\link{test_effect}} accept \code{boundary = "bound"}.
For Gaussian REML, \code{final_newton = TRUE} (the default) refines the
optimum after the outer optimizer stops. Each step uses the exact gradient and analytic outer
Hessian, holds boundary smoothing parameters fixed, replaces negative
eigenvalues by their absolute values so the step still descends, follows the
most negative eigenvector when the gradient is negligible but curvature is
negative, and is accepted only if the criterion decreases. It stops when the
predicted decrease falls below \code{newton_tolerance} (default \code{1e-3}
on the -2 REML scale) with a positive-definite Hessian, or after
\code{newton_max_steps} steps (default \code{5}). Each step costs about one
analytic Hessian, so the phase is skipped with a warning when that Hessian
exceeds the automatic memory budget unless \code{hessian = "analytic"}; the
warning appears only when \code{final_newton} is set explicitly. The
outcome is recorded in \code{fit$sparse$convergence$newton_refinement}. Reduced-model inference is conditional on the
recorded data-selected boundary reduction. Every applied reduction records
its source penalty, effective dimension, fit stage, and reduced-model
objective check.
\code{schur = "always"} enables the experimental
response-group Schur solver instead of the default \code{"never"}. It batches
block elimination, exact traces, and solves through BLAS; a threaded BLAS
implementation improves large-core runtime.
\code{crossprod_chunk_size} bounds the number of response rows used while
accumulating sparse normal equations, and \code{restarts} requests additional
deterministic smoothing-parameter starts.}
\item{rank_action}{Policy for nonidentifiability in a custom backend:
\code{"error"} (the default), \code{"minimum_norm"}, \code{"drop"}, or
\code{"penalize"}. Exact aliases among ordinary unpenalized parametric columns
are automatically removed and reported. \code{"drop"} will not arbitrarily
remove directions spanning smooth or IRF terms.}
\item{rank_tol}{Relative numerical rank tolerance; defaults to
\code{sqrt(.Machine$double.eps)}.}
\item{rank_penalty}{Relative fixed ridge used by
\code{rank_action="penalize"}; defaults to \code{1e-6}.}
\item{...}{Additional arguments passed to \code{mgcv::bam()} or
\code{mgcv::gam()}.}
\item{design}{A reusable \code{cdrgam_design} returned by
\code{prepare_cdrgam()}.}
}
\value{
A \code{cdrgam} object. The native backend also inherits from its selected
\pkg{mgcv} fitting class; custom backends provide corresponding methods.
}
\description{
\code{cdrgam()} compiles untiled impulse and response streams and fits the
expanded model. \code{cdrgam.fit()} fits an already prepared
\code{cdrgam_design} without rebuilding its stream histories.
}
