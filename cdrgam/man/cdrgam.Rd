% Generated manually from cdrgam/R/formula.R and compressed.R.
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
\item{drop.unused.levels}{Drop factor levels without training observations.
The default matches \code{mgcv::gam()}. Retained levels remain part of the
fitted factor vocabulary. For \code{cdrgam.fit()}, use \code{NULL} to inherit
the choice already compiled into \code{design}; a supplied value must match
it.}
\item{family}{A standard family object or a family name accepted by
\code{cdrgam_family()}. \code{gaulss} rejects non-unit prior weights because
mgcv's family accepts but does not use them.}
\item{method}{Smoothing-parameter estimation method.}
\item{engine}{Either \code{"bam"} or \code{"gam"}.}
\item{backend}{Fitting backend. \code{"mgcv"} uses native fitting,
\code{"block"} selects the dense reference solver (Gaussian REML,
generalized LAML, or joint Gaussian location--scale LAML), and
\code{"sparse"} uses streamed sparse Gaussian, supported generalized, or
joint Gaussian location--scale optimization. The distributional sparse solver
uses an exact outer score. Select the Gaussian trust-region
optimizer with
\code{sparse_control=list(gradient="exact", outer_optimizer="bfgs_trust")}.}
\item{checkpoint}{Optional checkpoint path for a custom backend. Checkpoints
are written atomically after periodic REML evaluations and contain validated
optimization and restart state. The trust optimizer checkpoints its dense
curvature approximation, trust radius, counters, parameters, score, and
criterion so interruption resumes the same optimization trajectory; older
checkpoints without that state remain parameter-only warm starts. Reusing the
path resumes an interrupted fit or skips an already completed outer
optimization. Distributional fits do not yet support checkpoints.}
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
\code{optimizer_gradient_tolerance}, and \code{optimizer_trust_radius}; its
exact gradient and safeguarded trust-region optimizer are fixed. The
remaining controls below apply to Gaussian sparse fits. The distributional
sparse solver accepts \code{crossprod_chunk_size}, \code{supernodal},
\code{optimizer_maxit}, \code{optimizer_gradient_tolerance},
\code{inner_tolerance}, and \code{inner_maxit}; its exact score and L-BFGS-B
optimizer are fixed.
\code{gradient} entry may be \code{"auto"} (the default), \code{"finite"},
\code{"exact"}, \code{"stochastic"}, or \code{"hybrid"}.
\code{gradient_probes} controls the fixed Rademacher trace probes used by the
latter two methods; \code{"hybrid"} uses stochastic scores for a warm start
and then refines with exact scores. \code{"auto"} compares predicted
exact-score time and memory with the cost of finite differences after the
first objective factorization. \code{gradient_cores} evaluates central
finite-difference directions in parallel on systems that support process
forking and falls back to one core on Windows; a single-threaded BLAS should
be used when it exceeds one.
\code{finite_difference_step} defaults to \code{1e-3}.
\code{outer_optimizer} may be \code{"auto"} (the default), \code{"lbfgsb"},
or the experimental safeguarded \code{"bfgs_trust"}. Automatic selection uses
trust-region BFGS for exact gradients and L-BFGS-B otherwise. Explicitly
requesting \code{"bfgs_trust"} with \code{gradient="auto"} forces an exact
gradient. \code{optimizer_maxit} defaults to \code{200}.
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
provide an explicit process limit. \code{hessian_step} controls the
finite-difference step and defaults
to \code{1e-2}.
\code{supernodal} optionally overrides the automatically selected CHOLMOD
factorization form. \code{trace_method} may be \code{"auto"} (the default),
\code{"solve"}, or \code{"schur_inverse"}; the latter contracts penalties
with block-local inverse entries and requires \code{schur="always"}.
\code{trace_chunk_size} bounds exact-score solves. \code{boundary_action}
controls empirically inactive
penalty subspaces: \code{"report"} (the default) records them,
\code{"reduce"} projects certified IRF curvature boundaries to their joint
null space and reoptimizes, and \code{"error"} turns optimizer nonconvergence
into an error. \code{boundary_log_sp} is the minimum log smoothing parameter
considered for reduction and defaults to \code{12}. Reduced-model inference
is conditional on the recorded data-selected boundary reduction. Full-rank
boundaries that would remove an entire term are reported but require user
confirmation.
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
