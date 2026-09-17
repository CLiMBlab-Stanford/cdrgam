% Generated manually from cdrgam/R/formula.R and compressed.R.
\name{fit_cdrgam}
\alias{fit_cdrgam}
\alias{fit_compressed_cdr_gam}
\title{Fit a continuous-time deconvolutional GAM}
\usage{
fit_cdrgam(
  formula,
  impulses = NULL,
  responses = NULL,
  series = character(),
  impulse_time = "time",
  response_time = "time",
  history = c("auto", "dense", "ragged"),
  history_length = Inf,
  chunk_size = 10000,
  family = stats::gaussian(),
  method = NULL,
  engine = c("bam", "gam"),
  backend = c("mgcv", "block", "sparse", "sparse_trust"),
  checkpoint = NULL,
  solver_trace = FALSE,
  sparse_control = list(),
  rank_action = c("error", "minimum_norm", "drop", "penalize"),
  rank_tol = NULL,
  rank_penalty = NULL,
  ...
)

fit_compressed_cdr_gam(
  y,
  terms,
  family = stats::gaussian(),
  method = NULL,
  engine = c("bam", "gam"),
  intercept = TRUE,
  ...
)
}
\arguments{
\item{formula}{An extended formula containing ordinary \pkg{mgcv} terms and
optional \code{irf()} terms, or a reusable \code{cdrgam_design}. An implicit
\code{irf(1)} is added unless suppressed with \code{- irf(1)}.}
\item{impulses}{Data frame with one row per impulse.}
\item{responses}{Data frame with one row per response.}
\item{series}{Character vector naming independent-series columns.}
\item{impulse_time}{Name of the impulse-time column.}
\item{response_time}{Name of the response-time column.}
\item{history}{History layout: \code{"auto"}, \code{"dense"}, or
\code{"ragged"}.}
\item{history_length}{Maximum number of most-recent impulses retained within
each IRF lag window; \code{Inf} retains all impulses.}
\item{chunk_size}{Maximum responses or impulse-response links transformed in
one chunk.}
\item{family}{An \pkg{mgcv} family.}
\item{method}{Smoothing-parameter estimation method.}
\item{engine}{Either \code{"bam"} or \code{"gam"}.}
\item{backend}{Fitting backend. \code{"mgcv"} uses native fitting,
\code{"block"} selects the dense Gaussian REML reference solver, and
\code{"sparse"} uses sparse penalized normal equations.
\code{"sparse_trust"} selects the sparse engine with its exact REML score
and safeguarded trust-region BFGS outer optimizer.}
\item{checkpoint}{Optional checkpoint path for a custom backend. Checkpoints
are written atomically after periodic REML evaluations and contain
validated optimization and restart state. Reusing the path resumes an
interrupted fit or skips an already completed outer optimization.}
\item{solver_trace}{Custom-backend progress reporting. \code{FALSE} or
\code{0} is silent; \code{TRUE} or \code{1} reports phases, improving
solutions, and trust-optimizer diagnostics; \code{2} reports every objective
evaluation; and \code{3} adds low-level chunk and factorization events. Trust
diagnostics include projected-gradient and tolerance ratios, step norms,
trust radius, predicted and actual improvement, acceptance statistics, and
the expected post-fit Hessian workload. A function receives the progress
events as named lists.}
\item{sparse_control}{Named control list for the sparse backend. The
\code{gradient} entry may be \code{"auto"} (the default), \code{"finite"},
\code{"exact"}, \code{"stochastic"}, or \code{"hybrid"}.
\code{gradient_probes} controls the fixed Rademacher trace probes used by the
latter two methods; \code{"hybrid"} always finishes with the exact
finite-difference REML objective. \code{gradient_cores} evaluates central
finite-difference directions in parallel on non-Windows systems; a
single-threaded BLAS should be used when it exceeds one.
\code{finite_difference_step} defaults to \code{1e-3}.
\code{outer_optimizer} may be \code{"lbfgsb"} (the default) or the
experimental safeguarded \code{"bfgs_trust"}; the latter requires
\code{gradient="exact"}. \code{optimizer_maxit},
\code{optimizer_gradient_tolerance}, and \code{optimizer_trust_radius}
control its iteration limit, projected-gradient tolerance, and initial trust
radius, respectively.
\code{hessian} selects the post-fit outer-Hessian calculation:
\code{"profiled"} (the default) reconstructs the unprofiled Gaussian REML
Hessian from the lower-dimensional profiled criterion, while
\code{"optimhess"} directly differences the unprofiled criterion;
\code{hessian_step} controls the finite-difference step for either method and
defaults to \code{1e-2}.
\code{supernodal} optionally overrides the automatically
selected CHOLMOD factorization form, and \code{trace_chunk_size} bounds
exact-score solves. \code{schur = "always"} enables the experimental
response-group Schur solver instead of the default \code{"never"}.
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
\item{y}{Numeric response vector for the low-level compressed interface.}
\item{terms}{A \code{cdrgam_term} or list of such terms.}
\item{intercept}{Whether to include a parametric intercept.}
}
\value{
A \code{cdrgam} object. The native backend also inherits from its selected
\pkg{mgcv} fitting class; custom backends provide corresponding methods.
}
\description{
The primary interface compiles untiled impulse and response streams and fits
the expanded model with a native or custom backend.
\code{fit_compressed_cdr_gam()} is the lower-level interface for already
compiled impulse-response terms.
}
