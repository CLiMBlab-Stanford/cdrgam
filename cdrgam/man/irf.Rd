% Generated manually from cdrgam/R/formula.R.
\name{irf}
\alias{irf}
\title{Specify an impulse-response term}
\usage{
irf(
  predictor,
  ...,
  window = NULL,
  k_l = 10,
  k_t = NULL,
  k_p = NULL,
  bs_l = "cr",
  bs_t = "cr",
  bs_p = "cr",
  by = NULL,
  group = NULL,
  k = NULL,
  bs = NULL,
  nonlinear = NULL,
  varying = NULL
)
}
\arguments{
\item{predictor}{Unquoted numeric impulse-stream column, or literal \code{1}
for the deconvolutional intercept.}
\item{...}{Additional unquoted numeric impulse-stream predictors. Multiple
predictors define one interaction IRF.}
\item{window}{Inclusive lag window as \code{c(minimum, maximum)}, where lag is
response time minus impulse time. Negative lags include impulses after the
response time. \code{NULL} inherits the model-level \code{window}; without
either, the default is \code{c(0, Inf)}.}
\item{k_l}{Lag-basis dimension. This axis is always smooth and cannot be
\code{NULL}.}
\item{k_t}{Optional response-time basis dimension. \code{NULL} makes the IRF
stationary; a number models nonstationarity against the configured response
time.}
\item{k_p}{Predictor-axis dimensions, one per predictor. \code{NULL} means
linear and a number gives a smooth marginal. Use a list such as
\code{list(NULL, 4)} for mixed axes because R reduces \code{c(NULL, 4)} to
\code{4}.}
\item{bs_l, bs_t}{Lag and response-time marginal basis names accepted by
\pkg{mgcv}.}
\item{bs_p}{Predictor marginal basis names; one value is recycled or supply
one per predictor.}
\item{k, bs, nonlinear, varying}{Compatibility controls for the original
single-predictor API. They cannot be mixed with explicitly supplied
axis-specific controls.}
\item{by}{Optional unquoted response-aligned numeric covariate multiplying the
complete IRF contribution.}
\item{group}{Optional unquoted response-stream factor defining grouped IRF
deviations. Add a corresponding ungrouped \code{irf()} term for a
population-plus-group-specific model.}
}
\value{An internal \code{cdrgam_irf_spec} used by the formula compiler.}
\description{
Use \code{irf()} as a top-level additive term inside a \code{cdrgam()}
formula. Ordinary terms retain their usual R and \pkg{mgcv} meanings. By
default each predictor enters linearly as a convolution weight. Numeric
\code{k_p} entries add smooth predictor marginals and \code{k_t} adds a
response-time marginal. Smooth non-lag marginals are centered, so a tensor
term represents the interaction beyond its lower-order terms.
Before constraints, a tensor term's coefficient count grows as the product of
its axis dimensions. Larger bases are supported, but each added axis therefore
increases preparation, fitting, and inference costs multiplicatively.
With \code{group}, the term is expanded into factor-specific deviations with
shared roughness penalties and an additional full-rank shrinkage penalty.
Unless explicitly removed with \code{- irf(1)}, the formula compiler inserts
an implicit \code{irf(1)}. An implicit rate term that has zero design after
the required centering constraint is removed with a warning; the same
situation is an error when \code{irf(1)} was explicitly requested.
Requested axis basis dimensions are upper bounds. The
stream compiler reduces a dimension when the linked data contain fewer
distinct values, subject to the marginal basis minimum, and
records the change in the prepared design.
}
