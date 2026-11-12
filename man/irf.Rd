% Generated manually from R/formula.R.
\name{irf}
\alias{irf}
\title{Specify an impulse-response term}
\usage{
irf(
  predictor,
  ...,
  window = NULL,
  knots_l = NULL,
  k_l = 10,
  k_t = NULL,
  k_p = NULL,
  bs_l = "cr",
  bs_t = "cr",
  bs_p = "cr",
  by = NULL,
  group = NULL,
  group_penalty = c("re", "t2", "fs"),
  center = NULL,
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
\item{knots_l}{Optional strictly increasing lag-basis construction points in
the original lag units. Supply exactly \code{k_l} values spanning the linked
training lags. This control is unavailable for \code{bs_l = "ps"}.}
\item{k_l}{Lag-basis dimension. This axis is always smooth and cannot be
\code{NULL}.}
\item{k_t}{Optional response-time basis dimension. \code{NULL} makes the IRF
stationary; a number models nonstationarity against the configured response
time.}
\item{k_p}{Predictor-axis dimensions. Supply one value to recycle it or one
entry per predictor. \code{NULL} means linear and a number gives a smooth
marginal. Use a list such as \code{list(NULL, 4)} for mixed axes because R
reduces \code{c(NULL, 4)} to \code{4}.}
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
\item{group_penalty}{Penalty construction for a grouped IRF. \code{"re"}
uses the marginal smoothness penalties plus one full-rank group shrinkage
penalty. \code{"t2"} matches \code{mgcv::t2(full = FALSE)}, and \code{"fs"}
matches \code{mgcv::t2(full = TRUE)} and a univariate factor smooth. The
\code{"t2"} and \code{"fs"} choices introduce progressively more smoothing
parameters as the number of smooth IRF axes increases.}
\item{center}{Centering of linear predictors. \code{NULL} inherits the
model-level \code{center_predictors} setting. \code{TRUE} or \code{FALSE}
centers all or none of this term's linear predictors. For interactions, a
logical vector with one entry per predictor, such as \code{c(FALSE, TRUE)}, or
a character vector of predictor names selects the predictors to center. A
centered predictor enters as its value minus its mean over the training
history links, so the term is zero for an event at that mean. Predictors with
a single linked value are not centered. Smooth predictor axes are always
centered by their own constraint. \code{irf(1)} has no predictors to
center.}
}
\value{An internal \code{cdrgam_irf_spec} used by the formula compiler.}
\description{
Use \code{irf()} as a top-level additive term inside a \code{cdrgam()}
formula. Ordinary terms retain their usual R and \pkg{mgcv} meanings. By
default each predictor enters linearly as a convolution weight. Numeric
\code{k_p} entries add smooth predictor marginals and \code{k_t} adds a
response-time marginal. Smooth non-lag marginals are centered, so a tensor
term represents the interaction beyond its lower-order terms. Linear
predictors are centered at their training means by default; see
\code{center} and the \code{center_predictors} argument of \code{cdrgam()}.
Before constraints, a tensor term's coefficient count grows as the product of
its axis dimensions. Larger bases are supported, but each added axis therefore
increases preparation, fitting, and inference costs multiplicatively.
With \code{group}, the term is expanded into factor-specific deviations. The
\code{group_penalty} argument selects its shared penalty decomposition.
Unless explicitly removed with \code{- irf(1)}, the formula compiler inserts
an implicit \code{irf(1)}. An implicit rate term that has zero design after
the required centering constraint is removed with a warning; the same
situation is an error when \code{irf(1)} was explicitly requested.
Requested axis basis dimensions are upper bounds. The
stream compiler reduces a dimension when the linked data contain fewer
distinct values, subject to the marginal basis minimum, and
records the change in the prepared design.
}
