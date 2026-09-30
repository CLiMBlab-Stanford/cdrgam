\name{estimate_effect}
\alias{estimate_effect}
\title{Evaluate fitted effects on declarative axis grids}
\usage{
estimate_effect(
  object,
  terms = NULL,
  axes = list(),
  impulses = NULL,
  responses = NULL,
  composition = c("term", "deviation", "total"),
  grouping = c("population", "deviation", "conditional"),
  groups = NULL,
  se = TRUE,
  unconditional = FALSE,
  level = 0.95,
  n = 200L
)
}
\arguments{
\item{object}{A fitted \code{cdrgam} model.}
\item{terms}{Term labels or a selector with \code{names}, \code{predictors},
\code{match}, and \code{grouped} fields.}
\item{axes}{Named axis requests. Predictor requests are nested below
\code{predictors}; \code{"*"} supplies a fallback.}
\item{impulses, responses}{Optional training streams used to recover axis
summaries from models fitted before summaries were stored.}
\item{composition}{Whether to return the selected term or add its fitted
marginal hierarchy.}
\item{grouping}{Whether grouped terms represent population effects,
deviations, or their conditional sum.}
\item{groups}{Optional grouping levels to evaluate.}
\item{se}{Whether to calculate pointwise standard errors.}
\item{unconditional}{Whether to include smoothing-parameter uncertainty when
the fitted backend supports it.}
\item{level}{Confidence level for \code{lower} and \code{upper}.}
\item{n}{Default number of points for fitted-domain grids.}
}
\value{
A \code{cdrgam_effect_grid} data frame containing evaluation axes, estimates,
standard errors, and pointwise confidence limits on the linear-predictor scale.
}
\description{
Evaluates selected impulse-response terms as joint linear functionals of the
fitted coefficients. Axes may vary over numeric values, quantiles, or the
fitted domain, or be fixed at a training-data summary with an optional standard
deviation offset.
}
