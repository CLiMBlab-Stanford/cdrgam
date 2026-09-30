\name{marginaleffects_view}
\alias{marginaleffects_view}
\title{Prepare a cdrgam model for marginaleffects}
\usage{
marginaleffects_view(
  object,
  impulses,
  responses,
  terms = NULL,
  component = c("response")
)
}
\arguments{
\item{object}{A fitted single-parameter \code{cdrgam} model.}

\item{impulses, responses}{Impulse and response streams on the source scale.}

\item{terms}{Optional impulse-response term labels or selector list. The
default includes every impulse-response term, including the rate term.}

\item{component}{Quantity represented by each evaluation row. The initial
implementation supports complete response predictions.}
}
\value{
A \code{cdrgam_marginaleffects} model view accepted by
\pkg{marginaleffects}. Its attached model data contain response, impulse, and
lag identifiers plus the scalar impulse predictors available for comparisons
and slopes.
}
\description{
Constructs an evaluation view with one row per response--impulse link.
Changing an impulse predictor in a row replaces that event's contribution
while holding the response's remaining impulse history and covariates fixed.
}
\details{
The returned object implements the \pkg{marginaleffects} coefficient,
covariance, and prediction extension methods. Calls such as
\code{marginaleffects::slopes()} therefore operate on scalar impulse predictors
without treating the fitted object as an \pkg{mgcv} \code{gam}.

Each row is one response--impulse incidence. An ungrouped average over rows is
incidence-weighted, so responses with longer histories receive more weight.
The internal \code{.cdrgam_response} and \code{.cdrgam_impulse} columns can be
used to define response- or event-level aggregation.

The view stores a response-level linear-predictor matrix for the supplied
streams. Its memory use therefore depends on the number of selected responses
and fitted coefficients. Covariance extraction is restricted to coefficients
that affect those responses, including selected sparse covariance solves for
the sparse backend. Distributional models are not yet supported.
}
\examples{
\dontrun{
view <- marginaleffects_view(fit, impulses, responses)
marginaleffects::slopes(view, variables = "surprisal")
marginaleffects::avg_comparisons(
  view,
  variables = list(surprisal = c(-1, 1))
)
}
}
