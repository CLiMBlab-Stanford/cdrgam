\name{response_scenarios}
\alias{response_scenarios}
\title{Construct response-history scenarios}
\usage{
response_scenarios(object, at = list(), rows = NULL)
}
\arguments{
\item{object}{A \code{cdrgam_marginaleffects} view.}

\item{at}{Named list of values for \code{lag} or impulse predictors.
\code{lag} is stored in the adapter column \code{.cdrgam_lag}.}

\item{rows}{Optional integer or logical subset of reference incidences.}
}
\value{
A data frame suitable for \code{predict()},
\code{estimate_response_effect()}, or direct use as \code{newdata} in
\pkg{marginaleffects}.
}
\description{
Expands response--impulse incidences from a \code{marginaleffects_view()} over
a Cartesian grid of lag and impulse-predictor values. Each output row retains
the identity of its reference incidence, so prediction replaces one event's
contribution while holding the rest of that response's history fixed.
}
\details{
Changing \code{lag} re-evaluates the selected impulse-response function
windows. A term stops contributing when the scenario lag lies outside its
window.
}
