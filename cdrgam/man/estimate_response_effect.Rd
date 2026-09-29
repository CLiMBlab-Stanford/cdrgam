\name{estimate_response_effect}
\alias{estimate_response_effect}
\title{Generate response-effect plot data}
\usage{
estimate_response_effect(
  object,
  newdata = NULL,
  estimand = c("prediction", "slope", "comparison"),
  variables = NULL,
  by = NULL,
  type = c("response", "link"),
  vcov = TRUE,
  level = 0.95,
  ...
)
}
\arguments{
\item{object}{A \code{cdrgam_marginaleffects} view.}

\item{newdata}{Scenario data, usually returned by
\code{response_scenarios()}. The default uses the unmodified incidences in the
view.}

\item{estimand}{Quantity to estimate: fitted response values, numeric slopes,
or finite comparisons.}

\item{variables}{Predictor specification passed to \pkg{marginaleffects}. It
is required for slopes and comparisons.}

\item{by}{Optional columns defining averages. \code{NULL} retains one result
per scenario row.}

\item{type}{Evaluate on the response or linear-predictor scale.}

\item{vcov}{Include fitted-coefficient uncertainty, or supply a covariance
matrix or function accepted by \code{marginaleffects_view()}.}

\item{level}{Confidence level for pointwise intervals.}

\item{...}{Additional arguments passed to the selected
\pkg{marginaleffects} function.}
}
\value{
A \code{cdrgam_response_effect_grid} data frame. Adapter identity columns and
scenario axes are retained when the requested aggregation permits it.
}
\description{
Evaluates predictions, slopes, or comparisons on response-history scenarios
through the \pkg{marginaleffects} extension API. This preserves each
response's remaining impulse history and returns point estimates with
pointwise uncertainty in a tidy table.
}
\details{
Grouping columns are passed to \pkg{marginaleffects}. Use
\code{.cdrgam_scenario} to average each scenario over its reference
incidences, and add \code{.cdrgam_response} or \code{.cdrgam_impulse} when a
response- or event-weighted result is required.
}
