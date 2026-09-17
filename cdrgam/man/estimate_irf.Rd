% Generated manually from cdrgam/R/diagnostics.R.
\name{estimate_irf}
\alias{estimate_irf}
\title{Evaluate fitted fixed-effect impulse-response functions}
\usage{
estimate_irf(
  object,
  term = NULL,
  lag = NULL,
  n = 200,
  predictor = NULL,
  n_predictor = 25,
  level = NULL,
  se = TRUE,
  unconditional = FALSE
)
}
\arguments{
\item{object}{A fitted \code{cdrgam} model.}
\item{term}{Term name or one-based term number; omit for every IRF.}
\item{lag}{Optional finite numeric lag grid.}
\item{n}{Number of default grid points.}
\item{predictor}{Optional predictor-value grid for nonlinear IRFs.}
\item{n_predictor}{Number of default predictor-value grid points.}
\item{level}{Optional grouping levels to evaluate for a random IRF. By default
all fitted levels are returned.}
\item{se}{Include pointwise standard errors when covariance is available.}
\item{unconditional}{Include the first-order smoothing-parameter uncertainty
correction when supported by the fitted backend.}
}
\value{A data frame containing term, grouping level, lag, predictor value,
estimate, and standard error.}
\description{Evaluates each fitted fixed-effect IRF using its stored marginal
basis, fitted coefficients, and covariance matrix.}
