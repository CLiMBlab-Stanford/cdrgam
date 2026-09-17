\name{predict_cdrgam}
\alias{predict_cdrgam}
\title{Predict from new impulse and response streams}
\usage{
predict_cdrgam(
  object,
  impulses,
  responses,
  type = c("response", "link", "lpmatrix"),
  se.fit = FALSE,
  allow_new_levels = FALSE,
  chunk_size = 10000,
  unconditional = FALSE
)
}
\arguments{
\item{object}{A fitted \code{cdrgam} model.}
\item{impulses}{A new untiled impulse stream.}
\item{responses}{New response-aligned covariates and timestamps.}
\item{type}{Prediction type: response, link, or linear-predictor matrix.}
\item{se.fit}{Return conditional prediction standard errors.}
\item{allow_new_levels}{Treat unseen grouping levels as zero deviations.}
\item{chunk_size}{Maximum history links transformed at once.}
\item{unconditional}{Include smoothing-parameter uncertainty in standard
errors.}
}
\value{A prediction vector, linear-predictor matrix, or a list with
\code{fit} and \code{se.fit}.}
\description{Constructs CDR prediction designs from new irregular impulse and
response streams using the fitted bases, constraints, knots, and grouping
levels. Standard \code{predict()} methods accept the same streams as a list
with \code{impulses} and \code{responses} entries.}
