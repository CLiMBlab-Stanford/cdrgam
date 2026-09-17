% Generated manually from cdrgam/R/formula.R.
\name{irf}
\alias{irf}
\title{Specify an impulse-response term}
\usage{
irf(
  predictor,
  window = c(0, Inf),
  k = 10,
  bs = "cr",
  nonlinear = FALSE,
  varying = NULL,
  by = NULL,
  group = NULL
)
}
\arguments{
\item{predictor}{Unquoted numeric impulse-stream column, or literal \code{1}
for the deconvolutional intercept.}
\item{window}{Inclusive lag window as \code{c(minimum, maximum)}.}
\item{k}{Lag-basis dimension. For a nonlinear IRF, one value or two values
ordered as lag and predictor dimensions.}
\item{bs}{Basis name; currently \code{"cr"}. For a nonlinear IRF, one value
or separate lag and predictor basis names.}
\item{nonlinear}{If true, fit a tensor surface over lag and predictor value.}
\item{varying}{Optional unquoted response-aligned numeric covariate defining
a centered tensor surface over lag and the covariate.}
\item{by}{Optional unquoted response-aligned numeric covariate multiplying the
complete IRF contribution.}
\item{group}{Optional unquoted response-stream factor defining random IRF
deviations. Add a corresponding ungrouped \code{irf()} term for a
population-plus-group-specific model.}
}
\value{An internal \code{cdrgam_irf_spec} used by the formula compiler.}
\description{
Use \code{irf()} as a top-level additive term inside a \code{fit_cdrgam()}
formula. Ordinary terms retain their usual R and \pkg{mgcv} meanings. By
default the predictor enters linearly as a convolution weight. With
\code{nonlinear = TRUE}, each impulse contributes a fitted bivariate function
of lag and predictor value.
With \code{group}, the term is expanded into factor-specific deviations with
shared roughness penalties and an additional full-rank shrinkage penalty.
Unless explicitly removed with \code{- irf(1)}, the formula compiler inserts
an implicit \code{irf(1)}. An implicit rate term that has zero design after
the required centering constraint is removed with a warning; the same
situation is an error when \code{irf(1)} was explicitly requested.
}
