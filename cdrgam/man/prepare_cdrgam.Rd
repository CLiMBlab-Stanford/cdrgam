% Generated manually from cdrgam/R/formula.R.
\name{prepare_cdrgam}
\alias{prepare_cdrgam}
\title{Compile untiled streams for CDR-GAM fitting}
\usage{
prepare_cdrgam(
  formula,
  impulses,
  responses,
  window = NULL,
  k_l = NULL,
  k_t = NULL,
  k_p = NULL,
  bs_l = NULL,
  bs_t = NULL,
  bs_p = NULL,
  series = character(),
  impulse_time = "time",
  response_time = "time",
  history = c("auto", "dense", "ragged"),
  chunk_size = 10000,
  rescale_predictors = FALSE,
  quiet = FALSE,
  drop.unused.levels = TRUE
)
}
\arguments{
\item{formula}{Extended formula containing ordinary terms and optional
\code{irf()} terms. An implicit \code{irf(1)} is added unless suppressed.}
\item{impulses}{Data frame with one row per impulse.}
\item{responses}{Data frame with one row per response.}
\item{window}{Default inclusive lag window inherited by each \code{irf()}
term that does not supply its own window. \code{NULL} retains the term-level
default \code{c(0, Inf)}. Every impulse in the applicable series and window is
included in the response-level design.}
\item{k_l, k_t, k_p}{Optional model-level defaults for the corresponding
\code{irf()} axis dimensions. An argument supplied by an individual term,
including an explicit \code{NULL}, overrides its model-level default.}
\item{bs_l, bs_t, bs_p}{Optional model-level defaults for the corresponding
\code{irf()} marginal basis names. Term-level arguments take precedence.}
\item{series}{Character vector naming independent-series columns.}
\item{impulse_time}{Name of the impulse-time column.}
\item{response_time}{Name of the response-time column.}
\item{history}{One of \code{"auto"}, \code{"dense"}, or \code{"ragged"}.}
\item{chunk_size}{Maximum responses or links transformed in a chunk.}
\item{rescale_predictors}{If \code{TRUE}, divide continuous numeric
predictors by their training-data standard deviations without centering.
Binary, constant, categorical, grouping, and series variables are unchanged.
Lag and explicit time predictors share the standard deviation of the untiled
impulse-time column, after history links have been constructed in native time
units.}
\item{drop.unused.levels}{Drop factor levels without training observations.
The default matches \code{mgcv::gam()}. Retained levels remain part of the
fitted factor vocabulary.}
\item{quiet}{Suppress preparation-plan messages.}
}
\value{A reusable \code{cdrgam_design} object. Its \code{simplifications}
data frame records automatic axis-basis reductions and
lower-dimensional representations of discrete predictor axes.}
\description{
Validates and aligns untiled streams, selects or applies a history layout, and
constructs response-level spline designs in bounded-memory chunks. Requested
axis basis dimensions are reduced when the linked data support fewer distinct
values and the marginal basis minimum remains.
Nonlinear IRFs with one or two linked predictor values are represented as
linear IRFs using their original numeric coding. Optional scaling is an
upstream data transformation and does not modify the user, normalized, or
effective model formula. Its training-data divisors and exclusions are stored
in the returned object's \code{scaling} component.
}
