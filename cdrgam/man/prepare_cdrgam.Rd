% Generated manually from cdrgam/R/formula.R.
\name{prepare_cdrgam}
\alias{prepare_cdrgam}
\title{Compile untiled streams for CDR-GAM fitting}
\usage{
prepare_cdrgam(
  formula,
  impulses,
  responses,
  series = character(),
  impulse_time = "time",
  response_time = "time",
  history = c("auto", "dense", "ragged"),
  history_length = Inf,
  chunk_size = 10000,
  quiet = FALSE
)
}
\arguments{
\item{formula}{Extended formula containing ordinary terms and optional
\code{irf()} terms. An implicit \code{irf(1)} is added unless suppressed.}
\item{impulses}{Data frame with one row per impulse.}
\item{responses}{Data frame with one row per response.}
\item{series}{Character vector naming independent-series columns.}
\item{impulse_time}{Name of the impulse-time column.}
\item{response_time}{Name of the response-time column.}
\item{history}{One of \code{"auto"}, \code{"dense"}, or \code{"ragged"}.}
\item{history_length}{Maximum number of most-recent impulses retained within
each IRF lag window; \code{Inf} retains all impulses.}
\item{chunk_size}{Maximum responses or links transformed in a chunk.}
\item{quiet}{Suppress preparation-plan messages.}
}
\value{A reusable \code{cdrgam_design} object.}
\description{
Validates and aligns untiled streams, selects or applies a history layout, and
constructs response-level spline designs in bounded-memory chunks.
}
