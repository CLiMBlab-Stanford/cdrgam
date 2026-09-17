% Generated manually from cdrgam/R/compressed.R.
\name{compress_cdr_smooth}
\alias{compress_cdr_smooth}
\title{Build a compressed stationary CDR impulse-response term}
\usage{
compress_cdr_smooth(
  t_delta,
  by,
  k = 10,
  bs = "cr",
  knots = NULL,
  chunk_size = 10000,
  name = NULL,
  response_multiplier = NULL,
  constraint_delays = NULL
)
}
\arguments{
\item{t_delta}{Numeric matrix of response-minus-impulse time differences.}
\item{by}{Numeric matrix of convolution weights with the same dimensions as
\code{t_delta}. Invalid or padded positions should have weight zero.}
\item{k}{Basis dimension.}
\item{bs}{Spline basis. Currently only \code{"cr"} is supported.}
\item{knots}{Optional numeric knot vector of length \code{k}.}
\item{chunk_size}{Number of responses transformed at once.}
\item{name}{Optional human-readable term name stored with the model.}
\item{response_multiplier}{Optional response-level multiplier used when
deciding whether linear-functional centering is required.}
\item{constraint_delays}{Optional delay values defining the centering measure.
The stream compiler uses actual, unpadded history links.}
}
\value{
An \code{cdrgam_term} object containing the
response-level design matrix, penalties, and underlying \pkg{mgcv} smooth.
}
\description{
Constructs the response-level design matrix for a stationary,
one-dimensional CDR impulse response without retaining the much larger basis
matrix over all response--impulse pairs. Constant weight sums receive the
same identifiability-centering constraint as an \pkg{mgcv} linear-functional
term. The initial implementation is limited to cubic regression splines.
}
