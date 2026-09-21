\name{formula.cdrgam}
\alias{formula.cdrgam}
\alias{formula.cdrgam_design}
\title{Extract a CDR-GAM formula}
\usage{
\method{formula}{cdrgam_design}(
  x,
  type = c("user", "normalized", "effective"),
  ...
)

\method{formula}{cdrgam}(
  x,
  type = c("user", "normalized", "effective", "mgcv"),
  ...
)
}
\arguments{
\item{x}{A fitted \code{cdrgam} object or prepared \code{cdrgam_design}.}
\item{type}{Formula representation to return. The \code{"mgcv"} translation
is available for fitted models.}
\item{...}{Unused.}
}
\value{The selected formula.}
\description{
Returns the user formula, the normalized CDR formula with structural defaults,
the effective formula after automatic simplification, or the translated
\pkg{mgcv} formula.
}
