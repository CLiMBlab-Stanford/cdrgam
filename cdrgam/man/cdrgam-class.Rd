% Generated manually from cdrgam/R/compressed.R.
\name{cdrgam-class}
\alias{is_cdrgam}
\alias{as_gam}
\title{CDR-GAM model class helpers}
\usage{
is_cdrgam(x)

as_gam(x)
}
\arguments{
\item{x}{An R object; for \code{as_gam()}, a fitted \code{cdrgam} object from
the native mgcv backend.}
}
\value{
\code{is_cdrgam()} returns a logical scalar. \code{as_gam()} returns the
same fitted model with the \code{cdrgam} subclass removed.
}
\description{
A fitted \code{cdrgam} object directly subclasses its native \pkg{mgcv}
\code{gam} or \code{bam} result, so standard model methods work without
unwrapping. These helpers test the class or remove only its CDR-GAM subclass.
}
