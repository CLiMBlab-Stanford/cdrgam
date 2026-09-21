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
Fits from the native backend subclass their \pkg{mgcv} \code{gam} or
\code{bam} result. Custom-backend fits implement the corresponding standard
model methods without inheriting from \code{gam}. These helpers test the common
CDR-GAM class or remove its subclass from a native-backend fit.
}
