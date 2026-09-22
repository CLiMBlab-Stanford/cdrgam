\name{cdrgam_family}
\alias{cdrgam_family}
\title{Construct a supported CDR-GAM response family}
\usage{
cdrgam_family(
  family = c("gaussian", "binomial", "poisson", "Gamma", "gaulss"),
  link = NULL
)
}
\arguments{
\item{family}{A supported family name or an existing family object.}

\item{link}{Optional link name. It may be supplied only when \code{family} is
a name.}
}
\value{
A standard R family object.
}
\description{
\code{cdrgam_family()} is the validated family registry used by the
command-line harness. Family objects may still be passed directly to
\code{cdrgam()} and \code{cdrgam.fit()}. The native \code{mgcv} backend
supports every combination returned here. \code{gaulss} selects mgcv's
two-predictor Gaussian location--scale family. The dense block backend
supports a joint location--scale reference solver plus
\code{gaussian(identity)}, \code{binomial(logit)}, \code{poisson(log)}, and
estimated-dispersion \code{Gamma(log)}. The sparse backend supports those
single-predictor combinations with streamed PIRLS and an exact Laplace score,
plus joint location--scale LAML with an exact outer score.
}
