\name{cdrgam_family}
\alias{cdrgam_family}
\title{Construct a supported CDR-GAM response family}
\usage{
cdrgam_family(
  family = c("gaussian", "binomial", "poisson", "Gamma"),
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
supports every combination returned here. The dense block backend supports
\code{gaussian(identity)} plus \code{binomial(logit)}, \code{poisson(log)},
and estimated-dispersion \code{Gamma(log)}. The sparse backend supports the
same combinations with streamed PIRLS and an exact Laplace score.
}
