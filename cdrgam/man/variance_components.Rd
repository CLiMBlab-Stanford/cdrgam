% Generated manually from cdrgam/R/diagnostics.R.
\name{variance_components}
\alias{variance_components}
\title{Extract smoothing variance components}
\usage{
variance_components(object, rescale = TRUE, conf.lev = 0.95)
}
\arguments{
\item{object}{A fitted \code{cdrgam} model.}

\item{rescale}{Apply the penalty rescaling used in the original smooth
specification, as in \code{mgcv::gam.vcomp()}.}

\item{conf.lev}{Confidence level for intervals when an outer REML Hessian is
available.}
}
\value{
The value returned by \code{mgcv::gam.vcomp()}: normally a matrix of standard
deviations and confidence intervals, or a named vector/list when intervals are
unavailable.
}
\description{
Reports the variance components implied by the smoothing parameters of the
expanded \pkg{mgcv} model. Native fits are passed directly to
\code{mgcv::gam.vcomp()}. For block-backend fits, the same calculation is
applied to the stored REML scale, penalty rescaling metadata, and outer
Hessian. The function does not estimate covariance parameters beyond those
present in the expanded model.
}
