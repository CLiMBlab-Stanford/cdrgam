\name{as_mgcv_spec}
\alias{as_mgcv_spec}
\title{Export an exact native mgcv fitting specification}
\usage{
as_mgcv_spec(x, family = stats::gaussian(), method = NULL, ...)
}
\arguments{
\item{x}{A \code{cdrgam_design} returned by \code{prepare_cdrgam()}.}
\item{family}{A family accepted by \code{mgcv::gam()}. Distributional designs
require their corresponding general family.}
\item{method}{Optional smoothing-parameter estimation method.}
\item{...}{Additional arguments to include in the \code{mgcv::gam()} call.}
}
\value{
A named argument list suitable for \code{do.call(mgcv::gam, spec)}.
}
\description{
Converts a prepared CDR-GAM design into native \pkg{mgcv} arguments.
Impulse-response functions become meaningfully named matrix terms. Scalar
models carry their exact penalties through \code{paraPen}. Distributional
models use \pkg{mgcv}'s prefit \code{G} interface because
\code{mgcv::gam()} does not accept \code{paraPen} with a formula list. The
export does not use CDR-GAM's private compressed smooth constructor.
}
\details{
The returned data materialize grouped impulse-response designs. This interface
supports interoperability and numerical validation; it is not the
memory-efficient fitting path used by the sparse backend. Boundary-reduced
fitted models are not supported.
}
