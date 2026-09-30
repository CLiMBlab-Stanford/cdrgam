\name{fit_report}
\alias{fit_metadata}
\alias{fit_report}
\title{Extract Stable Fit Metadata and Summary Tables}
\usage{
fit_metadata(object)

fit_report(object)
}
\arguments{
\item{object}{A fitted \code{cdrgam} model.}
}
\value{
\code{fit_metadata()} returns fitting metadata, formulas, distributional
parameter names, and default plotting metadata. \code{fit_report()} adds
printable summary text, parametric and smooth summary tables, and convergence
diagnostics.
}
\description{
\code{fit_report()} returns the fitted-model information needed by artifact
writers and other downstream tools without requiring them to inspect the
internal representation of a \code{cdrgam} object.
}
