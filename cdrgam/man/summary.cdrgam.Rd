\name{summary.cdrgam}
\alias{summary.cdrgam}
\alias{summary.cdrgam_block}
\alias{summary.cdrgam_sparse}
\alias{print.summary.cdrgam}
\title{Summarize a fitted CDR-GAM}
\usage{
\method{summary}{cdrgam}(
  object,
  dispersion = NULL,
  freq = FALSE,
  re.test = TRUE,
  all.coefficients = FALSE,
  ...
)

\method{summary}{cdrgam_block}(
  object,
  dispersion = NULL,
  freq = FALSE,
  re.test = TRUE,
  all.coefficients = FALSE,
  ...
)

\method{summary}{cdrgam_sparse}(
  object,
  dispersion = NULL,
  freq = FALSE,
  re.test = TRUE,
  all.coefficients = FALSE,
  ...
)

\method{print}{summary.cdrgam}(x, ...)
}
\arguments{
\item{object}{A fitted \code{cdrgam} model.}

\item{dispersion}{Optional known dispersion.}

\item{freq, re.test}{Controls passed to \code{mgcv::summary.gam()} for a
native \pkg{mgcv} fit. Custom backends retain these arguments for interface
compatibility.}

\item{all.coefficients}{Include the expanded smooth and random-effect
coefficient table. The default is false.}

\item{...}{Additional arguments passed to native
\code{mgcv::summary.gam()} or the print method.}

\item{x}{A CDR-GAM summary object.}
}
\value{
A \code{summary.gam}-style object containing parametric and smooth-term tables
plus the user, normalized, and effective CDR formulas. Native fits retain the
complete \code{summary.gam} result. Custom backends provide its commonly used
coefficient, smooth-term, fit-statistic, and formula fields.
}
\description{
The default output follows \code{mgcv::summary.gam()}: it reports parametric
coefficients and one row per smooth term instead of expanded smooth and
random-effect coefficients. Fully penalized random-effect and grouped IRF
terms report effective degrees of freedom without an approximate term-level
test. Set \code{all.coefficients = TRUE} to request expanded coefficients.
}
