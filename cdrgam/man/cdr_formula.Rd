% Generated manually from cdrgam/R/formula.R.
\name{cdr_formula}
\alias{cdr_formula}
\alias{mgcv_formula}
\title{Inspect CDR-GAM formulas}
\usage{
cdr_formula(object, type = c("user", "normalized", "effective"))

mgcv_formula(object)
}
\arguments{
\item{object}{A fitted \code{cdrgam} object; \code{cdr_formula()}
also accepts a prepared \code{cdrgam_design}.}
\item{type}{Return the user formula, normalized CDR formula including
structural defaults, or effective formula after automatic simplification.}
}
\value{The selected CDR formula or its expanded mgcv formula.}
\description{
Retrieve the transparent formula representations retained by the compiler.
}
