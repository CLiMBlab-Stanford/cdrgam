\name{effect_catalog}
\alias{effect_catalog}
\title{List estimable impulse-response effects}
\usage{
effect_catalog(object, impulses = NULL, responses = NULL)
}
\arguments{
\item{object}{A fitted \code{cdrgam} model.}
\item{impulses, responses}{Optional training streams used to recover axis
summaries from models fitted before summaries were stored.}
}
\value{
A data frame with one row per impulse-response term. The \code{predictors},
\code{axes}, and \code{axis_summaries} columns are list columns.
}
\description{
Returns the term and source-scale axis metadata needed to construct effect
evaluation grids.
}
