\name{predict.cdrgam}
\alias{predict.cdrgam}
\alias{predict.cdrgam_block}
\alias{predict.cdrgam_sparse}
\title{Predict from a fitted CDR-GAM}
\usage{
\method{predict}{cdrgam}(object, newdata = NULL, ...)

\method{predict}{cdrgam_block}(object, newdata = NULL, ...)

\method{predict}{cdrgam_sparse}(object, newdata = NULL, ...)
}
\arguments{
\item{object}{A fitted \code{cdrgam} model.}
\item{newdata}{A named list containing \code{impulses} and \code{responses}
data frames. Native \pkg{mgcv} fits also accept a response-side data frame.}
\item{...}{Prediction controls passed to the CDR stream predictor, including
\code{type}, \code{se.fit}, \code{chunk_size}, and \code{unconditional}.}
}
\value{A prediction vector, linear-predictor matrix, or a list containing
\code{fit} and \code{se.fit}.}
\description{
Rebuilds the response-level design from new untiled impulse and response
streams. Stored training-data scaling is applied before prediction. Unseen
random-effect and grouped-IRF levels produce a warning and have their
deviations set to zero.
}
