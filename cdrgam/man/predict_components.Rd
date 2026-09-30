\name{predict_components}
\alias{predict_components}
\title{Predict Response and Link Components}
\usage{
predict_components(object, newdata)
}
\arguments{
\item{object}{A fitted \code{cdrgam} model.}

\item{newdata}{Prediction data accepted by \code{predict.cdrgam()}.}
}
\value{
A named list containing primary response- and link-scale estimates and
standard errors. Distributional fits also include parameter names and the
complete parameter matrices on both scales.
}
\description{
\code{predict_components()} returns a stable representation of response- and
link-scale predictions for artifact writers and other downstream tools.
Ordinary models return vectors. Distributional models additionally return
matrices with one column per response-distribution parameter.
}
