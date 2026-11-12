\name{test_effect}
\alias{test_effect}
\title{Test fitted effects jointly over declarative axis grids}
\usage{
test_effect(
  object,
  terms = NULL,
  axes = list(),
  impulses = NULL,
  responses = NULL,
  composition = c("term", "total"),
  grouping = c("population", "deviation", "conditional"),
  groups = NULL,
  reference = list(),
  unconditional = FALSE,
  n = 50L,
  tolerance = 1e-08,
  boundary = c("estimate", "bound")
)
}
\arguments{
\item{object}{A fitted \code{cdrgam} model.}
\item{terms}{Term labels or a selector with \code{names}, \code{predictors},
\code{match}, and \code{grouped} fields.}
\item{axes}{Named axis requests, as in \code{\link{estimate_effect}}.}
\item{impulses, responses}{Optional training streams used to recover axis
summaries from models fitted before summaries were stored.}
\item{composition}{Whether to test the selected term or its fitted marginal
hierarchy.}
\item{grouping}{Whether grouped terms represent population effects,
deviations, or their conditional sum.}
\item{groups}{Optional grouping levels to test.}
\item{reference}{Optional named list of reference values for evaluation axes,
as in \code{\link{estimate_effect}}.}
\item{unconditional}{Whether to include smoothing-parameter uncertainty when
the fitted backend supports it. Smoothing parameters classified as on the
boundary of their range contribute no uncertainty, because the criterion has
no local curvature along them. Sparse fits stop with an error when the outer
Hessian is not positive definite after those parameters are excluded.}
\item{n}{Default number of points for fitted-domain grids.}
\item{tolerance}{Relative eigenvalue threshold for the rank of
\eqn{L V L^\top}.}

\item{boundary}{\code{"estimate"} (the default) uses the fitted smoothing
parameters. \code{"bound"} repeats the calculation with each boundary
smoothing parameter set to its profile bound, the largest variance component
the data still support, while the other smoothing parameters stay at their
estimates. Compare the two to check whether conclusions depend on treating
those variance components as zero. \code{"bound"} requires a sparse Gaussian
REML fit with \code{sparse_control$boundary_sensitivity = TRUE} and
\code{unconditional = FALSE}.}
}
\value{
A \code{cdrgam_effect_test} data frame with \code{term_id}, \code{term},
\code{group}, \code{composition}, \code{grouping}, \code{statistic},
\code{df}, \code{p_value}, and \code{points}. One row is returned per selected
term and grouping level. \code{df} is zero and \code{statistic} is missing when
the functional is identically zero.
}
\description{
Computes a Wald test of the null hypothesis that a functional returned by
\code{\link{estimate_effect}} is zero at every grid point.
}
\details{
The statistic is
\eqn{\hat\beta^\top L^\top (L V L^\top)^{-} L \hat\beta}, where \eqn{L}
stacks the grid rows of the functional and \eqn{V} is the coefficient
covariance. The pseudo-inverse keeps eigenvalues above \code{tolerance} times
the largest. The reference distribution is chi-squared with that rank as its
degrees of freedom. Use a grid dense enough to span the functional; the rank
stops growing once it does.

Penalized coefficients make this test approximate. Wood (2013) reports that
such tests can reject too rarely when the null holds on the boundary of the
smoothing-parameter space.
}
\references{
Wood, S. N. (2013). On p-values for smooth components of an extended
generalized additive model. \emph{Biometrika}, 100(1), 221--228.
}
\seealso{
\code{\link{estimate_effect}}
}
