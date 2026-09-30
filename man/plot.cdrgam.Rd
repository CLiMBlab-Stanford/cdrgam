\name{plot.cdrgam}
\alias{plot.cdrgam}
\title{Plot CDR effects and translated GAM components}
\usage{
\method{plot}{cdrgam}(
  x,
  view = c("auto", "irf", "predictor", "surface", "gam", "coef"),
  select = NULL,
  at = list(),
  component = c("term", "total", "population", "deviation", "conditional"),
  se = TRUE,
  unconditional = FALSE,
  ci_level = 0.95,
  n = 200,
  n_predictor = 50,
  surface = c("persp", "image", "contour"),
  pages = 0,
  ask = FALSE,
  draw = TRUE,
  coefficient_limit = 200L,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  col = NULL,
  lwd = 2,
  xlim = NULL,
  ylim = NULL,
  ...
)
}
\arguments{
\item{x}{A fitted \code{cdrgam} model.}
\item{view}{Plot type: automatic CDR view, lag IRF, predictor slice,
lag-by-predictor surface, translated GAM, or coefficients.}
\item{select}{IRF or coefficient-group labels or one-based indices. Native
GAM views pass this argument to \code{mgcv::plot.gam()}.}
\item{at}{Named conditioning values. Supported entries are \code{lag},
\code{predictor}, \code{group}, and \code{coefficient}.}
\item{component}{Effect composition. \code{"term"} shows the selected term.
Grouped terms also support population, deviation, and conditional effects.
\code{"total"} adds the corresponding baseline to tensor terms and adds the
population effect to grouped deviations.}
\item{se}{Draw or return pointwise standard errors.}
\item{unconditional}{Include smoothing-parameter uncertainty where supported.}
\item{ci_level}{Confidence level used for plotted intervals.}
\item{n}{Number of lag evaluation points.}
\item{n_predictor}{Number of predictor evaluation points.}
\item{surface}{Surface rendering style.}
\item{pages}{Number of pages. Zero places all panels on one page; a positive
value distributes panels across that many pages.}
\item{ask}{Pause before advancing graphical pages on interactive devices.}
\item{draw}{Draw the panels. Set to false to return evaluated data only.}
\item{coefficient_limit}{Maximum coefficients allowed in one panel unless
\code{at$coefficient} selects a smaller set. Oversized groups are omitted from
the default coefficient view and error when explicitly selected.}
\item{xlab, ylab, main, col, lwd, xlim, ylim}{Graphical controls.}
\item{...}{Additional arguments passed to surface plotting functions or to
\code{mgcv::plot.gam()} for \code{view="gam"}.}
}
\value{
Invisibly, a \code{cdrgam_plot_data} list containing evaluated panels.
}
\description{
Evaluate and plot impulse-response functions by lag, predictor slices at fixed
delays, lag-by-predictor surfaces, or coefficient groups. All CDR effects are
shown on the additive linear-predictor scale. Use \code{view="gam"} on a
native \code{mgcv} fit to plot the translated model. When predictor rescaling
was enabled, CDR and ordinary-smooth coordinates are restored to native source
units before display.
}
