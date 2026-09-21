\name{save_cdrgam_plots}
\alias{save_cdrgam_plots}
\title{Save CDR plots to PDF or raster images}
\usage{
save_cdrgam_plots(
  object,
  file,
  device = c("auto", "pdf", "png", "jpeg", "tiff"),
  width = 8,
  height = 6,
  dpi = 144,
  pages = 1,
  ...
)
}
\arguments{
\item{object}{A fitted \code{cdrgam} model.}
\item{file}{Output path.}
\item{device}{Output device. \code{"auto"} uses the file extension.}
\item{width, height}{Device dimensions in inches.}
\item{dpi}{Raster resolution.}
\item{pages}{Number of output pages passed to \code{plot.cdrgam()}.}
\item{...}{Additional arguments passed to \code{plot.cdrgam()}.}
}
\value{
Invisibly, the plot-data object returned by \code{plot.cdrgam()}.
}
\description{
Open a graphics device, call \code{plot.cdrgam()}, and close the device. PDF
output may contain multiple pages. Raster output inserts a numbered format
field before the extension when more than one page is requested.
}
