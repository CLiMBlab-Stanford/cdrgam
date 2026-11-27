# Save CDR plots to PDF or raster images

Open a graphics device, call
[`plot.cdrgam()`](https://climblab.org/cdrgam/dev/reference/plot.cdrgam.md),
and close the device. PDF output may contain multiple pages. Raster
output inserts a numbered format field before the extension when more
than one page is requested.

## Usage

``` r
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
```

## Arguments

- object:

  A fitted `cdrgam` model.

- file:

  Output path.

- device:

  Output device. `"auto"` uses the file extension.

- width, height:

  Device dimensions in inches.

- dpi:

  Raster resolution.

- pages:

  Number of output pages passed to
  [`plot.cdrgam()`](https://climblab.org/cdrgam/dev/reference/plot.cdrgam.md).

- ...:

  Additional arguments passed to
  [`plot.cdrgam()`](https://climblab.org/cdrgam/dev/reference/plot.cdrgam.md).

## Value

Invisibly, the plot-data object returned by
[`plot.cdrgam()`](https://climblab.org/cdrgam/dev/reference/plot.cdrgam.md).
