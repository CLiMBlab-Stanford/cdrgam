# Plot CDR effects and translated GAM components

Evaluate and plot impulse-response functions by lag, predictor slices at
fixed delays, lag-by-predictor surfaces, or coefficient groups. All CDR
effects are shown on the additive linear-predictor scale. Use
`view="gam"` on a native `mgcv` fit to plot the translated model. When
predictor rescaling was enabled, CDR and ordinary-smooth coordinates are
restored to native source units before display.

## Usage

``` r
# S3 method for class 'cdrgam'
plot(
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
```

## Arguments

- x:

  A fitted `cdrgam` model.

- view:

  Plot type: automatic CDR view, lag IRF, predictor slice,
  lag-by-predictor surface, translated GAM, or coefficients.

- select:

  IRF or coefficient-group labels or one-based indices. Native GAM views
  pass this argument to
  [`mgcv::plot.gam()`](https://rdrr.io/pkg/mgcv/man/plot.gam.html).

- at:

  Named conditioning values. Supported entries are `lag`, `predictor`,
  `group`, and `coefficient`.

- component:

  Effect composition. `"term"` shows the selected term. Grouped terms
  also support population, deviation, and conditional effects. `"total"`
  adds the corresponding baseline to tensor terms and adds the
  population effect to grouped deviations.

- se:

  Draw or return pointwise standard errors.

- unconditional:

  Include smoothing-parameter uncertainty where supported.

- ci_level:

  Confidence level used for plotted intervals.

- n:

  Number of lag evaluation points.

- n_predictor:

  Number of predictor evaluation points.

- surface:

  Surface rendering style.

- pages:

  Number of pages. Zero places all panels on one page; a positive value
  distributes panels across that many pages.

- ask:

  Pause before advancing graphical pages on interactive devices.

- draw:

  Draw the panels. Set to false to return evaluated data only.

- coefficient_limit:

  Maximum coefficients allowed in one panel unless `at$coefficient`
  selects a smaller set. Oversized groups are omitted from the default
  coefficient view and error when explicitly selected.

- xlab, ylab, main, col, lwd, xlim, ylim:

  Graphical controls.

- ...:

  Additional arguments passed to surface plotting functions or to
  [`mgcv::plot.gam()`](https://rdrr.io/pkg/mgcv/man/plot.gam.html) for
  `view="gam"`.

## Value

Invisibly, a `cdrgam_plot_data` list containing evaluated panels.
