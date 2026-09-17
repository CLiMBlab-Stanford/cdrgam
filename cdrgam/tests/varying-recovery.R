library(cdrgam)

base <- simulate_cdr(
    list(x=function(lag) exp(-1.7 * lag)),
    n_impulses=320,
    n_responses=420,
    duration=45,
    window=2,
    intercept=0,
    noise_sd=0,
    seed=8811
)
set.seed(8812)
base$responses$z <- as.numeric(scale(rnorm(nrow(base$responses))))

# Response-aligned by covariates multiply the complete convolved IRF.
base$responses$response <- 15 +
    base$responses$z * base$components[, 'x'] +
    rnorm(nrow(base$responses), sd=0.35)
by_design <- prepare_cdrgam(
    response ~ irf(x, window=c(0, 2), k=8, by=z) - irf(1),
    base$impulses,
    base$responses,
    history='ragged',
    quiet=TRUE
)
stopifnot(identical(names(by_design$terms), 'x:z'))
by_native <- fit_cdrgam(by_design, backend='mgcv', engine='gam', method='REML')
by_block <- fit_cdrgam(by_design, backend='block', method='REML')
by_sparse <- fit_cdrgam(by_design, backend='sparse', method='REML')
stopifnot(max(abs(fitted(by_native) - fitted(by_block))) < 5e-4)
stopifnot(max(abs(fitted(by_native) - fitted(by_sparse))) < 5e-4)

# A varying term is a centered lag-by-response-covariate tensor interaction.
links <- cdrgam:::.build_history_links(
    base$impulses,
    base$responses,
    series=character(),
    impulse_time='time',
    response_time='time',
    window=c(0, 2)
)
contribution <- base$impulses$x[links$impulse_index] *
    (0.7 + 0.4 * base$responses$z[links$response_index]) *
    exp(-1.7 * links$delay)
summed <- rowsum(contribution, links$response_index, reorder=FALSE)
signal <- numeric(nrow(base$responses))
signal[as.integer(rownames(summed))] <- summed[, 1L]
base$responses$response <- 20 + signal +
    rnorm(nrow(base$responses), sd=0.35)
varying_design <- prepare_cdrgam(
    response ~
        irf(x, window=c(0, 2), k=8) +
        irf(x, window=c(0, 2), k=c(8, 5), varying=z) - irf(1),
    base$impulses,
    base$responses,
    history='ragged',
    quiet=TRUE
)
stopifnot(identical(names(varying_design$terms), c('x', 'x~z')))
varying_native <- fit_cdrgam(
    varying_design,
    backend='mgcv',
    engine='gam',
    method='REML'
)
varying_block <- fit_cdrgam(varying_design, backend='block', method='REML')
varying_sparse <- fit_cdrgam(varying_design, backend='sparse', method='REML')
stopifnot(max(abs(fitted(varying_native) - fitted(varying_block))) < 8e-4)
stopifnot(max(abs(fitted(varying_native) - fitted(varying_sparse))) < 8e-4)

surface <- estimate_irf(
    varying_sparse,
    term='x~z',
    lag=seq(0, 2, length.out=31),
    predictor=c(-1, 0, 1)
)
stopifnot(nrow(surface) == 93L, all(is.finite(surface$estimate)))
prediction_data <- list(
    impulses=base$impulses,
    responses=base$responses[names(base$responses) != 'response']
)
stopifnot(max(abs(
    predict(varying_sparse, prediction_data) - fitted(varying_sparse)
)) < 1e-8)
