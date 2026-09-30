library(cdrgam)

simulation <- simulate_cdr(
    list(x=function(lag) exp(-lag)),
    n_impulses=100,
    n_responses=90,
    duration=20,
    window=2,
    seed=417
)

design <- prepare_cdrgam(
    response ~ irf(x, k_l=5, k_p=4),
    simulation$impulses,
    simulation$responses,
    window=c(0, 2),
    quiet=TRUE
)
fit <- cdrgam.fit(design, backend='mgcv', method='REML')

catalog <- effect_catalog(fit)
stopifnot(
    nrow(catalog) == 2L,
    all(c('term_id', 'label', 'predictors', 'axes', 'axis_summaries') %in%
        names(catalog)),
    any(vapply(catalog$predictors, identical, logical(1), 'x'))
)

curve <- estimate_effect(
    fit,
    terms=list(predictors='x'),
    axes=list(
        lag=list(grid='fitted', n=21L),
        predictors=list(x=list(at=list(summary='mean', `offset-sd`=1)))
    ),
    composition='total'
)
stopifnot(
    inherits(curve, 'cdrgam_effect_grid'),
    nrow(curve) == 21L,
    all(c('lag', 'x', 'estimate', 'se', 'lower', 'upper') %in% names(curve)),
    all(is.finite(curve$estimate)),
    all(curve$lower <= curve$estimate & curve$estimate <= curve$upper)
)

surface <- estimate_effect(
    fit,
    terms=list(predictors='x'),
    axes=list(
        lag=list(grid='fitted', n=9L),
        predictors=list(x=list(quantiles=c(0.1, 0.5, 0.9)))
    ),
    composition='term'
)
stopifnot(nrow(surface) == 27L, length(unique(surface$x)) == 3L)

legacy_fit <- fit
legacy_fit$cdrgam$terms <- lapply(legacy_fit$cdrgam$terms, function(term) {
    term$linear_predictor_summaries <- NULL
    if (!is.null(term$axis)) {
        term$axis <- lapply(term$axis, function(axis) {
            axis$summary <- NULL
            axis
        })
    }
    term
})
legacy_curve <- estimate_effect(
    legacy_fit,
    terms=list(predictors='x'),
    axes=list(
        lag=list(grid='fitted', n=11L),
        predictors=list(x=list(at=list(summary='mean', `offset-sd`=1)))
    ),
    impulses=simulation$impulses,
    responses=simulation$responses
)
stopifnot(nrow(legacy_curve) == 11L, all(is.finite(legacy_curve$x)))

varying_design <- prepare_cdrgam(
    response ~ irf(x, k_l=5, k_t=4, k_p=4),
    simulation$impulses,
    simulation$responses,
    window=c(0, 2),
    quiet=TRUE
)
varying_fit <- cdrgam.fit(varying_design, backend='mgcv', method='REML')
varying_surface <- estimate_effect(
    varying_fit,
    terms=list(predictors='x'),
    axes=list(
        lag=list(grid='fitted', n=7L),
        time=list(grid='fitted', n=6L),
        predictors=list(x=list(at=list(summary='median')))
    )
)
stopifnot(
    nrow(varying_surface) == 42L,
    all(c('lag', 'time', 'x') %in% names(varying_surface))
)

grouped_simulation <- simulation
grouped_simulation$responses$subject <- factor(rep(
    c('s1', 's2', 's3'), length.out=nrow(grouped_simulation$responses)
))
grouped_design <- prepare_cdrgam(
    response ~ irf(x, k_l=5) + irf(x, k_l=5, group=subject) - irf(1),
    grouped_simulation$impulses,
    grouped_simulation$responses,
    window=c(0, 2),
    quiet=TRUE
)
grouped_fit <- cdrgam.fit(grouped_design, backend='mgcv', method='REML')
conditional <- estimate_effect(
    grouped_fit,
    terms=list(predictors='x', grouped=TRUE),
    axes=list(lag=list(grid='fitted', n=8L)),
    grouping='conditional',
    groups=c('s1', 's2')
)
stopifnot(
    nrow(conditional) == 16L,
    identical(unique(conditional$group), c('s1', 's2')),
    all(is.finite(conditional$estimate))
)

legacy_linear_fit <- grouped_fit
legacy_linear_fit$cdrgam$terms[[1L]]$linear_predictors <- NULL
legacy_linear_fit$cdrgam$terms[[1L]]$linear_predictor_summaries <- NULL
legacy_specification <- legacy_linear_fit$cdrgam$preparation$specification[[1L]]
legacy_specification$predictor <- legacy_specification$predictors[[1L]]
legacy_specification$predictors <- NULL
legacy_specification$k_p <- NULL
legacy_specification$nonlinear <- FALSE
legacy_linear_fit$cdrgam$preparation$specification[[1L]] <- legacy_specification
legacy_linear <- estimate_effect(
    legacy_linear_fit,
    terms='x',
    axes=list(
        lag=list(grid='fitted', n=6L),
        predictors=list(x=list(at=list(summary='mean', `offset-sd`=1)))
    ),
    impulses=grouped_simulation$impulses,
    responses=grouped_simulation$responses
)
stopifnot(nrow(legacy_linear) == 6L, all(is.finite(legacy_linear$x)))

cat('Effect-grid tests passed\n')
