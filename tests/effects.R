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

# covariance = TRUE returns the rows' joint covariance: its diagonal is the
# squared standard errors, and a reference contrast's variance is
# C[i, i] + C[j, j] - 2 C[i, j].
joint <- estimate_effect(
    fit, terms=list(predictors='x'),
    axes=list(lag=list(grid='fitted', n=11L),
        predictors=list(x=list(at=list(summary='mean', `offset-sd`=1)))),
    composition='total', covariance=TRUE
)
joint_covariance <- attr(joint, 'covariance')
contrast <- estimate_effect(
    fit, terms=list(predictors='x'),
    axes=list(lag=list(values=joint$lag),
        predictors=list(x=list(at=list(summary='mean', `offset-sd`=1)))),
    reference=list(lag=joint$lag[[4L]]), composition='total'
)
stopifnot(
    is.matrix(joint_covariance), identical(dim(joint_covariance), c(11L, 11L)),
    isTRUE(all.equal(sqrt(diag(joint_covariance)), joint$se)),
    isTRUE(all.equal(
        contrast$se^2,
        diag(joint_covariance) + joint_covariance[4L, 4L] - 2 * joint_covariance[, 4L]
    )),
    is.null(attr(curve, 'covariance'))
)
every <- estimate_effect(fit, axes=list(lag=list(grid='fitted', n=5L)), covariance=TRUE)
stopifnot(
    identical(dim(attr(every, 'covariance')), rep(nrow(every), 2L)),
    isTRUE(all.equal(sqrt(diag(attr(every, 'covariance'))), every$se)),
    grepl('requires se', tryCatch(estimate_effect(
        fit, axes=list(lag=list(grid='fitted', n=5L)), se=FALSE, covariance=TRUE
    ), error=conditionMessage))
)

# The sparse backend's selected covariance gives the same joint covariance.
sparse_fit <- cdrgam.fit(design, backend='sparse')
sparse_joint <- estimate_effect(
    sparse_fit, axes=list(lag=list(grid='fitted', n=7L)), covariance=TRUE
)
stopifnot(
    is.matrix(attr(sparse_joint, 'covariance')),
    identical(dim(attr(sparse_joint, 'covariance')), rep(nrow(sparse_joint), 2L)),
    isTRUE(all.equal(sqrt(diag(attr(sparse_joint, 'covariance'))), sparse_joint$se))
)

# A term that boundary reduction removed leaves an extra specification in
# fits saved before the specifications were subset with the terms. Total
# composition still finds the hierarchy instead of falling back to the term.
removed <- fit
specifications <- removed$cdrgam$preparation$specification
removed$cdrgam$preparation$specification <- c(
    specifications[1L],
    list(utils::modifyList(specifications[[length(specifications)]], list(group='subject'))),
    specifications[-1L]
)
stopifnot(isTRUE(all.equal(
    estimate_effect(
        removed, terms=list(predictors='x'),
        axes=list(lag=list(grid='fitted', n=21L),
            predictors=list(x=list(at=list(summary='mean', `offset-sd`=1)))),
        composition='total'
    )$estimate,
    curve$estimate
)))

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

# simulate_effect() draws the estimate_effect() functional from the
# coefficient posterior: the draw mean and SD approach the estimate and its
# standard error, and a seed makes draws reproducible without changing the
# global random-number state.
request <- list(
    fit, terms=list(predictors='x'),
    axes=list(
        lag=list(grid='fitted', n=11L),
        predictors=list(x=list(at=list(summary='mean', `offset-sd`=1)))
    ),
    reference=list(x=list(at=list(summary='mean'))), composition='total'
)
target <- do.call(estimate_effect, request)
set.seed(3)
before <- .Random.seed
simulated <- do.call(simulate_effect, c(request, list(draws=4000L, seed=11)))
again <- do.call(simulate_effect, c(request, list(draws=4000L, seed=11)))
stopifnot(
    identical(.Random.seed, before),
    identical(simulated, again),
    is.matrix(simulated$draws),
    identical(dim(simulated$draws), c(nrow(target), 4000L)),
    isTRUE(all.equal(simulated$grid$estimate, target$estimate)),
    identical(simulated$grid$lag, target$lag),
    max(abs(rowMeans(simulated$draws) - target$estimate) / target$se) < 0.1,
    max(abs(apply(simulated$draws, 1L, stats::sd) / target$se - 1)) < 0.05
)
bad_draws <- tryCatch(
    do.call(simulate_effect, c(request, list(draws=0L))), error=conditionMessage
)
stopifnot(grepl('draws must be a positive integer', bad_draws, fixed=TRUE))

# Shared coefficient draws make separate simulations joint: the same draw
# of one functional, evaluated by two calls over overlapping rows, agrees.
# Their moments match the functional's, and a term the sample does not cover
# is an error.
shared <- sample_effect_coefficients(
    fit, terms=list(predictors='x'), composition='total', draws=4000L, seed=5
)
first <- do.call(simulate_effect, c(request, list(coefficients=shared)))
lagged <- simulate_effect(
    fit, terms=list(predictors='x'), composition='total', coefficients=shared,
    reference=request$reference,
    axes=list(lag=list(values=target$lag[2:3]), predictors=request$axes$predictors)
)
stopifnot(
    inherits(shared, 'cdrgam_coefficient_draws'),
    identical(dim(first$draws), c(nrow(target), 4000L)),
    isTRUE(all.equal(lagged$draws, first$draws[2:3, ])),
    max(abs(rowMeans(first$draws) - target$estimate) / target$se) < 0.1,
    max(abs(apply(first$draws, 1L, stats::sd) / target$se - 1)) < 0.05
)
narrow <- sample_effect_coefficients(
    fit, terms=list(predictors='x'), composition='term', draws=10L, seed=5
)
if (length(narrow$index) < length(shared$index)) {
    uncovered <- tryCatch(
        do.call(simulate_effect, c(request, list(coefficients=narrow))),
        error=conditionMessage
    )
    stopifnot(grepl('coefficients do not cover term', uncovered, fixed=TRUE))
}

# points evaluates paired axis values, one row per point, and matches the
# same rows of the crossed grid, with references and posterior draws.
grid <- estimate_effect(
    fit, terms=list(predictors='x'), composition='total',
    axes=list(lag=list(values=c(0, 0.5, 1)), predictors=list(x=list(values=c(-1, 1)))),
    reference=list(x=list(at=0))
)
path <- data.frame(lag=c(1, 0, 0.5), x=c(-1, 1, 1))
paired <- estimate_effect(
    fit, terms=list(predictors='x'), composition='total', points=path,
    reference=list(x=list(at=0))
)
rows <- match(paste(path$lag, path$x), paste(grid$lag, grid$x))
stopifnot(
    nrow(paired) == 3L,
    isTRUE(all.equal(paired$lag, path$lag)), isTRUE(all.equal(paired$x, path$x)),
    isTRUE(all.equal(paired$estimate, grid$estimate[rows])),
    isTRUE(all.equal(paired$se, grid$se[rows])),
    all(paired$x_reference == 0)
)
# Axes that are not columns take one value from axes.
held <- estimate_effect(
    fit, terms=list(predictors='x'), points=data.frame(lag=c(0, 1)),
    axes=list(predictors=list(x=list(at=1)))
)
crossed <- estimate_effect(
    fit, terms=list(predictors='x'),
    axes=list(lag=list(values=c(0, 1)), predictors=list(x=list(at=1)))
)
stopifnot(isTRUE(all.equal(held$estimate, crossed$estimate)))
drawn <- simulate_effect(
    fit, terms=list(predictors='x'), points=path, draws=50L, seed=3
)
stopifnot(identical(dim(drawn$draws), c(3L, 50L)), nrow(drawn$grid) == 3L)
tested <- test_effect(fit, terms=list(predictors='x'), points=path)
stopifnot(tested$points == 3L)
expect_error_text <- function(expression, pattern) {
    message <- tryCatch({
        force(expression)
        ''
    }, error=conditionMessage)
    stopifnot(grepl(pattern, message, fixed=TRUE))
}
expect_error_text(
    estimate_effect(fit, terms=list(predictors='x'), points=data.frame(x=1)),
    'Axis lag is not a points column, so it must resolve to one value'
)
expect_error_text(
    estimate_effect(fit, terms=list(predictors='x'), points=data.frame(lag=0, y=1)),
    'points columns are not evaluation axes of the term: y'
)
expect_error_text(
    estimate_effect(fit, terms=list(predictors='x'), points=data.frame(lag=NA_real_)),
    'points column lag must contain finite numbers'
)
expect_error_text(
    estimate_effect(fit, points=data.frame()),
    'points must be a data frame with at least one row'
)
