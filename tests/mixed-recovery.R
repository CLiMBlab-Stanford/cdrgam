library(cdrgam)

group_count <- 8
groups <- paste0('subject_', seq_len(group_count))
random_amplitude <- seq(-0.55, 0.55, length.out=group_count)
random_intercept <- seq(-1.4, 1.4, length.out=group_count)
population_irf <- function(lag) 0.9 * exp(-1.7 * lag)
deviation_shape <- function(lag) {
    exp(-((lag - 0.65) / 0.3)^2) - 0.25 * exp(-1.5 * lag)
}

impulse_streams <- vector('list', group_count)
response_streams <- vector('list', group_count)
for (i in seq_len(group_count)) {
    amplitude <- random_amplitude[[i]]
    total_irf <- function(lag) {
        population_irf(lag) + amplitude * deviation_shape(lag)
    }
    simulation <- simulate_cdr(
        list(x=total_irf),
        n_impulses=180,
        n_responses=220,
        duration=35,
        window=2,
        intercept=75 + random_intercept[[i]],
        noise_sd=0.45,
        seed=6600 + i
    )
    simulation$impulses$subject <- groups[[i]]
    simulation$responses$subject <- groups[[i]]
    impulse_streams[[i]] <- simulation$impulses
    response_streams[[i]] <- simulation$responses
}
impulses <- do.call(rbind, impulse_streams)
responses <- do.call(rbind, response_streams)
responses$subject <- factor(
    responses$subject,
    levels=c(groups, 'declared_but_unobserved')
)

formula <- response ~
    s(subject, bs='re') +
    irf(x, window=c(0, 2), k=9) +
    irf(x, window=c(0, 2), k=9, group=subject) - irf(1)
design <- prepare_cdrgam(
    formula,
    impulses,
    responses,
    series='subject',
    history='auto',
    chunk_size=250,
    quiet=TRUE
)
stopifnot(identical(names(design$terms), c('x', 'x|subject')))
stopifnot(identical(design$terms[[2]]$group_levels, groups))
retained_design <- prepare_cdrgam(
    formula,
    impulses,
    responses,
    series='subject',
    history='auto',
    chunk_size=250,
    drop.unused.levels=FALSE,
    quiet=TRUE
)
stopifnot(identical(
    retained_design$terms[[2]]$group_levels,
    c(groups, 'declared_but_unobserved')
))
# Grouped terms remain compact until a backend chooses dense or sparse
# materialization.
stopifnot(ncol(design$terms[[2]]$X) == 9)
stopifnot(design$terms[[2]]$expanded_dimension == 9 * group_count)
stopifnot(length(design$terms[[2]]$S) == 2L)

native <- cdrgam.fit(design, backend='mgcv', engine='gam', method='REML')
block <- cdrgam.fit(design, backend='block', method='REML')
sparse <- cdrgam.fit(
    design,
    backend='sparse',
    method='REML',
    sparse_control=list(schur='never')
)
stopifnot(is.null(sparse$sparse$schur_group))
stopifnot(sparse$sparse$streamed_grouped_terms == 1L)
stopifnot(max(abs(fitted(native) - fitted(block))) < 5e-4)
stopifnot(max(abs(fitted(native) - fitted(sparse))) < 5e-4)
stopifnot(isTRUE(sparse$converged))
stopifnot(isTRUE(sparse$sparse$convergence$hessian_positive_definite))
native_covariance <- vcov(native)
sparse_covariance <- vcov(sparse)
stopifnot(
    max(abs(native_covariance - sparse_covariance)) /
        max(abs(native_covariance)) < 1e-3
)
prediction_data <- list(
    impulses=impulses,
    responses=responses[names(responses) != 'response']
)
stopifnot(max(abs(predict(native, prediction_data) - fitted(native))) < 1e-8)
stopifnot(max(abs(predict(block, prediction_data) - fitted(block))) < 1e-8)
stopifnot(max(abs(predict(sparse, prediction_data) - fitted(sparse))) < 1e-8)

# Prediction uses the training vocabulary stored in the fit. Local factor
# codes in held-out data therefore do not affect the design after the fit is
# serialized and restored.
restored_path <- tempfile(fileext='.rds')
saveRDS(sparse, restored_path)
restored <- readRDS(restored_path)
unlink(restored_path)
reordered_impulses <- impulses
reordered_responses <- responses[names(responses) != 'response']
reordered_impulses$subject <- factor(
    reordered_impulses$subject,
    levels=rev(groups)
)
reordered_responses$subject <- factor(
    reordered_responses$subject,
    levels=rev(groups)
)
stopifnot(
    identical(
        restored$cdrgam$prediction$random_effects[[1L]]$levels$subject,
        groups
    ),
    identical(restored$cdrgam$terms[[2L]]$group_levels, groups),
    max(abs(predict(
        restored,
        newdata=list(
            impulses=reordered_impulses,
            responses=reordered_responses
        )
    ) - fitted(sparse))) < 1e-8
)

new_impulses <- impulses[impulses$subject == groups[[1L]], ]
new_responses <- responses[responses$subject == groups[[1L]], ]
new_impulses$subject <- 'new_subject'
new_responses$subject <- 'new_subject'
new_responses$response <- NULL
for (fit in list(native, block, sparse)) {
    prediction_warnings <- character()
    new_prediction <- withCallingHandlers(
        predict(
            fit,
            newdata=list(
                impulses=new_impulses,
                responses=new_responses
            )
        ),
        warning=function(warning) {
            prediction_warnings <<- c(
                prediction_warnings,
                conditionMessage(warning)
            )
            invokeRestart('muffleWarning')
        }
    )
    stopifnot(length(new_prediction) == nrow(new_responses))
    stopifnot(all(is.finite(new_prediction)))
    stopifnot(any(grepl('set to zero', prediction_warnings, fixed=TRUE)))
}
new_design <- suppressWarnings(predict(
    restored,
    newdata=list(impulses=new_impulses, responses=new_responses),
    type='lpmatrix'
))
random_columns <- restored$cdrgam$prediction$random_effects[[1L]]$coefficient_index
group_columns <- restored$cdrgam$terms[[2L]]$coefficient_index
deviation_columns <- unique(c(random_columns, group_columns))
stopifnot(
    all(new_design[, deviation_columns, drop=FALSE] == 0),
    any(new_design[, -deviation_columns, drop=FALSE] != 0)
)

# The opt-in compiled Schur factorization solves the same sparse REML model.
schur_sparse <- cdrgam.fit(
    design,
    backend='sparse',
    method='REML',
    sparse_control=list(schur='always')
)
stopifnot(identical(schur_sparse$sparse$schur_group, 'subject'))
stopifnot(max(abs(fitted(sparse) - fitted(schur_sparse))) < 5e-4)

# Both backends report only the variance components defined by the expanded
# mgcv penalties. This includes separate response-side random-intercept,
# population-IRF smoothness, and grouped-IRF penalties, but no extra
# correlations among them.
invisible(capture.output(native_vcomp <- variance_components(native)))
invisible(capture.output(block_vcomp <- variance_components(block)))
invisible(capture.output(sparse_vcomp <- variance_components(sparse)))
invisible(capture.output(schur_vcomp <- variance_components(schur_sparse)))
stopifnot(identical(rownames(native_vcomp), rownames(block_vcomp)))
stopifnot(nrow(native_vcomp) == length(native$sp) + 1L)
stopifnot(isTRUE(all.equal(
    unname(block_vcomp$vc),
    unname(native_vcomp$vc),
    tolerance=5e-3
)))
stopifnot(isTRUE(all.equal(
    unname(sparse_vcomp$vc),
    unname(native_vcomp$vc),
    tolerance=5e-3
)))
stopifnot(isTRUE(all.equal(
    unname(schur_vcomp$vc),
    unname(native_vcomp$vc),
    tolerance=5e-3
)))

lag <- seq(0, 2, length.out=101)
population <- estimate_irf(native, term='x', lag=lag)
deviations <- estimate_irf(
    native,
    term='x|subject',
    lag=lag,
    group=groups
)
block_deviations <- estimate_irf(
    block,
    term='x|subject',
    lag=lag,
    group=groups
)
sparse_deviations <- estimate_irf(
    sparse,
    term='x|subject',
    lag=lag,
    group=groups
)
stopifnot(max(abs(
    deviations$estimate - block_deviations$estimate
)) < 5e-4)
stopifnot(max(abs(
    deviations$estimate - sparse_deviations$estimate
)) < 5e-4)
stopifnot(max(abs(
    deviations$se - sparse_deviations$se
)) < 1e-5)

# Grouped plot data distinguish deviations from population and conditional
# effects, including their joint covariance.
plot_population <- plot(
    native,
    view='irf',
    select='x|subject',
    component='population',
    at=list(lag=lag),
    draw=FALSE
)$panels[[1L]]$data
plot_deviation <- plot(
    native,
    view='irf',
    select='x|subject',
    component='deviation',
    at=list(lag=lag, group=groups[[1L]]),
    draw=FALSE
)$panels[[1L]]$data
plot_conditional <- plot(
    native,
    view='irf',
    select='x|subject',
    component='conditional',
    at=list(lag=lag, group=groups[[1L]]),
    draw=FALSE
)$panels[[1L]]$data
sparse_conditional <- plot(
    sparse,
    view='irf',
    select='x|subject',
    component='conditional',
    at=list(lag=lag, group=groups[[1L]]),
    draw=FALSE
)$panels[[1L]]$data
stopifnot(
    max(abs(plot_population$estimate - population$estimate)) < 1e-10,
    max(abs(
        plot_conditional$estimate -
            (plot_population$estimate + plot_deviation$estimate)
    )) < 1e-10,
    max(abs(
        plot_conditional$estimate - sparse_conditional$estimate
    )) < 5e-4,
    all(is.finite(plot_conditional$se))
)

recovered <- numeric()
target <- numeric()
for (i in seq_along(groups)) {
    rows <- deviations$group == groups[[i]]
    recovered <- c(recovered, population$estimate + deviations$estimate[rows])
    target <- c(
        target,
        population_irf(lag) + random_amplitude[[i]] * deviation_shape(lag)
    )
}
stopifnot(sqrt(mean((recovered - target)^2)) < 0.09)
stopifnot(cor(recovered, target) > 0.97)

# The ordinary response-side random intercept is also recovered.
random_intercept_smooth <- which(vapply(
    native$smooth,
    function(smooth) identical(smooth$label, 's(subject)'),
    logical(1)
))
stopifnot(length(random_intercept_smooth) == 1L)
smooth <- native$smooth[[random_intercept_smooth]]
intercept_estimate <- coef(native)[smooth$first.para:smooth$last.para]
stopifnot(cor(intercept_estimate, random_intercept) > 0.95)

# Crossed response-aligned random intercepts are independently penalized. The
# sparse backend extracts simple bs="re" terms before mgcv setup, preventing
# dense n-by-level model matrices for these terms.
items <- paste0('item_', seq_len(15))
crossed_responses <- responses
crossed_responses$item <- factor(
    rep(items, length.out=nrow(crossed_responses)),
    levels=items
)
item_effect <- 0.7 * cos(seq(0, 2 * pi, length.out=length(items) + 1L))[
    seq_along(items)
]
crossed_responses$response <- crossed_responses$response +
    item_effect[as.integer(crossed_responses$item)]
crossed_design <- prepare_cdrgam(
    response ~ s(subject, bs='re') + s(item, bs='re') +
        irf(x, window=c(0, 2), k=9) - irf(1),
    impulses,
    crossed_responses,
    series='subject',
    quiet=TRUE
)
crossed_native <- cdrgam.fit(
    crossed_design,
    backend='mgcv',
    engine='gam',
    method='REML'
)
crossed_sparse <- cdrgam.fit(
    crossed_design,
    backend='sparse',
    method='REML',
    sparse_control=list(schur='never', hessian='optimhess')
)
crossed_sparse_profiled_hessian <- cdrgam.fit(
    crossed_design,
    backend='sparse',
    method='REML',
    sparse_control=list(schur='never', hessian='profiled')
)
crossed_sparse_analytic_hessian <- cdrgam.fit(
    crossed_design,
    backend='sparse',
    method='REML',
    sparse_control=list(schur='always', hessian='analytic')
)
crossed_hessian_scale <- max(1, max(abs(crossed_sparse$outer.info$hess)))
stopifnot(
    max(abs(
        crossed_sparse$outer.info$hess -
            crossed_sparse_profiled_hessian$outer.info$hess
    )) / crossed_hessian_scale < 3e-3,
    max(abs(
        crossed_sparse$outer.info$hess -
            crossed_sparse_analytic_hessian$outer.info$hess
    )) / crossed_hessian_scale < 3e-3,
    identical(
        crossed_sparse_analytic_hessian$sparse$factor_class,
        'cdrgam_schur_factor'
    ),
    max(abs(
        coef(crossed_sparse) - coef(crossed_sparse_profiled_hessian)
    )) < 1e-8
)
stopifnot(max(abs(
    fitted(crossed_native) - fitted(crossed_sparse)
)) < 5e-4)
stopifnot(identical(
    names(crossed_sparse$sp)[seq_len(2L)],
    c('s(subject)', 's(item)')
))
