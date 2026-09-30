library(cdrgam)
library(ggplot2)

arguments <- commandArgs(trailingOnly=TRUE)
output_directory <- if (length(arguments)) arguments[[1L]] else {
    file.path('validation', 'output', 'location-scale-criteria')
}
dir.create(output_directory, recursive=TRUE, showWarnings=FALSE)

seed <- 20260929L
window <- c(0, 2)
family <- cdrgam_family('gaulss')
minimum_sigma <- environment(family$ll)$b

location_irf <- function(lag) {
    0.72 * exp(-1.7 * lag) +
        0.30 * exp(-0.5 * ((lag - 0.55) / 0.20)^2) -
        0.10 * exp(-0.5 * ((lag - 1.25) / 0.32)^2)
}

scale_irf <- function(lag) {
    0.20 * exp(-0.5 * ((lag - 0.35) / 0.25)^2) -
        0.10 * exp(-1.1 * lag)
}

convolve_stream <- function(response_time, impulses, variable, kernel) {
    vapply(response_time, function(time) {
        lag <- time - impulses$time
        selected <- lag >= window[[1L]] & lag <= window[[2L]]
        sum(impulses[[variable]][selected] * kernel(lag[selected]))
    }, numeric(1))
}

simulate_streams <- function(seed, duration, n_impulses, n_responses) {
    set.seed(seed)
    impulses <- data.frame(
        time=sort(stats::runif(n_impulses, 0, duration)),
        location_signal=stats::rnorm(n_impulses),
        scale_signal=stats::rnorm(n_impulses)
    )
    response_time <- sort(stats::runif(
        n_responses, window[[2L]], duration
    ))
    location <- 0.25 + convolve_stream(
        response_time, impulses, 'location_signal', location_irf
    )
    scale_predictor <- log(0.80 - minimum_sigma) + convolve_stream(
        response_time, impulses, 'scale_signal', scale_irf
    )
    sigma <- minimum_sigma + exp(scale_predictor)
    responses <- data.frame(
        time=response_time,
        response=stats::rnorm(n_responses, location, sigma)
    )
    list(
        impulses=impulses,
        responses=responses,
        truth=data.frame(
            time=response_time,
            location=location,
            scale_predictor=scale_predictor,
            sigma=sigma,
            response=responses$response
        )
    )
}

training <- simulate_streams(
    seed=seed, duration=90, n_impulses=450L, n_responses=700L
)
test <- simulate_streams(
    seed=seed + 1L, duration=50, n_impulses=250L, n_responses=400L
)

formulas <- list(
    location=response ~
        irf(location_signal, window=c(0, 2), k_l=10) - irf(1),
    scale=~ irf(scale_signal, window=c(0, 2), k_l=10) - irf(1)
)
design <- prepare_cdrgam(
    formulas,
    training$impulses,
    training$responses,
    history='ragged',
    chunk_size=128L,
    quiet=TRUE
)

fit_specifications <- data.frame(
    id=c('mgcv_REML', 'mgcv_QNCV', 'cdrgam_REML', 'cdrgam_QNCV'),
    backend=c('mgcv', 'mgcv', 'sparse', 'sparse'),
    method=c('REML', 'QNCV', 'REML', 'QNCV'),
    stringsAsFactors=FALSE
)
fits <- vector('list', nrow(fit_specifications))
names(fits) <- fit_specifications$id
elapsed <- numeric(nrow(fit_specifications))

for (index in seq_len(nrow(fit_specifications))) {
    specification <- fit_specifications[index, ]
    message(
        'Fitting ', specification$id, ' (',
        index, '/', nrow(fit_specifications), ')'
    )
    timing <- system.time({
        fits[[specification$id]] <- if (identical(
                specification$backend, 'sparse'
            )) {
            cdrgam.fit(
                design,
                family=family,
                backend='sparse',
                method=specification$method,
                sparse_control=list(
                    gradient='exact',
                    optimizer_maxit=300L,
                    optimizer_gradient_tolerance=2e-3,
                    qncv_batch_size=128L
                )
            )
        } else {
            cdrgam.fit(
                design,
                family=family,
                backend='mgcv',
                engine='gam',
                method=specification$method
            )
        }
    })
    elapsed[[index]] <- unname(timing[['elapsed']])
}

test_data <- list(
    impulses=test$impulses,
    responses=test$responses['time']
)
prediction_rows <- list()
curve_rows <- list()
metric_rows <- list()

for (index in seq_len(nrow(fit_specifications))) {
    specification <- fit_specifications[index, ]
    id <- specification$id
    fit <- fits[[id]]
    link <- predict(fit, newdata=test_data, type='link')
    colnames(link) <- c('location', 'scale')
    sigma <- minimum_sigma + exp(link[, 'scale'])
    prediction_rows[[id]] <- rbind(
        data.frame(
            fit=id, parameter='location',
            truth=test$truth$location, estimate=link[, 'location']
        ),
        data.frame(
            fit=id, parameter='log sigma',
            truth=log(test$truth$sigma), estimate=log(sigma)
        )
    )
    location_curve <- estimate_irf(
        fit, term='location:location_signal', n=201L
    )
    scale_curve <- estimate_irf(
        fit, term='scale:scale_signal', n=201L
    )
    location_curve$fit <- id
    location_curve$parameter <- 'location'
    location_curve$truth <- location_irf(location_curve$lag)
    scale_curve$fit <- id
    scale_curve$parameter <- 'scale predictor'
    scale_curve$truth <- scale_irf(scale_curve$lag)
    curves <- rbind(location_curve, scale_curve)
    curves$lower <- curves$estimate - stats::qnorm(0.975) * curves$se
    curves$upper <- curves$estimate + stats::qnorm(0.975) * curves$se
    curve_rows[[id]] <- curves
    location_rows <- curves$parameter == 'location'
    scale_rows <- curves$parameter == 'scale predictor'
    nll <- mean(
        log(sigma) +
            0.5 * ((test$truth$response - link[, 'location']) / sigma)^2 +
            0.5 * log(2 * pi)
    )
    criterion <- if (length(fit$gcv.ubre)) {
        unname(fit$gcv.ubre[[1L]])
    } else if (length(fit$reml)) {
        unname(fit$reml[[1L]])
    } else NA_real_
    optimizer_gradient <- if (identical(specification$backend, 'sparse')) {
        gradient <- if (length(fit$optimizer$projected_gradient)) {
            fit$optimizer$projected_gradient
        } else fit$optimizer$gradient
        max(abs(gradient))
    } else if (length(fit$outer.info$grad)) {
        max(abs(fit$outer.info$grad))
    } else NA_real_
    converged <- if (identical(specification$backend, 'sparse')) {
        isTRUE(fit$converged)
    } else {
        identical(fit$outer.info$conv, 'full convergence')
    }
    metric_rows[[id]] <- data.frame(
        fit=id,
        backend=specification$backend,
        method=specification$method,
        converged=converged,
        elapsed_seconds=elapsed[[index]],
        criterion=criterion,
        outer_gradient=optimizer_gradient,
        location_irf_rmse=sqrt(mean(
            (curves$estimate[location_rows] - curves$truth[location_rows])^2
        )),
        location_irf_correlation=stats::cor(
            curves$estimate[location_rows], curves$truth[location_rows]
        ),
        location_irf_coverage=mean(
            curves$lower[location_rows] <= curves$truth[location_rows] &
                curves$truth[location_rows] <= curves$upper[location_rows]
        ),
        scale_irf_rmse=sqrt(mean(
            (curves$estimate[scale_rows] - curves$truth[scale_rows])^2
        )),
        scale_irf_correlation=stats::cor(
            curves$estimate[scale_rows], curves$truth[scale_rows]
        ),
        scale_irf_coverage=mean(
            curves$lower[scale_rows] <= curves$truth[scale_rows] &
                curves$truth[scale_rows] <= curves$upper[scale_rows]
        ),
        heldout_location_rmse=sqrt(mean(
            (link[, 'location'] - test$truth$location)^2
        )),
        heldout_log_sigma_rmse=sqrt(mean(
            (log(sigma) - log(test$truth$sigma))^2
        )),
        heldout_mean_nll=nll,
        heldout_95_coverage=mean(
            abs(test$truth$response - link[, 'location']) <=
                stats::qnorm(0.975) * sigma
        ),
        stringsAsFactors=FALSE
    )
}

curves <- do.call(rbind, curve_rows)
predictions <- do.call(rbind, prediction_rows)
metrics <- do.call(rbind, metric_rows)

backend_comparisons <- do.call(rbind, lapply(c('REML', 'QNCV'), function(
        method
) {
    native_id <- paste0('mgcv_', method)
    sparse_id <- paste0('cdrgam_', method)
    native_curve <- curves[curves$fit == native_id, ]
    sparse_curve <- curves[curves$fit == sparse_id, ]
    native_link <- predict(fits[[native_id]], newdata=test_data, type='link')
    sparse_link <- predict(fits[[sparse_id]], newdata=test_data, type='link')
    colnames(native_link) <- colnames(sparse_link) <- c('location', 'scale')
    data.frame(
        method=method,
        coefficient_max_abs_difference=max(abs(
            unname(stats::coef(fits[[native_id]])) -
                unname(stats::coef(fits[[sparse_id]]))
        )),
        irf_max_abs_difference=max(abs(
            native_curve$estimate - sparse_curve$estimate
        )),
        heldout_link_max_abs_difference=max(abs(native_link - sparse_link)),
        stringsAsFactors=FALSE
    )
}))

fit_labels <- c(
    mgcv_REML='mgcv / REML',
    mgcv_QNCV='mgcv / QNCV',
    cdrgam_REML='CDR-GAM / REML',
    cdrgam_QNCV='CDR-GAM / QNCV'
)
palette <- c(
    mgcv_REML='#0072B2',
    mgcv_QNCV='#56B4E9',
    cdrgam_REML='#D55E00',
    cdrgam_QNCV='#E69F00'
)
curves$fit <- factor(curves$fit, levels=names(fit_labels))
predictions$fit <- factor(predictions$fit, levels=names(fit_labels))

irf_plot <- ggplot(
    curves,
    aes(x=lag, y=estimate, color=fit, fill=fit)
) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.10, color=NA) +
    geom_line(linewidth=0.8) +
    geom_line(
        aes(y=truth),
        color='black',
        linewidth=1.0,
        linetype='dashed'
    ) +
    facet_grid(parameter ~ fit, scales='free_y') +
    scale_color_manual(values=palette, labels=fit_labels) +
    scale_fill_manual(values=palette, labels=fit_labels) +
    labs(
        x='Lag (s)', y='IRF estimate',
        color=NULL, fill=NULL,
        title='Location-scale IRF recovery',
        subtitle='Dashed black lines show the generating IRFs'
    ) +
    theme_bw(base_size=11) +
    theme(legend.position='none')

prediction_plot <- ggplot(
    predictions,
    aes(x=truth, y=estimate, color=fit)
) +
    geom_abline(slope=1, intercept=0, color='grey35', linetype='dashed') +
    geom_point(alpha=0.28, size=1.0) +
    facet_grid(parameter ~ fit, scales='free') +
    scale_color_manual(values=palette, labels=fit_labels, guide='none') +
    labs(
        x='True held-out value',
        y='Estimated held-out value',
        title='Held-out parameter recovery'
    ) +
    theme_bw(base_size=11)

metric_names <- c(
    'location_irf_rmse', 'scale_irf_rmse',
    'heldout_location_rmse', 'heldout_log_sigma_rmse',
    'heldout_mean_nll'
)
metric_long <- do.call(rbind, lapply(metric_names, function(name) {
    data.frame(
        fit=factor(metrics$fit, levels=names(fit_labels)),
        metric=name,
        value=metrics[[name]]
    )
}))
metric_plot <- ggplot(
    metric_long,
    aes(x=fit, y=value, fill=fit)
) +
    geom_col(width=0.72) +
    facet_wrap(~ metric, scales='free_y', ncol=2L) +
    scale_x_discrete(labels=fit_labels) +
    scale_fill_manual(values=palette, guide='none') +
    labs(x=NULL, y=NULL, title='Recovery and predictive metrics') +
    theme_bw(base_size=11) +
    theme(axis.text.x=element_text(angle=30, hjust=1))

utils::write.csv(
    metrics, file.path(output_directory, 'fit-metrics.csv'),
    row.names=FALSE
)
utils::write.csv(
    backend_comparisons,
    file.path(output_directory, 'backend-comparisons.csv'),
    row.names=FALSE
)
utils::write.csv(
    curves, file.path(output_directory, 'irf-estimates.csv'),
    row.names=FALSE
)
ggsave(
    file.path(output_directory, 'irf-recovery.png'),
    irf_plot, width=13, height=7, dpi=160
)
ggsave(
    file.path(output_directory, 'heldout-predictions.png'),
    prediction_plot, width=12, height=7, dpi=160
)
ggsave(
    file.path(output_directory, 'fit-metrics.png'),
    metric_plot, width=10, height=8, dpi=160
)
saveRDS(
    list(
        seed=seed,
        formulas=formulas,
        training=training,
        test=test,
        fits=fits,
        metrics=metrics,
        backend_comparisons=backend_comparisons
    ),
    file.path(output_directory, 'diagnostic.rds')
)

report <- c(
    'Location-scale criterion diagnostic',
    '',
    paste(capture.output(print(metrics, row.names=FALSE)), collapse='\n'),
    '',
    'Backend agreement',
    '',
    paste(
        capture.output(print(backend_comparisons, row.names=FALSE)),
        collapse='\n'
    ),
    '',
    paste(capture.output(sessionInfo()), collapse='\n')
)
writeLines(report, file.path(output_directory, 'report.txt'))
print(metrics, row.names=FALSE)
print(backend_comparisons, row.names=FALSE)
stopifnot(
    all(metrics$converged),
    all(metrics$location_irf_correlation > 0.98),
    all(metrics$scale_irf_correlation > 0.80),
    all(metrics$location_irf_rmse < 0.08),
    all(metrics$scale_irf_rmse < 0.06),
    all(backend_comparisons$coefficient_max_abs_difference < 5e-3),
    all(backend_comparisons$irf_max_abs_difference < 5e-3),
    all(backend_comparisons$heldout_link_max_abs_difference < 5e-3)
)
