# Randomized location-scale equivalence checks across all fitting backends.
# Run from the repository root with:
#   Rscript validation/location_scale_randomized.R

library(cdrgam)

scalar_or_na <- function(value, missing=NA_real_) {
    if (length(value)) value[[1L]] else missing
}

method <- toupper(Sys.getenv('CDRGAM_VALIDATION_METHOD', 'REML'))
if (!method %in% c('REML', 'QNCV')) {
    stop('CDRGAM_VALIDATION_METHOD must be REML or QNCV')
}
trace_solver <- identical(Sys.getenv('CDRGAM_VALIDATION_TRACE'), '1')

simulate_location_scale_streams <- function(seed, subjects=4L) {
    set.seed(seed)
    impulse_streams <- vector('list', subjects)
    response_streams <- vector('list', subjects)
    subject_offsets <- seq(-18, 18, length.out=subjects)
    for (subject_index in seq_len(subjects)) {
        subject <- paste0('s', subject_index)
        impulse_time <- sort(stats::runif(55, 0, 18))
        response_time <- sort(stats::runif(75, 0.5, 18))
        x <- stats::rnorm(length(impulse_time))
        z <- stats::rnorm(length(impulse_time))
        location <- numeric(length(response_time))
        log_sd <- numeric(length(response_time))
        for (row in seq_along(response_time)) {
            lag <- response_time[[row]] - impulse_time
            linked <- lag >= 0 & lag <= 1.5
            location[[row]] <- 320 + subject_offsets[[subject_index]] +
                sum(7 * x[linked] * exp(-1.4 * lag[linked]))
            log_sd[[row]] <- log(38) +
                sum(0.06 * z[linked] * exp(-2 * lag[linked]))
        }
        impulse_streams[[subject_index]] <- data.frame(
            subject=subject, time=impulse_time, x=x, z=z
        )
        response_streams[[subject_index]] <- data.frame(
            subject=subject,
            time=response_time,
            response=stats::rnorm(
                length(response_time), location, exp(log_sd)
            )
        )
    }
    impulses <- do.call(rbind, impulse_streams)
    responses <- do.call(rbind, response_streams)
    levels <- paste0('s', seq_len(subjects))
    impulses$subject <- factor(impulses$subject, levels=levels)
    responses$subject <- factor(responses$subject, levels=levels)
    list(impulses=impulses, responses=responses)
}

fit_replication <- function(seed) {
    streams <- simulate_location_scale_streams(seed)
    formulas <- list(
        location=response ~ s(subject, bs='re') +
            irf(x, window=c(0, 1.5), k_l=6) - irf(1),
        scale=~ s(subject, bs='re') +
            irf(z, window=c(0, 1.5), k_l=6) - irf(1)
    )
    design <- prepare_cdrgam(
        formulas,
        streams$impulses,
        streams$responses,
        series='subject',
        history='ragged',
        chunk_size=61,
        quiet=TRUE
    )
    elapsed <- numeric(2L)
    elapsed[[1L]] <- system.time({
        native <- cdrgam.fit(
            design, family='gaulss', backend='mgcv', engine='gam',
            method=method
        )
    })[['elapsed']]
    elapsed[[2L]] <- system.time({
        sparse <- cdrgam.fit(
            design, family='gaulss', backend='sparse', method=method,
            sparse_control=list(
                optimizer_maxit=500L,
                optimizer_gradient_tolerance=2e-3
            ),
            solver_trace=trace_solver
        )
    })[['elapsed']]
    fits <- list(mgcv=native, sparse=sparse)
    assemblies <- lapply(design$parameters, function(parameter_design) {
        getFromNamespace(
            '.cdrgam_distributional_sparse_setup', 'cdrgam'
        )(parameter_design, list(crossprod_chunk_size=61L))
    })
    native_state <- getFromNamespace(
        '.cdrgam_gaulss_sparse_fixed', 'cdrgam'
    )(
        assemblies,
        unname(native$sp),
        getFromNamespace('.cdrgam_gaulss_b', 'cdrgam')(
            cdrgam_family('gaulss')
        ),
        tolerance=1e-9,
        maxit=500L
    )
    initial_log_sp <- log(unlist(lapply(
        assemblies,
        getFromNamespace('.cdrgam_sparse_initial_sp', 'cdrgam')
    ), use.names=FALSE))
    native_state_criterion <- if (identical(method, 'QNCV')) {
        getFromNamespace('.cdrgam_gaulss_sparse_qncv', 'cdrgam')(
            assemblies, native_state
        )$criterion
    } else {
        counts <- vapply(
            assemblies,
            function(assembly) length(assembly$penalty_components),
            integer(1)
        )
        penalty_determinant <- 0
        offset <- 0L
        for (parameter in names(assemblies)) {
            indices <- offset + seq_len(counts[[parameter]])
            penalty_determinant <- penalty_determinant + getFromNamespace(
                '.sparse_penalty_logdet', 'cdrgam'
            )(
                assemblies[[parameter]]$blocks,
                unname(native$sp)[indices]
            )$value
            offset <- offset + counts[[parameter]]
        }
        2 * native_state$objective +
            getFromNamespace('.cdr_factor_logdet', 'cdrgam')(
                native_state$factor
            ) - penalty_determinant
    }
    sparse_criterion <- if (identical(method, 'QNCV')) {
        sparse$qncv
    } else sparse$reml
    newdata <- list(
        impulses=streams$impulses,
        responses=streams$responses[c('subject', 'time')]
    )
    predictions <- lapply(fits, predict, newdata=newdata, type='response')
    data.frame(
        seed=seed,
        converged=fits$sparse$converged,
        native_converged=identical(
            fits$mgcv$outer.info$conv, 'full convergence'
        ),
        native_seconds=elapsed[[1L]],
        sparse_seconds=elapsed[[2L]],
        inner_gradient=fits$sparse$distributional$gradient_norm,
        outer_gradient=max(abs(fits$sparse$optimizer$gradient)),
        native_outer_gradient=max(abs(fits$mgcv$outer.info$grad)),
        log_sp_max_difference=max(abs(
            log(unname(fits$sparse$sp)) - log(unname(fits$mgcv$sp))
        )),
        native_log_sp=paste(
            format(log(unname(fits$mgcv$sp)), digits=6), collapse=';'
        ),
        sparse_log_sp=paste(
            format(log(unname(fits$sparse$sp)), digits=6), collapse=';'
        ),
        initial_log_sp=paste(
            format(initial_log_sp, digits=6), collapse=';'
        ),
        fixed_state_coefficient_max_difference=max(abs(
            native_state$coefficients - unname(stats::coef(native))
        )),
        native_state_criterion_difference=
            native_state_criterion - sparse_criterion,
        condition_indicator=scalar_or_na(
            fits$sparse$sparse$condition_indicator
        ),
        coefficient_max_difference=max(
            abs(coef(fits$sparse) - coef(fits$mgcv))
        ),
        prediction_max_difference=max(
            abs(predictions$sparse - predictions$mgcv)
        ),
        log_likelihood_difference=abs(
            as.numeric(logLik(fits$sparse)) - as.numeric(logLik(fits$mgcv))
        )
    )
}

seed_text <- Sys.getenv('CDRGAM_VALIDATION_SEEDS', '1401,1402,1403')
seeds <- as.integer(strsplit(seed_text, ',', fixed=TRUE)[[1L]])
if (!length(seeds) || anyNA(seeds)) {
    stop('CDRGAM_VALIDATION_SEEDS must be comma-separated integers')
}
metrics <- do.call(rbind, lapply(seeds, fit_replication))
metrics$fixed_state_equivalent <-
    metrics$fixed_state_coefficient_max_difference < 2e-5
metrics$endpoint_agreement <-
    metrics$coefficient_max_difference < 2e-2 &
    metrics$prediction_max_difference < 3e-2 &
    metrics$log_likelihood_difference < 1e-2
metrics$classification <- ifelse(
    metrics$endpoint_agreement & metrics$converged &
        metrics$native_converged,
    'none',
    ifelse(
        metrics$fixed_state_equivalent & metrics$converged &
            metrics$native_converged,
        'optimizer reached a different stationary point',
        'unresolved'
    )
)
print(metrics, row.names=FALSE)
output <- Sys.getenv('CDRGAM_VALIDATION_OUTPUT')
if (nzchar(output)) {
    dir.create(dirname(output), recursive=TRUE, showWarnings=FALSE)
    utils::write.csv(metrics, output, row.names=FALSE)
}
stopifnot(
    all(metrics$native_converged),
    all(metrics$fixed_state_equivalent)
)
