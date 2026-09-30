library(cdrgam)

simulation <- simulate_cdr(
    list(signal=function(lag) 0.8 * exp(-2.2 * lag)),
    n_impulses=140,
    n_responses=180,
    duration=28,
    window=1.5,
    seed=9931
)
latent <- as.numeric(scale(simulation$responses$response))
simulation$responses$offset <- seq(
    -0.15,
    0.15,
    length.out=nrow(simulation$responses)
)
design_for <- function(response) {
    responses <- simulation$responses
    responses$response <- response
    prepare_cdrgam(
        response ~ offset(offset) +
            irf(signal, window=c(0, 1.5), k_l=7) - irf(1),
        simulation$impulses,
        responses,
        history='ragged',
        quiet=TRUE
    )
}

compare_family <- function(response, family, gamma=1) {
    design <- design_for(response)
    weights <- rep(c(1, 2), length.out=length(response))
    native <- cdrgam.fit(
        design,
        family=family,
        backend='mgcv',
        engine='gam',
        method='GCV.Cp',
        gamma=gamma,
        weights=weights
    )
    sparse <- cdrgam.fit(
        design,
        family=family,
        backend='sparse',
        method='GCV.Cp',
        gamma=gamma,
        weights=weights,
        sparse_control=list(
            crossprod_chunk_size=29L,
            optimizer_maxit=100L,
            optimizer_gradient_tolerance=2e-4,
            finite_difference_step=5e-4,
            cores=2L,
            score_workers=2L
        )
    )
    native_edf <- sum(native$edf)
    sparse_edf <- getFromNamespace('.sparse_effective_df', 'cdrgam')(sparse)
    relative_score_error <- abs(
        unname(sparse$gcv.ubre) - unname(native$gcv.ubre)
    ) / max(1, abs(unname(native$gcv.ubre)))
    predictor_error <- max(abs(
        sparse$linear.predictors - native$linear.predictors
    ))
    assembly <- getFromNamespace('.fit_sparse_gaussian', 'cdrgam')(
        design,
        family,
        method='REML',
        sparse_control=list(crossprod_chunk_size=29L),
        setup_only=TRUE,
        weights=weights
    )
    at_native <- getFromNamespace(
        '.cdrgam_streamed_sparse_prediction_error', 'cdrgam'
    )(
        assembly,
        family,
        log(native$sp),
        gamma=gamma,
        chunk_size=29L
    )
    gradient_point <- rep.int(
        log(0.8), length(assembly$penalty_components)
    )
    gradient_evaluate <- function(parameters) getFromNamespace(
        '.cdrgam_streamed_sparse_prediction_error', 'cdrgam'
    )(
        assembly, family, parameters, gamma=gamma,
        tolerance=1e-11, chunk_size=29L
    )
    gradient_base <- gradient_evaluate(gradient_point)
    analytic_gradient <- getFromNamespace(
        '.cdrgam_generalized_exact_prediction_error_gradient', 'cdrgam'
    )(
        assembly, family, gradient_base, gamma=gamma,
        chunk_size=29L, inverse_chunk_size=4L
    )
    gradient_step <- 1e-4
    numeric_gradient <- vapply(seq_along(gradient_point), function(index) {
        lower <- upper <- gradient_point
        lower[[index]] <- lower[[index]] - gradient_step
        upper[[index]] <- upper[[index]] + gradient_step
        (gradient_evaluate(upper)$criterion -
            gradient_evaluate(lower)$criterion) / (2 * gradient_step)
    }, numeric(1))
    stopifnot(
        isTRUE(sparse$converged),
        identical(sparse$method, 'GCV.Cp'),
        identical(sparse$sparse$control$gradient, 'exact'),
        identical(sparse$sparse$control$gradient_requested, 'auto'),
        identical(sparse$sparse$control$outer_optimizer, 'lbfgsb'),
        sparse$sparse$control$score_workers >= 1L,
        sparse$sparse$control$score_blas_threads >= 1L,
        is.logical(sparse$sparse$control$score_memory_limited),
        is.finite(sparse$sparse$control$score_worker_bytes),
        identical(
            sparse$sparse$score_plan$workers,
            sparse$sparse$control$score_workers
        ),
        identical(sparse$sparse$gradient_selection$method, 'exact'),
        identical(
            sparse$sparse$control$criterion,
            if (identical(family$family, 'Gamma')) 'GCV' else 'UBRE'
        ),
        relative_score_error < 2e-4,
        abs(at_native$criterion - unname(native$gcv.ubre)) < 2e-6,
        max(abs(analytic_gradient - numeric_gradient)) < 1e-6,
        abs(sparse_edf - native_edf) < 0.02,
        predictor_error < 0.02,
        identical(
            summary(sparse)$method,
            if (identical(family$family, 'Gamma')) 'GCV' else 'UBRE'
        )
    )
    if (identical(family$family, 'Gamma')) {
        stopifnot(abs(log(sparse$scale / native$scale)) < 0.02)
    } else stopifnot(identical(sparse$scale, 1))
    invisible(sparse)
}

set.seed(9932)
compare_family(
    rpois(length(latent), exp(0.5 + 0.3 * latent)),
    poisson(link='log')
)

set.seed(9933)
compare_family(
    rbinom(length(latent), 1, plogis(-0.1 + 0.7 * latent)),
    binomial(link='logit'),
    gamma=1.2
)

set.seed(9934)
mean <- exp(0.8 + 0.25 * latent)
compare_family(
    rgamma(length(latent), shape=6, scale=mean / 6),
    Gamma(link='log')
)

set.seed(9935)
multi_impulses <- simulation$impulses
multi_impulses$signal2 <- stats::rnorm(nrow(multi_impulses))
multi_responses <- simulation$responses
multi_responses$response <- stats::rpois(
    nrow(multi_responses), exp(0.4 + 0.25 * latent)
)
multi_design <- prepare_cdrgam(
    response ~ offset(offset) +
        irf(signal, window=c(0, 1.5), k_l=7) +
        irf(signal2, window=c(0, 1.5), k_l=6) - irf(1),
    multi_impulses, multi_responses, history='ragged', quiet=TRUE
)
multi_assembly <- getFromNamespace('.fit_sparse_gaussian', 'cdrgam')(
    multi_design, poisson(link='log'), method='REML',
    sparse_control=list(crossprod_chunk_size=29L), setup_only=TRUE
)
multi_parameters <- log(seq(
    0.7, 1.4, length.out=length(multi_assembly$penalty_components)
))
multi_evaluate <- function(parameters) getFromNamespace(
    '.cdrgam_streamed_sparse_prediction_error', 'cdrgam'
)(
    multi_assembly, poisson(link='log'), parameters, gamma=1.1,
    tolerance=1e-11, chunk_size=29L
)
multi_base <- multi_evaluate(multi_parameters)
multi_analytic <- getFromNamespace(
    '.cdrgam_generalized_exact_prediction_error_gradient', 'cdrgam'
)(
    multi_assembly, poisson(link='log'), multi_base, gamma=1.1,
    chunk_size=29L, batch_size=1L, inverse_chunk_size=4L
)
gradient_step <- 1e-4
multi_numeric <- vapply(seq_along(multi_parameters), function(index) {
    lower <- upper <- multi_parameters
    lower[[index]] <- lower[[index]] - gradient_step
    upper[[index]] <- upper[[index]] + gradient_step
    (multi_evaluate(upper)$criterion - multi_evaluate(lower)$criterion) /
        (2 * gradient_step)
}, numeric(1))
selection <- getFromNamespace(
    '.cdrgam_generalized_gcv_gradient_selection', 'cdrgam'
)
low_memory <- selection(
    'auto', multi_assembly,
    memory=list(
        source='synthetic-low', limit_bytes=1, used_bytes=0,
        available_bytes=1
    )
)
high_memory <- selection(
    'auto', multi_assembly,
    memory=list(
        source='synthetic-high', limit_bytes=1e12, used_bytes=0,
        available_bytes=1e12
    )
)
stopifnot(
    length(multi_parameters) > 1L,
    max(abs(multi_analytic - multi_numeric)) < 1e-6,
    identical(low_memory$method, 'finite'),
    identical(high_memory$method, 'exact')
)

finite_fit <- cdrgam.fit(
    multi_design,
    family=poisson(link='log'),
    backend='sparse',
    method='GCV.Cp',
    gamma=1.1,
    sparse_control=list(
        gradient='finite', crossprod_chunk_size=29L,
        optimizer_maxit=100L, optimizer_gradient_tolerance=2e-4,
        finite_difference_step=5e-4, cores=2L, score_workers=2L
    )
)
stopifnot(
    isTRUE(finite_fit$converged),
    identical(finite_fit$sparse$control$gradient, 'finite'),
    identical(finite_fit$sparse$control$gradient_requested, 'finite')
)
