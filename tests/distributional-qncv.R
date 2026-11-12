library(cdrgam)

set.seed(9461)
n_impulses <- 100L
n_responses <- 120L
impulses <- data.frame(
    time=sort(stats::runif(n_impulses, 0, 30)),
    location_signal=stats::rnorm(n_impulses),
    scale_signal=stats::rnorm(n_impulses)
)
responses <- data.frame(time=sort(stats::runif(n_responses, 0.5, 30)))
history_effect <- function(variable) vapply(responses$time, function(time) {
    lag <- time - impulses$time
    rows <- lag >= 0 & lag <= 1.5
    sum(impulses[[variable]][rows] * exp(-2 * lag[rows]))
}, numeric(1))
location_effect <- history_effect('location_signal')
scale_effect <- history_effect('scale_signal')
responses$response <- stats::rnorm(
    n_responses,
    0.7 * location_effect,
    exp(-0.3 + 0.18 * scale_effect)
)
formulas <- list(
    location=response ~
        irf(location_signal, window=c(0, 1.5), k_l=5) - irf(1),
    scale=~ irf(scale_signal, window=c(0, 1.5), k_l=5) - irf(1)
)
design <- prepare_cdrgam(
    formulas,
    impulses,
    responses,
    history='ragged',
    chunk_size=37L,
    quiet=TRUE
)
family <- cdrgam_family('gaulss')
gamma <- 1.1
native <- cdrgam.fit(
    design,
    family=family,
    backend='mgcv',
    method='QNCV',
    gamma=gamma
)
assemblies <- lapply(design$parameters, function(parameter_design) {
    getFromNamespace('.cdrgam_distributional_sparse_setup', 'cdrgam')(
        parameter_design,
        list(crossprod_chunk_size=37L)
    )
})
direction_plan <- getFromNamespace(
    '.cdrgam_gaulss_qncv_direction_plan', 'cdrgam'
)
constrained_plan <- direction_plan(
    assemblies,
    direction_count=2L,
    cores=4L,
    workers=4L,
    memory=list(source='synthetic-low', available_bytes=1)
)
roomy_plan <- direction_plan(
    assemblies,
    direction_count=2L,
    cores=4L,
    workers=4L,
    memory=list(source='synthetic-high', available_bytes=1e12)
)
stopifnot(
    constrained_plan$workers == 1L,
    constrained_plan$blas_threads == 4L,
    isTRUE(constrained_plan$memory_limited) ==
        (.Platform$OS.type != 'windows'),
    !isTRUE(constrained_plan$cache_parallel),
    roomy_plan$workers == if (.Platform$OS.type == 'windows') 1L else 2L,
    roomy_plan$blas_threads == if (
        .Platform$OS.type == 'windows'
    ) 4L else 2L
)
fit_fixed <- getFromNamespace('.cdrgam_gaulss_sparse_fixed', 'cdrgam')
evaluate_qncv <- getFromNamespace(
    '.cdrgam_gaulss_sparse_qncv', 'cdrgam'
)
b <- getFromNamespace('.cdrgam_gaulss_b', 'cdrgam')(family)
native_solution <- fit_fixed(
    assemblies,
    unname(native$sp),
    b,
    tolerance=1e-9,
    maxit=500L
)
evaluate_moments <- getFromNamespace(
    '.cdrgam_gaulss_sparse_moments', 'cdrgam'
)
moment_arguments <- list(
    assemblies=assemblies,
    coefficients=native_solution$parameter_coefficients,
    sigma=native_solution$sigma,
    q=exp(native_solution$linear_predictors[, 'scale']),
    residual=native_solution$residuals,
    weights=native_solution$prior_weights,
    include_observed=TRUE
)
serial_moments <- do.call(evaluate_moments, c(
    moment_arguments, list(workers=1L)
))
parallel_moments <- do.call(evaluate_moments, c(
    moment_arguments, list(workers=2L, blas_threads=1L)
))
stopifnot(
    max(abs(
        serial_moments$gradients$location -
            parallel_moments$gradients$location
    )) < 1e-10,
    max(abs(
        serial_moments$gradients$scale -
            parallel_moments$gradients$scale
    )) < 1e-10,
    max(abs(
        serial_moments$information$location -
            parallel_moments$information$location
    )) < 1e-10,
    max(abs(
        serial_moments$information$scale -
            parallel_moments$information$scale
    )) < 1e-10,
    max(abs(
        serial_moments$information$cross -
            parallel_moments$information$cross
    )) < 1e-10
)
native_criterion <- evaluate_qncv(
    assemblies,
    native_solution,
    gamma=gamma,
    batch_size=23L
)
small_batch <- evaluate_qncv(
    assemblies,
    native_solution,
    gamma=gamma,
    batch_size=7L
)
parallel_batch <- evaluate_qncv(
    assemblies,
    native_solution,
    gamma=gamma,
    batch_size=7L,
    workers=2L,
    blas_threads=1L
)
partial_qncv <- NULL
qncv_interrupted <- tryCatch(
    evaluate_qncv(
        assemblies,
        native_solution,
        gamma=gamma,
        batch_size=7L,
        checkpoint=function(state) {
            partial_qncv <<- state
            stop('simulated QNCV response interruption')
        }
    ),
    error=function(error) conditionMessage(error)
)
qncv_resume_updates <- 0L
resumed_batch <- evaluate_qncv(
    assemblies,
    native_solution,
    gamma=gamma,
    batch_size=7L,
    resume=partial_qncv,
    checkpoint=function(state) {
        qncv_resume_updates <<- qncv_resume_updates + 1L
    }
)
criterion_scale <- max(
    1,
    abs(native_criterion$criterion),
    abs(unname(native$gcv.ubre))
)
stopifnot(
    identical(qncv_interrupted, 'simulated QNCV response interruption'),
    length(partial_qncv$results) > 0L,
    length(unique(unlist(lapply(
        partial_qncv$results, `[[`, 'starts'
    )))) < length(seq.int(1L, n_responses, by=7L)),
    qncv_resume_updates > 0L,
    abs(resumed_batch$criterion - small_batch$criterion) < 1e-10,
    max(abs(
        resumed_batch$predictor_change - small_batch$predictor_change
    )) < 1e-12,
    abs(native_criterion$criterion - unname(native$gcv.ubre)) /
        criterion_scale < 1e-6,
    abs(native_criterion$criterion - small_batch$criterion) < 1e-10,
    abs(native_criterion$criterion - parallel_batch$criterion) < 1e-10,
    max(abs(
        native_criterion$predictor_change -
            parallel_batch$predictor_change
    )) < 1e-12,
    max(abs(
        unname(native_solution$coefficients) - unname(stats::coef(native))
    )) < 2e-5,
    native_criterion$minimum_determinant > 0
)

interior_sp <- unname(native$sp) * exp(c(-2, -10))
interior_solution <- fit_fixed(
    assemblies,
    interior_sp,
    b,
    tolerance=1e-10,
    maxit=500L
)
score_qncv <- getFromNamespace(
    '.cdrgam_gaulss_sparse_qncv_score', 'cdrgam'
)
analytic <- score_qncv(
    assemblies,
    interior_solution,
    interior_sp,
    family,
    gamma=gamma,
    batch_size=23L
)
analytic_parallel <- score_qncv(
    assemblies,
    interior_solution,
    interior_sp,
    family,
    gamma=gamma,
    batch_size=23L,
    cores=2L,
    workers=2L
)
partial_score <- NULL
score_interrupted <- tryCatch(
    score_qncv(
        assemblies,
        interior_solution,
        interior_sp,
        family,
        gamma=gamma,
        batch_size=7L,
        cores=1L,
        workers=1L,
        score_batch_size=1L,
        checkpoint=function(state) {
            partial_score <<- state
            if (is.list(state$active) &&
                    length(state$active$qncv$results) > 0L) {
                stop('simulated exact-score interruption')
            }
        }
    ),
    error=function(error) conditionMessage(error)
)
resumed_score_updates <- 0L
resumed_score <- score_qncv(
    assemblies,
    interior_solution,
    interior_sp,
    family,
    gamma=gamma,
    batch_size=7L,
    cores=1L,
    workers=1L,
    score_batch_size=1L,
    resume=partial_score,
    checkpoint=function(state) {
        resumed_score_updates <<- resumed_score_updates + 1L
    }
)
step <- 1e-3
numeric <- vapply(seq_along(interior_sp), function(index) {
    lower <- upper <- log(interior_sp)
    lower[[index]] <- lower[[index]] - step
    upper[[index]] <- upper[[index]] + step
    criterion <- function(log_sp) {
        solution <- fit_fixed(
            assemblies,
            exp(log_sp),
            b,
            tolerance=1e-10,
            maxit=500L
        )
        evaluate_qncv(
            assemblies,
            solution,
            gamma=gamma,
            batch_size=23L
        )$criterion
    }
    (criterion(upper) - criterion(lower)) / (2 * step)
}, numeric(1))
gradient_scale <- max(1, abs(analytic), abs(numeric))
stopifnot(
    grepl(
        'simulated exact-score interruption',
        score_interrupted,
        fixed=TRUE
    ),
    is.list(partial_score$active$qncv),
    length(partial_score$active$qncv$results) > 0L,
    resumed_score_updates > 0L,
    max(abs(resumed_score - analytic)) / gradient_scale < 1e-10,
    max(abs(analytic - numeric)) / gradient_scale < 5e-4,
    max(abs(analytic_parallel - analytic)) / gradient_scale < 1e-10
)

progress <- list()
checkpoint <- tempfile(fileext='.rds')
sparse <- cdrgam.fit(
    design,
    family=family,
    backend='sparse',
    method='QNCV',
    gamma=gamma,
    checkpoint=checkpoint,
    solver_trace=function(record) {
        progress[[length(progress) + 1L]] <<- record
    },
    sparse_control=list(
        crossprod_chunk_size=37L,
        qncv_batch_size=23L,
        optimizer_gradient_tolerance=2e-3,
        cores=2L
    )
)
resumed <- cdrgam.fit(
    design,
    family=family,
    backend='sparse',
    method='QNCV',
    gamma=gamma,
    checkpoint=checkpoint,
    sparse_control=list(
        crossprod_chunk_size=37L,
        qncv_batch_size=23L,
        optimizer_gradient_tolerance=2e-3,
        cores=2L
    )
)
unlink(checkpoint)
sparse_summary <- summary(sparse)
sparse_diagnostics <- fit_diagnostics(sparse)
continued <- cdrgam.fit(
    design,
    family=family,
    backend='sparse',
    method='QNCV',
    gamma=gamma,
    sparse_control=list(
        crossprod_chunk_size=37L,
        qncv_batch_size=23L,
        inner_maxit=1L,
        inner_adaptive_maxit=500L,
        optimizer_maxit=1L,
        optimizer_gradient_tolerance=1e6,
        cores=1L
    )
)
resume_coefficient_difference <- max(abs(
    stats::coef(resumed) - stats::coef(sparse)
))
resume_criterion_difference <- abs(resumed$qncv - sparse$qncv)
stopifnot(
    isTRUE(sparse$converged),
    identical(sparse$method, 'QNCV'),
    identical(sparse_summary$method, 'QNCV'),
    abs(sparse_summary$sp.criterion - sparse$qncv) < 1e-12,
    identical(sparse$sparse$control$criterion, 'QNCV'),
    identical(sparse$sparse$control$gradient, 'exact'),
    identical(sparse$sparse$control$gamma, gamma),
    identical(sparse$sparse$control$qncv_batch_size, 23L),
    sparse$distributional$stochastic_score_evaluations == 0L,
    sparse$distributional$score_evaluations > 0L,
    continued$distributional$inner_continuations > 0L,
    sparse$distributional$qncv$minimum_determinant > 0,
    resume_coefficient_difference < 1e-4,
    resume_criterion_difference < 1e-5,
    identical(
        sparse$distributional$score_plan$response_batch_size,
        23L
    ),
    sparse$distributional$score_plan$workers >= 1L,
    sparse$distributional$score_plan$blas_threads >= 1L,
    is.logical(sparse$distributional$score_plan$memory_limited),
    abs(sparse$qncv - unname(native$gcv.ubre)) < 2e-3,
    max(abs(unname(stats::coef(sparse)) - unname(stats::coef(native)))) <
        5e-4,
    isTRUE(sparse_diagnostics$converged),
    sparse_diagnostics$gradient_norm <= 2e-3,
    any(vapply(progress, function(record) {
        identical(record$event, 'outer exact score complete') &&
            is.finite(suppressWarnings(as.numeric(record$criterion)))
    }, logical(1)))
)
