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
    isTRUE(constrained_plan$memory_limited),
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
stopifnot(
    abs(native_criterion$criterion - unname(native$gcv.ubre)) < 2e-5,
    abs(native_criterion$criterion - small_batch$criterion) < 1e-10,
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
analytic <- getFromNamespace(
    '.cdrgam_gaulss_sparse_qncv_score', 'cdrgam'
)(
    assemblies,
    interior_solution,
    interior_sp,
    family,
    gamma=gamma,
    batch_size=23L
)
step <- 1e-4
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
stopifnot(max(abs(analytic - numeric)) < 1e-3)

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
