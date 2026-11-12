library(cdrgam)

# Generalized REML fits difference their exact score to obtain an outer
# Hessian, which supports boundary classification, optimizer checks, and
# unconditional inference. `noise` has no effect, so its penalty is expected
# to reach the boundary while the `signal` penalty stays interior.
simulation <- simulate_cdr(
    list(signal=function(lag) exp(-2 * lag)),
    n_impulses=400,
    n_responses=360,
    duration=60,
    window=1.5,
    seed=4207
)
set.seed(4208)
simulation$impulses$noise <- stats::rnorm(nrow(simulation$impulses))
latent <- as.numeric(scale(simulation$responses$response))
simulation$responses$count <- stats::rpois(
    nrow(simulation$responses), exp(0.3 + 0.6 * latent)
)
design <- prepare_cdrgam(
    count ~ irf(signal, window=c(0, 1.5), k_l=6) +
        irf(noise, window=c(0, 1.5), k_l=5) - irf(1),
    impulses=simulation$impulses,
    responses=simulation$responses,
    quiet=TRUE
)
fit <- suppressWarnings(cdrgam.fit(design, family=poisson, backend='sparse',
    method='REML'))
convergence <- fit$sparse$convergence
diagnostics <- fit_diagnostics(fit)
hessian <- fit$outer.info$hess
stopifnot(
    is.matrix(hessian),
    nrow(hessian) == length(fit$sp),
    identical(convergence$hessian_method, 'gradient'),
    is.logical(convergence$hessian_positive_definite),
    is.data.frame(convergence$hessian_boundary),
    is.finite(convergence$predicted_decrease),
    convergence$predicted_decrease < 0.1,
    !('optimizer_stopped_short' %in%
        vapply(diagnostics$reasons, `[[`, character(1), 'code'))
)
if (isTRUE(convergence$hessian_positive_definite)) {
    unconditional <- vcov(fit, unconditional=TRUE)
    conditional <- vcov(fit, unconditional=FALSE)
    stopifnot(
        all(is.finite(unconditional)),
        all(diag(unconditional) >= diag(conditional) - 1e-10)
    )
}
boundary <- convergence$hessian_boundary
stopifnot(
    isTRUE(convergence$hessian_positive_definite),
    nrow(boundary) == 1L,
    identical(boundary$smoothing_parameter, 's(cdr_term_2)'),
    identical(boundary$status, 'bounded'),
    boundary$bound_log_sp < boundary$log_sp,
    identical(vapply(diagnostics$notes, `[[`, character(1), 'code'),
        'boundary_variance_components')
)

# hessian = "none" skips the curvature work.
skipped <- suppressWarnings(cdrgam.fit(design, family=poisson,
    backend='sparse', method='REML', sparse_control=list(hessian='none')))
stopifnot(
    is.null(skipped$outer.info$hess),
    is.na(skipped$sparse$convergence$hessian_positive_definite)
)
