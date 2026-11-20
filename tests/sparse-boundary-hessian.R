library(cdrgam)

classify <- cdrgam:::.sparse_outer_hessian_boundary
profile <- cdrgam:::.sparse_boundary_profile

# A parameter is on the boundary only when its whole outer-Hessian row is
# negligible at an extreme log smoothing parameter.
hessian <- diag(c(2, 3, 1e-4, 1))
stopifnot(
    identical(
        classify(hessian, c(1, 2, 15), curvature_tolerance=0.05,
            log_sp_threshold=8),
        3L
    ),
    identical(
        classify(hessian, c(1, 2, 5), curvature_tolerance=0.05,
            log_sp_threshold=8),
        integer()
    )
)

# Two parameters trading off along an interior ridge make the Hessian
# singular, but their rows are large, so they stay unclassified.
ridge <- matrix(c(4, 4, 0, 4, 4, 0, 0, 0, 1), 3, 3)
stopifnot(
    min(eigen(ridge, only.values=TRUE)$values) < 1e-12,
    identical(
        classify(ridge, c(15, 15), curvature_tolerance=0.05,
            log_sp_threshold=8),
        integer()
    )
)

# The profile search finds where the criterion rises by the requested amount.
# For 10 * exp(-rho) the exact crossing is log(10 / increase).
increase <- stats::qchisq(0.9, 1)
flat_tail <- function(log_sp) 10 * exp(-log_sp[[1L]])
bounded <- profile(flat_tail, c(20, 0), 1L, base_value=flat_tail(20),
    increase=increase)
exact <- log(10 / (increase + flat_tail(20)))
stopifnot(
    identical(bounded$status, 'bounded'),
    bounded$bound_log_sp >= exact,
    bounded$bound_log_sp - exact < 0.3,
    bounded$criterion_increase <= increase
)
unbounded <- profile(function(log_sp) 0, c(20, 0), 1L, base_value=0,
    increase=increase)
stopifnot(identical(unbounded$status, 'unbounded'), is.na(unbounded$bound_log_sp))

# End to end: subjects deviate from a curved population IRF only by lines in
# lag, so REML drives the grouped lag-smoothness penalty to its upper
# boundary while the grouped ridge stays interior.
set.seed(4021)
subjects <- 10L
per_subject <- 240L
levels <- sprintf('s%02d', seq_len(subjects))
impulses <- do.call(rbind, lapply(seq_len(subjects), function(index) {
    data.frame(
        subject=factor(levels[[index]], levels=levels),
        time=sort(cumsum(stats::rexp(per_subject, 8))),
        x=stats::rnorm(per_subject)
    )
}))
intercepts <- stats::rnorm(subjects, 0, 0.6)
slopes <- stats::rnorm(subjects, 0, 0.8)
responses <- do.call(rbind, lapply(split(impulses, impulses$subject), function(rows) {
    index <- as.integer(rows$subject[[1L]])
    y <- vapply(seq_len(nrow(rows)), function(j) {
        lag <- rows$time[[j]] - rows$time
        keep <- lag >= 0 & lag <= 0.9
        sum(rows$x[keep] * (sin(pi * lag[keep] / 0.9) +
            intercepts[[index]] + slopes[[index]] * lag[keep]))
    }, numeric(1))
    data.frame(subject=rows$subject, time=rows$time,
        y=y + stats::rnorm(nrow(rows), 0, 0.3))
}))
design <- prepare_cdrgam(
    y ~ irf(x, window=c(0, 0.9), k_l=6) +
        irf(x, window=c(0, 0.9), k_l=6, group=subject) - irf(1),
    impulses=impulses,
    responses=responses,
    series='subject',
    quiet=TRUE
)
control <- list(hessian='analytic', boundary_action='report')
fit <- cdrgam.fit(design, backend='sparse',
    sparse_control=c(control, list(boundary_sensitivity=TRUE)))
convergence <- fit$sparse$convergence
boundary <- convergence$hessian_boundary
diagnostics <- fit_diagnostics(fit)
stopifnot(
    nrow(boundary) == 1L,
    identical(boundary$smoothing_parameter, 's(cdr_term_2)1'),
    identical(boundary$side, 'upper'),
    identical(boundary$status, 'bounded'),
    boundary$bound_log_sp < boundary$log_sp,
    boundary$criterion_increase <= increase,
    isTRUE(convergence$hessian_positive_definite),
    convergence$hessian_full_min_eigenvalue < convergence$hessian_min_eigenvalue,
    identical(diagnostics$outcome, 'success'),
    fit$sparse$condition_indicator > diagnostics$condition_tolerance,
    fit$sparse$condition_indicator_full <= fit$sparse$condition_indicator,
    identical(
        vapply(diagnostics$notes, `[[`, character(1), 'code'),
        'boundary_variance_components'
    )
)

# Unconditional covariance inverts only the interior Hessian block.
unconditional <- vcov(fit, unconditional=TRUE)
stopifnot(all(is.finite(unconditional)))

# An indefinite interior block is refused instead of being inverted.
indefinite <- fit
interior_index <- match('s(cdr_term_1)', names(fit$sp))
indefinite$outer.info$hess[interior_index, interior_index] <- -1
indefinite_error <- tryCatch(
    vcov(indefinite, unconditional=TRUE),
    error=conditionMessage
)
stopifnot(
    is.character(indefinite_error),
    grepl('positive-definite outer Hessian', indefinite_error, fixed=TRUE),
    all(is.finite(vcov(indefinite, unconditional=FALSE)))
)

# boundary = "bound" reuses the stored fit at the profile bound. Its system
# differs from the fitted one by exactly the change in the boundary penalty.
axes <- list(
    lag=list(grid='fitted', n=11L),
    predictors=list(x=list(at=list(summary='mean', `offset-sd`=1)))
)
deviation <- function(boundary) estimate_effect(
    fit, terms='x|subject', axes=axes, grouping='deviation',
    groups=levels[[1L]], boundary=boundary
)
at_estimate <- deviation('estimate')
at_bound <- deviation('bound')
bound_object <- cdrgam:::.sparse_boundary_bound_object(fit)
indices <- seq_along(coef(fit))
system_at <- function(object) solve(
    cdrgam:::.sparse_selected_vcov(object, indices) / object$scale
)
penalty_at <- function(sp) as.matrix(Reduce(
    `+`, Map(`*`, sp, fit$sparse$penalty_components)
))
system_change <- system_at(fit) - system_at(bound_object)
penalty_change <- penalty_at(fit$sp) - penalty_at(bound_object$sp)
stopifnot(
    all(is.finite(at_bound$se)),
    !isTRUE(all.equal(at_bound$estimate, at_estimate$estimate)),
    max(abs(system_change - penalty_change)) <
        1e-4 * max(abs(penalty_change)),
    isTRUE(all.equal(
        log(bound_object$sp[['s(cdr_term_2)1']]),
        boundary$bound_log_sp
    )),
    is.data.frame(test_effect(fit, terms=list(predictors='x'),
        axes=list(lag=list(grid='fitted', n=7L)), boundary='bound'))
)

# Without the stored sensitivity fit, or with unconditional inference, the
# bound mode stops with an explanation.
plain <- cdrgam.fit(design, backend='sparse', sparse_control=control)
plain_error <- tryCatch(
    estimate_effect(plain, terms=list(predictors='x'), axes=axes,
        boundary='bound'),
    error=conditionMessage
)
unconditional_error <- tryCatch(
    estimate_effect(fit, terms=list(predictors='x'), axes=axes,
        boundary='bound', unconditional=TRUE),
    error=conditionMessage
)
stopifnot(
    is.null(plain$sparse$boundary_sensitivity),
    nrow(plain$sparse$convergence$hessian_boundary) == 1L,
    grepl('boundary_sensitivity = TRUE', plain_error, fixed=TRUE),
    grepl('unconditional = FALSE', unconditional_error, fixed=TRUE)
)
