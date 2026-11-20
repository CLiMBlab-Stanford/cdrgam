library(cdrgam)

refine <- cdrgam:::.sparse_newton_refinement
box <- list(lower=c(-25, -25), upper=c(25, 25))

# A convex quadratic is solved by one Newton step.
target <- c(1.5, -2)
curvature_matrix <- matrix(c(3, 0.5, 0.5, 1), 2)
quadratic <- function(x) {
    difference <- x - target
    drop(crossprod(difference, curvature_matrix %*% difference))
}
quadratic_curvature <- function(x) list(
    gradient=2 * drop(curvature_matrix %*% (x - target)),
    hessian=2 * curvature_matrix
)
solved <- refine(quadratic, quadratic_curvature, c(0, 0), quadratic(c(0, 0)),
    box$lower, box$upper)
stopifnot(
    identical(solved$status, 'converged'),
    solved$improved,
    max(abs(solved$log_sp - target)) < 1e-8,
    identical(solved$steps$step[[1L]], 'newton')
)

# At a saddle with zero gradient the phase follows negative curvature to a
# minimum of x^2 - y^2 + y^4, where y^2 = 1/2 and the value is -1/4.
saddle <- function(x) x[[1L]]^2 - x[[2L]]^2 + x[[2L]]^4
saddle_curvature <- function(x) list(
    gradient=c(2 * x[[1L]], -2 * x[[2L]] + 4 * x[[2L]]^3),
    hessian=diag(c(2, -2 + 12 * x[[2L]]^2))
)
escaped <- refine(saddle, saddle_curvature, c(0, 0), 0, box$lower, box$upper,
    max_steps=20L, tolerance=1e-12)
stopifnot(
    identical(escaped$status, 'converged'),
    identical(escaped$steps$step[[1L]], 'negative_curvature'),
    abs(escaped$value + 0.25) < 1e-10,
    abs(escaped$log_sp[[2L]]^2 - 0.5) < 1e-6
)

# Parameters reported by `boundary` stay fixed.
frozen <- refine(quadratic, quadratic_curvature, c(0, 0), quadratic(c(0, 0)),
    box$lower, box$upper, boundary=function(hessian, x) 2L)
stopifnot(
    frozen$log_sp[[2L]] == 0,
    abs(frozen$log_sp[[1L]] - (target[[1L]] + 0.5 * target[[2L]] / 3)) < 1e-8
)

# End to end: L-BFGS-B with finite-difference gradients is stopped after two
# iterations, and the Newton phase then reaches a positive-definite optimum
# with a lower criterion.
set.seed(1209)
subjects <- 8L
levels <- sprintf('s%02d', seq_len(subjects))
impulses <- do.call(rbind, lapply(seq_len(subjects), function(index) {
    data.frame(
        subject=factor(levels[[index]], levels=levels),
        time=sort(cumsum(stats::rexp(200, 8))),
        x=stats::rnorm(200)
    )
}))
amplitudes <- stats::rnorm(subjects, 0, 0.5)
responses <- do.call(rbind, lapply(split(impulses, impulses$subject), function(rows) {
    index <- as.integer(rows$subject[[1L]])
    y <- vapply(seq_len(nrow(rows)), function(j) {
        lag <- rows$time[[j]] - rows$time
        keep <- lag >= 0 & lag <= 0.9
        sum(rows$x[keep] * (1 + amplitudes[[index]]) * sin(pi * lag[keep] / 0.9))
    }, numeric(1))
    data.frame(subject=rows$subject, time=rows$time,
        y=y + stats::rnorm(nrow(rows), 0, 0.3))
}))
design <- prepare_cdrgam(
    y ~ irf(x, window=c(0, 0.9), k_l=6) +
        irf(x, window=c(0, 0.9), k_l=5, group=subject) - irf(1),
    impulses=impulses, responses=responses, series='subject', quiet=TRUE
)
control <- list(gradient='finite', outer_optimizer='lbfgsb',
    optimizer_maxit=2L, hessian='analytic', boundary_action='report')
plain <- cdrgam.fit(design, backend='sparse',
    sparse_control=c(control, list(final_newton=FALSE)))
refined <- cdrgam.fit(design, backend='sparse', sparse_control=control)
newton <- refined$sparse$convergence$newton_refinement
final_step <- utils::tail(newton$steps, 1L)
stopifnot(
    is.null(plain$sparse$convergence$newton_refinement),
    identical(newton$status, 'converged'),
    isTRUE(newton$improved),
    refined$reml < plain$reml - 1e-6,
    final_step$predicted_decrease < 1e-3,
    final_step$min_eigenvalue > 0,
    isTRUE(refined$sparse$convergence$hessian_positive_definite)
)

# The unrefined fit is diagnosed as stopped short; the refined one, which uses
# the default, is not.
reason_codes <- function(fit) {
    vapply(fit_diagnostics(fit)$reasons, `[[`, character(1), 'code')
}
stopifnot(
    plain$sparse$convergence$predicted_decrease > 0.1,
    'optimizer_stopped_short' %in% reason_codes(plain),
    refined$sparse$convergence$predicted_decrease < 1e-3,
    !('optimizer_stopped_short' %in% reason_codes(refined))
)

# A fit reports the Newton step size and the resolvable criterion change
# beside the predicted decrease.
refined_diagnostics <- fit_diagnostics(refined)
stopifnot(
    is.finite(refined_diagnostics$newton_step_max),
    refined_diagnostics$newton_step_max >= 0,
    is.finite(refined_diagnostics$predicted_decrease_floor),
    refined_diagnostics$predicted_decrease_floor >= 1e-8
)

# A large gradient along a well-curved direction gives a negligible Newton
# step and a predicted decrease below the floor.
newton <- cdrgam:::.sparse_newton_summary(
    diag(c(45, 1)), c(0.033, 0), 3.57e7
)
stopifnot(
    isTRUE(all.equal(newton$predicted_decrease, 0.5 * 0.033^2 / 45)),
    isTRUE(all.equal(newton$newton_step_max, 0.033 / 45)),
    isTRUE(all.equal(newton$predicted_decrease_floor, 1e-10 * (1 + 3.57e7))),
    newton$predicted_decrease < newton$predicted_decrease_floor
)
