library(cdrgam)

healthy <- structure(list(
    converged=TRUE,
    sparse=list(
        convergence=list(
            converged=TRUE,
            code=0L,
            message='full convergence',
            boundary=numeric(),
            hessian_positive_definite=TRUE
        ),
        condition_indicator=0.5
    ),
    cdrgam=list(rank=list(resolution='none', tolerance=1e-8))
), class=c('cdrgam_sparse', 'cdrgam'))

healthy_diagnostics <- fit_diagnostics(healthy)
stopifnot(
    identical(healthy_diagnostics$outcome, 'success'),
    length(healthy_diagnostics$reasons) == 0L,
    identical(healthy_diagnostics$converged, TRUE)
)

indefinite <- healthy
indefinite$sparse$convergence$hessian_positive_definite <- FALSE
indefinite_diagnostics <- fit_diagnostics(indefinite)
stopifnot(
    identical(indefinite_diagnostics$outcome, 'warning'),
    identical(
        vapply(indefinite_diagnostics$reasons, `[[`, character(1), 'code'),
        'indefinite_outer_hessian'
    )
)

boundary <- healthy
boundary$sparse$convergence$boundary <- c(term='upper')
boundary_diagnostics <- fit_diagnostics(boundary)
stopifnot(
    identical(boundary_diagnostics$outcome, 'warning'),
    'smoothing_parameter_boundary' %in% vapply(
        boundary_diagnostics$reasons, `[[`, character(1), 'code'
    )
)

ill_conditioned <- healthy
ill_conditioned$sparse$condition_indicator <- 1e-20
ill_conditioned_diagnostics <- fit_diagnostics(ill_conditioned)
stopifnot(
    identical(ill_conditioned_diagnostics$outcome, 'warning'),
    isTRUE(all.equal(
        ill_conditioned_diagnostics$condition_tolerance, 1e-16
    )),
    'ill_conditioned_penalized_system' %in% vapply(
        ill_conditioned_diagnostics$reasons, `[[`, character(1), 'code'
    )
)

regularized <- healthy
regularized$cdrgam$rank$resolution <- 'minimum-norm-ridge-approximation'
regularized_diagnostics <- fit_diagnostics(regularized)
stopifnot(
    identical(regularized_diagnostics$outcome, 'warning'),
    'rank_deficiency_regularized' %in% vapply(
        regularized_diagnostics$reasons, `[[`, character(1), 'code'
    )
)

nonconverged <- indefinite
nonconverged$converged <- FALSE
nonconverged$sparse$convergence$converged <- FALSE
nonconverged$sparse$convergence$code <- 1L
nonconverged$sparse$convergence$message <- 'iteration limit reached'
nonconverged_diagnostics <- fit_diagnostics(nonconverged)
stopifnot(
    identical(nonconverged_diagnostics$outcome, 'nonconverged'),
    'optimizer_nonconvergence' %in% vapply(
        nonconverged_diagnostics$reasons, `[[`, character(1), 'code'
    )
)

distributional <- structure(list(
    converged=TRUE,
    optimizer=list(
        convergence=0L, message='full convergence',
        projected_gradient=c(location=0, scale=0), gradient=numeric()
    ),
    distributional=list(evaluations=4L, gradient_norm=0),
    sparse=list(condition_indicator=0.5),
    cdrgam=list(rank=list(resolution='identified', tolerance=1e-8))
), class=c('cdrgam_distributional_sparse', 'cdrgam_sparse', 'cdrgam'))
stopifnot(identical(fit_diagnostics(distributional)$outcome, 'success'))

block <- structure(list(
    optimizer=list(convergence=0L, message='full convergence'),
    outer.info=list(hess=diag(2)),
    cdrgam=list(
        identifiability=list(global=list(condition_indicator=0.5)),
        rank=list(resolution='none', tolerance=1e-8)
    )
), class=c('cdrgam_block', 'cdrgam'))
stopifnot(identical(fit_diagnostics(block)$outcome, 'success'))

native <- structure(list(
    converged=TRUE,
    outer.info=list(conv='full convergence', hess=diag(2)),
    cdrgam=list(
        identifiability=list(global=list(condition_indicator=0.5)),
        rank=list(resolution='identified', tolerance=1e-8)
    )
), class='cdrgam')
stopifnot(identical(fit_diagnostics(native)$outcome, 'success'))
