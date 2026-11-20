library(cdrgam)

# The exact sparse gradient agrees with central differences of the defining
# GCV score for overlapping penalties.
set.seed(9820)
X <- matrix(stats::rnorm(30 * 8), 30, 8)
y <- stats::rnorm(30)
XtX <- Matrix::Matrix(crossprod(X), sparse=TRUE)
Xty <- as.numeric(crossprod(X, y))
components <- list(
    Matrix::sparseMatrix(
        i=1:5,
        j=1:5,
        x=c(1, 2, 1.5, 0.7, 1.2),
        dims=c(8, 8)
    ),
    {
        local <- tcrossprod(matrix(stats::rnorm(5 * 3), 5, 3))
        output <- Matrix::Matrix(0, 8, 8, sparse=TRUE)
        output[4:8, 4:8] <- local
        output
    }
)
supports <- cdrgam:::.sparse_penalty_supports(components)
gcv_evaluate <- function(log_sp, retain=FALSE) {
    sp <- exp(log_sp)
    penalty <- Reduce(`+`, Map(`*`, sp, components))
    system <- Matrix::forceSymmetric(XtX + penalty, uplo='U')
    factor <- Matrix::Cholesky(system, LDL=FALSE, perm=TRUE)
    coefficients <- as.numeric(cdrgam:::.cdr_factor_solve(factor, Xty))
    residual_rss <- sum((y - as.numeric(X %*% coefficients))^2)
    effective_df <- sum(diag(as.matrix(
        cdrgam:::.cdr_factor_solve(factor, XtX)
    )))
    score <- length(y) * residual_rss /
        (length(y) - effective_df)^2
    if (!retain) return(score)
    list(
        score=score,
        factor=factor,
        sp=sp,
        coefficients=coefficients,
        residual_rss=residual_rss,
        effective_df=effective_df
    )
}
gradient_point <- c(0.3, -0.5)
gradient_solution <- gcv_evaluate(gradient_point, retain=TRUE)
exact_gradient <- cdrgam:::.sparse_exact_gcv_gradient(
    gradient_solution$factor,
    gradient_solution$sp,
    components,
    supports,
    gradient_solution$coefficients,
    gradient_solution$residual_rss,
    gradient_solution$effective_df,
    length(y)
)
step <- 1e-5
finite_gradient <- vapply(seq_along(gradient_point), function(index) {
    lower <- upper <- gradient_point
    lower[[index]] <- lower[[index]] - step
    upper[[index]] <- upper[[index]] + step
    (gcv_evaluate(upper) - gcv_evaluate(lower)) / (2 * step)
}, numeric(1))
stopifnot(max(abs(exact_gradient - finite_gradient)) < 1e-7)

gradient_function <- function(log_sp) {
    retained <- gcv_evaluate(log_sp, retain=TRUE)
    cdrgam:::.sparse_exact_gcv_gradient(
        retained$factor,
        retained$sp,
        components,
        supports,
        retained$coefficients,
        retained$residual_rss,
        retained$effective_df,
        length(y)
    )
}
gradient_hessian <- cdrgam:::.central_difference_jacobian(
    gradient_function,
    gradient_point,
    step=1e-4
)$hessian
objective_hessian <- cdrgam:::.central_difference_hessian(
    gcv_evaluate,
    gradient_point,
    step=1e-4
)$hessian
stopifnot(
    max(abs(gradient_hessian - objective_hessian)) /
        max(1, max(abs(objective_hessian))) < 1e-5
)

# The sparse Gaussian GCV objective must reproduce native mgcv on the same
# prepared design, including its exact effective degrees of freedom and scale.
simulation <- simulate_cdr(
    list(
        x1=function(lag) 1.1 * exp(-1.8 * lag),
        x2=function(lag) -0.7 * exp(-((lag - 0.55) / 0.28)^2)
    ),
    n_impulses=180,
    n_responses=240,
    duration=32,
    window=2,
    noise_sd=0.7,
    seed=9821
)
design <- prepare_cdrgam(
    response ~ irf(x1, window=c(0, 2), k=7) +
        irf(x2, window=c(0, 2), k=7) - irf(1),
    simulation$impulses,
    simulation$responses,
    history='ragged',
    quiet=TRUE
)

native <- cdrgam.fit(
    design,
    backend='mgcv',
    engine='gam',
    method='GCV.Cp'
)
sparse <- cdrgam.fit(
    design,
    backend='sparse',
    method='GCV.Cp',
    sparse_control=list(
        gradient='finite',
        outer_optimizer='lbfgsb',
        hessian='none',
        optimizer_maxit=200,
        restarts=0
    )
)

native_edf <- sum(native$edf)
sparse_edf <- cdrgam:::.sparse_effective_df(sparse)
stopifnot(
    isTRUE(sparse$converged),
    identical(sparse$method, 'GCV.Cp'),
    identical(sparse$sparse$control$gradient, 'finite'),
    sparse$sparse$control$gradient_workers >= 1L,
    sparse$sparse$control$gradient_blas_threads >= 1L,
    is.logical(sparse$sparse$control$gradient_memory_limited),
    is.finite(sparse$sparse$control$gradient_worker_bytes),
    identical(names(sparse$gcv.ubre), 'GCV.Cp'),
    abs(unname(sparse$gcv.ubre - native$gcv.ubre)) /
        unname(native$gcv.ubre) < 1e-5,
    sqrt(mean((fitted(sparse) - fitted(native))^2)) < 1e-3,
    abs(native_edf - sparse_edf) < 1e-3,
    abs(native$sig2 - sparse$sig2) < 1e-3,
    identical(summary(sparse)$method, 'GCV'),
    isTRUE(sparse$scale.estimated)
)

sparse_newton <- cdrgam.fit(
    design,
    backend='sparse',
    method='GCV.Cp',
    sparse_control=list(
        gradient='exact',
        outer_optimizer='newton',
        hessian='gradient',
        optimizer_maxit=100,
        restarts=0
    )
)
gcv_unconditional_error <- tryCatch(
    vcov(sparse_newton, unconditional=TRUE),
    error=conditionMessage
)
stopifnot(
    isTRUE(sparse_newton$converged),
    identical(sparse_newton$sparse$outer_optimizer, 'newton'),
    identical(sparse_newton$sparse$hessian, 'gradient'),
    identical(dim(sparse_newton$outer.info$hess), c(2L, 2L)),
    isTRUE(sparse_newton$sparse$convergence$hessian_positive_definite),
    isTRUE(all.equal(
        sparse_newton$outer.info$hess,
        t(sparse_newton$outer.info$hess),
        tolerance=1e-10
    )),
    sparse_newton$sparse$newton_hessian_gradient_evaluations > 0L,
    grepl('does not define smoothing-parameter covariance',
        gcv_unconditional_error, fixed=TRUE),
    abs(unname(sparse_newton$gcv.ubre - native$gcv.ubre)) /
        unname(native$gcv.ubre) < 1e-5,
    sqrt(mean((fitted(sparse_newton) - fitted(native))^2)) < 1e-3
)

# Gamma must enter the GCV denominator exactly as it does in mgcv.
native_gamma <- cdrgam.fit(
    design,
    backend='mgcv',
    engine='gam',
    method='GCV.Cp',
    gamma=1.4
)
previous_memory_limit <- getOption('cdrgam.memory_limit_bytes')
process_memory <- max(c(
    cdrgam:::.cdrgam_proc_memory('VmRSS'),
    cdrgam:::.cdrgam_r_memory()
), na.rm=TRUE)
options(cdrgam.memory_limit_bytes=process_memory + 8 * 1024^3)
sparse_gamma <- cdrgam.fit(
    design,
    backend='sparse',
    method='GCV.Cp',
    gamma=1.4,
    # Automatic selection times both gradients, so request the exact one.
    sparse_control=list(
        gradient='exact',
        hessian='none',
        optimizer_maxit=200,
        restarts=0
    )
)
options(cdrgam.memory_limit_bytes=previous_memory_limit)
stopifnot(
    isTRUE(sparse_gamma$converged),
    identical(sparse_gamma$sparse$control$gradient, 'exact'),
    identical(sparse_gamma$sparse$outer_optimizer, 'newton'),
    abs(unname(sparse_gamma$gcv.ubre - native_gamma$gcv.ubre)) /
        unname(native_gamma$gcv.ubre) < 1e-5,
    abs(sum(native_gamma$edf) -
        cdrgam:::.sparse_effective_df(sparse_gamma)) < 1e-3
)

stopifnot(
    identical(sparse_gamma$sparse$control$boundary_action, 'reduce'),
    identical(sparse_gamma$sparse$boundary_action, 'reduce'),
    is.numeric(sparse_gamma$optimizer$gradient),
    length(sparse_gamma$optimizer$gradient) == length(sparse_gamma$sp)
)
