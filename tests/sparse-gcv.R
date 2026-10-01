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
    sparse_control=list(
        hessian='none',
        optimizer_maxit=200,
        restarts=0
    )
)
options(cdrgam.memory_limit_bytes=previous_memory_limit)
stopifnot(
    isTRUE(sparse_gamma$converged),
    identical(sparse_gamma$sparse$control$gradient, 'exact'),
    abs(unname(sparse_gamma$gcv.ubre - native_gamma$gcv.ubre)) /
        unname(native_gamma$gcv.ubre) < 1e-5,
    abs(sum(native_gamma$edf) -
        cdrgam:::.sparse_effective_df(sparse_gamma)) < 1e-3
)

gcv_reduction_error <- tryCatch(
    cdrgam.fit(
        design,
        backend='sparse',
        method='GCV.Cp',
        sparse_control=list(boundary_action='reduce')
    ),
    error=conditionMessage
)
stopifnot(
    is.character(gcv_reduction_error),
    grepl('does not yet support boundary_action', gcv_reduction_error)
)
