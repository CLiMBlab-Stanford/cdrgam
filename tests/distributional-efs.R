library(cdrgam)

efs <- cdrgam:::.cdrgam_efs_optimize

# The ratio a_j / lambda_j has its fixed point at log(a_j); the criterion is
# minimized there, so the updates converge to it.
target <- c(1.5, -2)
quadratic <- function(log_sp) list(
    criterion=sum((log_sp - target)^2), solution=list()
)
converged <- efs(
    quadratic,
    function(evaluation, log_sp) exp(target) / exp(log_sp),
    log_sp=c(0, 0), lower=c(-25, -25), upper=c(25, 25)
)
stopifnot(
    identical(converged$status, 'converged'),
    max(abs(converged$log_sp - target)) < 0.01,
    length(converged$drifting) == 0L
)

# A parameter whose ratio stays above one keeps increasing; once beyond the
# boundary threshold it is reported as drifting and stops blocking
# convergence of the others.
drifting <- efs(
    function(log_sp) list(
        criterion=(log_sp[[1L]] - 1)^2 + exp(-log_sp[[2L]]), solution=list()
    ),
    function(evaluation, log_sp) c(exp(1) / exp(log_sp[[1L]]), exp(3)),
    log_sp=c(0, 0), lower=c(-25, -25), upper=c(25, 25), maxit=20L
)
stopifnot(
    identical(drifting$status, 'converged'),
    identical(drifting$drifting, 2L),
    drifting$log_sp[[2L]] >= 8,
    abs(drifting$log_sp[[1L]] - 1) < 0.01
)

# End to end: an EFS warm-up on a small location-scale model reaches a REML
# criterion at least as low as the exact optimizer alone.
set.seed(9341)
impulses <- data.frame(
    time=sort(stats::runif(180L, 0, 45)),
    location_signal=stats::rnorm(180L),
    scale_signal=stats::rnorm(180L)
)
responses <- data.frame(time=sort(stats::runif(260L, 0.5, 45)))
responses$response <- stats::rnorm(260L)
design <- prepare_cdrgam(
    list(
        location=response ~ irf(location_signal, window=c(0, 1.5), k_l=6) -
            irf(1),
        scale=~ irf(scale_signal, window=c(0, 1.5), k_l=6) - irf(1)
    ),
    impulses, responses, history='ragged', quiet=TRUE
)
exact <- cdrgam.fit(design, family='gaulss', backend='sparse', method='REML')
warmed <- cdrgam.fit(design, family='gaulss', backend='sparse', method='REML',
    sparse_control=list(warmup='efs'))
summary <- warmed$distributional$efs_warmup
stopifnot(
    is.null(exact$distributional$efs_warmup),
    is.list(summary),
    summary$iterations >= 1L,
    isTRUE(warmed$converged),
    is.finite(warmed$reml),
    is.finite(exact$reml),
    warmed$reml <= exact$reml + 1e-3 * max(1, abs(exact$reml))
)
