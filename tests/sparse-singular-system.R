library(cdrgam)

# A penalized system that is singular to rounding is factored after a small
# relative ridge on its diagonal, and the ridge used is reported.
singular <- Matrix::forceSymmetric(Matrix::Matrix(
    c(1, 1, 0,
      1, 1, 0,
      0, 0, 2),
    3, 3, sparse=TRUE
), uplo='U')
singular <- methods::as(singular, 'CsparseMatrix')
diagonal_positions <- which(
    methods::slot(singular, 'i') ==
        rep.int(0:2, diff(methods::slot(singular, 'p')))
)
factor_values <- function(values) {
    system <- singular
    methods::slot(system, 'x') <- values
    factor <- tryCatch(
        Matrix::Cholesky(system, LDL=FALSE, perm=TRUE, super=FALSE),
        warning=function(w) NULL,
        error=function(e) NULL
    )
    list(system=system, factor=factor)
}
values <- methods::slot(singular, 'x')
values[diagonal_positions[[2L]]] <- 1 - 1e-12
stopifnot(is.null(factor_values(values)$factor))
ridged <- cdrgam:::.sparse_ridged_factor(
    values, diagonal_positions, factor_values, c(1e-14, 1e-10, 1e-8)
)
stopifnot(
    !is.null(ridged$factor),
    ridged$stabilization == 1e-10,
    is.null(cdrgam:::.sparse_ridged_factor(
        values, diagonal_positions, factor_values, 1e-14
    ))
)

# The gradient guard computes each point once and skips points whose
# objective failed, returning the last valid gradient instead.
calls <- 0L
skips <- 0L
failed_points <- character()
guarded <- cdrgam:::.sparse_guarded_gradient(
    function(x) {
        calls <<- calls + 1L
        2 * x
    },
    key=function(x) paste(x, collapse=','),
    failed=function(key) key %in% failed_points,
    on_skip=function() skips <<- skips + 1L
)
stopifnot(identical(guarded(c(0, 0)), c(0, 0)))
failed_points <- '25,-25'
stopifnot(identical(guarded(c(25, -25)), c(0, 0)), skips == 1L)
stopifnot(identical(guarded(c(1, 2)), c(2, 4)))
stopifnot(identical(guarded(c(1, 2)), c(2, 4)), calls == 2L)
stopifnot(identical(guarded(c(25, -25)), c(2, 4)), skips == 2L)

cat('Singular-system checks passed.\n')

# A fit logs its full named smoothing-parameter vector at the end, as one
# field that the reporter does not shorten.
simulation <- simulate_cdr(
    list(x=function(lag) exp(-lag)),
    n_impulses=60, n_responses=80, duration=10, window=2, noise_sd=1, seed=9
)
design <- prepare_cdrgam(
    response ~ irf(x, window=c(0, 2), k=5) + s(time, k=8),
    simulation$impulses, simulation$responses, quiet=TRUE
)
log_lines <- c(
    utils::capture.output(
        messages <- utils::capture.output(
            fit <- cdrgam:::.fit_sparse_gaussian(
                design, family=gaussian(), method='REML', trace=1,
                sparse_control=list(final_newton=FALSE)
            ),
            type='message'
        )
    ),
    messages
)
final_line <- grep('final smoothing parameters', log_lines, value=TRUE)
stopifnot(length(final_line) == 1L)
logged <- strsplit(sub('.*log_sp=', '', final_line), ',', fixed=TRUE)[[1L]]
stopifnot(
    length(logged) == length(fit$sp),
    !any(logged == '...'),
    all(grepl('^s\\(', logged))
)

cat('Smoothing-parameter logging checks passed.\n')
