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
