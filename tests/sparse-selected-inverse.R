library(cdrgam)

# The streamed selected-inverse traces used by the analytic outer Hessian and
# the exact GCV gradient agree with dense formulas. The penalties include a
# shared grouped support, an overlapping support, and an empty support.
set.seed(4417)
dimension <- 14L
X <- matrix(stats::rnorm(60 * dimension), 60, dimension)
XtX <- Matrix::Matrix(crossprod(X), sparse=TRUE)
coefficients <- stats::rnorm(dimension)

grouped <- 1:9
group_block <- function(local) {
    output <- Matrix::Matrix(0, dimension, dimension, sparse=TRUE)
    for (start in c(1L, 4L, 7L)) {
        rows <- start:(start + 2L)
        output[rows, rows] <- local
    }
    output
}
components <- list(
    group_block(tcrossprod(matrix(stats::rnorm(6), 3, 2))),
    group_block(diag(3)),
    {
        local <- tcrossprod(matrix(stats::rnorm(7 * 4), 7, 4))
        output <- Matrix::Matrix(0, dimension, dimension, sparse=TRUE)
        output[8:14, 8:14] <- local
        output
    },
    Matrix::Matrix(0, dimension, dimension, sparse=TRUE)
)
supports <- cdrgam:::.sparse_penalty_supports(components)
sp <- c(0.7, 1.9, 0.4, 2.0)
penalty <- Reduce(`+`, Map(`*`, sp, components))
system <- Matrix::forceSymmetric(XtX + penalty, uplo='U')
factor <- Matrix::Cholesky(system, LDL=FALSE, perm=TRUE)

inverse <- solve(as.matrix(system))
derivatives <- lapply(seq_along(sp), function(i) {
    sp[[i]] * as.matrix(components[[i]])
})
reference_system <- vapply(derivatives, function(d) {
    sum(diag(inverse %*% d))
}, numeric(1))
reference_cross <- outer(seq_along(sp), seq_along(sp), Vectorize(
    function(i, j) sum(diag(
        inverse %*% derivatives[[i]] %*% inverse %*% derivatives[[j]]
    ))
))
limit_message <- list(invalid='invalid', prefix='limit ', fallback='fallback')

groups <- cdrgam:::.sparse_trace_batches(components, supports, chunk_size=3L)
grouped_batches <- groups[[1L]]$batches
stopifnot(
    length(groups) == 2L,
    identical(groups[[1L]]$penalties, 1:2),
    identical(lengths(grouped_batches), rep(3L, 3L)),
    identical(sort(unlist(grouped_batches)), grouped),
    identical(lengths(groups[[2L]]$batches), 7L)
)

for (chunk_size in c(1L, 3L, 64L)) {
    traces <- cdrgam:::.sparse_selected_inverse_traces(
        factor, sp, components, supports, coefficients,
        chunk_size=chunk_size,
        element_limit=1e6,
        limit_message=limit_message
    )
    union <- traces$union_support
    reference_directions <- lapply(derivatives, function(d) {
        as.numeric(inverse %*% d %*% coefficients)[union]
    })
    stopifnot(
        max(abs(traces$system_scores - reference_system)) < 1e-10,
        max(abs(traces$trace_cross - reference_cross)) < 1e-10,
        max(abs(unlist(traces$coefficient_derivatives) -
            unlist(reference_directions))) < 1e-10,
        isSymmetric(traces$trace_cross)
    )
}

wide <- cdrgam:::.sparse_selected_inverse_traces(
    factor, sp, components, supports, coefficients,
    chunk_size=64L, element_limit=1e6, limit_message=limit_message
)
stopifnot(wide$working_elements == length(wide$union_support) * 9)
limited <- tryCatch(
    cdrgam:::.sparse_selected_inverse_traces(
        factor, sp, components, supports, coefficients,
        chunk_size=3L, element_limit=1, limit_message=limit_message
    ),
    error=function(error) conditionMessage(error)
)
stopifnot(is.character(limited), startsWith(limited, 'limit '))

# The memory estimate uses the batch width when penalty components are known
# and falls back to the widest support otherwise.
estimate <- cdrgam:::.sparse_analytic_hessian_memory(
    supports, dimension, chunk_size=3L, components=components
)
fallback <- cdrgam:::.sparse_analytic_hessian_memory(
    supports, dimension, chunk_size=3L
)
stopifnot(
    estimate$retained_elements == 14 * 7,
    fallback$retained_elements == 14 * 9
)
