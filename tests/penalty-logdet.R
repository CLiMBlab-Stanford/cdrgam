library(cdrgam)

logdet <- cdrgam:::.sparse_penalty_logdet
derivatives <- cdrgam:::.sparse_penalty_logdet_derivatives
score_only <- cdrgam:::.sparse_penalty_logdet_score

block_of <- function(penalties) list(
    S=penalties, sp_index=seq_along(penalties), repetitions=2L,
    dimension=nrow(penalties[[1L]])
)
naive <- function(penalties, sp) {
    values <- eigen(Reduce(`+`, Map(`*`, sp, penalties)), symmetric=TRUE,
        only.values=TRUE)$values
    values <- values[values > max(values) * 1e-12]
    list(value=2 * sum(log(values)), rank=2L * length(values))
}

set.seed(31415)
dimension <- 6L
rotation <- qr.Q(qr(matrix(stats::rnorm(dimension^2), dimension)))
A <- matrix(stats::rnorm(dimension * 2), dimension, 2)
overlapping <- list(tcrossprod(A), diag(dimension))
mixed <- list(
    tcrossprod(A),
    crossprod(matrix(stats::rnorm(3 * dimension), 3, dimension)),
    diag(dimension)
)

# Moderate smoothing-parameter ratios agree with a plain eigen decomposition.
for (sp in list(c(0.3, 2), c(5, 0.1, 1.7))) {
    penalties <- if (length(sp) == 2L) overlapping else mixed
    stable <- logdet(list(block_of(penalties)), sp)
    reference <- naive(penalties, sp)
    stopifnot(
        abs(stable$value - reference$value) < 1e-9,
        stable$rank == reference$rank
    )
}

# Extreme ratios: lambda_1 A A' + lambda_2 I has eigenvalues
# lambda_1 mu_i + lambda_2 on range(A) and lambda_2 elsewhere.
mu <- eigen(tcrossprod(A), symmetric=TRUE, only.values=TRUE)$values[1:2]
for (log_sp in list(c(30, -8), c(-8, 30), c(40, -15))) {
    sp <- exp(log_sp)
    big <- sp[[1L]] * mu
    total <- big + sp[[2L]]
    exact_value <- 2 * (sum(log(total)) + (dimension - 2) * log(sp[[2L]]))
    exact_score <- 2 * c(
        sum(big / total),
        sum(sp[[2L]] / total) + dimension - 2
    )
    curvature <- 2 * sum(big * sp[[2L]] / total^2)
    exact_hessian <- matrix(c(curvature, -curvature, -curvature, curvature), 2)
    block <- list(block_of(overlapping))
    stable <- logdet(block, sp)
    second <- derivatives(block, sp, 2L)
    stopifnot(
        stable$rank == 2L * dimension,
        abs(stable$value - exact_value) < 1e-8 * max(1, abs(exact_value)),
        max(abs(second$score - exact_score)) < 1e-8,
        max(abs(second$hessian - exact_hessian)) < 1e-8,
        max(abs(score_only(block, sp, 2L) - exact_score)) < 1e-8
    )
}

# Penalties on orthogonal subspaces keep their exact rank and log
# determinant even when their scales differ by about 10^20.
orthogonal <- list(
    rotation %*% diag(c(1, 1, 0, 0, 0, 0)) %*% t(rotation),
    rotation %*% diag(c(0, 0, 1, 1, 0, 0)) %*% t(rotation)
)
sp <- exp(c(35, -12))
stable <- logdet(list(block_of(orthogonal)), sp)
second <- derivatives(list(block_of(orthogonal)), sp, 2L)
stopifnot(
    stable$rank == 8L,
    abs(stable$value - 2 * (2 * 35 + 2 * -12)) < 1e-8,
    max(abs(second$score - c(4, 4))) < 1e-8,
    max(abs(second$hessian)) < 1e-8,
    naive(orthogonal, sp)$rank < 8L
)

# Analytic derivatives match central differences of the stable value at a
# ratio where the plain eigen decomposition already loses the small penalty.
log_sp <- c(18, -9, 2)
block <- list(block_of(mixed))
second <- derivatives(block, exp(log_sp), 3L)
step <- 1e-4
finite_score <- vapply(seq_along(log_sp), function(j) {
    upper <- lower <- log_sp
    upper[[j]] <- upper[[j]] + step
    lower[[j]] <- lower[[j]] - step
    (logdet(block, exp(upper))$value - logdet(block, exp(lower))$value) /
        (2 * step)
}, numeric(1))
finite_hessian <- sapply(seq_along(log_sp), function(j) {
    upper <- lower <- log_sp
    upper[[j]] <- upper[[j]] + step
    lower[[j]] <- lower[[j]] - step
    (derivatives(block, exp(upper), 3L)$score -
        derivatives(block, exp(lower), 3L)$score) / (2 * step)
})
stopifnot(
    max(abs(second$score - finite_score)) < 1e-6,
    max(abs(second$hessian - finite_hessian)) < 1e-6,
    isSymmetric(second$hessian)
)

# A block without penalties contributes nothing.
empty <- list(S=list(), sp_index=integer(), repetitions=1L, dimension=3L)
stopifnot(
    identical(logdet(list(empty), numeric())$rank, 0L),
    logdet(list(empty), numeric())$value == 0,
    length(derivatives(list(empty), numeric(), 0L)$score) == 0L
)
