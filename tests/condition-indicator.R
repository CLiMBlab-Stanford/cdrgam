library(cdrgam)

indicator <- cdrgam:::.cdr_factor_condition_indicator

# Dense reference: eliminate in `order` and take min L_ii^2 / A_ii.
reference <- function(system, order) {
    dense <- as.matrix(system)[order, order]
    min(diag(chol(dense))^2 / diag(dense))
}
cholmod <- function(system) {
    Matrix::Cholesky(Matrix::forceSymmetric(system, uplo='U'), LDL=FALSE,
        perm=TRUE)
}
cholmod_order <- function(factor) {
    as.integer(Matrix::expand(factor)$P %*% seq_len(factor@Dim[[1L]]))
}
layout <- list(core=1:2, blocks=list(3:4, 5:6))
schur <- function(system) {
    cdrgam:::.factor_schur_system(
        Matrix::forceSymmetric(Matrix::Matrix(system, sparse=TRUE), uplo='U'),
        layout,
        supernodal=FALSE
    )
}
schur_order <- c(3:6, 1:2)

set.seed(2718)
X <- matrix(stats::rnorm(60 * 6), 60, 6)
base <- crossprod(X) + diag(0.5, 6)
base[1:2, 5:6] <- base[5:6, 1:2] <- 0
base[3:4, 5:6] <- base[5:6, 3:4] <- 0
system <- Matrix::Matrix(base, sparse=TRUE)

# Both factor types match the dense reference in their elimination order.
factor <- cholmod(system)
stopifnot(
    abs(indicator(factor, system=system) -
        reference(system, cholmod_order(factor))) < 1e-10,
    abs(indicator(schur(system), system=system) -
        reference(system, schur_order)) < 1e-10,
    # Without the system, a CHOLMOD factor still recovers A_ii exactly.
    abs(indicator(factor) - indicator(factor, system=system)) < 1e-10
)

# A very large penalty on one coefficient leaves the indicator essentially
# unchanged; the old ratio of extreme pivots collapsed instead.
penalized <- system
penalized[6, 6] <- penalized[6, 6] + 1e12
penalized_factor <- cholmod(penalized)
pivots <- Matrix::diag(Matrix::expand(penalized_factor)$L)
stopifnot(
    (min(pivots) / max(pivots))^2 < 1e-10,
    indicator(penalized_factor, system=penalized) >
        0.5 * indicator(factor, system=system),
    indicator(schur(penalized), system=penalized) >
        0.5 * indicator(schur(system), system=system)
)

# A coefficient that is nearly a combination of others still drives the
# indicator toward zero.
collinear_X <- X
collinear_X[, 4] <- collinear_X[, 3] + 1e-7 * stats::rnorm(60)
collinear <- crossprod(collinear_X)
collinear[1:2, 5:6] <- collinear[5:6, 1:2] <- 0
collinear[3:4, 5:6] <- collinear[5:6, 3:4] <- 0
collinear <- Matrix::Matrix(collinear, sparse=TRUE)
stopifnot(
    indicator(cholmod(collinear), system=collinear) < 1e-10,
    indicator(schur(collinear), system=collinear) < 1e-10
)

# A very large rank-deficient penalty makes its own coefficients nearly
# dependent in the metric of A. Excluding them restores the indicator for the
# rest of the system.
lag_penalty <- matrix(c(1, -1, -1, 1), 2, 2)
saturated <- system
saturated[5:6, 5:6] <- saturated[5:6, 5:6] + 1e12 * lag_penalty
saturated_factor <- cholmod(saturated)
stopifnot(
    indicator(saturated_factor, system=saturated) < 1e-8,
    indicator(saturated_factor, system=saturated, exclude=5:6) >
        0.5 * indicator(cholmod(system[1:4, 1:4]), system=system[1:4, 1:4])
)
