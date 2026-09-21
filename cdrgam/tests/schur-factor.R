library(cdrgam)

# The experimental Schur factor must be algebraically identical to a direct
# factorization for an arrowhead system with two independent grouped blocks.
system <- Matrix::Matrix(matrix(c(
    7.0, 0.5, 0.8, 0.1, 0.3, 0.0,
    0.5, 6.0, 0.2, 0.6, 0.0, 0.4,
    0.8, 0.2, 4.0, 0.3, 0.0, 0.0,
    0.1, 0.6, 0.3, 3.5, 0.0, 0.0,
    0.3, 0.0, 0.0, 0.0, 3.0, 0.2,
    0.0, 0.4, 0.0, 0.0, 0.2, 2.5
), nrow=6, byrow=TRUE), sparse=TRUE)
system <- Matrix::forceSymmetric(system, uplo='U')
layout <- list(core=1:2, blocks=list(3:4, 5:6))
schur <- cdrgam:::.factor_schur_system(system, layout, supernodal=FALSE)
stopifnot(inherits(schur, 'cdrgam_schur_factor'))

rhs <- matrix(seq_len(18) / 7, nrow=6)
direct_solution <- as.matrix(Matrix::solve(system, rhs))
schur_solution <- cdrgam:::.cdr_factor_solve(schur, rhs)
stopifnot(max(abs(direct_solution - schur_solution)) < 1e-11)

legacy_schur <- schur
legacy_schur$cross <- lapply(layout$blocks, function(indices) {
    as.matrix(system[layout$core, indices, drop=FALSE])
})
legacy_schur$transfer <- NULL
legacy_solution <- cdrgam:::.cdr_factor_solve(legacy_schur, rhs)
stopifnot(max(abs(direct_solution - legacy_solution)) < 1e-11)

direct_logdet <- as.numeric(determinant(as.matrix(system), logarithm=TRUE)$modulus)
schur_logdet <- cdrgam:::.cdr_factor_logdet(schur)
stopifnot(abs(direct_logdet - schur_logdet) < 1e-11)
direct_diagonal <- diag(solve(as.matrix(system)))
schur_diagonal <- cdrgam:::.sparse_schur_inverse_diagonal(
    schur,
    chunk_size=2L
)
legacy_diagonal <- cdrgam:::.sparse_schur_inverse_diagonal(
    legacy_schur,
    chunk_size=2L
)
stopifnot(
    max(abs(direct_diagonal - schur_diagonal)) < 1e-11,
    max(abs(direct_diagonal - legacy_diagonal)) < 1e-11
)

# Exact score traces need only inverse entries within the core and independent
# blocks when the penalty respects the same partition.
components <- list(
    Matrix::forceSymmetric(Matrix::sparseMatrix(
        i=c(1, 1, 2), j=c(1, 2, 2), x=c(1.5, -0.2, 0.8), dims=c(6, 6)
    ), uplo='U'),
    Matrix::Diagonal(6, x=c(0, 0, 1, 2, 3, 4))
)
supports <- cdrgam:::.sparse_penalty_supports(components)
plan <- cdrgam:::.sparse_schur_trace_plan(layout, components)
stopifnot(!is.null(plan))
sp <- c(0.7, 1.3)
solve_scores <- cdrgam:::.sparse_logdet_scores(
    schur, components, sp, supports=supports, chunk_size=2L
)
schur_scores <- cdrgam:::.sparse_logdet_scores(
    schur, components, sp, supports=supports, chunk_size=2L,
    schur_plan=plan
)
stopifnot(max(abs(solve_scores - schur_scores)) < 1e-11)
legacy_scores <- cdrgam:::.sparse_logdet_scores(
    legacy_schur, components, sp, supports=supports, chunk_size=2L,
    schur_plan=plan
)
stopifnot(max(abs(solve_scores - legacy_scores)) < 1e-11)

# A failed Schur update may fall back to CHOLMOD for one objective evaluation.
# The score calculation must then fall back to inverse-column solves as well.
direct_factor <- Matrix::Cholesky(system, LDL=FALSE, perm=TRUE)
fallback_scores <- cdrgam:::.sparse_logdet_scores(
    direct_factor, components, sp, supports=supports, chunk_size=2L,
    schur_plan=plan
)
stopifnot(max(abs(solve_scores - fallback_scores)) < 1e-11)

cross_block <- Matrix::sparseMatrix(
    i=c(3, 5), j=c(5, 3), x=1, dims=c(6, 6)
)
stopifnot(is.null(cdrgam:::.sparse_schur_trace_plan(
    layout, list(cross_block)
)))

automatic_schur <- cdrgam:::.sparse_optimizer_selection(
    gradient_requested='auto',
    outer_optimizer_requested='auto',
    trace_method='schur_inverse',
    factor=schur,
    objective_seconds=1,
    exact_trace_rhs=6L,
    penalty_count=2L,
    gradient_cores=1L,
    dimension=6L
)
stopifnot(
    identical(automatic_schur$gradient, 'exact'),
    identical(automatic_schur$outer_optimizer, 'bfgs_trust'),
    automatic_schur$schur_trace
)

automatic_solve <- cdrgam:::.sparse_optimizer_selection(
    gradient_requested='auto',
    outer_optimizer_requested='auto',
    trace_method='solve',
    factor=direct_factor,
    objective_seconds=1,
    exact_trace_rhs=10000L,
    penalty_count=2L,
    gradient_cores=1L,
    dimension=6L
)
stopifnot(
    identical(automatic_solve$gradient, 'finite'),
    identical(automatic_solve$outer_optimizer, 'lbfgsb')
)

explicit_trust <- cdrgam:::.sparse_optimizer_selection(
    gradient_requested='auto',
    outer_optimizer_requested='bfgs_trust',
    trace_method='solve',
    factor=direct_factor,
    objective_seconds=1,
    exact_trace_rhs=10000L,
    penalty_count=2L,
    gradient_cores=1L,
    dimension=6L
)
stopifnot(
    identical(explicit_trust$gradient, 'exact'),
    identical(explicit_trust$outer_optimizer, 'bfgs_trust')
)

saved_resolution <- cdrgam:::.sparse_optimizer_selection(
    gradient_requested='auto',
    outer_optimizer_requested='auto',
    trace_method='schur_inverse',
    factor=schur,
    objective_seconds=1,
    exact_trace_rhs=6L,
    penalty_count=2L,
    gradient_cores=1L,
    dimension=6L,
    saved=list(policy_version=1L, gradient='finite')
)
stopifnot(
    identical(saved_resolution$gradient, 'finite'),
    identical(saved_resolution$outer_optimizer, 'lbfgsb'),
    identical(saved_resolution$reason, 'checkpoint resolution reused')
)
