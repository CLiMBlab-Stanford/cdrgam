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

direct_logdet <- as.numeric(determinant(as.matrix(system), logarithm=TRUE)$modulus)
schur_logdet <- cdrgam:::.cdr_factor_logdet(schur)
stopifnot(abs(direct_logdet - schur_logdet) < 1e-11)
