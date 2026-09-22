# Compare the block-at-a-time and batched Schur kernels on the same system.

library(cdrgam)

core_dimension <- as.integer(Sys.getenv('CDRGAM_SCHUR_CORE', '800'))
block_count <- as.integer(Sys.getenv('CDRGAM_SCHUR_BLOCKS', '2000'))
cross_density <- as.numeric(Sys.getenv('CDRGAM_SCHUR_DENSITY', '0.25'))
right_hand_sides <- as.integer(Sys.getenv('CDRGAM_SCHUR_RHS', '16'))
seed <- as.integer(Sys.getenv('CDRGAM_SCHUR_SEED', '20260918'))
stopifnot(
    core_dimension > 0L,
    block_count > 0L,
    cross_density > 0,
    cross_density <= 1,
    right_hand_sides > 0L
)

set.seed(seed)
cross <- matrix(
    stats::rnorm(core_dimension * block_count, sd=0.002),
    core_dimension
)
cross[matrix(stats::runif(length(cross)), core_dimension) > cross_density] <- 0
core <- tcrossprod(cross) + diag(core_dimension) * 2
blocks <- rep(list(matrix(1, 1L, 1L)), block_count)
crosses <- lapply(seq_len(block_count), function(index) {
    cross[, index, drop=FALSE]
})
dimension <- core_dimension + block_count
system <- Matrix::bdiag(
    Matrix::Matrix(core, sparse=TRUE),
    Matrix::Diagonal(block_count, 1)
)
system[seq_len(core_dimension), core_dimension + seq_len(block_count)] <- cross
system <- Matrix::forceSymmetric(system, uplo='U')
layout <- list(
    core=seq_len(core_dimension),
    blocks=lapply(core_dimension + seq_len(block_count), as.integer)
)

legacy_seconds <- system.time({
    raw_legacy <- .Call(
        'cdrgam_schur_factor', core, blocks, crosses, PACKAGE='cdrgam'
    )
})[['elapsed']]
legacy <- list(
    core=layout$core,
    blocks=layout$blocks,
    cross=crosses,
    block_cholesky=raw_legacy[[2L]],
    core_cholesky=raw_legacy[[1L]],
    core_logdet=2 * sum(log(diag(raw_legacy[[1L]]))),
    dimension=dimension
)
class(legacy) <- 'cdrgam_schur_factor'

batched_seconds <- system.time({
    batched <- cdrgam:::.factor_schur_system(system, layout, FALSE)
})[['elapsed']]
full_sparse_seconds <- system.time({
    full_sparse <- Matrix::Cholesky(system, LDL=FALSE, perm=TRUE)
})[['elapsed']]

rhs <- matrix(stats::rnorm(dimension * right_hand_sides), dimension)
legacy_solve_seconds <- system.time({
    legacy_solution <- cdrgam:::.cdr_factor_solve(legacy, rhs)
})[['elapsed']]
batched_solve_seconds <- system.time({
    batched_solution <- cdrgam:::.cdr_factor_solve(batched, rhs)
})[['elapsed']]
full_sparse_solve_seconds <- system.time({
    full_sparse_solution <- Matrix::solve(full_sparse, rhs)
})[['elapsed']]

component <- Matrix::Diagonal(
    dimension,
    x=c(rep.int(0, core_dimension), rep.int(1, block_count))
)
plan <- cdrgam:::.sparse_schur_trace_plan(layout, list(component))
legacy_score_seconds <- system.time({
    legacy_score <- cdrgam:::.sparse_logdet_scores(
        legacy, list(component), 1, schur_plan=plan
    )
})[['elapsed']]
batched_score_seconds <- system.time({
    batched_score <- cdrgam:::.sparse_logdet_scores(
        batched, list(component), 1, schur_plan=plan
    )
})[['elapsed']]

result <- data.frame(
    core_dimension=core_dimension,
    block_count=block_count,
    cross_density=cross_density,
    right_hand_sides=right_hand_sides,
    legacy_factor_seconds=legacy_seconds,
    batched_factor_seconds=batched_seconds,
    factor_speedup=legacy_seconds / batched_seconds,
    full_sparse_factor_seconds=full_sparse_seconds,
    legacy_solve_seconds=legacy_solve_seconds,
    batched_solve_seconds=batched_solve_seconds,
    solve_speedup=legacy_solve_seconds / batched_solve_seconds,
    full_sparse_solve_seconds=full_sparse_solve_seconds,
    legacy_score_seconds=legacy_score_seconds,
    batched_score_seconds=batched_score_seconds,
    score_speedup=legacy_score_seconds / batched_score_seconds,
    maximum_solve_difference=max(abs(legacy_solution - batched_solution)),
    full_sparse_solve_difference=max(abs(
        as.matrix(full_sparse_solution) - batched_solution
    )),
    batched_factor_mib=as.numeric(object.size(batched)) / 1024^2,
    full_sparse_factor_mib=as.numeric(object.size(full_sparse)) / 1024^2,
    score_difference=abs(legacy_score - batched_score)
)
print(result, row.names=FALSE)
