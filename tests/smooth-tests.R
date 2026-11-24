library(cdrgam)

# The port of Wood's (2013) test gives mgcv's p-value and reference rank for
# integer and fractional ranks, with and without an estimated scale. The
# reported statistic may be the second of the two statistics averaged for a
# fractional rank, because which one is reported depends on the sign
# convention of the square root of X'WX.
set.seed(11)
for (case in seq_len(200L)) {
    k <- sample(2:8, 1L)
    X <- matrix(stats::rnorm(40L * k), 40L, k)
    V <- crossprod(matrix(stats::rnorm(k * k), k)) / k
    p <- stats::rnorm(k) * stats::runif(1L)
    rank <- if (case %% 4L == 0L) sample(seq_len(k), 1L) else stats::runif(1L, 0.3, k)
    residual_df <- if (case %% 2L) -1 else 30
    expected <- mgcv:::testStat(p, X, V, rank, type=0, res.df=residual_df)
    observed <- cdrgam:::.cdrgam_wood_test(p, crossprod(X), V, rank, residual_df)
    stopifnot(
        abs(expected$pval - observed$p_value) < 1e-10,
        abs(expected$rank - observed$rank) < 1e-12
    )
}

simulation <- simulate_cdr(
    list(x=function(lag) exp(-lag)),
    n_impulses=100, n_responses=90, duration=20, window=2, seed=417
)
design <- prepare_cdrgam(
    response ~ irf(x, k_l=5, k_p=4),
    simulation$impulses, simulation$responses, window=c(0, 2), quiet=TRUE
)

# The block backend fits the same model as native mgcv, so its smooth table
# now reproduces summary.gam(): Wood's statistic with edf1 as reference rank.
native <- summary(cdrgam.fit(design, backend='mgcv', method='REML'))$s.table
block <- summary(cdrgam.fit(design, backend='block', method='REML'))
stopifnot(
    all(attr(block$s.table, 'test.status') == 'Wood (2013)'),
    max(abs(block$s.table[, 'Ref.df'] - native[, 'Ref.df'])) < 1e-3,
    max(abs(block$s.table[, 'F'] - native[, 'F']) / native[, 'F']) < 1e-3,
    max(abs(block$s.table[, 'p-value'] - native[, 'p-value'])) < 1e-4
)

# The sparse backend recovers X'WX = H - S and edf1 from its factor. Both
# match a dense reconstruction from the coefficient covariance.
sparse <- cdrgam.fit(design, backend='sparse', method='REML')
structure <- cdrgam:::.cdrgam_sparse_test_structure(sparse)
dimension <- length(coef(sparse))
precision <- solve(cdrgam:::.sparse_selected_vcov(sparse, seq_len(dimension)) /
    sparse$scale)
penalty <- as.matrix(Reduce(`+`, Map(`*`, sparse$sp, sparse$sparse$penalty_components)))
gram <- precision - penalty
influence <- solve(precision, gram)
edf1 <- 2 * diag(influence) - rowSums(influence * t(influence))
for (i in seq_along(sparse$smooth)) {
    indices <- sparse$smooth[[i]]$first.para:sparse$smooth[[i]]$last.para
    stopifnot(
        max(abs(structure$gram(indices) - gram[indices, indices])) <
            1e-8 * max(1, abs(gram)),
        abs(structure$edf1[[i]] - sum(edf1[indices])) < 1e-8
    )
}
sparse_table <- summary(sparse)$s.table
stopifnot(
    all(attr(sparse_table, 'test.status') == 'Wood (2013)'),
    all(sparse_table[, 'Ref.df'] >= sparse_table[, 'edf'] - 1e-8)
)

# A Poisson block fit uses the chi-squared reference and agrees with mgcv.
counts <- simulation$responses
counts$response <- stats::rpois(nrow(counts), exp(0.5 + 0.3 * counts$response))
poisson_design <- prepare_cdrgam(
    response ~ irf(x, k_l=5, k_p=4),
    simulation$impulses, counts, window=c(0, 2), quiet=TRUE
)
poisson_native <- summary(cdrgam.fit(
    poisson_design, backend='mgcv', method='REML', family=poisson()
))$s.table
poisson_block <- summary(cdrgam.fit(
    poisson_design, backend='block', method='REML', family=poisson()
))$s.table
stopifnot(
    'Chi.sq' %in% colnames(poisson_block),
    max(abs(poisson_block[, 'Ref.df'] - poisson_native[, 'Ref.df'])) < 1e-2,
    max(abs(poisson_block[, 'p-value'] - poisson_native[, 'p-value'])) < 1e-3
)
