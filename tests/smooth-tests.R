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

# Columns of the factored system H are recovered from both factor types.
system <- Matrix::forceSymmetric(Matrix::Matrix(matrix(c(
    7.0, 0.5, 0.8, 0.1, 0.3, 0.0,
    0.5, 6.0, 0.2, 0.6, 0.0, 0.4,
    0.8, 0.2, 4.0, 0.3, 0.0, 0.0,
    0.1, 0.6, 0.3, 3.5, 0.0, 0.0,
    0.3, 0.0, 0.0, 0.0, 3.0, 0.2,
    0.0, 0.4, 0.0, 0.0, 0.2, 2.5
), nrow=6, byrow=TRUE), sparse=TRUE), uplo='U')
schur <- cdrgam:::.factor_schur_system(
    system, list(core=1:2, blocks=list(3:4, 5:6)), supernodal=FALSE
)
cholmod <- Matrix::Cholesky(system, LDL=FALSE)
for (columns in list(1:6, c(5L, 1L, 4L))) {
    stopifnot(
        max(abs(cdrgam:::.cdr_factor_system_columns(schur, columns) -
            as.matrix(system)[, columns])) < 1e-12,
        max(abs(cdrgam:::.cdr_factor_system_columns(cholmod, columns) -
            as.matrix(system)[, columns])) < 1e-12
    )
}

# Fully penalized terms get Wood's (2013b) random-effect test, with other
# random effects integrated out, as in mgcv's reTest(). Block fits and
# sparse fits with either factorization reproduce the native table.
grouped_impulses <- list()
grouped_responses <- list()
for (i in 1:4) {
    amplitude <- c(-0.5, -0.2, 0.2, 0.5)[[i]]
    grouped <- simulate_cdr(
        list(x=function(lag) {
            0.9 * exp(-1.7 * lag) + amplitude * exp(-((lag - 0.65) / 0.3)^2)
        }),
        n_impulses=120, n_responses=150, duration=30, window=2,
        intercept=75 + i / 2, noise_sd=0.45, seed=6600 + i
    )
    grouped$impulses$subject <- grouped$responses$subject <- paste0('s', i)
    grouped_impulses[[i]] <- grouped$impulses
    grouped_responses[[i]] <- grouped$responses
}
grouped_responses <- do.call(rbind, grouped_responses)
grouped_responses$subject <- factor(grouped_responses$subject)
grouped_design <- prepare_cdrgam(
    response ~ s(subject, bs='re') + irf(x, window=c(0, 2), k=6) +
        irf(x, window=c(0, 2), k=6, group=subject) - irf(1),
    do.call(rbind, grouped_impulses), grouped_responses,
    series='subject', history='auto', quiet=TRUE
)
grouped_native <- summary(cdrgam.fit(
    grouped_design, backend='mgcv', method='REML'
))$s.table
grouped_tables <- list(
    block=summary(cdrgam.fit(grouped_design, backend='block', method='REML')),
    sparse=summary(cdrgam.fit(grouped_design, backend='sparse', method='REML')),
    schur=summary(cdrgam.fit(
        grouped_design, backend='sparse', method='REML',
        sparse_control=list(schur='always')
    ))
)
for (grouped_summary in grouped_tables) {
    table <- grouped_summary$s.table
    stopifnot(
        identical(unname(attr(table, 'test.status')),
            c('Wood (2013b) random effect', 'Wood (2013)',
                'Wood (2013b) random effect')),
        max(abs(table[, 'Ref.df'] - grouped_native[, 'Ref.df'])) < 1e-3,
        max(abs(table[, 'F'] - grouped_native[, 'F']) / grouped_native[, 'F']) < 1e-3
    )
}

# The random-effect port gives mgcv's p-value when it is not negligible:
# refitting with a weak subject effect keeps the statistic moderate.
weak <- grouped_responses
weak$response <- weak$response - (as.integer(weak$subject) / 2) +
    stats::rnorm(4L, sd=0.02)[as.integer(weak$subject)]
weak_design <- prepare_cdrgam(
    response ~ s(subject, bs='re') + irf(x, window=c(0, 2), k=6),
    do.call(rbind, grouped_impulses), weak,
    series='subject', history='auto', quiet=TRUE
)
weak_native <- summary(cdrgam.fit(weak_design, backend='mgcv', method='REML'))$s.table
weak_block <- summary(cdrgam.fit(weak_design, backend='block', method='REML'))$s.table
stopifnot(abs(weak_block[1L, 'p-value'] - weak_native[1L, 'p-value']) < 1e-3)

# The whole factored system is also recovered as a sparse matrix.
stopifnot(
    max(abs(as.matrix(cdrgam:::.cdr_factor_system(schur)) - as.matrix(system))) < 1e-12,
    max(abs(as.matrix(cdrgam:::.cdr_factor_system(cholmod)) - as.matrix(system))) < 1e-12
)

# Schur factors saved before the transfer was stored hold H_cb as `cross`,
# from which the same tests are computed.
legacy <- cdrgam.fit(
    grouped_design, backend='sparse', method='REML',
    sparse_control=list(schur='always')
)
legacy_system <- cdrgam:::.cdr_factor_system(legacy$sparse$factor)
legacy$sparse$factor$cross <- lapply(legacy$sparse$factor$blocks, function(indices) {
    as.matrix(legacy_system[legacy$sparse$factor$core, indices, drop=FALSE])
})
legacy$sparse$factor$transfer <- NULL
stopifnot(max(abs(
    summary(legacy)$s.table[, 'F'] - grouped_tables$schur$s.table[, 'F']
)) < 1e-6)

# Sparse fits saved before smooths carried a random flag are recognized
# from their random-effect prediction entries.
unflagged <- cdrgam.fit(grouped_design, backend='sparse', method='REML')
flags <- cdrgam:::.cdrgam_random_smooths(unflagged)
unflagged$smooth <- lapply(unflagged$smooth, function(smooth) {
    smooth$random <- NULL
    smooth
})
stopifnot(
    identical(flags, c(TRUE, FALSE, FALSE)),
    identical(cdrgam:::.cdrgam_random_smooths(unflagged), flags),
    max(abs(summary(unflagged)$s.table[, 'F'] - grouped_tables$sparse$s.table[, 'F'])) <
        1e-8
)
