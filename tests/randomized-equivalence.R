library(cdrgam)

# Seeded differential tests exercise independently randomized signals,
# irregular sampling patterns, history layouts, and basis dimensions. The
# three fitters solve the same Gaussian REML problem by different numerical
# routes, so agreement here detects drift in data compilation as well as in
# fitting, prediction, and IRF evaluation.
for (case in seq_len(3L)) {
    seed <- 9100L + case
    set.seed(seed)
    amplitude <- runif(2L, 0.35, 1.1) * c(1, -1)
    rate <- runif(2L, 0.9, 2.2)
    peak <- runif(1L, 0.35, 0.85)
    width <- runif(1L, 0.2, 0.4)
    truth <- list(
        x1=local({
            a <- amplitude[[1L]]
            r <- rate[[1L]]
            function(lag) a * exp(-r * lag)
        }),
        x2=local({
            a <- amplitude[[2L]]
            p <- peak
            w <- width
            function(lag) a * exp(-((lag - p) / w)^2)
        })
    )
    simulation <- simulate_cdr(
        truth,
        n_impulses=150L + 10L * case,
        n_responses=190L + 10L * case,
        duration=25 + case,
        window=1.6,
        intercept=runif(1L, 5, 15),
        noise_sd=0.45,
        seed=seed
    )
    basis_dimension <- 5L + case %% 2L
    formula <- response ~
        irf(x1, window=c(0, 1.6), k=basis_dimension) +
        irf(x2, window=c(0, 1.6), k=basis_dimension) - irf(1)
    dense <- prepare_cdrgam(
        formula,
        simulation$impulses,
        simulation$responses,
        history='dense',
        chunk_size=17L + case,
        quiet=TRUE
    )
    ragged <- prepare_cdrgam(
        formula,
        simulation$impulses,
        simulation$responses,
        history='ragged',
        chunk_size=29L - case,
        quiet=TRUE
    )
    for (term in seq_along(dense$terms)) {
        stopifnot(isTRUE(all.equal(
            dense$terms[[term]]$X,
            ragged$terms[[term]]$X,
            tolerance=2e-12
        )))
    }

    native <- cdrgam.fit(dense, backend='mgcv', engine='gam', method='REML')
    block <- cdrgam.fit(ragged, backend='block', method='REML')
    sparse <- cdrgam.fit(
        ragged,
        backend='sparse',
        method='REML',
        sparse_control=list(
            gradient='exact',
            crossprod_chunk_size=31L + case
        )
    )
    stopifnot(isTRUE(sparse$converged))
    stopifnot(max(abs(fitted(native) - fitted(block))) < 5e-4)
    stopifnot(max(abs(fitted(native) - fitted(sparse))) < 5e-4)
    stopifnot(max(abs(coef(native) - coef(block))) < 5e-4)
    stopifnot(max(abs(coef(native) - coef(sparse))) < 5e-4)

    prediction_data <- list(
        impulses=simulation$impulses,
        responses=simulation$responses[names(simulation$responses) != 'response']
    )
    stopifnot(max(abs(predict(native, prediction_data) - fitted(native))) < 1e-8)
    stopifnot(max(abs(predict(block, prediction_data) - fitted(block))) < 1e-8)
    stopifnot(max(abs(predict(sparse, prediction_data) - fitted(sparse))) < 1e-8)

    evaluation_lags <- seq(0, 1.6, length.out=37L)
    native_irf <- estimate_irf(native, lag=evaluation_lags)
    block_irf <- estimate_irf(block, lag=evaluation_lags)
    sparse_irf <- estimate_irf(sparse, lag=evaluation_lags)
    stopifnot(max(abs(native_irf$estimate - block_irf$estimate)) < 5e-4)
    stopifnot(max(abs(native_irf$estimate - sparse_irf$estimate)) < 5e-4)
}
