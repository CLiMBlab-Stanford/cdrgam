library(cdrgam)

truth <- list(
    decay=function(lag) 1.2 * exp(-2.5 * lag),
    biphasic=function(lag) {
        0.9 * exp(-((lag - 0.35) / 0.2)^2) -
            0.5 * exp(-((lag - 1.05) / 0.3)^2)
    }
)
simulation <- simulate_cdr(
    truth,
    n_impulses=500,
    n_responses=650,
    duration=60,
    window=2,
    intercept=100,
    noise_sd=0.8,
    seed=4404
)

predictors <- as.matrix(simulation$impulses[names(truth)])
predictor_crossproduct <- crossprod(predictors) / nrow(predictors)
stopifnot(max(abs(predictor_crossproduct - diag(length(truth)))) < 1e-12)
stopifnot(any(diff(simulation$impulses$time) !=
    diff(simulation$impulses$time)[[1L]]))
stopifnot(any(diff(simulation$responses$time) !=
    diff(simulation$responses$time)[[1L]]))

formula <- response ~
    irf(decay, window=c(0, 2), k=10) +
    irf(biphasic, window=c(0, 2), k=10) - irf(1)
design <- prepare_cdrgam(
    formula,
    simulation$impulses,
    simulation$responses,
    history='auto',
    chunk_size=100,
    quiet=TRUE
)
native <- cdrgam.fit(design, backend='mgcv', engine='gam', method='REML')
block <- cdrgam.fit(design, backend='block', method='REML')
sparse <- cdrgam.fit(
    design,
    backend='sparse',
    method='REML',
    sparse_control=list(gradient='finite', crossprod_chunk_size=37)
)
sparse_exact <- cdrgam.fit(
    design,
    backend='sparse',
    method='REML',
    sparse_control=list(gradient='exact')
)
prediction_data <- list(
    impulses=simulation$impulses,
    responses=simulation$responses[names(simulation$responses) != 'response']
)
stopifnot(max(abs(predict(native, prediction_data) - fitted(native))) < 1e-8)
stopifnot(max(abs(predict(block, prediction_data) - fitted(block))) < 1e-8)
stopifnot(max(abs(predict(sparse, prediction_data) - fitted(sparse))) < 1e-8)
stopifnot(
    is.numeric(fitted(sparse)),
    is.numeric(residuals(sparse)),
    is.numeric(predict(sparse, prediction_data))
)
prediction_with_se <- predict(
    sparse,
    prediction_data,
    se.fit=TRUE
)
stopifnot(all(is.finite(prediction_with_se$se.fit)))
stopifnot(all(prediction_with_se$se.fit >= 0))
stopifnot(isTRUE(sparse$converged))
stopifnot(isTRUE(sparse$sparse$convergence$hessian_positive_definite))
stopifnot(length(sparse$sparse$convergence$boundary) == 0L)
stopifnot(isTRUE(all.equal(predict(sparse), fitted(sparse))))
sparse_summary <- summary(sparse)
stopifnot(inherits(sparse_summary, 'summary.cdrgam_sparse'))
stopifnot(abs(sum(sparse_summary$edf) - sum(summary(native)$edf)) < 1e-4)
stopifnot(
    is.null(sparse_summary$coefficients),
    nrow(sparse_summary$p.table) < length(coef(sparse)),
    nrow(sparse_summary$s.table) == length(sparse$smooth),
    identical(rownames(sparse_summary$s.table), sparse$cdrgam$term_labels),
    identical(names(sparse_summary$s.pv), sparse$cdrgam$term_labels),
    identical(
        sparse_summary$formula_strings$user,
        paste(deparse(formula(sparse, type='user')), collapse=' ')
    )
)
native_summary <- summary(native)
stopifnot(
    inherits(native_summary, 'summary.cdrgam'),
    identical(rownames(native_summary$s.table), native$cdrgam$term_labels),
    identical(
        native_summary$formula_strings$effective,
        paste(deparse(formula(native, type='effective')), collapse=' ')
    )
)
expanded_summary <- summary(sparse, all.coefficients=TRUE)
stopifnot(nrow(expanded_summary$coefficients) == length(coef(sparse)))
stopifnot(abs(
    sparse_summary$df.residual - native$df.residual
) < 1e-4)
stopifnot(is.finite(deviance(sparse)), is.finite(as.numeric(logLik(sparse))))
stopifnot(is.finite(AIC(sparse)), nobs(sparse) == nrow(simulation$responses))
native_covariance <- vcov(native)
sparse_covariance <- vcov(sparse)
stopifnot(max(abs(native_covariance - sparse_covariance)) < 1e-6)
stopifnot(
    max(abs(native_covariance - sparse_covariance)) /
        max(abs(native_covariance)) < 1e-4
)
native_unconditional <- vcov(native, unconditional=TRUE)
sparse_unconditional <- vcov(sparse, unconditional=TRUE)
stopifnot(
    max(abs(native_unconditional - sparse_unconditional)) /
        max(abs(native_unconditional)) < 0.03
)
stopifnot(all(diag(sparse_unconditional) >= diag(sparse_covariance) - 1e-12))

lag <- seq(0, 2, length.out=201)
native_irfs <- estimate_irf(native, lag=lag)
block_irfs <- estimate_irf(block, lag=lag)
sparse_irfs <- estimate_irf(sparse, lag=lag)
sparse_exact_irfs <- estimate_irf(sparse_exact, lag=lag)
for (term in names(truth)) {
    native_rows <- native_irfs$term == term
    block_rows <- block_irfs$term == term
    sparse_rows <- sparse_irfs$term == term
    sparse_exact_rows <- sparse_exact_irfs$term == term
    target <- truth[[term]](lag)
    native_rmse <- sqrt(mean((native_irfs$estimate[native_rows] - target)^2))
    block_rmse <- sqrt(mean((block_irfs$estimate[block_rows] - target)^2))
    native_correlation <- cor(native_irfs$estimate[native_rows], target)
    stopifnot(native_rmse < 0.18)
    stopifnot(block_rmse < 0.18)
    stopifnot(native_correlation > 0.9)
    stopifnot(max(abs(
        native_irfs$estimate[native_rows] -
            block_irfs$estimate[block_rows]
    )) < 2e-4)
    stopifnot(max(abs(
        native_irfs$estimate[native_rows] -
            sparse_irfs$estimate[sparse_rows]
    )) < 2e-4)
    stopifnot(max(abs(
        sparse_irfs$estimate[sparse_rows] -
            sparse_exact_irfs$estimate[sparse_exact_rows]
    )) < 2e-4)
    stopifnot(max(abs(
        native_irfs$se[native_rows] - sparse_irfs$se[sparse_rows]
    )) < 1e-5)
}

# Exercise the visual recovery path without retaining check artifacts.
plot_path <- tempfile(fileext='.png')
grDevices::png(plot_path, width=900, height=700)
graphics::par(mfrow=c(2, 1))
for (term in names(truth)) {
    rows <- native_irfs$term == term
    graphics::plot(
        lag,
        truth[[term]](lag),
        type='l',
        lwd=2,
        xlab='Lag',
        ylab='IRF',
        main=term
    )
    graphics::lines(lag, native_irfs$estimate[rows], col='blue', lwd=2)
}
grDevices::dev.off()
stopifnot(file.exists(plot_path), file.info(plot_path)$size > 0)
unlink(plot_path)
