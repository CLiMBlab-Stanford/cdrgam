library(cdrgam)

set.seed(9341)
n_impulses <- 180L
n_responses <- 260L
impulses <- data.frame(
    time=sort(stats::runif(n_impulses, 0, 45)),
    location_signal=stats::rnorm(n_impulses),
    scale_signal=stats::rnorm(n_impulses)
)
responses <- data.frame(time=sort(stats::runif(n_responses, 0.5, 45)))
responses$response <- stats::rnorm(n_responses)

formulas <- list(
    location=response ~
        irf(location_signal, window=c(0, 1.5), k_l=6) - irf(1),
    scale=~ irf(scale_signal, window=c(0, 1.5), k_l=6) - irf(1)
)
design <- prepare_cdrgam(
    formulas,
    impulses,
    responses,
    history='ragged',
    chunk_size=47,
    quiet=TRUE
)
stopifnot(
    inherits(design, 'cdrgam_distributional_design'),
    identical(names(design$parameters), c('location', 'scale')),
    identical(design$response_name, 'response')
)

fit <- cdrgam.fit(
    design,
    family='gaulss',
    backend='mgcv',
    method='REML'
)
one_step <- cdrgam(
    formulas,
    impulses,
    responses,
    family='gaulss',
    backend='mgcv',
    method='REML',
    history='ragged',
    chunk_size=47
)
stopifnot(
    inherits(fit, 'cdrgam'),
    inherits(fit, 'gam'),
    isTRUE(fit$cdrgam$distributional),
    identical(fit$cdrgam$engine, 'gam'),
    identical(one_step$cdrgam$engine, 'gam'),
    max(abs(stats::coef(one_step) - stats::coef(fit))) < 1e-9,
    identical(
        fit$cdrgam$term_labels,
        c('location:location_signal', 'scale:scale_signal')
    )
)

block <- cdrgam.fit(
    design,
    family='gaulss',
    backend='block',
    method='REML'
)
sparse <- cdrgam.fit(
    design,
    family='gaulss',
    backend='sparse',
    method='REML'
)
stopifnot(
    inherits(block, 'cdrgam_distributional_block'),
    isTRUE(block$converged),
    identical(block$cdrgam$parameter_names, c('location', 'scale')),
    max(abs(unname(stats::coef(block)) - unname(stats::coef(fit)))) < 2e-4,
    max(abs(unname(stats::vcov(block)) - unname(stats::vcov(fit)))) < 2e-4,
    abs(as.numeric(stats::logLik(block)) - as.numeric(stats::logLik(fit))) <
        2e-3,
    abs(stats::deviance(block) - stats::deviance(fit)) < 1e-4,
    inherits(sparse, 'cdrgam_distributional_sparse'),
    isTRUE(sparse$converged),
    max(abs(unname(stats::coef(sparse)) - unname(stats::coef(block)))) < 2e-4,
    abs(as.numeric(stats::logLik(sparse)) - as.numeric(stats::logLik(block))) <
        2e-3
)

training_link <- predict(fit, type='link')
training_response <- predict(fit, type='response')
stream_data <- list(impulses=impulses, responses=responses['time'])
stream_link <- predict(fit, newdata=stream_data, type='link')
stream_response <- predict(fit, newdata=stream_data, type='response')
block_link <- predict(block, newdata=stream_data, type='link')
block_response <- predict(block, newdata=stream_data, type='response')
sparse_response <- predict(sparse, newdata=stream_data, type='response')
stopifnot(
    identical(dim(training_link), c(n_responses, 2L)),
    identical(dim(training_response), c(n_responses, 2L)),
    max(abs(training_link - stream_link)) < 1e-9,
    max(abs(training_response - stream_response)) < 1e-9,
    max(abs(block_link - stream_link)) < 2e-4,
    max(abs(block_response - stream_response)) < 2e-4,
    max(abs(sparse_response - block_response)) < 2e-4,
    all(stream_response[, 2L] > 0)
)

with_se <- predict(fit, newdata=stream_data, type='response', se.fit=TRUE)
stopifnot(
    identical(dim(with_se$fit), c(n_responses, 2L)),
    identical(dim(with_se$se.fit), c(n_responses, 2L)),
    all(is.finite(with_se$se.fit)),
    all(with_se$se.fit >= 0)
)
block_with_se <- predict(
    block, newdata=stream_data, type='response', se.fit=TRUE
)
sparse_with_se <- predict(
    sparse, newdata=stream_data, type='response', se.fit=TRUE
)
stopifnot(
    max(abs(block_with_se$fit - with_se$fit)) < 2e-4,
    max(abs(block_with_se$se.fit - with_se$se.fit)) < 2e-4,
    max(abs(sparse_with_se$fit - block_with_se$fit)) < 2e-4,
    max(abs(sparse_with_se$se.fit - block_with_se$se.fit)) < 2e-4
)

reference <- mgcv::gam(
    unname(fit$cdrgam$formula$mgcv),
    family=mgcv::gaulss(),
    data=as.list(fit$model),
    method='REML'
)
stopifnot(
    max(abs(stats::coef(fit) - stats::coef(reference))) < 1e-9,
    max(abs(stats::vcov(fit) - stats::vcov(reference))) < 1e-9
)

location_irf <- estimate_irf(
    fit, term='location:location_signal', n=11
)
scale_irf <- estimate_irf(fit, term='scale:scale_signal', n=11)
block_location_irf <- estimate_irf(
    block, term='location:location_signal', n=11
)
block_scale_irf <- estimate_irf(
    block, term='scale:scale_signal', n=11
)
sparse_location_irf <- estimate_irf(
    sparse, term='location:location_signal', n=11
)
sparse_scale_irf <- estimate_irf(
    sparse, term='scale:scale_signal', n=11
)
stopifnot(
    nrow(location_irf) == 11L,
    nrow(scale_irf) == 11L,
    all(is.finite(location_irf$estimate)),
    all(is.finite(scale_irf$estimate)),
    max(abs(block_location_irf$estimate - location_irf$estimate)) < 2e-4,
    max(abs(block_scale_irf$estimate - scale_irf$estimate)) < 2e-4,
    max(abs(sparse_location_irf$estimate - block_location_irf$estimate)) <
        2e-4,
    max(abs(sparse_scale_irf$estimate - block_scale_irf$estimate)) < 2e-4
)

fit_summary <- summary(fit)
block_summary <- summary(block)
sparse_summary <- summary(sparse)
stopifnot(
    inherits(fit_summary, 'summary.cdrgam'),
    identical(rownames(sparse_summary$s.table), rownames(fit_summary$s.table)),
    identical(
        names(fit_summary$formula_strings$user),
        c('location', 'scale')
    ),
    identical(
        rownames(fit_summary$s.table),
        c('location:location_signal', 'scale:scale_signal')
    ),
    identical(rownames(block_summary$s.table), rownames(fit_summary$s.table)),
    identical(
        rownames(block_summary$p.table),
        c('location:(Intercept)', 'scale:(Intercept)')
    )
)

smooth_formulas <- list(
    location=response ~ stats::offset(location_offset) +
        s(time, k=5) +
        irf(location_signal, window=c(0, 1.5), k_l=6) - irf(1),
    scale=~ stats::offset(scale_offset) + s(time, k=5) +
        irf(scale_signal, window=c(0, 1.5), k_l=6) - irf(1)
)
responses$location_offset <- seq(-0.1, 0.1, length.out=n_responses)
responses$scale_offset <- seq(0.05, -0.05, length.out=n_responses)
smooth_design <- prepare_cdrgam(
    smooth_formulas,
    impulses,
    responses,
    history='ragged',
    chunk_size=47,
    quiet=TRUE
)
smooth_native <- cdrgam.fit(
    smooth_design,
    family='gaulss',
    backend='mgcv',
    engine='gam',
    method='REML'
)
smooth_block <- cdrgam.fit(
    smooth_design,
    family='gaulss',
    backend='block',
    method='REML'
)
smooth_sparse <- cdrgam.fit(
    smooth_design,
    family='gaulss',
    backend='sparse',
    method='REML'
)
smooth_streams <- list(
    impulses=impulses,
    responses=responses[c('time', 'location_offset', 'scale_offset')]
)
smooth_block_se <- predict(
    smooth_block, smooth_streams, type='response', se.fit=TRUE
)$se.fit
smooth_native_se <- predict(
    smooth_native, smooth_streams, type='response', se.fit=TRUE
)$se.fit
smooth_sparse_prediction <- predict(
    smooth_sparse, smooth_streams, type='response', se.fit=TRUE
)
stopifnot(
    isTRUE(smooth_block$converged),
    max(abs(
        unname(stats::coef(smooth_block)) -
            unname(stats::coef(smooth_native))
    )) < 6e-3,
    max(abs(
        unname(stats::vcov(smooth_block)) -
            unname(stats::vcov(smooth_native))
    )) < 5e-2,
    max(abs(
        predict(smooth_block, smooth_streams, type='response') -
            predict(smooth_native, smooth_streams, type='response')
    )) < 6e-3,
    max(abs(smooth_block_se - smooth_native_se)) < 6e-3,
    abs(
        as.numeric(stats::logLik(smooth_block)) -
            as.numeric(stats::logLik(smooth_native))
    ) < 2e-3,
    isTRUE(smooth_sparse$converged),
    max(abs(
        smooth_sparse_prediction$fit -
            predict(smooth_block, smooth_streams, type='response')
    )) < 6e-3,
    max(abs(smooth_sparse_prediction$se.fit - smooth_block_se)) < 6e-3,
    identical(
        rownames(summary(smooth_block)$s.table),
        c(
            'location:s(time)', 'location:location_signal',
            'scale:s(time)', 'scale:scale_signal'
        )
    )
)

expect_error <- function(expression, pattern) {
    message <- tryCatch(
        {
            force(expression)
            NA_character_
        },
        error=function(error) conditionMessage(error)
    )
    stopifnot(!is.na(message), grepl(pattern, message, fixed=TRUE))
}
expect_error(
    prepare_cdrgam(
        list(mean=formulas$location, scale=formulas$scale),
        impulses,
        responses,
        quiet=TRUE
    ),
    'requires formulas named location and scale'
)
expect_error(
    cdrgam.fit(design, family='gaussian', engine='gam'),
    'requires family="gaulss"'
)
expect_error(
    prepare_cdrgam(
        formulas,
        impulses,
        responses,
        rescale_predictors=TRUE,
        quiet=TRUE
    ),
    'rescaling is not yet supported'
)
expect_error(
    cdrgam.fit(
        design,
        family='gaulss',
        backend='mgcv',
        weights=rep(c(1, 2), length.out=n_responses)
    ),
    'does not support non-unit prior weights'
)

cat('Native distributional regression checks passed.\n')
