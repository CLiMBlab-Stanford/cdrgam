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

training_link <- predict(fit, type='link')
training_response <- predict(fit, type='response')
stream_data <- list(impulses=impulses, responses=responses['time'])
stream_link <- predict(fit, newdata=stream_data, type='link')
stream_response <- predict(fit, newdata=stream_data, type='response')
stopifnot(
    identical(dim(training_link), c(n_responses, 2L)),
    identical(dim(training_response), c(n_responses, 2L)),
    max(abs(training_link - stream_link)) < 1e-9,
    max(abs(training_response - stream_response)) < 1e-9,
    all(stream_response[, 2L] > 0)
)

with_se <- predict(fit, newdata=stream_data, type='response', se.fit=TRUE)
stopifnot(
    identical(dim(with_se$fit), c(n_responses, 2L)),
    identical(dim(with_se$se.fit), c(n_responses, 2L)),
    all(is.finite(with_se$se.fit)),
    all(with_se$se.fit >= 0)
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
stopifnot(
    nrow(location_irf) == 11L,
    nrow(scale_irf) == 11L,
    all(is.finite(location_irf$estimate)),
    all(is.finite(scale_irf$estimate))
)

fit_summary <- summary(fit)
stopifnot(
    inherits(fit_summary, 'summary.cdrgam'),
    identical(
        names(fit_summary$formula_strings$user),
        c('location', 'scale')
    ),
    identical(
        rownames(fit_summary$s.table),
        c('location:location_signal', 'scale:scale_signal')
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
    cdrgam.fit(design, family='gaulss', backend='sparse'),
    'currently require backend="mgcv"'
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

cat('Native distributional regression checks passed.\n')
