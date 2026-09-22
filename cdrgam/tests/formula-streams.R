library(cdrgam)

exports <- getNamespaceExports('cdrgam')
stopifnot(
    all(c('cdrgam', 'cdrgam.fit', 'cdrgam_family') %in% exports),
    !any(c(
        'fit_cdrgam', 'fit_compressed_cdr_gam', 'predict_cdrgam',
        'cdr_formula', 'mgcv_formula', 'compress_cdr_smooth'
    ) %in% exports)
)

set.seed(20260917)
impulses <- data.frame(
    document=rep(c('a', 'b'), each=30),
    time=c(seq(0, 5.8, length.out=30), seq(0, 5.8, length.out=30)),
    surprisal=rnorm(60)
)
responses <- data.frame(
    document=rep(c('a', 'b'), each=20),
    time=c(seq(0.35, 6, length.out=20), seq(0.42, 6.07, length.out=20)),
    trial=rep(seq(0, 1, length.out=20), 2)
)
responses$rt <- 400 + 25 * responses$trial + rnorm(nrow(responses), sd=5)

warning_responses <- responses
warning_responses$document <- factor(warning_responses$document)
warning_responses$position <- factor(ifelse(
    warning_responses$document == 'a',
    rep(c('first', 'second'), length.out=nrow(warning_responses)),
    'first'
))
interaction_warning <- character()
invisible(withCallingHandlers(
    prepare_cdrgam(
        rt ~ s(document, position, bs='re') +
            irf(surprisal, window=c(0, 1.25), k=6),
        impulses,
        warning_responses,
        series='document',
        quiet=TRUE
    ),
    warning=function(condition) {
        interaction_warning <<- c(
            interaction_warning, conditionMessage(condition)
        )
        invokeRestart('muffleWarning')
    }
))
stopifnot(any(grepl(
    'unattested factor combination', interaction_warning, fixed=TRUE
)))

user_formula <- rt ~ s(trial, k=5) +
    irf(surprisal, window=c(0, 1.25), k=6)

dense <- prepare_cdrgam(
    user_formula,
    impulses,
    responses,
    series='document',
    history='dense',
    chunk_size=7,
    quiet=TRUE
)
ragged <- prepare_cdrgam(
    user_formula,
    impulses,
    responses,
    series='document',
    history='ragged',
    chunk_size=13,
    quiet=TRUE
)
automatic <- prepare_cdrgam(
    user_formula,
    impulses,
    responses,
    series='document',
    history='auto',
    quiet=TRUE
)
stopifnot(inherits(dense, 'cdrgam_design'))
stopifnot(dense$plan[[1]]$selected == 'dense')
stopifnot(ragged$plan[[1]]$selected == 'ragged')
stopifnot(automatic$plan[[1]]$selected %in% c('dense', 'ragged'))
expected_links <- sum(vapply(seq_len(nrow(responses)), function(i) {
    same_series <- impulses$document == responses$document[[i]]
    delay <- responses$time[[i]] - impulses$time
    sum(same_series & delay >= 0 & delay <= 1.25)
}, integer(1)))
stopifnot(
    identical(dense$plan[[1]]$links, expected_links),
    !('history_length' %in% names(formals(prepare_cdrgam))),
    !('history_length' %in% names(formals(cdrgam)))
)
stopifnot(!('history_length' %in% names(dense$stream)))
stopifnot(isTRUE(all.equal(
    dense$terms[[1]]$X,
    ragged$terms[[1]]$X,
    tolerance=1e-11
)))
stopifnot(isTRUE(all.equal(
    dense$terms[[1]]$S,
    ragged$terms[[1]]$S,
    tolerance=1e-11
)))

future_formula <- rt ~ irf(surprisal, window=c(-0.8, -0.05), k=5) - irf(1)
future_dense <- prepare_cdrgam(
    future_formula,
    impulses,
    responses,
    series='document',
    history='dense',
    quiet=TRUE
)
future_ragged <- prepare_cdrgam(
    future_formula,
    impulses,
    responses,
    series='document',
    history='ragged',
    quiet=TRUE
)
stopifnot(
    all(future_dense$terms[[1L]]$knots < 0),
    isTRUE(all.equal(
        future_dense$terms[[1L]]$X,
        future_ragged$terms[[1L]]$X,
        tolerance=1e-11
    ))
)
future_fit <- cdrgam.fit(future_dense, engine='gam', method='REML')
future_prediction <- predict(
    future_fit,
    list(
        impulses=impulses,
        responses=responses[names(responses) != 'rt']
    )
)
stopifnot(max(abs(future_prediction - fitted(future_fit))) < 1e-8)

fit <- cdrgam.fit(
    dense,
    engine='gam',
    method='REML'
)
stopifnot(is_cdrgam(fit))
stopifnot(inherits(fit, 'gam'))
stopifnot(identical(fit$call[[1L]], as.name('cdrgam.fit')))
stopifnot(identical(formula(fit), user_formula))
stopifnot(grepl(
    'irf\\(1',
    paste(deparse(formula(fit, type='normalized')), collapse='')
))
stopifnot(grepl(
    'irf\\(1',
    paste(deparse(formula(fit, type='effective')), collapse='')
))
stopifnot(inherits(formula(fit, type='mgcv'), 'formula'))
stopifnot(grepl(
    'cdr_term_1',
    paste(deparse(formula(fit, type='mgcv')), collapse='')
))
stopifnot(grepl(
    's\\(trial',
    paste(deparse(formula(fit, type='mgcv')), collapse='')
))
stopifnot(length(fitted(fit)) == nrow(responses))
stopifnot(all(is.finite(vcov(fit))))
stopifnot(identical(fit$cdrgam$identifiability$rate$status, 'retained'))
stopifnot(identical(fit$cdrgam$identifiability$global$resolution, 'identified'))
stopifnot(identical(fit$cdrgam$term_labels[[1L]], 'irf(1)'))
stream_prediction <- predict(
    fit,
    list(
        impulses=impulses,
        responses=responses[names(responses) != 'rt']
    )
)
stopifnot(max(abs(stream_prediction - fitted(fit))) < 1e-8)

fit_from_streams <- cdrgam(
    user_formula,
    impulses=impulses,
    responses=responses,
    series='document',
    history='ragged',
    chunk_size=13,
    engine='gam',
    method='REML'
)
stopifnot(isTRUE(all.equal(
    fitted(fit_from_streams),
    fitted(fit),
    tolerance=1e-8
)))
stopifnot(identical(fit_from_streams$call[[1L]], as.name('cdrgam')))

block_fit <- cdrgam.fit(
    dense,
    backend='block',
    method='REML'
)
stopifnot(is_cdrgam(block_fit))
stopifnot(inherits(block_fit, 'cdrgam_block'))
stopifnot(identical(summary(block_fit)$call, block_fit$call))
stopifnot(!inherits(block_fit, 'gam'))
stopifnot(isTRUE(all.equal(
    fitted(block_fit),
    fitted(fit),
    tolerance=2e-4
)))
stopifnot(all(is.finite(coef(block_fit))))
stopifnot(all(is.finite(vcov(block_fit))))
stopifnot(inherits(summary(block_fit), 'summary.cdrgam_block'))

noncanonical_error <- tryCatch(
    {
        cdrgam.fit(dense, backend='block', family=binomial(link='probit'))
        NA_character_
    },
    error=function(e) conditionMessage(e)
)
stopifnot(
    !is.na(noncanonical_error),
    grepl('currently supports', noncanonical_error, fixed=TRUE)
)
