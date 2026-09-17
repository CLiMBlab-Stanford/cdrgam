library(cdrgam)

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
limited <- prepare_cdrgam(
    user_formula,
    impulses,
    responses,
    series='document',
    history='ragged',
    history_length=3,
    quiet=TRUE
)

stopifnot(inherits(dense, 'cdrgam_design'))
stopifnot(dense$plan[[1]]$selected == 'dense')
stopifnot(ragged$plan[[1]]$selected == 'ragged')
stopifnot(automatic$plan[[1]]$selected %in% c('dense', 'ragged'))
stopifnot(all(vapply(limited$plan, `[[`, numeric(1), 'maximum_history') <= 3))
stopifnot(identical(limited$stream$history_length, 3))
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

fit <- fit_cdrgam(
    dense,
    engine='gam',
    method='REML'
)
stopifnot(is_cdrgam(fit))
stopifnot(inherits(fit, 'gam'))
stopifnot(identical(cdr_formula(fit), user_formula))
stopifnot(grepl(
    'irf\\(1',
    paste(deparse(cdr_formula(fit, type='normalized')), collapse='')
))
stopifnot(grepl(
    'irf\\(1',
    paste(deparse(cdr_formula(fit, type='effective')), collapse='')
))
stopifnot(inherits(mgcv_formula(fit), 'formula'))
stopifnot(grepl('cdr_term_1', paste(deparse(mgcv_formula(fit)), collapse='')))
stopifnot(grepl('s\\(trial', paste(deparse(mgcv_formula(fit)), collapse='')))
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

fit_from_streams <- fit_cdrgam(
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

block_fit <- fit_cdrgam(
    dense,
    backend='block',
    method='REML'
)
stopifnot(is_cdrgam(block_fit))
stopifnot(inherits(block_fit, 'cdrgam_block'))
stopifnot(!inherits(block_fit, 'gam'))
stopifnot(isTRUE(all.equal(
    fitted(block_fit),
    fitted(fit),
    tolerance=2e-4
)))
stopifnot(all(is.finite(coef(block_fit))))
stopifnot(all(is.finite(vcov(block_fit))))
stopifnot(inherits(summary(block_fit), 'summary.cdrgam_block'))

poisson_error <- tryCatch(
    {
        fit_cdrgam(dense, backend='block', family=poisson())
        NA_character_
    },
    error=function(e) conditionMessage(e)
)
stopifnot(!is.na(poisson_error), grepl('only gaussian', poisson_error))
