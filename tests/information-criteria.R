library(cdrgam)

# AIC and BIC of block and sparse fits use mgcv's log-likelihood and its
# smoothing-parameter-corrected effective degrees of freedom (edf2), so they
# reproduce AIC() and BIC() of the native fit at the same optimum.
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
responses <- do.call(rbind, grouped_responses)
responses$subject <- factor(responses$subject)
design <- prepare_cdrgam(
    response ~ s(subject, bs='re') + irf(x, window=c(0, 2), k=6) +
        irf(x, window=c(0, 2), k=6, group=subject) - irf(1),
    do.call(rbind, grouped_impulses), responses,
    series='subject', history='auto', quiet=TRUE
)
native <- as_gam(cdrgam.fit(design, backend='mgcv', method='REML', engine='gam'))
fits <- list(
    block=cdrgam.fit(design, backend='block', method='REML'),
    sparse=cdrgam.fit(design, backend='sparse', method='REML')
)
for (fit in fits) {
    edf2 <- cdrgam:::.cdrgam_edf2(fit)
    stopifnot(
        max(abs(edf2 - native$edf2)) < 1e-4,
        abs(as.numeric(logLik(fit)) - as.numeric(logLik(native))) < 1e-3,
        abs(AIC(fit) - AIC(native)) < 1e-2,
        abs(BIC(fit) - BIC(native)) < 1e-2,
        # The correction adds degrees of freedom to the plain edf.
        sum(edf2) > sum(fit$edf)
    )
}

# Several fits give a table, as stats::AIC() does.
table <- AIC(fits$block, fits$sparse)
stopifnot(
    identical(dim(table), c(2L, 2L)), identical(names(table), c('df', 'AIC')),
    abs(table$AIC[[1L]] - AIC(fits$block)) < 1e-10
)

# The correction runs or fails with its reason: above the size limit, and
# for GCV fits, it is an error rather than an uncorrected value.
failure <- function(expression, pattern) {
    message <- tryCatch({
        force(expression)
        NULL
    }, error=conditionMessage)
    !is.null(message) && grepl(pattern, message, fixed=TRUE)
}
old <- options(cdrgam.edf2_max_coefficients=5L)
stopifnot(failure(AIC(fits$sparse), 'cdrgam.edf2_max_coefficients = 5'))
options(old)
gcv <- cdrgam.fit(design, backend='sparse', method='GCV.Cp')
stopifnot(failure(AIC(gcv), 'needs a REML or ML fit'))
