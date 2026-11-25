library(cdrgam)

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
fit <- cdrgam.fit(design, backend='sparse', method='REML')
lags <- data.frame(lag=seq(0, 2, length.out=21))
grouped_terms <- list(predictors='x', grouped=TRUE)
population <- list(terms=list(predictors='x', grouped=FALSE), points=lags)
s1 <- list(terms=grouped_terms, grouping='conditional', groups='s1', points=lags)
s4 <- list(terms=grouped_terms, grouping='conditional', groups='s4', points=lags)

# With no estimand, every effect row is an estimand, and the estimates and
# standard errors are those of estimate_effect().
identity <- effect_inference(fit, list(population=population))
reference <- estimate_effect(fit, terms=population$terms, points=lags)
stopifnot(
    nrow(identity) == 21L, all(identity$method == 'wald'),
    max(abs(identity$estimate - reference$estimate)) < 1e-10,
    max(abs(identity$std_error - reference$se)) < 1e-10,
    max(abs(identity$lower - reference$lower)) < 1e-10,
    all(identity$simultaneous_lower <= identity$lower + 1e-12),
    all(identity$adjusted_p_value >= identity$p_value - 1e-12)
)

# The joint Wald test of all rows is the grid Wald test.
grid_test <- test_effect(fit, terms=population$terms, points=lags, method='wald')
joint <- attr(identity, 'joint_test')
stopifnot(
    joint$df == grid_test$df,
    abs(joint$statistic - grid_test$statistic) < 1e-6 * grid_test$statistic
)

# A contrast between two groups uses their joint covariance, which
# estimate_effect(covariance = TRUE) also gives.
contrast <- effect_inference(
    fit, list(s1=s1, s4=s4),
    estimand=list(difference_at_peak=list(s1=c(rep(0, 7), 1, rep(0, 13)),
        s4=c(rep(0, 7), -1, rep(0, 13))))
)
both <- estimate_effect(fit, terms=grouped_terms, grouping='conditional',
    groups=c('s1', 's4'), points=lags, covariance=TRUE)
joint_covariance <- attr(both, 'covariance')
stopifnot(
    identical(contrast$estimand, 'difference_at_peak'),
    abs(contrast$estimate - (both$estimate[[8]] - both$estimate[[29]])) < 1e-10,
    abs(contrast$std_error^2 - (joint_covariance[8, 8] + joint_covariance[29, 29] -
        2 * joint_covariance[8, 29])) < 1e-10,
    contrast$std_error < sqrt(joint_covariance[8, 8] + joint_covariance[29, 29])
)

# A function estimand that happens to be linear: the delta method recovers
# the Wald result, and simulation agrees within Monte Carlo error.
average <- function(values, grids) c(area=mean(values$population))
wald_area <- effect_inference(fit, list(population=population),
    estimand=list(area=list(population=1 / 21)))
delta_area <- effect_inference(fit, list(population=population),
    estimand=average, method='delta')
simulated_area <- effect_inference(fit, list(population=population),
    estimand=average, draws=4000L, seed=3)
stopifnot(
    abs(delta_area$estimate - wald_area$estimate) < 1e-10,
    abs(delta_area$std_error - wald_area$std_error) < 1e-6 * wald_area$std_error,
    identical(simulated_area$method, 'simulation'),
    abs(simulated_area$std_error / wald_area$std_error - 1) < 0.06,
    abs(simulated_area$lower - wald_area$lower) < 0.1 * wald_area$std_error * 4
)

# A nonlinear estimand: the lag at which each group's response peaks and
# their difference, all from the same draws.
peaks <- function(values, grids) {
    lag <- grids$s1$lag
    c(s1=lag[[which.max(values$s1)]], s4=lag[[which.max(values$s4)]],
        difference=lag[[which.max(values$s1)]] - lag[[which.max(values$s4)]])
}
peak <- effect_inference(fit, list(s1=s1, s4=s4), estimand=peaks, seed=5,
    draws=500L, keep_draws=TRUE)
stopifnot(
    identical(peak$estimand, c('s1', 's4', 'difference')),
    all(peak$lower <= peak$upper),
    identical(dim(attr(peak, 'draws')), c(3L, 500L)),
    all(abs(attr(peak, 'draws')['difference', ] - (attr(peak, 'draws')['s1', ] -
        attr(peak, 'draws')['s4', ])) < 1e-12)
)

# A seed makes simulation reproducible and leaves the random state alone.
set.seed(99)
before <- .Random.seed
again <- effect_inference(fit, list(s1=s1, s4=s4), estimand=peaks, seed=5,
    draws=500L)
stopifnot(identical(.Random.seed, before),
    identical(again$std_error, peak$std_error))

# Mismatched estimands and methods fail.
failure <- function(expression) {
    inherits(tryCatch(expression, error=function(error) error), 'error')
}
stopifnot(
    failure(effect_inference(fit, list(population=population), method='simulation')),
    failure(effect_inference(fit, list(population=population), estimand=average,
        method='wald')),
    failure(effect_inference(fit, list(population=population), null=c(0, 1))),
    failure(effect_inference(fit, list(population)))
)
