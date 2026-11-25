library(cdrgam)

# Wood whole-functional tests of grouped IRF effects. Each test of a single
# term reproduces that term's row of the smooth table, on every backend.
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
impulses <- do.call(rbind, grouped_impulses)
responses <- do.call(rbind, grouped_responses)
responses$subject <- factor(responses$subject)
design <- prepare_cdrgam(
    response ~ s(subject, bs='re') + irf(x, window=c(0, 2), k=6) +
        irf(x, window=c(0, 2), k=6, group=subject),
    impulses, responses, series='subject', history='auto', quiet=TRUE
)
fits <- list(
    native=cdrgam.fit(design, backend='mgcv', method='REML'),
    block=cdrgam.fit(design, backend='block', method='REML'),
    sparse=cdrgam.fit(design, backend='sparse', method='REML')
)
population_terms <- list(predictors='x', grouped=FALSE)
grouped_terms <- list(predictors='x', grouped=TRUE)
for (fit in fits) {
    table <- summary(fit)$s.table
    statistic <- colnames(table)[[3L]]
    population <- test_effect(fit, terms=population_terms, method='wood')
    deviation <- test_effect(
        fit, terms=grouped_terms, grouping='deviation', method='wood'
    )
    conditional <- test_effect(
        fit, terms=grouped_terms, grouping='conditional', method='wood'
    )
    total <- test_effect(
        fit, terms=population_terms, composition='total', method='wood'
    )
    matches <- function(test, label) {
        abs(test$df - table[label, 'Ref.df']) < 1e-8 &&
            abs(test$statistic - table[label, statistic]) <
                1e-8 * max(1, table[label, statistic]) &&
            abs(test$p_value - table[label, 'p-value']) < 1e-10
    }
    ranges <- lapply(fit$smooth, function(smooth) smooth$first.para:smooth$last.para)
    names(ranges) <- rownames(table)
    stopifnot(
        identical(population$method, 'Wood (2013)'), matches(population, 'x'),
        identical(deviation$method, 'Wood (2013b) random effect'),
        is.na(deviation$group), matches(deviation, 'x|subject'),
        identical(conditional$group, levels(responses$subject)),
        all(conditional$method == 'Wood (2013)'),
        all(conditional$coefficients == length(ranges$x) +
            length(ranges$`x|subject`) / 4L),
        total$coefficients == length(ranges$x) + length(ranges$`irf(1)`),
        is.na(total$points), all(is.na(population$points))
    )
}

# A total composition is tested over the union of its terms' coefficients
# with the summed edf1 as the rank, which mgcv's testStat() reproduces from
# the native fit's R factor.
native <- fits$native
union <- c(ranges$`irf(1)`, ranges$x)
expected <- mgcv:::testStat(
    coef(native)[union], native$R[, union, drop=FALSE], native$Vp[union, union],
    min(length(union), sum(native$edf1[union])), type=0,
    res.df=native$df.residual
)
total <- test_effect(native, terms=population_terms, composition='total',
    method='wood')
stopifnot(
    abs(total$p_value - expected$pval) < 1e-10,
    abs(total$df - expected$rank) < 1e-10,
    abs(total$residual_df - native$df.residual) < 1e-10
)

# The sparse fit reaches the native smoothing parameters and gives the same
# total test. (The block fit stops at a smoother irf(1) term.)
sparse_total <- test_effect(fits$sparse, terms=population_terms,
    composition='total', method='wood')
stopifnot(
    abs(sparse_total$df - total$df) < 1e-3,
    abs(sparse_total$statistic - total$statistic) < 1e-3 * total$statistic
)

# The test runs or fails explicitly; it never falls back.
failure <- function(expression) {
    inherits(tryCatch(expression, error=function(error) error), 'error')
}
stopifnot(
    failure(test_effect(native, terms=population_terms, method='wood',
        axes=list(lag=list(grid='fitted', n=5L)))),
    failure(test_effect(native, terms=population_terms, method='wood',
        reference=list(lag=0))),
    failure(test_effect(native, terms=grouped_terms, grouping='deviation',
        groups='s1', method='wood')),
    failure(test_effect(native, terms=grouped_terms, grouping='deviation',
        method='wood', unconditional=TRUE))
)

# Wood's test is the default; arguments that define a grid select the Wald
# test unless a method is given.
default <- test_effect(native, terms=population_terms)
wald <- test_effect(native, terms=population_terms, method='wald')
gridded <- test_effect(native, terms=population_terms, axes=list(lag=list(grid='fitted', n=20L)))
stopifnot(
    identical(default$method, 'Wood (2013)'),
    identical(wald$method, 'Wald'), wald$points == 50L, is.na(wald$residual_df),
    identical(gridded$method, 'Wald'), gridded$points == 20L
)
