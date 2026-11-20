library(cdrgam)

# Hierarchical composition, references, and joint tests of effect functionals.
simulation <- simulate_cdr(
    list(x=function(lag) exp(-lag), z=function(lag) 0.5 * exp(-2 * lag)),
    n_impulses=300,
    n_responses=250,
    duration=40,
    window=2,
    seed=1103
)
design <- prepare_cdrgam(
    response ~ irf(x, k_l=5) + irf(z, k_l=5, k_p=3) +
        irf(x, z, k_l=5, k_p=list(NULL, 3)),
    simulation$impulses,
    simulation$responses,
    window=c(0, 2),
    quiet=TRUE
)
fit <- cdrgam.fit(design, backend='mgcv', method='REML')
labels <- fit$cdrgam$term_labels
interaction <- which(labels == 'x:z')
stopifnot(length(interaction) == 1L)
stopifnot(setequal(
    labels[cdrgam:::.cdr_effect_hierarchy(fit, interaction)],
    c('irf(1)', 'x', 'z')
))

lags <- seq(0.05, 1.95, length.out=9L)
z_values <- c(-1, 0, 1)
axes <- list(lag=lags, predictors=list(x=1.5, z=z_values))
evaluate <- function(term, composition='term', reference=list(), se=TRUE) {
    estimate_effect(
        fit, terms=term, axes=axes, composition=composition,
        reference=reference, se=se
    )
}
total <- evaluate('x:z', 'total')
stopifnot(nrow(total) == length(lags) * length(z_values))

# The total effect is the sum of the hierarchy evaluated on the same grid.
pieces <- lapply(c('x:z', 'x', 'z', 'irf(1)'), function(term) {
    value <- evaluate(term, se=FALSE)
    key <- paste(value$lag, if ('z' %in% names(value)) value$z else NA)
    stats::setNames(value$estimate, key)
})
key <- paste(total$lag, total$z)
summed <- pieces[[1L]][key] + pieces[[2L]][paste(total$lag, NA)] +
    pieces[[3L]][key] + pieces[[4L]][paste(total$lag, NA)]
stopifnot(isTRUE(all.equal(unname(summed), total$estimate, tolerance=1e-10)))

# Standard errors use the joint covariance of all summed coefficients.
terms <- match(c('x:z', 'x', 'z', 'irf(1)'), labels)
grid <- total[c('lag', 'x', 'z')]
columns <- sort(unique(unlist(lapply(
    terms, function(term) fit$cdrgam$terms[[term]]$coefficient_index
))))
L <- matrix(0, nrow(grid), length(columns))
for (term in terms) {
    info <- fit$cdrgam$terms[[term]]
    destination <- match(info$coefficient_index, columns)
    L[, destination] <- L[, destination] + cdrgam:::.cdr_effect_basis(info, grid)
}
V <- stats::vcov(fit)[columns, columns]
stopifnot(isTRUE(all.equal(
    sqrt(rowSums((L %*% V) * L)), total$se, tolerance=1e-8
)))

# A reference value of x cancels terms that do not involve x.
contrast <- evaluate('x:z', 'total', reference=list(x=0.5))
stopifnot(
    'x_reference' %in% names(contrast),
    all(contrast$x_reference == 0.5),
    all(contrast$se > 0)
)
partial <- function(value) {
    axes$predictors$x <- value
    a <- estimate_effect(fit, terms='x:z', axes=axes, se=FALSE)$estimate
    b <- estimate_effect(fit, terms='x', axes=axes, se=FALSE)
    a + rep(b$estimate, times=length(z_values))
}
stopifnot(isTRUE(all.equal(
    contrast$estimate, partial(1.5) - partial(0.5), tolerance=1e-10
)))

# A reference value of z gives differences between interaction lines.
level_difference <- evaluate('x:z', 'total', reference=list(z=-1))
stopifnot(
    all(abs(level_difference$estimate[level_difference$z == -1]) < 1e-10),
    all(level_difference$se[level_difference$z == -1] < 1e-10)
)
reference <- total[total$z == -1, 'estimate']
stopifnot(isTRUE(all.equal(
    level_difference$estimate[level_difference$z == 1],
    total$estimate[total$z == 1] - reference,
    tolerance=1e-10
)))

# A predictor left at the model default drops every term that involves it from
# the reference. With every predictor at its default, the reference is irf(1).
rate <- evaluate('irf(1)', se=FALSE)$estimate
all_default <- evaluate('x:z', 'total', reference=list(x=NA, z='default'))
stopifnot(
    all(is.na(all_default$x_reference)), all(is.na(all_default$z_reference)),
    isTRUE(all.equal(
        all_default$estimate, total$estimate - rep(rate, times=length(z_values)),
        tolerance=1e-10
    ))
)
# A centered linear predictor's default equals its center.
center <- fit$cdrgam$terms[[match('x', labels)]]$linear_predictor_centers[['x']]
at_center <- evaluate('x:z', 'total', reference=list(x=center))
x_default <- evaluate('x:z', 'total', reference=list(x=NA))
stopifnot(isTRUE(all.equal(x_default$estimate, at_center$estimate, tolerance=1e-8)))
tryCatch({
    evaluate('x:z', 'total', reference=list(lag=NA))
    stop('expected a lag default error')
}, error=function(error) stopifnot(grepl('no model default', conditionMessage(error))))

# Joint Wald tests over the grid.
tested <- test_effect(
    fit, terms='x:z', axes=list(lag=list(grid='fitted', n=20L),
        predictors=list(x=1.5, z=z_values)),
    composition='total', reference=list(x=0.5)
)
stopifnot(
    inherits(tested, 'cdrgam_effect_test'),
    nrow(tested) == 1L,
    tested$df >= 1L,
    is.finite(tested$statistic),
    tested$p_value >= 0 && tested$p_value <= 1
)
null_test <- test_effect(
    fit, terms='x:z', axes=list(lag=lags, predictors=list(x=1.5, z=0)),
    composition='total', reference=list(x=1.5)
)
stopifnot(null_test$df == 0L, is.na(null_test$statistic))
tryCatch({
    estimate_effect(fit, terms='x', axes=axes, reference=list(z=0))
    stop('expected an unknown reference axis error')
}, error=function(error) stopifnot(grepl('not evaluation axes', conditionMessage(error))))

cat('Effect composition tests passed\n')
