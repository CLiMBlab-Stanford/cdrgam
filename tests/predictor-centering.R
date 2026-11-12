library(cdrgam)

set.seed(20261107)
n_impulse <- 240L
n_response <- 200L
impulses <- data.frame(
    series=rep(c('a', 'b'), each=n_impulse / 2L),
    time=c(sort(runif(n_impulse / 2L, 0, 60)), sort(runif(n_impulse / 2L, 0, 60))),
    x=round(runif(n_impulse, 1, 12)),
    z=rnorm(n_impulse, 3, 2),
    binary=rbinom(n_impulse, 1L, 0.3)
)
responses <- data.frame(
    series=rep(c('a', 'b'), each=n_response / 2L),
    time=c(seq(2, 60, length.out=n_response / 2L), seq(2, 60, length.out=n_response / 2L))
)
responses$response <- rnorm(n_response)

fit_model <- function(formula, ..., backend='mgcv') {
    suppressMessages(cdrgam(
        formula, impulses=impulses, responses=responses, window=c(0, 3),
        series='series', backend=backend, method='REML', engine='gam', ...
    ))
}
term <- function(fit, label) {
    labels <- fit$cdrgam$term_labels
    fit$cdrgam$terms[[match(label, labels)]]
}
linked_mean <- function(variable) {
    links <- getFromNamespace('.build_history_links', 'cdrgam')(
        impulses, responses, 'series', 'time', 'time', c(0, 3)
    )
    mean(impulses[[variable]][links$impulse_index])
}

# Linear predictors, including binary ones, are centered by default at their
# means over the training history links, in source units.
formula <- response ~ irf(x, k_l=5) + irf(binary, k_l=5) + irf(1, k_l=5)
centered <- fit_model(formula)
stopifnot(
    isTRUE(all.equal(unname(term(centered, 'x')$linear_predictor_centers), linked_mean('x'))),
    isTRUE(all.equal(unname(term(centered, 'binary')$linear_predictor_centers), linked_mean('binary'))),
    !length(term(centered, 'irf(1)')$linear_predictor_centers)
)

# Prediction on the training streams reproduces the fitted values.
predicted <- predict(centered, newdata=list(impulses=impulses, responses=responses))
stopifnot(isTRUE(all.equal(as.numeric(predicted), as.numeric(fitted(centered)), tolerance=1e-6)))

# A term-alone effect is zero at the center and is the deviation from an
# event at the center elsewhere.
at_center <- estimate_effect(
    centered, terms='x', axes=list(x=list(at=linked_mean('x'))), n=7L
)
stopifnot(max(abs(at_center$estimate)) < 1e-10)
catalog <- effect_catalog(centered)
stopifnot(
    isTRUE(all.equal(unname(catalog$centers[[match('x', catalog$label)]]), linked_mean('x'))),
    !length(catalog$centers[[match('irf(1)', catalog$label)]])
)

# Without centering, irf(1) is the response at zero. The two fits describe
# nearly the same total response, so the centered baseline approximates the
# uncentered baseline plus the predictor term at the centers.
uncentered <- fit_model(formula, center_predictors=FALSE)
stopifnot(
    !length(term(uncentered, 'x')$linear_predictor_centers),
    cor(as.numeric(fitted(centered)), as.numeric(fitted(uncentered))) > 0.99
)
lag_grid <- list(lag=list(values=seq(0, 3, length.out=9L)))
baseline <- function(fit) estimate_effect(fit, terms='irf(1)', axes=lag_grid)$estimate
term_at <- function(fit, variable) estimate_effect(
    fit, terms=variable,
    axes=c(lag_grid, stats::setNames(list(list(at=linked_mean(variable))), variable))
)$estimate
reconstructed <- baseline(uncentered) + term_at(uncentered, 'x') + term_at(uncentered, 'binary')
stopifnot(cor(baseline(centered), reconstructed) > 0.95)

# Term-level settings override the model-level default in both directions,
# and a character vector selects predictors within an interaction.
overridden <- fit_model(
    response ~ irf(x, k_l=5, center=FALSE) + irf(x, z, k_l=5, center='z') +
        irf(binary, k_l=5)
)
stopifnot(
    !length(term(overridden, 'x')$linear_predictor_centers),
    identical(names(term(overridden, 'x:z')$linear_predictor_centers), 'z'),
    identical(names(term(overridden, 'binary')$linear_predictor_centers), 'binary')
)
opted_in <- fit_model(
    response ~ irf(x, k_l=5, center=TRUE) + irf(binary, k_l=5),
    center_predictors=FALSE
)
logical_vector <- fit_model(
    response ~ irf(x, k_l=5) + irf(x, z, k_l=5, center=c(FALSE, TRUE))
)
stopifnot(identical(names(term(logical_vector, 'x:z')$linear_predictor_centers), 'z'))
stopifnot(
    identical(names(term(opted_in, 'x')$linear_predictor_centers), 'x'),
    !length(term(opted_in, 'binary')$linear_predictor_centers)
)
interaction <- estimate_effect(
    overridden, terms='x:z',
    axes=list(x=list(at=5), z=list(at=linked_mean('z'))), n=5L
)
stopifnot(max(abs(interaction$estimate)) < 1e-10)
predicted <- predict(overridden, newdata=list(impulses=impulses, responses=responses))
stopifnot(isTRUE(all.equal(as.numeric(predicted), as.numeric(fitted(overridden)), tolerance=1e-6)))

# Smooth predictor axes are unaffected; only their linear margins are centered.
mixed <- fit_model(
    response ~ irf(x, k_l=4, k_p=4) + irf(x, z, k_l=4, k_p=list(4, NULL))
)
stopifnot(identical(names(term(mixed, 'x:z')$linear_predictor_centers), 'z'))

# Centering does not depend on which other terms are present, so removing
# irf(1) or a lower-order term removes that component at average predictors.
prepare <- function(formula, ...) prepare_cdrgam(
    formula, impulses, responses, window=c(0, 3), series='series', quiet=TRUE, ...
)
no_rate <- prepare(response ~ irf(x, k_l=5) - irf(1))
interaction_only <- prepare(response ~ irf(x, z, k_l=5))
stopifnot(
    isTRUE(all.equal(unname(no_rate$terms$x$linear_predictor_centers), linked_mean('x'))),
    all(interaction_only$centering$centered),
    isTRUE(all.equal(
        interaction_only$centering$center[interaction_only$centering$predictor == 'z'],
        linked_mean('z')
    ))
)
constant_impulses <- impulses
constant_impulses$constant <- 2
single_value <- prepare_cdrgam(
    response ~ irf(constant, k_l=5) - irf(1), constant_impulses, responses,
    window=c(0, 3), series='series', quiet=TRUE
)
stopifnot(
    !single_value$centering$centered,
    identical(single_value$centering$reason, 'single linked value')
)
uncentered_design <- prepare(response ~ irf(x, k_l=5), center_predictors=FALSE)
stopifnot(identical(uncentered_design$centering$reason, 'center_predictors = FALSE'))

# Centers are stored in source units when predictors are rescaled.
rescaled <- fit_model(formula, rescale_predictors=TRUE)
stopifnot(
    isTRUE(all.equal(unname(term(rescaled, 'x')$linear_predictor_centers), linked_mean('x'))),
    isTRUE(all.equal(
        as.numeric(predict(rescaled, newdata=list(impulses=impulses, responses=responses))),
        as.numeric(fitted(rescaled)), tolerance=1e-6
    ))
)

# The custom backends store the same centers and predict consistently.
for (backend in c('block', 'sparse')) {
    fit <- fit_model(formula, backend=backend)
    stopifnot(
        isTRUE(all.equal(unname(term(fit, 'x')$linear_predictor_centers), linked_mean('x'))),
        isTRUE(all.equal(
            as.numeric(predict(fit, newdata=list(impulses=impulses, responses=responses))),
            as.numeric(fitted(fit)), tolerance=1e-6
        ))
    )
}

# The normalized formula records term-level settings.
formula_text <- paste(deparse(formula(overridden, type='normalized')), collapse='')
stopifnot(
    grepl('center = FALSE', formula_text, fixed=TRUE),
    grepl('center = c("z")', formula_text, fixed=TRUE)
)

expect_error <- function(expression, pattern) {
    message <- tryCatch({ expression; NA_character_ }, error=conditionMessage)
    stopifnot(grepl(pattern, message, fixed=TRUE))
}
expect_error(irf(1, center=TRUE), 'irf(1) has no predictors to center')
expect_error(irf(x, center='w'), 'center names predictors that are not in this IRF: w')
expect_error(irf(x, center=NA), 'center must be NULL, TRUE, FALSE, one logical value per predictor')
expect_error(irf(x, z, center=c(TRUE, FALSE, TRUE)), 'one logical value per predictor')
expect_error(fit_model(formula, center_predictors=NA), 'center_predictors must be TRUE or FALSE')

cat('predictor centering: ok\n')
