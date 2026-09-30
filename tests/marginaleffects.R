library(cdrgam)

simulation <- simulate_cdr(
    list(x=function(lag, value) {
        (0.4 + 0.25 * value + 0.12 * value^2) * exp(-lag)
    }),
    n_impulses=80,
    n_responses=70,
    duration=15,
    window=1.5,
    intercept=2,
    noise_sd=0.2,
    seed=922
)

design <- prepare_cdrgam(
    response ~ irf(x, k_l=5, k_p=4),
    simulation$impulses,
    simulation$responses,
    window=c(0, 1.5),
    quiet=TRUE
)
fit <- cdrgam.fit(design, backend='mgcv', method='REML')

view <- marginaleffects_view(
    fit,
    simulation$impulses,
    simulation$responses
)
evaluation <- attr(view, 'marginaleffects_modeldata')
stopifnot(
    inherits(view, 'cdrgam_marginaleffects'),
    nrow(evaluation) > nrow(simulation$responses),
    nrow(response_scenarios(view)) == nrow(evaluation),
    all(c(
        '.cdrgam_pair', '.cdrgam_response', '.cdrgam_impulse',
        '.cdrgam_lag', 'x'
    ) %in% names(evaluation)),
    identical(all.vars(formula(view)), c('.cdrgam_outcome', 'x'))
)

baseline <- predict(view, type='link')
stream_baseline <- predict(
    fit,
    newdata=list(
        impulses=simulation$impulses,
        responses=simulation$responses
    ),
    type='link'
)
stopifnot(max(abs(
    baseline - stream_baseline[evaluation$.cdrgam_response]
)) < 1e-9)

row <- which(evaluation$x != 0)[[1L]]
changed <- evaluation[row, , drop=FALSE]
changed$x <- changed$x + 0.15
adapter_prediction <- predict(view, newdata=changed, type='link')

changed_impulses <- simulation$impulses
impulse <- changed$.cdrgam_impulse[[1L]]
changed_impulses$x[[impulse]] <- changed$x[[1L]]
direct_prediction <- predict(
    fit,
    newdata=list(
        impulses=changed_impulses,
        responses=simulation$responses
    ),
    type='link'
)[changed$.cdrgam_response[[1L]]]
stopifnot(abs(adapter_prediction - direct_prediction) < 1e-8)

# Repeated incidences form independent scenario rows, including counterfactual
# lag values. Prediction and uncertainty use the same scenario design.
scenarios <- response_scenarios(
    view,
    at=list(lag=c(0.2, 0.8), x=c(-0.5, 0.5)),
    rows=row
)
stopifnot(
    nrow(scenarios) == 4L,
    identical(scenarios$.cdrgam_scenario, 1:4),
    length(unique(scenarios$.cdrgam_pair)) == 1L
)
scenario_prediction <- predict(
    view, newdata=scenarios, type='link', se.fit=TRUE
)
stopifnot(
    length(scenario_prediction$fit) == 4L,
    length(scenario_prediction$se.fit) == 4L,
    all(is.finite(scenario_prediction$fit)),
    all(is.finite(scenario_prediction$se.fit)),
    all(scenario_prediction$se.fit >= 0)
)
scenario_design <- predict(view, newdata=scenarios, type='lpmatrix')
scenario_link <- as.numeric(
    scenario_design %*% coef(fit) +
        view$context$offset[view$context$pair_response[row]]
)
stopifnot(max(abs(scenario_prediction$fit - scenario_link)) < 1e-9)
for (scenario in seq_len(nrow(scenarios))) {
    current <- scenarios[scenario, , drop=FALSE]
    changed_impulses <- simulation$impulses
    changed_impulses$x[[current$.cdrgam_impulse]] <- current$x
    response_row <- current$.cdrgam_response
    response_time <- simulation$responses$time[[response_row]]
    changed_impulses$time[[current$.cdrgam_impulse]] <-
        response_time - current$.cdrgam_lag
    direct <- predict(
        fit,
        newdata=list(
            impulses=changed_impulses,
            responses=simulation$responses
        ),
        type='link'
    )[[response_row]]
    stopifnot(abs(scenario_prediction$fit[[scenario]] - direct) < 1e-8)
}

backend_fits <- list(
    block=cdrgam.fit(design, backend='block', method='REML'),
    sparse=cdrgam.fit(
        design,
        backend='sparse',
        method='REML',
        sparse_control=list(hessian='optimhess')
    )
)
for (backend_fit in backend_fits) {
    backend_view <- marginaleffects_view(
        backend_fit,
        simulation$impulses,
        simulation$responses
    )
    backend_data <- attr(backend_view, 'marginaleffects_modeldata')
    backend_changed <- backend_data[row, , drop=FALSE]
    backend_changed$x <- backend_changed$x + 0.15
    backend_impulses <- simulation$impulses
    backend_impulses$x[[backend_changed$.cdrgam_impulse]] <- backend_changed$x
    stopifnot(abs(
        predict(backend_view, newdata=backend_changed, type='link') -
        predict(
            backend_fit,
            newdata=list(
                impulses=backend_impulses,
                responses=simulation$responses
            ),
            type='link'
        )[[backend_changed$.cdrgam_response]]
    ) < 1e-8)
    backend_scenarios <- response_scenarios(
        backend_view, at=list(x=c(-0.2, 0.2)), rows=row
    )
    backend_scenario_prediction <- predict(
        backend_view, newdata=backend_scenarios,
        type='link', se.fit=TRUE
    )
    stopifnot(
        length(backend_scenario_prediction$fit) == 2L,
        all(is.finite(backend_scenario_prediction$se.fit))
    )
}

grouped_responses <- simulation$responses
grouped_responses$subject <- factor(rep(
    c('s1', 's2', 's3'), length.out=nrow(grouped_responses)
))
grouped_fit <- cdrgam(
    response ~
        irf(x, k_l=5) +
        irf(x, k_l=5, group=subject) -
        irf(1),
    simulation$impulses,
    grouped_responses,
    window=c(0, 1.5),
    backend='mgcv',
    method='REML',
    quiet=TRUE
)
grouped_view <- marginaleffects_view(
    grouped_fit,
    simulation$impulses,
    grouped_responses
)
grouped_data <- attr(grouped_view, 'marginaleffects_modeldata')
grouped_changed <- grouped_data[row, , drop=FALSE]
grouped_changed$x <- grouped_changed$x + 0.2
grouped_impulses <- simulation$impulses
grouped_impulses$x[[grouped_changed$.cdrgam_impulse]] <- grouped_changed$x
stopifnot(abs(
    predict(grouped_view, newdata=grouped_changed, type='link') -
    predict(
        grouped_fit,
        newdata=list(
            impulses=grouped_impulses,
            responses=grouped_responses
        ),
        type='link'
    )[[grouped_changed$.cdrgam_response]]
) < 1e-8)

varying_responses <- simulation$responses
varying_responses$exposure <- seq(
    0.5, 1.5, length.out=nrow(varying_responses)
)
varying_fit <- cdrgam(
    response ~ irf(
        x,
        k_l=4,
        k_t=3,
        k_p=3,
        by=exposure
    ),
    simulation$impulses,
    varying_responses,
    window=c(0, 1.5),
    backend='mgcv',
    method='REML',
    quiet=TRUE
)
varying_view <- marginaleffects_view(
    varying_fit,
    simulation$impulses,
    varying_responses
)
varying_data <- attr(varying_view, 'marginaleffects_modeldata')
varying_changed <- varying_data[row, , drop=FALSE]
varying_changed$x <- varying_changed$x - 0.2
varying_impulses <- simulation$impulses
varying_impulses$x[[varying_changed$.cdrgam_impulse]] <- varying_changed$x
stopifnot(abs(
    predict(varying_view, newdata=varying_changed, type='link') -
    predict(
        varying_fit,
        newdata=list(
            impulses=varying_impulses,
            responses=varying_responses
        ),
        type='link'
    )[[varying_changed$.cdrgam_response]]
) < 1e-8)

# A predictor change also updates every selected IRF containing that predictor.
simulation$impulses$z <- sin(simulation$impulses$time)
two_term_fit <- cdrgam(
    response ~
        irf(x, k_l=5) +
        irf(x, z, k_l=5, k_p=list(NULL, NULL)),
    simulation$impulses,
    simulation$responses,
    window=c(0, 1.5),
    backend='mgcv',
    method='REML',
    quiet=TRUE
)
two_term_view <- marginaleffects_view(
    two_term_fit,
    simulation$impulses,
    simulation$responses
)
two_term_data <- attr(two_term_view, 'marginaleffects_modeldata')
changed <- two_term_data[row, , drop=FALSE]
changed$x <- changed$x - 0.1
changed_impulses <- simulation$impulses
changed_impulses$x[[changed$.cdrgam_impulse]] <- changed$x
two_term_adapter <- predict(two_term_view, newdata=changed, type='link')
two_term_direct <- predict(
        two_term_fit,
        newdata=list(
            impulses=changed_impulses,
            responses=simulation$responses
        ),
        type='link'
    )[[changed$.cdrgam_response]]
stopifnot(abs(two_term_adapter - two_term_direct) < 1e-8)

if (requireNamespace('marginaleffects', quietly=TRUE)) {
    row_plot_data <- estimate_response_effect(
        view, newdata=scenarios, type='link', vcov=FALSE
    )
    stopifnot(nrow(row_plot_data) == nrow(scenarios))
    plot_data <- estimate_response_effect(
        view,
        newdata=scenarios,
        by=c('.cdrgam_scenario', '.cdrgam_lag', 'x'),
        type='link'
    )
    stopifnot(
        inherits(plot_data, 'cdrgam_response_effect_grid'),
        nrow(plot_data) == 4L,
        all(c('estimate', 'std.error', 'conf.low', 'conf.high') %in%
            names(plot_data)),
        identical(attr(plot_data, 'estimand'), 'prediction'),
        identical(attr(plot_data, 'scale'), 'link')
    )
    slope_plot_data <- estimate_response_effect(
        view,
        newdata=scenarios,
        estimand='slope',
        variables='x',
        by=c('.cdrgam_scenario', '.cdrgam_lag', 'x'),
        type='link'
    )
    stopifnot(
        nrow(slope_plot_data) == 4L,
        all(is.finite(slope_plot_data$estimate))
    )
    inferred <- marginaleffects::slopes(
        view,
        newdata=evaluation[row, , drop=FALSE],
        vcov=FALSE,
        type='link'
    )
    stopifnot(identical(inferred$term, 'x'))
    slope <- marginaleffects::slopes(
        view,
        variables='x',
        newdata=evaluation[row, , drop=FALSE],
        vcov=FALSE,
        type='link'
    )
    epsilon <- 1e-5
    high <- evaluation[row, , drop=FALSE]
    low <- high
    high$x <- high$x + epsilon / 2
    low$x <- low$x - epsilon / 2
    expected <- (
        predict(view, newdata=high, type='link') -
        predict(view, newdata=low, type='link')
    ) / epsilon
    stopifnot(abs(slope$estimate - expected) < 1e-5)
    slope_with_uncertainty <- marginaleffects::slopes(
        view,
        variables='x',
        newdata=evaluation[row, , drop=FALSE],
        vcov=TRUE,
        type='link'
    )
    stopifnot(
        is.finite(slope_with_uncertainty$std.error),
        slope_with_uncertainty$std.error >= 0
    )
    comparison_values <- evaluation$x[[row]] + c(-0.1, 0.1)
    comparison <- marginaleffects::comparisons(
        view,
        variables=list(x=comparison_values),
        newdata=evaluation[row, , drop=FALSE],
        vcov=FALSE,
        type='link'
    )
    low <- high <- evaluation[row, , drop=FALSE]
    low$x <- comparison_values[[1L]]
    high$x <- comparison_values[[2L]]
    expected_comparison <-
        predict(view, newdata=high, type='link') -
        predict(view, newdata=low, type='link')
    stopifnot(abs(comparison$estimate - expected_comparison) < 1e-8)
    sparse_view <- marginaleffects_view(
        backend_fits$sparse,
        simulation$impulses,
        simulation$responses
    )
    sparse_data <- attr(sparse_view, 'marginaleffects_modeldata')
    sparse_slope <- marginaleffects::slopes(
        sparse_view,
        variables='x',
        newdata=sparse_data[row, , drop=FALSE],
        vcov=TRUE,
        type='link'
    )
    stopifnot(is.finite(sparse_slope$std.error))
}

cat('marginaleffects adapter tests passed\n')
