library(cdrgam)

set.seed(20260917)
n_impulse <- 72L
n_response <- 60L
origin <- 1.7e12
impulses <- data.frame(
    series=rep(c('a', 'b'), each=n_impulse / 2L),
    time=origin + c(seq(0, 70, length.out=n_impulse / 2L),
        seq(0, 70, length.out=n_impulse / 2L)),
    x=runif(n_impulse, -30, 45),
    binary=rep(c(0, 1), length.out=n_impulse)
)
responses <- data.frame(
    series=rep(c('a', 'b'), each=n_response / 2L),
    time=origin + c(seq(5, 70, length.out=n_response / 2L),
        seq(5, 70, length.out=n_response / 2L)),
    z=seq(100, 900, length.out=n_response),
    binary_response=rep(c(0, 1), length.out=n_response),
    group=factor(rep(letters[1:3], length.out=n_response)),
    response=rnorm(n_response)
)

model_formula <- response ~ z + binary_response + s(time, k=5) +
    irf(x, window=c(0, 8), k=c(6, 5), nonlinear=TRUE) +
    irf(binary, window=c(0, 8), k=c(6, 5), nonlinear=TRUE) - irf(1)

native <- prepare_cdrgam(
    model_formula,
    impulses,
    responses,
    series='series',
    history='ragged',
    quiet=TRUE
)
scaled <- prepare_cdrgam(
    model_formula,
    impulses,
    responses,
    series='series',
    history='ragged',
    rescale_predictors=TRUE,
    quiet=TRUE
)

metadata <- scaled$scaling$variables
divisor <- function(stream, variable) {
    row <- metadata$stream == stream & metadata$variable == variable
    stopifnot(sum(row) == 1L)
    metadata$divisor[row]
}
applied <- function(stream, variable) {
    row <- metadata$stream == stream & metadata$variable == variable
    stopifnot(sum(row) == 1L)
    metadata$applied[row]
}

stopifnot(
    isTRUE(scaled$scaling$enabled),
    identical(native$scaling$enabled, FALSE),
    isTRUE(all.equal(scaled$scaling$time_divisor, sd(impulses$time))),
    isTRUE(all.equal(divisor('impulses', 'x'), sd(impulses$x))),
    isTRUE(all.equal(divisor('responses', 'z'), sd(responses$z))),
    isTRUE(all.equal(
        divisor('responses', 'time'),
        scaled$scaling$time_divisor
    )),
    applied('impulses', 'x'),
    applied('responses', 'z'),
    applied('responses', 'time'),
    !applied('impulses', 'binary'),
    !applied('responses', 'binary_response'),
    !applied('impulses', 'series'),
    !applied('responses', 'series'),
    identical(scaled$responses$binary_response, responses$binary_response),
    identical(scaled$responses$group, responses$group),
    isTRUE(all.equal(
        scaled$responses$z,
        responses$z / sd(responses$z)
    )),
    identical(
        paste(deparse(formula(native, 'user')), collapse=''),
        paste(deparse(formula(scaled, 'user')), collapse='')
    ),
    identical(
        paste(deparse(formula(native, 'normalized')), collapse=''),
        paste(deparse(formula(scaled, 'normalized')), collapse='')
    ),
    identical(
        paste(deparse(formula(native, 'effective')), collapse=''),
        paste(deparse(formula(scaled, 'effective')), collapse='')
    ),
    identical(native$plan, scaled$plan)
)

# The binary nonlinear request is simplified using its original coding, and
# scaling does not turn it into a continuous predictor.
binary_index <- which(vapply(
    scaled$specification,
    function(specification) identical(specification$predictor, 'binary'),
    logical(1)
))
stopifnot(
    length(binary_index) == 1L,
    !isTRUE(scaled$specification[[binary_index]]$nonlinear),
    identical(scaled$terms[[binary_index]]$amplitude_scale, 1)
)

# Link membership is computed in native time, then the resulting lag and any
# explicit time predictor use one common divisor.
x_index <- which(vapply(
    scaled$specification,
    function(specification) identical(specification$predictor, 'x'),
    logical(1)
))
stopifnot(
    length(x_index) == 1L,
    isTRUE(all.equal(
        scaled$terms[[x_index]]$lag_scale,
        scaled$scaling$time_divisor
    )),
    isTRUE(all.equal(
        scaled$terms[[x_index]]$predictor_scale,
        sd(impulses$x)
    ))
)

# Prediction accepts native-unit streams and reuses the training divisors.
scaled_fit <- cdrgam.fit(scaled, backend='mgcv', engine='gam', method='REML')
native_fit <- cdrgam.fit(native, backend='mgcv', engine='gam', method='REML')
prediction_data <- responses[names(responses) != 'response']
scaled_prediction <- predict(
    scaled_fit,
    newdata=list(impulses=impulses, responses=prediction_data)
)
stopifnot(max(abs(scaled_prediction - fitted(scaled_fit))) < 1e-7)
stopifnot(max(abs(fitted(native_fit) - fitted(scaled_fit))) < 2e-4)

# Parametric summary coefficients and their uncertainty are converted back to
# source units; test all fitting paths that expose coefficient tables.
native_table <- summary(native_fit)$p.table
scaled_table <- summary(scaled_fit)$p.table
stopifnot(
    abs(native_table['z', 'Estimate'] - scaled_table['z', 'Estimate']) < 1e-7,
    abs(native_table['z', 'Std. Error'] -
        scaled_table['z', 'Std. Error']) < 1e-7
)
for (backend in c('block', 'sparse')) {
    control <- if (identical(backend, 'sparse')) {
        list(
            outer_optimizer='lbfgsb', optimizer_maxit=20, hessian='none'
        )
    } else list()
    custom_fit <- cdrgam.fit(
        scaled,
        backend=backend,
        rank_action='minimum_norm',
        sparse_control=control
    )
    custom_prediction <- predict(
        custom_fit,
        newdata=list(impulses=impulses, responses=prediction_data)
    )
    custom_table <- summary(custom_fit)$p.table
    stopifnot(
        max(abs(custom_prediction - fitted(custom_fit))) < 1e-7,
        abs(custom_table['z', 'Estimate'] -
            native_table['z', 'Estimate']) < 1e-7
    )
}

# Public IRF coordinates are restored to source units.
x_irf <- estimate_irf(scaled_fit, term=x_index, n=11, n_predictor=7, se=FALSE)
native_x_irf <- estimate_irf(
    native_fit,
    term=x_index,
    n=11,
    n_predictor=7,
    se=FALSE
)
stopifnot(
    min(x_irf$lag) >= -1e-10,
    max(x_irf$lag) <= 8 + 1e-8,
    min(x_irf$predictor) >= min(impulses$x) - 1e-8,
    max(x_irf$predictor) <= max(impulses$x) + 1e-8,
    max(abs(x_irf$lag - native_x_irf$lag)) < 1e-10,
    max(abs(x_irf$predictor - native_x_irf$predictor)) < 1e-10,
    max(abs(x_irf$estimate - native_x_irf$estimate)) < 1e-4
)

# The translated GAM view restores ordinary smooth axes while evaluating the
# fitted smooth on its internal coordinates.
plot_gam <- cdrgam:::.cdr_source_scale_gam_plot(scaled_fit)
time_smooth <- which(vapply(
    plot_gam$smooth,
    function(smooth) identical(smooth$label, 's(time)'),
    logical(1)
))
time_panel <- mgcv:::plot.mgcv.smooth(
    plot_gam$smooth[[time_smooth]],
    data=plot_gam$model,
    n=11
)
stopifnot(
    length(time_smooth) == 1L,
    max(abs(range(time_panel$x) - range(responses$time))) < 1e-8
)

# Invalid training-time scales fail early; scaling remains opt-in.
constant_time <- impulses
constant_time$time <- 1
time_error <- try(prepare_cdrgam(
    response ~ irf(x, window=c(0, 8), k=5) - irf(1),
    constant_time,
    responses,
    series='series',
    rescale_predictors=TRUE,
    quiet=TRUE
), silent=TRUE)
stopifnot(
    inherits(time_error, 'try-error'),
    grepl('zero or non-finite standard deviation', time_error, fixed=TRUE)
)
