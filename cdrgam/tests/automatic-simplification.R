library(cdrgam)

set.seed(20260917)
impulses <- data.frame(
    time=0:19,
    binary=rep(c(0, 1), 10),
    constant=rep(2, 20),
    continuous=seq(-1, 1, length.out=20)
)
responses <- data.frame(
    time=2:21,
    response=rnorm(20)
)
window <- c(0, 4)

binary <- prepare_cdrgam(
    response ~ irf(
        binary,
        window=window,
        k=c(40, 8),
        nonlinear=TRUE
    ) - irf(1),
    impulses,
    responses,
    history='ragged',
    quiet=TRUE
)
binary_specification <- binary$specification[[1L]]
binary_links <- cdrgam:::.build_history_links(
    impulses,
    responses,
    character(),
    'time',
    'time',
    window
)
delay_count <- length(unique(binary_links$delay))

manual <- prepare_cdrgam(
    stats::as.formula(paste0(
        'response ~ irf(binary, window=c(0, 4), k=',
        binary_specification$k, ') - irf(1)'
    )),
    impulses,
    responses,
    history='ragged',
    quiet=TRUE
)

stopifnot(
    !isTRUE(binary_specification$nonlinear),
    identical(binary_specification$k, as.integer(delay_count)),
    identical(binary_specification$bs, 'cr'),
    is.null(binary_specification$predictor_center),
    identical(binary$terms[[1L]]$type, 'linear'),
    isTRUE(all.equal(
        binary$terms[[1L]]$X,
        manual$terms[[1L]]$X,
        tolerance=1e-12
    )),
    all(c('lag', 'predictor') %in% binary$simplifications$axis),
    all(c('basis_reduced', 'nonlinear_to_linear') %in%
        binary$simplifications$action),
    grepl(
        'k_p *= *list\\(8\\)',
        paste(deparse(formula(binary, 'normalized')), collapse='')
    ),
    grepl(
        'k_p *= *list\\(NULL\\)',
        paste(deparse(formula(binary, 'effective')), collapse='')
    )
)

continuous <- prepare_cdrgam(
    response ~ irf(
        continuous,
        window=window,
        k=c(40, 40),
        nonlinear=TRUE
    ) - irf(1),
    impulses,
    responses,
    history='ragged',
    quiet=TRUE
)
continuous_specification <- continuous$specification[[1L]]
linked_continuous <- impulses$continuous[binary_links$impulse_index]
stopifnot(
    isTRUE(continuous_specification$nonlinear),
    identical(
        continuous_specification$k,
        as.integer(c(
            length(unique(binary_links$delay)),
            length(unique(linked_continuous))
        ))
    ),
    nrow(continuous$simplifications) == 2L,
    all(continuous$simplifications$action == 'basis_reduced')
)

constant <- prepare_cdrgam(
    response ~ irf(
        constant,
        window=window,
        k=c(40, 8),
        nonlinear=TRUE
    ) - irf(1),
    impulses,
    responses,
    history='ragged',
    quiet=TRUE
)
stopifnot(
    !isTRUE(constant$specification[[1L]]$nonlinear),
    is.null(constant$specification[[1L]]$predictor_center),
    identical(constant$terms[[1L]]$type, 'linear'),
    any(grepl(
        '1 distinct linked predictor value',
        constant$simplifications$reason,
        fixed=TRUE
    ))
)

constant_with_rate <- suppressWarnings(prepare_cdrgam(
    response ~ irf(
        constant,
        window=window,
        k=c(40, 8),
        nonlinear=TRUE
    ),
    impulses,
    responses,
    history='ragged',
    quiet=TRUE
))
stopifnot(
    identical(names(constant_with_rate$terms), 'constant'),
    identical(constant_with_rate$identifiability$rate$status, 'removed'),
    any(constant_with_rate$simplifications$action == 'term_removed'),
    any(grepl(
        'duplicates explicit term constant',
        constant_with_rate$simplifications$reason,
        fixed=TRUE
    ))
)

# Prediction uses the original binary coding without a hidden transformation.
binary_fit <- cdrgam.fit(binary, backend='mgcv', engine='gam', method='REML')
binary_prediction <- predict(
    binary_fit,
    newdata=list(
        impulses=impulses,
        responses=responses[names(responses) != 'response']
    )
)
stopifnot(max(abs(binary_prediction - fitted(binary_fit))) < 1e-8)
stopifnot(identical(
    binary_fit$cdrgam$preparation$simplifications,
    binary$simplifications
))
