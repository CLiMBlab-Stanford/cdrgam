library(cdrgam)

simulation <- simulate_cdr(
    list(A=function(lag) exp(-lag)),
    n_impulses=140,
    n_responses=120,
    duration=24,
    window=2,
    seed=812
)
simulation$impulses$B <- sin(seq_len(nrow(simulation$impulses)) / 9)

formula <- response ~
    irf(1, k_l=5) +
    irf(A, k_l=5) +
    irf(
        A, B,
        k_l=5,
        k_t=4,
        k_p=list(NULL, 4),
        bs_l='tp',
        bs_t='ps',
        bs_p=list('cr', 'cs')
    )
design <- prepare_cdrgam(
    formula,
    simulation$impulses,
    simulation$responses,
    window=c(-0.25, 2),
    history='ragged',
    quiet=TRUE
)

stopifnot(
    all(vapply(design$specification, function(specification) {
        identical(specification$window, c(-0.25, 2))
    }, logical(1))),
    identical(design$specification[[3L]]$predictors, c('A', 'B')),
    identical(design$specification[[3L]]$k_p, list(NULL, 4L)),
    identical(design$specification[[3L]]$time, 'time'),
    identical(vapply(design$terms[[3L]]$axis, `[[`, character(1), 'role'),
        c('lag', 'time', 'predictor')),
    all(c('deviation-centering:time', 'deviation-centering:B') %in%
        design$terms[[3L]]$constraints),
    ncol(design$terms[[3L]]$X) == 5L * 3L * 3L
)

normalized <- gsub(
    '[[:space:]]+', ' ',
    paste(deparse(design$normalized_formula), collapse=' ')
)
stopifnot(
    grepl('k_l = 5', normalized, fixed=TRUE),
    grepl('k_t = 4', normalized, fixed=TRUE),
    grepl('k_p = list(NULL, 4)', normalized, fixed=TRUE)
)

native <- cdrgam.fit(design, backend='mgcv', engine='gam', method='REML')
block <- cdrgam.fit(design, backend='block', method='REML')
sparse <- cdrgam.fit(
    design,
    backend='sparse',
    method='REML',
    sparse_control=list(gradient='finite')
)
stopifnot(
    max(abs(fitted(native) - fitted(block))) < 1e-5,
    max(abs(fitted(native) - fitted(sparse))) < 1e-5
)

newdata <- list(
    impulses=simulation$impulses,
    responses=simulation$responses[setdiff(
        names(simulation$responses), 'response'
    )]
)
stopifnot(max(abs(predict(sparse, newdata) - fitted(sparse))) < 1e-7)

surface <- estimate_irf(
    sparse,
    term=3,
    lag=c(0, 1),
    predictor=c(-0.5, 0.5),
    at=list(B=c(-0.5, 0.5), time=12),
    se=FALSE
)
stopifnot(
    all(c('lag', 'predictor', 'time', 'B') %in% names(surface)),
    nrow(surface) == 4L
)

plot_path <- tempfile(fileext='.pdf')
grDevices::pdf(plot_path)
plot(sparse, view='surface', select=3, se=FALSE)
grDevices::dev.off()
stopifnot(file.info(plot_path)$size > 0)
unlink(plot_path)

for (basis in c('cr', 'cs', 'cc', 'tp', 'ts', 'ps')) {
    candidate <- prepare_cdrgam(
        stats::as.formula(paste0(
            'response ~ irf(A, k_l=6, bs_l="', basis, '") - irf(1)'
        )),
        simulation$impulses,
        simulation$responses,
        window=c(0, 2),
        history='ragged',
        quiet=TRUE
    )
    stopifnot(ncol(candidate$terms[[1L]]$X) > 0L)
}

focused_knots <- c(0, 0.05, 0.15, 0.4, 1, 2)
custom_dense <- prepare_cdrgam(
    response ~ irf(A, k_l=6, knots_l=focused_knots) - irf(1),
    simulation$impulses,
    simulation$responses,
    window=c(0, 2),
    history='dense',
    rescale_predictors=TRUE,
    quiet=TRUE
)
custom_ragged <- prepare_cdrgam(
    response ~ irf(A, k_l=6, knots_l=focused_knots) - irf(1),
    simulation$impulses,
    simulation$responses,
    window=c(0, 2),
    history='ragged',
    rescale_predictors=TRUE,
    quiet=TRUE
)
stopifnot(
    identical(custom_dense$specification[[1L]]$knots_l, focused_knots),
    max(abs(custom_dense$terms[[1L]]$knots *
        custom_dense$terms[[1L]]$lag_scale - focused_knots)) < 1e-12,
    max(abs(custom_dense$terms[[1L]]$X - custom_ragged$terms[[1L]]$X)) < 1e-10,
    max(abs(custom_dense$terms[[1L]]$S[[1L]] -
        custom_ragged$terms[[1L]]$S[[1L]])) < 1e-10
)
custom_fit <- cdrgam.fit(custom_ragged, backend='mgcv', engine='gam', method='REML')
custom_block <- cdrgam.fit(custom_ragged, backend='block', method='REML')
custom_sparse <- cdrgam.fit(
    custom_ragged,
    backend='sparse',
    method='REML',
    sparse_control=list(gradient='finite')
)
custom_estimate <- estimate_irf(custom_fit, term=1, n=17, se=FALSE)
stopifnot(
    max(abs(fitted(custom_fit) - fitted(custom_block))) < 2e-5,
    max(abs(fitted(custom_fit) - fitted(custom_sparse))) < 2e-5,
    min(custom_estimate$lag) == min(focused_knots),
    max(custom_estimate$lag) == max(focused_knots),
    grepl(
        'knots_l = c(0, 0.05, 0.15, 0.4, 1, 2)',
        gsub('[[:space:]]+', ' ', paste(
            deparse(custom_ragged$normalized_formula), collapse=' '
        )),
        fixed=TRUE
    )
)

for (basis in c('cr', 'cs', 'cc', 'tp', 'ts')) {
    custom_basis <- prepare_cdrgam(
        stats::as.formula(paste0(
            'response ~ irf(A, k_l=6, knots_l=focused_knots, bs_l="',
            basis, '") - irf(1)'
        )),
        simulation$impulses,
        simulation$responses,
        window=c(0, 2),
        history='ragged',
        quiet=TRUE
    )
    stopifnot(max(abs(custom_basis$terms[[1L]]$knots - focused_knots)) < 1e-12)
}

invalid_ps_knots <- tryCatch({
    irf(A, k_l=6, bs_l='ps', knots_l=focused_knots)
    NULL
}, error=identity)
stopifnot(
    inherits(invalid_ps_knots, 'error'),
    grepl('not supported', conditionMessage(invalid_ps_knots), fixed=TRUE)
)
stopifnot(
    inherits(try(
        irf(A, k_l=6, knots_l=focused_knots[-1L]),
        silent=TRUE
    ), 'try-error'),
    inherits(try(
        irf(A, k_l=6, knots_l=rev(focused_knots)),
        silent=TRUE
    ), 'try-error'),
    inherits(try(
        prepare_cdrgam(
            response ~ irf(A, k_l=6, knots_l=c(0.1, 0.2, 0.4, 0.8, 1.2, 1.8)) -
                irf(1),
            simulation$impulses,
            simulation$responses,
            window=c(0, 2),
            quiet=TRUE
        ),
        silent=TRUE
    ), 'try-error')
)

default_knots <- prepare_cdrgam(
    response ~ irf(A) - irf(1),
    simulation$impulses,
    simulation$responses,
    window=c(0, 2),
    knots_l=focused_knots,
    k_l=6,
    history='ragged',
    quiet=TRUE
)
stopifnot(
    identical(default_knots$specification[[1L]]$knots_l, focused_knots),
    identical(default_knots$configuration$knots_l, focused_knots)
)

override <- prepare_cdrgam(
    response ~ irf(A, window=c(0, 1), k_l=5) - irf(1),
    simulation$impulses,
    simulation$responses,
    window=c(-0.25, 2),
    history='ragged',
    quiet=TRUE
)
stopifnot(identical(override$specification[[1L]]$window, c(0, 1)))

defaults <- prepare_cdrgam(
    response ~ irf(A, B) +
        irf(A, k_l=5, k_t=NULL, k_p=NULL, bs_l='cr', bs_p='tp') -
        irf(1),
    simulation$impulses,
    simulation$responses,
    window=c(0, 2),
    k_l=7,
    k_t=4,
    k_p=3,
    bs_l='tp',
    bs_t='ps',
    bs_p='cs',
    history='ragged',
    quiet=TRUE
)
stopifnot(
    identical(defaults$specification[[1L]]$k_l, 7L),
    identical(defaults$specification[[1L]]$k_t, 4L),
    identical(defaults$specification[[1L]]$k_p, list(3L, 3L)),
    identical(defaults$specification[[1L]]$bs_l, 'tp'),
    identical(defaults$specification[[1L]]$bs_t, 'ps'),
    identical(defaults$specification[[1L]]$bs_p, list('cs', 'cs')),
    identical(defaults$specification[[2L]]$k_l, 5L),
    is.null(defaults$specification[[2L]]$k_t),
    identical(defaults$specification[[2L]]$k_p, list(NULL)),
    identical(defaults$specification[[2L]]$bs_l, 'cr'),
    identical(defaults$specification[[2L]]$bs_p, list('tp')),
    identical(defaults$configuration$k_l, 7),
    identical(defaults$configuration$k_t, 4),
    identical(defaults$configuration$k_p, 3)
)
implicit_rate_defaults <- prepare_cdrgam(
    response ~ irf(A),
    simulation$impulses,
    simulation$responses,
    window=c(0, 2),
    k_l=7,
    k_t=4,
    history='ragged',
    quiet=TRUE
)
stopifnot(
    isTRUE(implicit_rate_defaults$specification[[1L]]$constant),
    is.null(implicit_rate_defaults$specification[[1L]]$k_t),
    identical(implicit_rate_defaults$specification[[1L]]$k_l, 7L),
    identical(implicit_rate_defaults$specification[[2L]]$k_t, 4L)
)

# Centering removes the stationary lag-only subspace from a nonstationary
# tensor, so the two terms remain jointly identifiable.
hierarchy <- prepare_cdrgam(
    response ~ irf(1, k_l=5) + irf(1, k_l=5, k_t=4),
    simulation$impulses,
    simulation$responses,
    window=c(0, 2),
    history='ragged',
    quiet=TRUE
)
hierarchy_matrix <- do.call(cbind, lapply(hierarchy$terms, `[[`, 'X'))
stopifnot(
    qr(hierarchy_matrix)$rank == ncol(hierarchy_matrix),
    'deviation-centering:time' %in% hierarchy$terms[[2L]]$constraints
)
