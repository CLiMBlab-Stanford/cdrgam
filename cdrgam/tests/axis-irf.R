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

override <- prepare_cdrgam(
    response ~ irf(A, window=c(0, 1), k_l=5) - irf(1),
    simulation$impulses,
    simulation$responses,
    window=c(-0.25, 2),
    history='ragged',
    quiet=TRUE
)
stopifnot(identical(override$specification[[1L]]$window, c(0, 1)))

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
