library(cdrgam)

truth <- list(
    linear_decay=function(lag) 0.7 * exp(-1.8 * lag),
    nonlinear_peak=function(lag, value) {
        0.65 * (value^2 - 1) * exp(-((lag - 0.65) / 0.32)^2)
    }
)
simulation <- simulate_cdr(
    truth,
    n_impulses=700,
    n_responses=850,
    duration=70,
    window=2,
    intercept=50,
    noise_sd=0.7,
    seed=5505
)
stopifnot(identical(
    unname(simulation$nonlinear),
    c(FALSE, TRUE)
))

formula <- response ~
    irf(linear_decay, window=c(0, 2), k=10) +
    irf(
        nonlinear_peak,
        window=c(0, 2),
        nonlinear=TRUE,
        k=c(10, 7)
    ) - irf(1)
design <- prepare_cdrgam(
    formula,
    simulation$impulses,
    simulation$responses,
    history='ragged',
    chunk_size=250,
    quiet=TRUE
)
stopifnot(identical(design$terms[[1]]$type, 'linear'))
stopifnot(identical(design$terms[[2]]$type, 'nonlinear'))
stopifnot(ncol(design$terms[[2]]$X) == 10 * (7 - 1))

native <- fit_cdrgam(design, backend='mgcv', engine='gam', method='REML')
block <- fit_cdrgam(design, backend='block', method='REML')
sparse <- fit_cdrgam(
    design,
    backend='sparse',
    method='REML',
    sparse_control=list(hessian='optimhess')
)
sparse_profiled_hessian <- fit_cdrgam(
    design,
    backend='sparse',
    method='REML',
    sparse_control=list(hessian='profiled')
)
hessian_scale <- max(1, max(abs(sparse$outer.info$hess)))
stopifnot(
    max(abs(
        sparse$outer.info$hess - sparse_profiled_hessian$outer.info$hess
    )) / hessian_scale < 3e-3,
    max(abs(coef(sparse) - coef(sparse_profiled_hessian))) < 1e-8
)
lag <- seq(0, 2, length.out=81)
values <- c(-1.5, -0.75, 0, 0.75, 1.5)
native_surface <- estimate_irf(
    native,
    term='nonlinear_peak',
    lag=lag,
    predictor=values
)
block_surface <- estimate_irf(
    block,
    term='nonlinear_peak',
    lag=lag,
    predictor=values
)
sparse_surface <- estimate_irf(
    sparse,
    term='nonlinear_peak',
    lag=lag,
    predictor=values
)
target <- truth$nonlinear_peak(
    native_surface$lag,
    native_surface$predictor
)
rmse <- sqrt(mean((native_surface$estimate - target)^2))
correlation <- cor(native_surface$estimate, target)
stopifnot(rmse < 0.11)
stopifnot(correlation > 0.97)
stopifnot(max(abs(
    native_surface$estimate - block_surface$estimate
)) < 3e-4)
stopifnot(max(abs(
    native_surface$estimate - sparse_surface$estimate
)) < 3e-4)

# Surface extraction uses the full Cartesian grid and provides uncertainty.
stopifnot(nrow(native_surface) == length(lag) * length(values))
stopifnot(all(is.finite(native_surface$se)))

# Render predictor-conditioned IRF slices as a visual smoke test.
plot_path <- tempfile(fileext='.png')
grDevices::png(plot_path, width=900, height=650)
colors <- grDevices::hcl.colors(length(values), 'Blue-Red 3')
graphics::plot(
    range(lag),
    range(target, native_surface$estimate),
    type='n',
    xlab='Lag',
    ylab='IRF',
    main='Nonlinear IRF recovery by predictor value'
)
for (i in seq_along(values)) {
    rows <- native_surface$predictor == values[[i]]
    graphics::lines(
        lag,
        truth$nonlinear_peak(lag, values[[i]]),
        col=colors[[i]],
        lwd=3
    )
    graphics::lines(
        lag,
        native_surface$estimate[rows],
        col=colors[[i]],
        lwd=2,
        lty=2
    )
}
grDevices::dev.off()
stopifnot(file.exists(plot_path), file.info(plot_path)$size > 0)
unlink(plot_path)
