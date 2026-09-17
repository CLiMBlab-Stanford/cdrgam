# Synthetic fixed-IRF recovery suite.
# Run from the repository root with:
#   Rscript validation/synthetic_recovery.R
# An optional first argument selects the output directory.

library(cdrgam)

synthetic_cases <- list(
    two_distinct_irfs=list(
        seed=1101,
        noise_sd=0.75,
        k=12,
        irfs=list(
            excitation=function(lag) 1.25 * exp(-2.3 * lag),
            biphasic=function(lag) {
                0.95 * exp(-((lag - 0.35) / 0.18)^2) -
                    0.55 * exp(-((lag - 1.05) / 0.28)^2)
            }
        )
    ),
    three_simultaneous_irfs=list(
        seed=2202,
        noise_sd=1,
        k=14,
        irfs=list(
            early_positive=function(lag) 1.1 * exp(-((lag - 0.2) / 0.16)^2),
            late_negative=function(lag) -0.9 * exp(-((lag - 1.15) / 0.3)^2),
            damped_wave=function(lag) {
                0.8 * exp(-1.25 * lag) * sin(2.4 * pi * lag)
            }
        )
    ),
    unequal_scales_noisy=list(
        seed=3303,
        noise_sd=1.5,
        k=12,
        irfs=list(
            broad=function(lag) 1.4 * exp(-((lag - 0.7) / 0.55)^2),
            subtle=function(lag) -0.4 * exp(-((lag - 0.45) / 0.2)^2)
        )
    ),
    nonlinear_tensor_irfs=list(
        seed=4404,
        noise_sd=0.9,
        k=c(10, 7),
        irfs=list(
            saturating_decay=function(lag, value) {
                1.05 * tanh(1.1 * value) * exp(-1.7 * lag)
            },
            quadratic_peak=function(lag, value) {
                0.55 * (value^2 - 1) *
                    exp(-((lag - 0.7) / 0.32)^2)
            }
        )
    )
)

run_synthetic_case <- function(name, case, output_dir) {
    window <- 2
    simulation <- simulate_cdr(
        case$irfs,
        n_impulses=900,
        n_responses=1100,
        duration=80,
        window=window,
        intercept=300,
        noise_sd=case$noise_sd,
        seed=case$seed
    )
    nonlinear <- simulation$nonlinear
    irf_terms <- vapply(names(case$irfs), function(term) {
        if (nonlinear[[term]]) {
            sprintf(
                'irf(%s, window=c(0, %s), nonlinear=TRUE, k=c(%s, %s))',
                term,
                window,
                case$k[[1L]],
                case$k[[2L]]
            )
        } else {
            sprintf(
                'irf(%s, window=c(0, %s), k=%s)',
                term,
                window,
                case$k[[1L]]
            )
        }
    }, character(1))
    formula <- stats::as.formula(paste(
        'response ~',
        paste(irf_terms, collapse=' + ')
    ))
    design <- prepare_cdrgam(
        formula,
        simulation$impulses,
        simulation$responses,
        history='auto',
        chunk_size=250,
        quiet=FALSE
    )
    fits <- list(
        mgcv=fit_cdrgam(design, backend='mgcv', engine='gam', method='REML'),
        block=fit_cdrgam(design, backend='block', method='REML')
    )
    lag <- seq(0, window, length.out=301)
    predictor_slices <- c(-1.5, -0.75, 0, 0.75, 1.5)
    estimates <- lapply(
        fits,
        estimate_irf,
        lag=lag,
        predictor=predictor_slices
    )
    metrics <- list()
    for (backend in names(estimates)) {
        estimate <- estimates[[backend]]
        for (term in names(case$irfs)) {
            rows <- estimate$term == term
            truth <- if (nonlinear[[term]]) {
                case$irfs[[term]](
                    estimate$lag[rows],
                    estimate$predictor[rows]
                )
            } else {
                case$irfs[[term]](estimate$lag[rows])
            }
            error <- estimate$estimate[rows] - truth
            metrics[[length(metrics) + 1L]] <- data.frame(
                case=name,
                backend=backend,
                term=term,
                rmse=sqrt(mean(error^2)),
                correlation=stats::cor(estimate$estimate[rows], truth),
                max_abs_error=max(abs(error)),
                coverage_95=mean(
                    truth >= estimate$estimate[rows] - 1.96 * estimate$se[rows] &
                    truth <= estimate$estimate[rows] + 1.96 * estimate$se[rows]
                )
            )
        }
    }
    metrics <- do.call(rbind, metrics)

    grDevices::png(
        file.path(output_dir, paste0(name, '.png')),
        width=1400,
        height=420 * length(case$irfs),
        res=140
    )
    graphics::par(mfrow=c(length(case$irfs), 1), mar=c(4, 4, 3, 1))
    mgcv_estimate <- estimates$mgcv
    block_estimate <- estimates$block
    marker_index <- seq.int(1L, length(lag), by=25L)
    for (term in names(case$irfs)) {
        mgcv_rows <- mgcv_estimate$term == term
        block_rows <- block_estimate$term == term
        if (nonlinear[[term]]) {
            estimate_term <- mgcv_estimate[mgcv_rows, , drop=FALSE]
            block_term <- block_estimate[block_rows, , drop=FALSE]
            truth_values <- case$irfs[[term]](
                estimate_term$lag,
                estimate_term$predictor
            )
            limits <- range(
                truth_values,
                estimate_term$estimate,
                block_term$estimate,
                finite=TRUE
            )
            graphics::plot(
                range(lag),
                limits,
                type='n',
                xlab='Lag',
                ylab='Impulse response',
                main=paste(name, '-', term, '(predictor-conditioned slices)')
            )
            colors <- grDevices::hcl.colors(
                length(predictor_slices),
                'Blue-Red 3'
            )
            for (slice in seq_along(predictor_slices)) {
                value <- predictor_slices[[slice]]
                rows <- estimate_term$predictor == value
                block_slice <- block_term$predictor == value
                graphics::lines(
                    lag,
                    case$irfs[[term]](lag, value),
                    col=colors[[slice]],
                    lwd=3
                )
                graphics::lines(
                    lag,
                    estimate_term$estimate[rows],
                    col=colors[[slice]],
                    lwd=2,
                    lty=2
                )
                graphics::lines(
                    lag,
                    block_term$estimate[block_slice],
                    col=colors[[slice]],
                    lwd=1,
                    lty=3
                )
                graphics::points(
                    lag[marker_index],
                    block_term$estimate[block_slice][marker_index],
                    col=colors[[slice]],
                    bg='white',
                    pch=21,
                    cex=0.7,
                    lwd=1.2
                )
            }
            graphics::abline(h=0, col='grey70', lty=3)
            graphics::legend(
                'topright',
                legend=paste('x =', predictor_slices),
                col=colors,
                lwd=3,
                bty='n',
                ncol=2
            )
            graphics::legend(
                'bottomright',
                legend=c('truth', 'mgcv', 'block'),
                col='grey20',
                lty=c(1, 2, 3),
                lwd=c(3, 2, 1),
                pch=c(NA, NA, 21),
                pt.bg='white',
                pt.cex=0.8,
                bty='n'
            )
            next
        }
        truth <- case$irfs[[term]](lag)
        limits <- range(
            truth,
            mgcv_estimate$estimate[mgcv_rows] -
                1.96 * mgcv_estimate$se[mgcv_rows],
            mgcv_estimate$estimate[mgcv_rows] +
                1.96 * mgcv_estimate$se[mgcv_rows],
            block_estimate$estimate[block_rows],
            finite=TRUE
        )
        graphics::plot(
            lag,
            truth,
            type='n',
            ylim=limits,
            xlab='Lag',
            ylab='Impulse response',
            main=paste(name, '-', term)
        )
        graphics::polygon(
            c(lag, rev(lag)),
            c(
                mgcv_estimate$estimate[mgcv_rows] -
                    1.96 * mgcv_estimate$se[mgcv_rows],
                rev(mgcv_estimate$estimate[mgcv_rows] +
                    1.96 * mgcv_estimate$se[mgcv_rows])
            ),
            col=grDevices::adjustcolor('#377eb8', alpha.f=0.16),
            border=NA
        )
        graphics::lines(lag, truth, lwd=3, col='black')
        graphics::lines(
            lag,
            mgcv_estimate$estimate[mgcv_rows],
            lwd=2,
            col='#377eb8'
        )
        graphics::lines(
            lag,
            block_estimate$estimate[block_rows],
            lwd=2,
            lty=2,
            col='#e41a1c'
        )
        graphics::points(
            lag[marker_index],
            block_estimate$estimate[block_rows][marker_index],
            col='#e41a1c',
            bg='white',
            pch=21,
            cex=0.7,
            lwd=1.2
        )
        graphics::abline(h=0, col='grey70', lty=3)
        graphics::legend(
            'topright',
            legend=c('truth', 'mgcv', 'block', 'mgcv 95% interval'),
            col=c('black', '#377eb8', '#e41a1c', '#377eb8'),
            lty=c(1, 1, 2, NA),
            lwd=c(3, 2, 2, NA),
            pch=c(NA, NA, 21, 15),
            pt.bg=c(NA, NA, 'white', '#377eb8'),
            pt.cex=c(NA, NA, 0.8, 1.5),
            bty='n'
        )
    }
    grDevices::dev.off()
    metrics
}

run_synthetic_suite <- function(output_dir='validation/output') {
    dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
    metrics <- lapply(names(synthetic_cases), function(name) {
        message('Running synthetic case: ', name)
        run_synthetic_case(name, synthetic_cases[[name]], output_dir)
    })
    metrics <- do.call(rbind, metrics)
    utils::write.csv(
        metrics,
        file.path(output_dir, 'recovery_metrics.csv'),
        row.names=FALSE
    )
    print(metrics, row.names=FALSE)
    invisible(metrics)
}

if (sys.nframe() == 0L) {
    arguments <- commandArgs(trailingOnly=TRUE)
    output_dir <- if (length(arguments)) arguments[[1L]] else 'validation/output'
    run_synthetic_suite(output_dir)
}
