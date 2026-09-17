# Synthetic random-intercept and random-IRF recovery validation.
# Run from the repository root with:
#   Rscript validation/mixed_recovery.R

library(cdrgam)

run_mixed_recovery <- function(output_dir='validation/output') {
    dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
    group_count <- 12
    groups <- paste0('subject_', seq_len(group_count))
    random_amplitude <- seq(-0.65, 0.65, length.out=group_count)
    random_intercept <- 1.8 * sin(seq(0, 2 * pi, length.out=group_count + 1))[
        seq_len(group_count)
    ]
    population_irf <- function(lag) 0.95 * exp(-1.6 * lag)
    deviation_shape <- function(lag) {
        exp(-((lag - 0.7) / 0.3)^2) - 0.22 * exp(-1.4 * lag)
    }

    impulse_streams <- vector('list', group_count)
    response_streams <- vector('list', group_count)
    for (i in seq_len(group_count)) {
        amplitude <- random_amplitude[[i]]
        total_irf <- function(lag) {
            population_irf(lag) + amplitude * deviation_shape(lag)
        }
        simulation <- simulate_cdr(
            list(x=total_irf),
            n_impulses=240,
            n_responses=280,
            duration=40,
            window=2,
            intercept=100 + random_intercept[[i]],
            noise_sd=0.55,
            seed=7700 + i
        )
        simulation$impulses$subject <- groups[[i]]
        simulation$responses$subject <- groups[[i]]
        impulse_streams[[i]] <- simulation$impulses
        response_streams[[i]] <- simulation$responses
    }
    impulses <- do.call(rbind, impulse_streams)
    responses <- do.call(rbind, response_streams)
    responses$subject <- factor(responses$subject, levels=groups)

    formula <- response ~
        s(subject, bs='re') +
        irf(x, window=c(0, 2), k=10) +
        irf(x, window=c(0, 2), k=10, group=subject)
    design <- prepare_cdrgam(
        formula,
        impulses,
        responses,
        series='subject',
        history='auto',
        chunk_size=300
    )
    fits <- list(
        mgcv=fit_cdrgam(design, backend='mgcv', engine='gam', method='REML'),
        block=fit_cdrgam(design, backend='block', method='REML')
    )
    lag <- seq(0, 2, length.out=301)
    selected_groups <- groups[c(1, 3, 5, 8, 10, 12)]
    estimates <- lapply(fits, function(fit) {
        list(
            population=estimate_irf(fit, term='x', lag=lag),
            deviation=estimate_irf(
                fit,
                term='x|subject',
                lag=lag,
                level=selected_groups
            )
        )
    })

    random_smooth <- which(vapply(
        fits$mgcv$smooth,
        function(smooth) identical(smooth$label, 's(subject)'),
        logical(1)
    ))
    random_smooth <- fits$mgcv$smooth[[random_smooth]]
    random_index <- random_smooth$first.para:random_smooth$last.para
    intercept_estimates <- list(
        mgcv=coef(fits$mgcv)[random_index],
        block=coef(fits$block)[random_index]
    )

    metrics <- list()
    for (backend in names(fits)) {
        recovered <- numeric()
        target <- numeric()
        for (group in selected_groups) {
            group_number <- match(group, groups)
            rows <- estimates[[backend]]$deviation$group == group
            recovered <- c(
                recovered,
                estimates[[backend]]$population$estimate +
                    estimates[[backend]]$deviation$estimate[rows]
            )
            target <- c(
                target,
                population_irf(lag) +
                    random_amplitude[[group_number]] * deviation_shape(lag)
            )
        }
        metrics[[backend]] <- data.frame(
            backend=backend,
            group_irf_rmse=sqrt(mean((recovered - target)^2)),
            group_irf_correlation=stats::cor(recovered, target),
            random_intercept_correlation=stats::cor(
                intercept_estimates[[backend]],
                random_intercept
            ),
            fitted_difference_from_mgcv=max(abs(
                fitted(fits[[backend]]) - fitted(fits$mgcv)
            ))
        )
    }
    metrics <- do.call(rbind, metrics)
    utils::write.csv(
        metrics,
        file.path(output_dir, 'mixed_recovery_metrics.csv'),
        row.names=FALSE
    )

    grDevices::png(
        file.path(output_dir, 'mixed_random_irfs.png'),
        width=1400,
        height=1050,
        res=140
    )
    graphics::par(mfrow=c(2, 1), mar=c(4, 4, 3, 1))
    colors <- grDevices::hcl.colors(length(selected_groups), 'Dark 3')
    all_truth <- unlist(lapply(selected_groups, function(group) {
        i <- match(group, groups)
        population_irf(lag) + random_amplitude[[i]] * deviation_shape(lag)
    }))
    graphics::plot(
        range(lag),
        range(all_truth),
        type='n',
        xlab='Lag',
        ylab='Total impulse response',
        main='Population + subject-specific random IRFs'
    )
    marker_index <- seq.int(1L, length(lag), by=25L)
    for (i in seq_along(selected_groups)) {
        group <- selected_groups[[i]]
        group_number <- match(group, groups)
        mgcv_rows <- estimates$mgcv$deviation$group == group
        block_rows <- estimates$block$deviation$group == group
        mgcv_curve <- estimates$mgcv$population$estimate +
            estimates$mgcv$deviation$estimate[mgcv_rows]
        block_curve <- estimates$block$population$estimate +
            estimates$block$deviation$estimate[block_rows]
        graphics::lines(
            lag,
            population_irf(lag) +
                random_amplitude[[group_number]] * deviation_shape(lag),
            col=colors[[i]],
            lwd=3
        )
        graphics::lines(lag, mgcv_curve, col=colors[[i]], lwd=2, lty=2)
        graphics::lines(lag, block_curve, col=colors[[i]], lwd=1, lty=3)
        graphics::points(
            lag[marker_index],
            block_curve[marker_index],
            col=colors[[i]],
            bg='white',
            pch=21,
            cex=0.7
        )
    }
    graphics::legend(
        'topright',
        legend=selected_groups,
        col=colors,
        lwd=3,
        ncol=2,
        bty='n'
    )
    graphics::legend(
        'bottomright',
        legend=c('truth', 'mgcv', 'block'),
        col='grey20',
        lty=c(1, 2, 3),
        pch=c(NA, NA, 21),
        pt.bg='white',
        lwd=c(3, 2, 1),
        bty='n'
    )

    limits <- range(random_intercept, unlist(intercept_estimates))
    graphics::plot(
        random_intercept,
        intercept_estimates$mgcv,
        xlim=limits,
        ylim=limits,
        pch=19,
        col='#377eb8',
        xlab='Ground-truth random intercept',
        ylab='Recovered random intercept',
        main='Subject random-intercept recovery'
    )
    graphics::points(
        random_intercept,
        intercept_estimates$block,
        pch=1,
        cex=1.3,
        col='#e41a1c'
    )
    graphics::abline(a=0, b=1, lty=2)
    graphics::legend(
        'topleft',
        legend=c('mgcv', 'block'),
        col=c('#377eb8', '#e41a1c'),
        pch=c(19, 1),
        bty='n'
    )
    grDevices::dev.off()
    print(metrics, row.names=FALSE)
    invisible(metrics)
}

if (sys.nframe() == 0L) {
    arguments <- commandArgs(trailingOnly=TRUE)
    output_dir <- if (length(arguments)) arguments[[1L]] else 'validation/output'
    run_mixed_recovery(output_dir)
}
