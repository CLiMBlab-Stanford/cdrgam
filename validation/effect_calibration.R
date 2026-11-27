# Calibration of whole-term tests and effect inference on simulated data.
#
# Run from the repository root with:
#   Rscript validation/effect_calibration.R [replicates] [first_seed] [backend]
#   Rscript validation/effect_calibration.R summarize
#
# Each replicate simulates fresh data, fits it, and records, per scenario:
#   null_term    p-values for a term whose true IRF is zero (type-I error of
#                the Wood and grid Wald tests) and for a nonzero term (power);
#   coverage     whether intervals from effect_inference() cover the truth
#                for pointwise values, a simultaneous band, an area, a
#                contrast between lags, and the peak lag, and the p-value of
#                a test of the true area;
#   null_random  p-values of random-effect tests for a subject intercept and
#                grouped IRF deviations that are zero.
# Results are appended to validation/output/effect_calibration_<backend>.csv.
# `summarize` reports rejection rates and coverage for every results file.

library(cdrgam)

lags <- data.frame(lag=seq(0, 2, length.out=21))
zero <- function(lag) 0 * lag
decay <- function(lag) 0.9 * exp(-1.5 * lag)
biphasic <- function(lag) {
    0.8 * exp(-((lag - 0.4) / 0.2)^2) - 0.45 * exp(-((lag - 1.1) / 0.3)^2)
}

.calibration_fit <- function(design, backend) {
    if (identical(backend, 'mgcv')) {
        return(cdrgam.fit(design, backend='mgcv', engine='gam', method='REML'))
    }
    cdrgam.fit(design, backend=backend, method='REML')
}

.calibration_null_term <- function(seed, backend) {
    simulation <- simulate_cdr(
        list(x=zero, z=decay), n_impulses=300, n_responses=320, duration=60,
        window=2, noise_sd=1, seed=seed
    )
    design <- prepare_cdrgam(
        response ~ irf(x, window=c(0, 2), k=8) + irf(z, window=c(0, 2), k=8),
        simulation$impulses, simulation$responses, quiet=TRUE
    )
    fit <- .calibration_fit(design, backend)
    rows <- list()
    for (term in c('x', 'z')) {
        selector <- list(predictors=term)
        wood <- test_effect(fit, terms=selector)
        points <- lags
        points[[term]] <- 1
        wald <- test_effect(fit, terms=selector, points=points, method='wald')
        rows[[length(rows) + 1L]] <- data.frame(
            scenario='null_term', seed=seed, quantity=paste0(term, '_wood'),
            truth_null=identical(term, 'x'), p_value=wood$p_value, covered=NA
        )
        rows[[length(rows) + 1L]] <- data.frame(
            scenario='null_term', seed=seed, quantity=paste0(term, '_wald'),
            truth_null=identical(term, 'x'), p_value=wald$p_value, covered=NA
        )
    }
    do.call(rbind, rows)
}

.calibration_coverage <- function(seed, backend) {
    simulation <- simulate_cdr(
        list(x=biphasic, z=decay), n_impulses=300, n_responses=320,
        duration=60, window=2, noise_sd=1, seed=seed
    )
    design <- prepare_cdrgam(
        response ~ irf(x, window=c(0, 2), k=8) + irf(z, window=c(0, 2), k=8),
        simulation$impulses, simulation$responses, quiet=TRUE
    )
    fit <- .calibration_fit(design, backend)
    truth <- biphasic(lags$lag)
    # The truth is the response to a unit impulse.
    effect <- list(x=list(terms=list(predictors='x'), points=transform(lags, x=1)))
    pointwise <- effect_inference(fit, effect, draws=4000L, seed=seed)
    early <- which.min(abs(lags$lag - 0.4))
    late <- which.min(abs(lags$lag - 1.1))
    contrast_weights <- numeric(nrow(lags))
    contrast_weights[c(early, late)] <- c(1, -1)
    true_area <- mean(truth)
    linear <- effect_inference(
        fit, effect,
        estimand=list(area=list(x=1 / nrow(lags)), contrast=list(x=contrast_weights)),
        null=c(true_area, truth[[early]] - truth[[late]])
    )
    peak <- effect_inference(
        fit, effect, estimand=function(values, grids) {
            c(peak_lag=grids$x$lag[[which.max(values$x)]])
        }, draws=500L, seed=seed
    )
    true_peak <- lags$lag[[which.max(truth)]]
    data.frame(
        scenario='coverage', seed=seed,
        quantity=c('pointwise', 'simultaneous', 'area', 'contrast', 'peak_lag',
            'area_test', 'contrast_test'),
        truth_null=c(NA, NA, NA, NA, NA, TRUE, TRUE),
        p_value=c(NA, NA, NA, NA, NA, linear$p_value),
        covered=as.numeric(c(
            mean(truth >= pointwise$lower & truth <= pointwise$upper),
            all(truth >= pointwise$simultaneous_lower &
                truth <= pointwise$simultaneous_upper),
            linear$lower[[1L]] <= true_area && true_area <= linear$upper[[1L]],
            linear$lower[[2L]] <= truth[[early]] - truth[[late]] &&
                truth[[early]] - truth[[late]] <= linear$upper[[2L]],
            peak$lower <= true_peak && true_peak <= peak$upper,
            NA, NA
        ))
    )
}

.calibration_null_random <- function(seed, backend) {
    impulses <- list()
    responses <- list()
    for (i in 1:6) {
        simulation <- simulate_cdr(
            list(x=decay), n_impulses=120, n_responses=130, duration=30,
            window=2, noise_sd=1, seed=seed * 10L + i
        )
        simulation$impulses$subject <- simulation$responses$subject <- paste0('s', i)
        impulses[[i]] <- simulation$impulses
        responses[[i]] <- simulation$responses
    }
    responses <- do.call(rbind, responses)
    responses$subject <- factor(responses$subject)
    design <- prepare_cdrgam(
        response ~ s(subject, bs='re') + irf(x, window=c(0, 2), k=6) +
            irf(x, window=c(0, 2), k=6, group=subject),
        do.call(rbind, impulses), responses, series='subject', history='auto',
        quiet=TRUE
    )
    fit <- .calibration_fit(design, backend)
    table <- summary(fit)$s.table
    # Boundary reduction removes a deviation whose variance is estimated at
    # zero. That is a non-rejection, recorded as p = 1 and counted under
    # irf_deviation_removed.
    removed <- !any(effect_catalog(fit)$grouped)
    deviation_p <- if (removed) 1 else test_effect(
        fit, terms=list(predictors='x', grouped=TRUE), grouping='deviation'
    )$p_value
    data.frame(
        scenario='null_random', seed=seed,
        quantity=c('subject_intercept', 'irf_deviation', 'irf_deviation_removed'),
        truth_null=c(TRUE, TRUE, NA),
        p_value=c(table[grep('subject', rownames(table))[[1L]], 'p-value'],
            deviation_p, NA),
        covered=c(NA, NA, as.numeric(removed))
    )
}

run_effect_calibration <- function(
        replicates=100L, first_seed=1L, backend='sparse',
        output_dir='validation/output'
) {
    dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
    path <- file.path(output_dir, paste0('effect_calibration_', backend, '.csv'))
    for (seed in first_seed + seq_len(replicates) - 1L) {
        for (scenario in list(.calibration_null_term, .calibration_coverage,
                .calibration_null_random)) {
            rows <- tryCatch(scenario(seed, backend), error=function(error) {
                message('seed ', seed, ': ', conditionMessage(error))
                NULL
            })
            if (is.null(rows)) next
            rows$backend <- backend
            utils::write.table(
                rows, path, sep=',', row.names=FALSE,
                col.names=!file.exists(path), append=file.exists(path)
            )
        }
        message('seed ', seed, ' done')
    }
    invisible(path)
}

summarize_effect_calibration <- function(output_dir='validation/output') {
    files <- list.files(output_dir, '^effect_calibration_.*[.]csv$', full.names=TRUE)
    results <- do.call(rbind, lapply(files, utils::read.csv))
    groups <- split(results, list(results$backend, results$scenario,
        results$quantity), drop=TRUE)
    summary <- do.call(rbind, lapply(groups, function(rows) {
        p <- rows$p_value
        data.frame(
            backend=rows$backend[[1L]], scenario=rows$scenario[[1L]],
            quantity=rows$quantity[[1L]], replicates=nrow(rows),
            truth_null=rows$truth_null[[1L]],
            reject_05=if (all(is.na(p))) NA else mean(p < 0.05, na.rm=TRUE),
            reject_01=if (all(is.na(p))) NA else mean(p < 0.01, na.rm=TRUE),
            coverage=if (all(is.na(rows$covered))) NA else
                mean(as.numeric(rows$covered), na.rm=TRUE)
        )
    }))
    rownames(summary) <- NULL
    summary
}

if (sys.nframe() == 0L) {
    arguments <- commandArgs(trailingOnly=TRUE)
    if (length(arguments) && identical(arguments[[1L]], 'summarize')) {
        print(summarize_effect_calibration(), digits=3)
    } else {
        run_effect_calibration(
            replicates=if (length(arguments) >= 1L) as.integer(arguments[[1L]]) else 100L,
            first_seed=if (length(arguments) >= 2L) as.integer(arguments[[2L]]) else 1L,
            backend=if (length(arguments) >= 3L) arguments[[3L]] else 'sparse'
        )
    }
}
