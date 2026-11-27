# Diagnoses the null calibration of the random-effect test of grouped IRF
# deviations by fitting the same simulations with native mgcv and with
# cdrgam's own backends. Run from the repository root with:
#   Rscript validation/deviation_test_diagnostic.R [replicates] [first_seed]
# Results go to validation/output/deviation_test_diagnostic.csv.

library(cdrgam)

source('validation/effect_calibration.R')

.deviation_design <- function(seed) {
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
    prepare_cdrgam(
        response ~ s(subject, bs='re') + irf(x, window=c(0, 2), k=6) +
            irf(x, window=c(0, 2), k=6, group=subject),
        do.call(rbind, impulses), responses, series='subject', history='auto',
        quiet=TRUE
    )
}

.deviation_row <- function(fit, backend, seed) {
    table <- summary(fit)$s.table
    status <- attr(table, 'test.status')
    row <- grep('|', rownames(table), fixed=TRUE)
    removed <- !length(row)
    data.frame(
        seed=seed, backend=backend, removed=removed,
        edf=if (removed) NA else table[row, 'edf'],
        ref_df=if (removed) NA else table[row, 'Ref.df'],
        statistic=if (removed) NA else table[row, 3L],
        p_value=if (removed) 1 else table[row, 'p-value'],
        status=if (removed) 'removed' else if (is.null(status)) 'mgcv' else
            unname(status[[row]]),
        subject_p=table[grep('subject', rownames(table))[[1L]], 'p-value'],
        log_sp=if (removed) NA else paste(round(log(fit$sp), 3), collapse=' ')
    )
}

run_deviation_diagnostic <- function(
        replicates=200L, first_seed=1L, output_dir='validation/output'
) {
    dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
    path <- file.path(output_dir, 'deviation_test_diagnostic.csv')
    for (seed in first_seed + seq_len(replicates) - 1L) {
        design <- .deviation_design(seed)
        for (backend in c('mgcv', 'sparse', 'block')) {
            row <- tryCatch(
                .deviation_row(.calibration_fit(design, backend), backend, seed),
                error=function(error) {
                    message('seed ', seed, ' ', backend, ': ', conditionMessage(error))
                    NULL
                }
            )
            if (is.null(row)) next
            utils::write.table(row, path, sep=',', row.names=FALSE,
                col.names=!file.exists(path), append=file.exists(path))
        }
        message('seed ', seed, ' done')
    }
    invisible(path)
}

if (sys.nframe() == 0L) {
    arguments <- commandArgs(trailingOnly=TRUE)
    run_deviation_diagnostic(
        replicates=if (length(arguments) >= 1L) as.integer(arguments[[1L]]) else 200L,
        first_seed=if (length(arguments) >= 2L) as.integer(arguments[[2L]]) else 1L
    )
}
