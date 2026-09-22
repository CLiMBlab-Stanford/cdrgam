# Benchmark the profiled Gaussian REML Hessian reconstruction against direct
# finite differencing of the unprofiled criterion. Both fits use the same
# deterministic data, model, starts, optimizer, and BLAS configuration.

library(cdrgam)

subject_count <- as.integer(Sys.getenv('CDRGAM_BENCH_SUBJECTS', '20'))
item_count <- as.integer(Sys.getenv('CDRGAM_BENCH_ITEMS', '100'))
predictor_count <- as.integer(Sys.getenv('CDRGAM_BENCH_PREDICTORS', '4'))
response_count <- as.integer(Sys.getenv('CDRGAM_BENCH_RESPONSES', '1500'))
replicates <- as.integer(Sys.getenv('CDRGAM_BENCH_REPLICATES', '2'))
gradient <- Sys.getenv('CDRGAM_BENCH_GRADIENT', 'exact')
hessian_step <- as.numeric(Sys.getenv('CDRGAM_BENCH_HESSIAN_STEP', '1e-2'))
hessian_tolerance <- as.numeric(Sys.getenv(
    'CDRGAM_BENCH_HESSIAN_TOLERANCE',
    '3e-3'
))
output_path <- Sys.getenv(
    'CDRGAM_BENCH_OUTPUT',
    file.path('validation', 'output', 'hessian_benchmark.csv')
)

stopifnot(
    subject_count >= 2L,
    item_count >= 2L,
    predictor_count >= 1L,
    response_count >= 100L,
    replicates >= 1L,
    gradient %in% c('finite', 'exact'),
    is.finite(hessian_step),
    hessian_step > 0,
    is.finite(hessian_tolerance),
    hessian_tolerance > 0
)

truth <- stats::setNames(lapply(seq_len(predictor_count), function(i) {
    force(i)
    function(lag) (-1)^(i + 1L) * (0.3 + 0.1 * i) *
        exp(-(1 + i / 4) * lag)
}), paste0('x', seq_len(predictor_count)))
simulation <- simulate_cdr(
    truth,
    n_impulses=response_count,
    n_responses=response_count,
    duration=100,
    window=2,
    intercept=20,
    noise_sd=1,
    seed=9201
)
set.seed(9202)
simulation$responses$subject <- factor(sample(
    sprintf('s%03d', seq_len(subject_count)),
    response_count,
    replace=TRUE
))
simulation$responses$item <- factor(sample(
    sprintf('i%05d', seq_len(item_count)),
    response_count,
    replace=TRUE
))

# Full IRFs by subject, random intercept only by item.
groups <- c(NA_character_, 'subject')
irf_terms <- unlist(lapply(c('1', names(truth)), function(predictor) {
    vapply(groups, function(group) {
        paste0(
            'irf(', predictor, ', window=c(0,2), k=6',
            if (is.na(group)) '' else paste0(', group=', group),
            ')'
        )
    }, character(1))
}), use.names=FALSE)
formula <- stats::as.formula(paste(
    'response ~ s(subject, bs="re") + s(item, bs="re") +',
    paste(irf_terms, collapse=' + '),
    '- irf(1)'
))
design <- prepare_cdrgam(
    formula,
    simulation$impulses,
    simulation$responses,
    history='auto',
    history_length=8,
    quiet=TRUE
)
stopifnot(
    !any(vapply(
        design$terms,
        function(term) identical(term$group, 'item'),
        logical(1)
    )),
    any(vapply(
        design$terms,
        function(term) identical(term$group, 'subject'),
        logical(1)
    ))
)

run_fit <- function(method, replicate, order) {
    events <- list()
    elapsed <- system.time({
        fit <- cdrgam.fit(
            design,
            backend='sparse',
            method='REML',
            solver_trace=function(event) {
                events[[length(events) + 1L]] <<- event
            },
            sparse_control=list(
                gradient=gradient,
                gradient_cores=1L,
                hessian=method,
                hessian_step=hessian_step,
                schur='always'
            )
        )
    })[['elapsed']]
    phases <- vapply(events, `[[`, character(1), 'phase')
    event_names <- vapply(events, `[[`, character(1), 'event')
    hessian_start <- events[[which(
        phases == 'outer Hessian' & event_names == 'phase'
    )[[1L]]]]$elapsed_seconds
    finalization_start <- events[[which(
        phases == 'finalization' & event_names == 'phase'
    )[[1L]]]]$elapsed_seconds
    status <- readLines('/proc/self/status', warn=FALSE)
    peak_rss <- status[grepl('^VmHWM:', status)]
    peak_rss_mib <- if (length(peak_rss)) {
        as.numeric(sub('^VmHWM:\\s*([0-9]+)\\s+kB.*$', '\\1', peak_rss)) /
            1024
    } else {
        NA_real_
    }
    row <- data.frame(
        method=method,
        gradient=gradient,
        hessian_step=if (method == 'profiled') hessian_step else NA_real_,
        replicate=replicate,
        order=order,
        observations=nrow(simulation$responses),
        coefficients=length(coef(fit)),
        smoothing_parameters=length(fit$sp),
        optimization_seconds=hessian_start,
        hessian_seconds=finalization_start - hessian_start,
        total_seconds=unname(elapsed),
        reml=fit$reml,
        objective_evaluations=
            fit$sparse$convergence$total_objective_evaluations,
        hessian_evaluations=fit$sparse$hessian_evaluations,
        factor_updates=fit$sparse$numeric_updates,
        peak_rss_mib=peak_rss_mib,
        converged=fit$converged
    )
    list(fit=fit, row=row)
}

rows <- list()
comparisons <- list()
row_index <- 0L
for (replicate in seq_len(replicates)) {
    methods <- if (replicate %% 2L) {
        c('optimhess', 'profiled')
    } else {
        c('profiled', 'optimhess')
    }
    fits <- list()
    for (order in seq_along(methods)) {
        method <- methods[[order]]
        message('Benchmark replicate ', replicate, ', method ', method)
        result <- run_fit(method, replicate, order)
        row_index <- row_index + 1L
        rows[[row_index]] <- result$row
        fits[[method]] <- result$fit
    }
    hessian_scale <- max(1, max(abs(fits$optimhess$outer.info$hess)))
    comparisons[[replicate]] <- data.frame(
        replicate=replicate,
        reml_difference=abs(fits$profiled$reml - fits$optimhess$reml),
        coefficient_max_difference=max(abs(
            coef(fits$profiled) - coef(fits$optimhess)
        )),
        sp_max_relative_difference=max(abs(
            fits$profiled$sp - fits$optimhess$sp
        ) / pmax(1, abs(fits$optimhess$sp))),
        hessian_max_relative_difference=max(abs(
            fits$profiled$outer.info$hess - fits$optimhess$outer.info$hess
        )) / hessian_scale
    )
}

results <- do.call(rbind, rows)
comparison <- do.call(rbind, comparisons)
dir.create(dirname(output_path), recursive=TRUE, showWarnings=FALSE)
utils::write.csv(results, output_path, row.names=FALSE)
utils::write.csv(
    comparison,
    sub('\\.csv$', '_equivalence.csv', output_path),
    row.names=FALSE
)
print(results, row.names=FALSE)
print(comparison, row.names=FALSE)

stopifnot(
    all(results$converged),
    max(comparison$reml_difference) < 1e-8,
    max(comparison$coefficient_max_difference) < 1e-8,
    max(comparison$sp_max_relative_difference) < 1e-8,
    max(comparison$hessian_max_relative_difference) < hessian_tolerance
)
