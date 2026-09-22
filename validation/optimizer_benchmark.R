# Compare sparse smoothing-parameter strategies on a deterministic crossed
# random-effects problem. Environment overrides keep this useful from a
# laptop to a Slurm node.

library(cdrgam)

subject_count <- as.integer(Sys.getenv('CDRGAM_OPT_SUBJECTS', '20'))
item_count <- as.integer(Sys.getenv('CDRGAM_OPT_ITEMS', '100'))
predictor_count <- as.integer(Sys.getenv('CDRGAM_OPT_PREDICTORS', '4'))
response_count <- as.integer(Sys.getenv('CDRGAM_OPT_RESPONSES', '1500'))
gradient <- Sys.getenv('CDRGAM_OPT_GRADIENT', 'finite')
probes <- as.integer(Sys.getenv('CDRGAM_OPT_PROBES', '12'))
gradient_cores <- as.integer(Sys.getenv('CDRGAM_OPT_GRADIENT_CORES', '1'))
hessian <- Sys.getenv('CDRGAM_OPT_HESSIAN', 'profiled')
outer_optimizer <- Sys.getenv('CDRGAM_OPT_OUTER', 'auto')
trace_method <- Sys.getenv('CDRGAM_OPT_TRACE_METHOD', 'auto')
simulation_seed <- as.integer(Sys.getenv('CDRGAM_OPT_SEED', '9201'))
restarts <- as.integer(Sys.getenv('CDRGAM_OPT_RESTARTS', '0'))
optimizer_maxit <- as.integer(Sys.getenv('CDRGAM_OPT_MAXIT', '100'))

truth <- stats::setNames(lapply(seq_len(predictor_count), function(i) {
    force(i)
    function(lag) (-1)^(i + 1L) * (0.3 + 0.1 * i) * exp(-(1 + i / 4) * lag)
}), paste0('x', seq_len(predictor_count)))
simulation <- simulate_cdr(
    truth,
    n_impulses=response_count,
    n_responses=response_count,
    duration=100,
    window=2,
    intercept=20,
    noise_sd=1,
    seed=simulation_seed
)
set.seed(simulation_seed + 1L)
simulation$responses$subject <- factor(sample(
    sprintf('s%03d', seq_len(subject_count)),
    response_count,
    replace=TRUE
))
simulation$responses$item <- factor(sample(
    sprintf('i%04d', seq_len(item_count)),
    response_count,
    replace=TRUE
))

# Match the Brown model: complete subject-specific IRFs, but only an ordinary
# random intercept for item.
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

elapsed <- system.time({
    fit <- cdrgam.fit(
        design,
        backend='sparse',
        method='REML',
        sparse_control=list(
            gradient=gradient,
            gradient_probes=probes,
            gradient_cores=gradient_cores,
            hessian=hessian,
            outer_optimizer=outer_optimizer,
            optimizer_maxit=optimizer_maxit,
            trace_method=trace_method,
            restarts=restarts,
            schur='always'
        )
    )
})[['elapsed']]

result <- data.frame(
    gradient_requested=fit$sparse$gradient_requested,
    gradient=fit$sparse$gradient,
    probes=if (fit$sparse$gradient %in% c(
        'stochastic', 'hybrid'
    )) probes else 0L,
    gradient_cores=gradient_cores,
    hessian=hessian,
    outer_optimizer_requested=fit$sparse$outer_optimizer_requested,
    outer_optimizer=fit$sparse$outer_optimizer,
    trace_method=fit$sparse$trace_method,
    selection_reason=fit$sparse$optimizer_selection$reason,
    probe_objective_seconds=
        fit$sparse$optimizer_selection$objective_seconds,
    predicted_exact_seconds=
        fit$sparse$optimizer_selection$predicted_exact_seconds,
    predicted_finite_seconds=
        fit$sparse$optimizer_selection$predicted_finite_seconds,
    median_exact_seconds=
        fit$sparse$optimizer_selection$median_exact_seconds,
    exact_core_mib=
        fit$sparse$optimizer_selection$exact_core_bytes / 1024^2,
    simulation_seed=simulation_seed,
    restarts=restarts,
    observations=nrow(simulation$responses),
    coefficients=length(coef(fit)),
    smoothing_parameters=length(fit$sp),
    reml=fit$reml,
    elapsed_seconds=unname(elapsed),
    objective_evaluations=
        fit$sparse$convergence$total_objective_evaluations,
    optimizer_convergence=fit$optimizer$convergence,
    optimizer_iterations=if (is.null(fit$optimizer$iterations)) {
        NA_integer_
    } else {
        fit$optimizer$iterations
    },
    optimizer_rejected_steps=if (is.null(fit$optimizer$rejected_steps)) {
        NA_integer_
    } else {
        fit$optimizer$rejected_steps
    },
    optimizer_gradient_max=if (is.null(fit$optimizer$gradient)) {
        NA_real_
    } else {
        max(abs(fit$optimizer$gradient))
    },
    optimizer_message=if (is.null(fit$optimizer$message)) {
        ''
    } else {
        fit$optimizer$message
    },
    hessian_evaluations=fit$sparse$hessian_evaluations,
    factor_updates=fit$sparse$numeric_updates,
    factor_nonzeros=fit$sparse$factor_nonzeros,
    converged=fit$converged
)
print(result, row.names=FALSE)
