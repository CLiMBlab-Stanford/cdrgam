library(cdrgam)

set.seed(5291)
impulses <- data.frame(
    time=sort(stats::runif(100, 0, 30)),
    location_signal=stats::rnorm(100),
    scale_signal=stats::rnorm(100)
)
responses <- data.frame(
    time=sort(stats::runif(140, 0.5, 30)),
    response=stats::rnorm(140)
)
impulses$subject <- factor(ifelse(impulses$time < 15, 's1', 's2'))
responses$subject <- factor(ifelse(responses$time < 15, 's1', 's2'))
design <- prepare_cdrgam(
    list(
        location=response ~ s(subject, bs='re') +
            irf(location_signal, window=c(0, 1.5), k_l=5) - irf(1),
        scale=~ s(subject, bs='re') +
            irf(scale_signal, window=c(0, 1.5), k_l=5) - irf(1)
    ),
    impulses,
    responses,
    series='subject',
    history='ragged',
    chunk_size=31,
    quiet=TRUE
)
assemblies <- lapply(design$parameters, function(parameter_design) {
    cdrgam:::.cdrgam_distributional_sparse_setup(
        parameter_design,
        list(crossprod_chunk_size=31)
    )
})
family <- cdrgam_family('gaulss')

cached <- cdrgam:::.cdrgam_distributional_cache_designs(assemblies)
stopifnot(
    is.logical(cached$enabled),
    length(cached$enabled) == 1L,
    is.numeric(cached$projected_bytes),
    is.numeric(cached$bytes)
)
if (cached$enabled) {
    rows <- seq_len(min(25L, assemblies$location$observation_count))
    for (parameter in names(assemblies)) stopifnot(isTRUE(all.equal(
        cdrgam:::.cdrgam_sparse_design_chunk(
            assemblies[[parameter]], rows
        ),
        cdrgam:::.cdrgam_sparse_design_chunk(
            cached$assemblies[[parameter]], rows
        )
    )))
}

singular_system <- Matrix::Matrix(
    matrix(c(1, 1, 1, 1), 2L, 2L), sparse=TRUE
)
damped_factor <- cdrgam:::.cdrgam_gaulss_step_factor(
    singular_system, supernodal=FALSE
)
stopifnot(
    damped_factor$damping > 0,
    all(is.finite(cdrgam:::.cdr_factor_solve(
        damped_factor$factor, c(1, 0)
    )))
)

numerically_converged <- cdrgam:::.cdrgam_gaulss_numerically_converged
stopifnot(
    numerically_converged(-6.2743e-9, 1.0594e-8, 299647.5249),
    !numerically_converged(6.2743e-9, 1.0594e-8, 299647.5249),
    !numerically_converged(-1, 1.0594e-8, 299647.5249),
    !numerically_converged(-6.2743e-9, 1, 299647.5249)
)

valid_cache <- list(criterion=12, solution=list(converged=TRUE))
invalid_cache <- list(
    criterion=1e50,
    score=c(0, 0),
    solution=NULL,
    invalid_reason='test failure'
)
stopifnot(
    identical(
        cdrgam:::.cdrgam_gaulss_cache_candidate(
            c(1, 2), invalid_cache, c(1, 2), valid_cache
        ),
        valid_cache
    ),
    identical(
        cdrgam:::.cdrgam_gaulss_cache_candidate(
            c(1, 3), invalid_cache, c(1, 2), valid_cache
        ),
        invalid_cache
    )
)
log_sp <- log(unlist(lapply(
    assemblies,
    cdrgam:::.cdrgam_sparse_initial_sp
), use.names=FALSE))

criterion <- function(parameters, retain=FALSE) {
    sp <- exp(parameters)
    solution <- cdrgam:::.cdrgam_gaulss_sparse_fixed(
        assemblies,
        sp,
        cdrgam:::.cdrgam_gaulss_b(family)
    )
    counts <- vapply(
        assemblies,
        function(assembly) length(assembly$penalty_components),
        integer(1)
    )
    penalty_determinant <- 0
    offset <- 0L
    for (parameter in names(assemblies)) {
        count <- counts[[parameter]]
        indices <- offset + seq_len(count)
        penalty_determinant <- penalty_determinant +
            cdrgam:::.sparse_penalty_logdet(
                assemblies[[parameter]]$blocks,
                sp[indices]
            )$value
        offset <- offset + count
    }
    value <- 2 * solution$objective +
        cdrgam:::.cdr_factor_logdet(solution$factor) - penalty_determinant
    if (retain) list(value=value, solution=solution) else value
}

step <- 1e-4
points <- list(
    log_sp,
    log_sp + c(-3, 2, -1, 1),
    log_sp + c(2, -3, 1, -1)
)
differences <- vapply(points, function(point) {
    retained <- criterion(point, retain=TRUE)
    analytic <- cdrgam:::.cdrgam_gaulss_sparse_score(
        assemblies,
        retained$solution,
        exp(point),
        family
    )
    finite <- vapply(seq_along(point), function(i) {
        plus <- minus <- point
        plus[[i]] <- plus[[i]] + step
        minus[[i]] <- minus[[i]] - step
        (criterion(plus) - criterion(minus)) / (2 * step)
    }, numeric(1))
    max(abs(analytic - finite))
}, numeric(1))
stopifnot(max(differences) < 2e-4)

if (.Platform$OS.type != 'windows') {
    retained <- criterion(log_sp, retain=TRUE)
    reference <- cdrgam:::.cdrgam_gaulss_sparse_score_reference(
        assemblies, retained$solution, exp(log_sp), family
    )
    serial <- cdrgam:::.cdrgam_gaulss_sparse_score(
        assemblies, retained$solution, exp(log_sp), family
    )
    previous_memory_limit <- getOption('cdrgam.memory_limit_bytes')
    process_rss <- max(c(
        cdrgam:::.cdrgam_proc_memory('VmRSS'),
        cdrgam:::.cdrgam_r_memory()
    ), na.rm=TRUE)
    options(cdrgam.memory_limit_bytes=process_rss + 1024^3)
    concurrent <- cdrgam:::.cdrgam_gaulss_sparse_score(
        assemblies, retained$solution, exp(log_sp), family, workers=2L
    )
    plan <- attr(concurrent, 'score_plan')
    options(cdrgam.memory_limit_bytes=process_rss + 100 * 1024)
    streamed <- cdrgam:::.cdrgam_gaulss_sparse_score(
        assemblies, retained$solution, exp(log_sp), family, workers=2L
    )
    options(cdrgam.memory_limit_bytes=previous_memory_limit)
    streamed_timing <- attr(streamed, 'score_plan')$batch_timings[[1L]]
    stopifnot(
        isTRUE(all.equal(
            as.numeric(reference), as.numeric(serial), tolerance=1e-10
        )),
        isTRUE(all.equal(
            as.numeric(serial), as.numeric(concurrent), tolerance=1e-10
        )),
        isTRUE(all.equal(
            as.numeric(serial), as.numeric(streamed), tolerance=1e-10
        )),
        identical(
            plan$parallel_axis,
            'likelihood directions and inverse columns'
        ),
        identical(plan$workers, 2L),
        length(plan$batch_timings) == plan$batches,
        all(c(
            'rhs', 'coefficient_solve', 'likelihood_derivatives',
            'likelihood_direction_workers',
            'likelihood_design_cache_mb',
            'selected_inverse_trace', 'quadratic'
        ) %in% names(plan$batch_timings[[1L]])),
        plan$batch_timings[[1L]][['likelihood_direction_workers']] == 2,
        plan$batch_timings[[1L]][['likelihood_design_cache_mb']] > 0,
        streamed_timing[['likelihood_direction_workers']] == 1,
        streamed_timing[['likelihood_design_cache_mb']] == 0
    )
}

retained <- criterion(log_sp, retain=TRUE)
probe_count <- 17L
probes <- cdrgam:::.deterministic_rademacher(
    length(retained$solution$coefficients), probe_count
)
stochastic <- cdrgam:::.cdrgam_gaulss_sparse_stochastic_score(
    assemblies, retained$solution, exp(log_sp), family, probes
)
dimensions <- vapply(assemblies, `[[`, integer(1), 'dimension')
ranges <- split(seq_len(sum(dimensions)), rep(names(dimensions), dimensions))
counts <- vapply(
    assemblies,
    function(assembly) length(assembly$penalty_components),
    integer(1)
)
tasks <- list()
penalty_scores <- numeric()
offset <- 0L
for (parameter in names(assemblies)) {
    count <- counts[[parameter]]
    indices <- offset + seq_len(count)
    local_scores <- cdrgam:::.sparse_penalty_logdet_score(
        assemblies[[parameter]]$blocks, exp(log_sp)[indices], count
    )
    for (local_index in seq_len(count)) {
        global_index <- offset + local_index
        tasks[[global_index]] <- list(
            parameter=parameter,
            local_index=local_index,
            global_index=global_index
        )
        penalty_scores[[global_index]] <- local_scores[[local_index]]
    }
    offset <- offset + count
}
penalties <- lapply(tasks, function(task) {
    cdrgam:::.cdrgam_distributional_global_penalty(
        assemblies, ranges, task$parameter, task$local_index,
        exp(log_sp)[[task$global_index]]
    )
})
right_hand_sides <- do.call(cbind, lapply(penalties, function(penalty) {
    as.numeric(penalty %*% retained$solution$coefficients)
}))
coefficient_derivatives <- -as.matrix(cdrgam:::.cdr_factor_solve(
    retained$solution$factor, right_hand_sides
))
likelihood_derivatives <-
    cdrgam:::.cdrgam_gaulss_sparse_hessian_directions(
        assemblies, retained$solution, coefficient_derivatives,
        cdrgam:::.cdrgam_gaulss_b(family)
    )
system_derivatives <- Map(function(likelihood, penalty) {
    Matrix::forceSymmetric(likelihood + penalty, uplo='U')
}, likelihood_derivatives, penalties)
inverse_probes <- as.matrix(cdrgam:::.cdr_factor_solve(
    retained$solution$factor, probes
))
explicit_trace <- vapply(system_derivatives, function(derivative) {
    mean(colSums(inverse_probes * as.matrix(derivative %*% probes)))
}, numeric(1))
quadratic <- vapply(seq_along(penalties), function(index) {
    as.numeric(Matrix::crossprod(
        retained$solution$coefficients,
        penalties[[index]] %*% retained$solution$coefficients
    ))
}, numeric(1))
explicit_stochastic <- quadratic + explicit_trace - penalty_scores
stochastic_plan <- attr(stochastic, 'score_plan')
stopifnot(
    isTRUE(all.equal(
        as.numeric(stochastic), explicit_stochastic, tolerance=1e-10
    )),
    identical(stochastic_plan$method, 'stochastic matrix-free'),
    identical(stochastic_plan$probes, probe_count)
)

perturbed_sp <- exp(log_sp + c(0.05, -0.03, 0.02, -0.04))
cold <- cdrgam:::.cdrgam_gaulss_sparse_fixed(
    assemblies,
    perturbed_sp,
    cdrgam:::.cdrgam_gaulss_b(family)
)
warm <- cdrgam:::.cdrgam_gaulss_sparse_fixed(
    assemblies,
    perturbed_sp,
    cdrgam:::.cdrgam_gaulss_b(family),
    initial=retained$solution
)
stopifnot(
    isTRUE(warm$warm_started),
    isTRUE(cold$converged),
    isTRUE(warm$converged),
    warm$step_mode %in% c('joint', 'location-only', 'scale-only'),
    is.finite(warm$maximum_damping),
    abs(cold$objective - warm$objective) < 1e-7,
    max(abs(cold$coefficients - warm$coefficients)) < 2e-5,
    warm$iterations <= cold$iterations
)

cat('Sparse distributional score checks passed.\n')
