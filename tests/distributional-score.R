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

indefinite_system <- Matrix::Diagonal(x=c(1, -0.005))
indefinite_factor <- cdrgam:::.cdrgam_gaulss_step_factor(
    indefinite_system, supernodal=FALSE
)
stopifnot(
    indefinite_factor$damping >= 0.01,
    all(is.finite(cdrgam:::.cdr_factor_solve(
        indefinite_factor$factor, c(1, 0)
    )))
)

precondition_location <- Matrix::Diagonal(x=c(2, 3))
precondition_scale <- Matrix::Diagonal(x=c(4, 5))
precondition_cross <- matrix(c(0.2, -0.1, 0.3, 0.4), nrow=2L)
precondition_factors <- list(
    location=cdrgam:::.cdrgam_gaulss_step_factor(
        precondition_location, supernodal=FALSE
    )$factor,
    scale=cdrgam:::.cdrgam_gaulss_step_factor(
        precondition_scale, supernodal=FALSE
    )$factor
)
precondition_split <- function(value) list(
    location=value[1:2], scale=value[3:4]
)
precondition_rhs <- c(1, -2, 3, 4)
precondition_triangular <- cdrgam:::.cdrgam_gaulss_block_precondition(
    precondition_rhs,
    precondition_factors,
    precondition_split,
    mode='triangular',
    cross_scale=function(location) {
        as.numeric(Matrix::t(precondition_cross) %*% location)
    }
)
precondition_system <- rbind(
    cbind(as.matrix(precondition_location), matrix(0, 2, 2)),
    cbind(t(precondition_cross), as.matrix(precondition_scale))
)
stopifnot(isTRUE(all.equal(
    precondition_triangular,
    as.numeric(solve(precondition_system, precondition_rhs)),
    tolerance=1e-12
)))

schur_systems <- list(
    location=Matrix::Diagonal(x=c(4, 5, 6)),
    scale=Matrix::Diagonal(x=c(3, 7))
)
schur_cross <- matrix(c(0.4, -0.2, 0.1, 0.3, 0.2, -0.1), nrow=3L)
schur_factors <- Map(function(system) {
    cdrgam:::.cdrgam_gaulss_step_factor(
        system, supernodal=FALSE
    )$factor
}, schur_systems)
schur_gradients <- list(location=c(1, -2, 0.5), scale=c(-1, 3))
schur_reference_system <- rbind(
    cbind(as.matrix(schur_systems$location), schur_cross),
    cbind(t(schur_cross), as.matrix(schur_systems$scale))
)
schur_reference <- as.numeric(solve(
    schur_reference_system,
    -unlist(schur_gradients, use.names=FALSE)
))
schur_auto <- cdrgam:::.cdrgam_gaulss_schur_solve(
    schur_systems,
    schur_factors,
    schur_gradients,
    cross_location=function(scale) as.numeric(schur_cross %*% scale),
    cross_scale=function(location) {
        as.numeric(t(schur_cross) %*% location)
    },
    tolerance=1e-10,
    maxit=10L
)
schur_reverse <- cdrgam:::.cdrgam_gaulss_schur_solve(
    schur_systems,
    schur_factors,
    schur_gradients,
    cross_location=function(scale) as.numeric(schur_cross %*% scale),
    cross_scale=function(location) {
        as.numeric(t(schur_cross) %*% location)
    },
    tolerance=1e-10,
    maxit=10L,
    eliminate='scale'
)
schur_inexact <- cdrgam:::.cdrgam_gaulss_schur_solve(
    schur_systems,
    schur_factors,
    schur_gradients,
    cross_location=function(scale) as.numeric(schur_cross %*% scale),
    cross_scale=function(location) {
        as.numeric(t(schur_cross) %*% location)
    },
    tolerance=1e-14,
    maxit=1L,
    preconditioner_rank=0L
)
batch_applications <- 0L
batched_preconditioner <- cdrgam:::.cdrgam_gaulss_schur_preconditioner(
    diagonal_solve=function(values) values / 2,
    correction=function(values) {
        batch_applications <<- batch_applications + 1L
        values / 4
    },
    right_hand_side=c(1, 2, 3, 4),
    rank=3L
)
stopifnot(
    isTRUE(schur_auto$converged),
    schur_auto$schur_preconditioner_rank > 0L,
    identical(schur_auto$eliminate, 'location'),
    identical(schur_auto$remaining, 'scale'),
    isTRUE(all.equal(
        schur_auto$solution, schur_reference, tolerance=1e-8
    )),
    isTRUE(schur_reverse$converged),
    schur_reverse$schur_preconditioner_rank > 0L,
    identical(schur_reverse$eliminate, 'scale'),
    identical(schur_reverse$remaining, 'location'),
    isTRUE(all.equal(
        schur_reverse$solution, schur_reference, tolerance=1e-8
    )),
    !isTRUE(schur_inexact$converged),
    all(is.finite(schur_inexact$solution)),
    sum(unlist(schur_gradients, use.names=FALSE) *
        schur_inexact$solution) < 0,
    identical(batch_applications, 1L),
    identical(batched_preconditioner$rank, 3L),
    all(is.finite(batched_preconditioner$apply(c(1, 2, 3, 4))))
)

numerically_converged <- cdrgam:::.cdrgam_gaulss_numerically_converged
inner_algorithm <- cdrgam:::.cdrgam_gaulss_inner_algorithm
forcing_tolerance <- cdrgam:::.cdrgam_gaulss_forcing_tolerance
relaxed_forcing_tolerance <-
    cdrgam:::.cdrgam_gaulss_relaxed_forcing_tolerance
relative_convergence <- cdrgam:::.cdrgam_gaulss_relative_convergence
anderson_candidate <- cdrgam:::.cdrgam_anderson_candidate
anderson_map <- function(value) c(0.8 * value[[1L]] + 1, 0.5 * value[[2L]] - 1)
anderson_states <- anderson_residuals <- list()
anderson_iterate <- c(0, 0)
for (anderson_index in seq_len(3L)) {
    mapped <- anderson_map(anderson_iterate)
    anderson_states[[anderson_index]] <- anderson_iterate
    anderson_residuals[[anderson_index]] <- mapped - anderson_iterate
    anderson_iterate <- mapped
}
anderson_accelerated <- anderson_candidate(
    anderson_states, anderson_residuals, c(location=1L, scale=1L)
)
anderson_fixed_point <- c(5, -2)
krylov_events <- character()
krylov_system <- matrix(c(4, 2, 1, 3), nrow=2L)
krylov_solution <- cdrgam:::.cdrgam_bicgstab(
    c(1, 2),
    product=function(value) as.numeric(krylov_system %*% value),
    preconditioner=identity,
    tolerance=1e-10,
    maxit=10L,
    progress=function(event, ...) {
        krylov_events <<- c(krylov_events, event)
    }
)
krylov_warm <- cdrgam:::.cdrgam_bicgstab(
    c(1, 2),
    product=function(value) as.numeric(krylov_system %*% value),
    preconditioner=identity,
    tolerance=1e-3,
    maxit=10L,
    initial_solution=krylov_solution$solution + c(1e-5, -1e-5)
)
krylov_rejected <- cdrgam:::.cdrgam_bicgstab(
    c(1, 2),
    product=function(value) as.numeric(krylov_system %*% value),
    preconditioner=identity,
    tolerance=1e-10,
    maxit=10L,
    initial_solution=c(100, -100)
)
fgmres_solution <- cdrgam:::.cdrgam_fgmres(
    c(1, 2),
    product=function(value) as.numeric(krylov_system %*% value),
    preconditioner=identity,
    tolerance=1e-10,
    maxit=10L
)
stopifnot(
    numerically_converged(-6.2743e-9, 1.0594e-8, 299647.5249),
    !numerically_converged(6.2743e-9, 1.0594e-8, 299647.5249),
    !numerically_converged(-1, 1.0594e-8, 299647.5249),
    !numerically_converged(-6.2743e-9, 1, 299647.5249),
    identical(inner_algorithm('auto', 255L), 'joint'),
    identical(inner_algorithm('auto', 256L), 'mixed'),
    identical(inner_algorithm('schur', 1000L), 'schur'),
    identical(inner_algorithm('joint', 1000L), 'joint'),
    identical(forcing_tolerance(100), 0.1),
    identical(forcing_tolerance(50, 100), 0.1),
    isTRUE(all.equal(forcing_tolerance(0.01, 1), 9e-4)),
    identical(forcing_tolerance(0, 100), 1e-6),
    isTRUE(all.equal(relaxed_forcing_tolerance(0.1), sqrt(0.1))),
    identical(relaxed_forcing_tolerance(0.36), 0.5),
    isTRUE(all.equal(relaxed_forcing_tolerance(1e-4), 0.01)),
    isTRUE(relative_convergence(
        1e5, 1e5 - 1e-4, -5e-4, tolerance=1e-8
    )$converged),
    !isTRUE(relative_convergence(
        1e5, 1e5 - 1e-4, -5, tolerance=1e-8
    )$converged),
    !isTRUE(relative_convergence(
        1e5, 1e5 + 1e-4, -5e-4, tolerance=1e-8
    )$converged),
    all(is.finite(anderson_accelerated)),
    sum((anderson_accelerated - anderson_fixed_point)^2) <
        sum((anderson_iterate - anderson_fixed_point)^2),
    isTRUE(krylov_solution$converged),
    isTRUE(all.equal(
        krylov_solution$solution,
        as.numeric(solve(krylov_system, c(1, 2))),
        tolerance=1e-8
    )),
    krylov_solution$operator_applications > 0L,
    krylov_solution$preconditioner_applications > 0L,
    isTRUE(krylov_warm$converged),
    isTRUE(krylov_warm$warm_started),
    identical(krylov_warm$iterations, 0L),
    krylov_warm$operator_applications <
        krylov_solution$operator_applications,
    isTRUE(krylov_rejected$converged),
    !isTRUE(krylov_rejected$warm_started),
    krylov_rejected$initial_residual > 1,
    isTRUE(fgmres_solution$converged),
    isTRUE(all.equal(
        fgmres_solution$solution,
        as.numeric(solve(krylov_system, c(1, 2))),
        tolerance=1e-8
    )),
    fgmres_solution$operator_applications <=
        fgmres_solution$iterations + 1L,
    identical(
        fgmres_solution$operator_applications,
        fgmres_solution$preconditioner_applications
    ),
    all(c(
        'Krylov solve started',
        'Krylov operator started',
        'Krylov operator complete',
        'Krylov solve complete'
    ) %in% krylov_events)
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
    cdrgam:::.cdrgam_gaulss_b(family),
    algorithm='alternating'
)
joint <- cdrgam:::.cdrgam_gaulss_sparse_fixed(
    assemblies,
    perturbed_sp,
    cdrgam:::.cdrgam_gaulss_b(family),
    algorithm='joint'
)
fixed_progress_events <- character()
mixed <- cdrgam:::.cdrgam_gaulss_sparse_fixed(
    assemblies,
    perturbed_sp,
    cdrgam:::.cdrgam_gaulss_b(family),
    algorithm='mixed',
    progress=function(event, ...) {
        fixed_progress_events <<- c(fixed_progress_events, event)
    }
)
schur <- cdrgam:::.cdrgam_gaulss_sparse_fixed(
    assemblies,
    perturbed_sp,
    cdrgam:::.cdrgam_gaulss_b(family),
    algorithm='schur'
)
mixed_partial <- cdrgam:::.cdrgam_gaulss_sparse_fixed(
    assemblies,
    perturbed_sp,
    cdrgam:::.cdrgam_gaulss_b(family),
    maxit=2L,
    algorithm='mixed'
)
mixed_continued <- cdrgam:::.cdrgam_gaulss_sparse_fixed(
    assemblies,
    perturbed_sp,
    cdrgam:::.cdrgam_gaulss_b(family),
    initial=mixed_partial,
    algorithm='mixed'
)
q <- exp(mixed$linear_predictors[, 'scale'])
observed_moments <- cdrgam:::.cdrgam_gaulss_sparse_moments(
    assemblies,
    mixed$parameter_coefficients,
    mixed$sigma,
    q,
    mixed$residuals,
    mixed$prior_weights,
    include_observed=TRUE
)
cached_assemblies <- lapply(assemblies, function(assembly) {
    assembly$cached_design <- cdrgam:::.cdrgam_sparse_design_chunk(
        assembly, seq_len(assembly$observation_count)
    )
    assembly$crossprod_chunk_size <- max(
        1L, assembly$observation_count %/% 5L
    )
    assembly
})
cached_observed_moments <- cdrgam:::.cdrgam_gaulss_sparse_moments(
    cached_assemblies,
    mixed$parameter_coefficients,
    mixed$sigma,
    q,
    mixed$residuals,
    mixed$prior_weights,
    include_observed=TRUE,
    workers=2L
)
scale_gradient_only <- cdrgam:::.cdrgam_gaulss_sparse_moments(
    assemblies,
    mixed$parameter_coefficients,
    mixed$sigma,
    q,
    mixed$residuals,
    mixed$prior_weights,
    include_information=FALSE,
    parameters='scale'
)
cross_vectors <- lapply(assemblies, function(assembly) {
    stats::rnorm(assembly$dimension)
})
cross_weight <- mixed$prior_weights * 2 * q * mixed$residuals / mixed$sigma^3
matrix_free_cross <- cdrgam:::.cdrgam_gaulss_cross_product(
    assemblies, cross_weight, cross_vectors
)
matrix_free_scale <- cdrgam:::.cdrgam_gaulss_cross_product(
    assemblies, cross_weight, cross_vectors, outputs='scale'
)
cross_batches <- lapply(cross_vectors, function(vector) {
    cbind(vector, stats::rnorm(length(vector)), -0.5 * vector)
})
matrix_free_batch <- cdrgam:::.cdrgam_gaulss_cross_product(
    assemblies, cross_weight, cross_batches
)
matrix_free_columns <- lapply(seq_len(ncol(cross_batches$location)), function(i) {
    cdrgam:::.cdrgam_gaulss_cross_product(
        assemblies, cross_weight,
        lapply(cross_batches, function(values) values[, i])
    )
})
stopifnot(
    cached_observed_moments$batches == 2L,
    cached_observed_moments$batches < cached_observed_moments$chunks,
    isTRUE(all.equal(
        cached_observed_moments$gradients,
        observed_moments$gradients,
        tolerance=1e-9
    )),
    isTRUE(all.equal(
        cached_observed_moments$information,
        observed_moments$information,
        tolerance=1e-9
    )),
    is.null(scale_gradient_only$information$scale),
    isTRUE(all.equal(
        scale_gradient_only$gradients$scale,
        observed_moments$gradients$scale,
        tolerance=1e-10
    )),
    isTRUE(all.equal(
        matrix_free_cross$location,
        as.numeric(observed_moments$information$cross %*% cross_vectors$scale),
        tolerance=1e-10
    )),
    isTRUE(all.equal(
        matrix_free_cross$scale,
        as.numeric(Matrix::t(observed_moments$information$cross) %*%
            cross_vectors$location),
        tolerance=1e-10
    )),
    is.null(matrix_free_scale$location),
    isTRUE(all.equal(
        matrix_free_scale$scale,
        matrix_free_cross$scale,
        tolerance=1e-10
    )),
    isTRUE(all.equal(
        unname(matrix_free_batch$location),
        unname(do.call(cbind, lapply(
            matrix_free_columns, `[[`, 'location'
        ))),
        tolerance=1e-10
    )),
    isTRUE(all.equal(
        unname(matrix_free_batch$scale),
        unname(do.call(cbind, lapply(
            matrix_free_columns, `[[`, 'scale'
        ))),
        tolerance=1e-10
    ))
)
warm <- cdrgam:::.cdrgam_gaulss_sparse_fixed(
    assemblies,
    perturbed_sp,
    cdrgam:::.cdrgam_gaulss_b(family),
    initial=retained$solution,
    algorithm='alternating'
)
retry_delay <- cdrgam:::.cdrgam_gaulss_retry_delay
stopifnot(
    identical(retry_delay(0L), 0L),
    identical(retry_delay(1L), 1L),
    identical(retry_delay(2L), 2L),
    identical(retry_delay(3L), 4L),
    identical(retry_delay(4L), 8L),
    identical(retry_delay(8L), 8L),
    isTRUE(warm$warm_started),
    isTRUE(cold$converged),
    isTRUE(warm$converged),
    isTRUE(joint$converged),
    isTRUE(mixed$converged),
    isTRUE(schur$converged),
    !isTRUE(mixed_partial$converged),
    isTRUE(mixed_continued$converged),
    identical(cold$algorithm, 'alternating'),
    identical(joint$algorithm, 'joint'),
    identical(mixed$algorithm, 'mixed'),
    identical(schur$algorithm, 'schur'),
    is.logical(mixed$coupled_disabled),
    length(mixed$coupled_disabled) == 1L,
    is.integer(mixed$coupled_failures),
    length(mixed$coupled_failures) == 1L,
    mixed$coupled_failures >= 0L,
    is.integer(mixed$coupled_consecutive_failures),
    length(mixed$coupled_consecutive_failures) == 1L,
    mixed$coupled_consecutive_failures >= 0L,
    is.integer(mixed$coupled_retry_delay),
    length(mixed$coupled_retry_delay) == 1L,
    mixed$coupled_retry_delay >= 0L,
    all(c(
        'solve started', 'iteration started', 'moments started',
        'moments complete', 'solve complete'
    ) %in% fixed_progress_events),
    mixed$linear_iterations >= 0L,
    is.logical(mixed$linear_warm_started),
    length(mixed$linear_warm_started) == 1L,
    is.numeric(mixed$linear_initial_residual),
    length(mixed$linear_initial_residual) == 1L,
    is.character(mixed$linear_preconditioner),
    length(mixed$linear_preconditioner) == 1L,
    is.na(mixed$linear_preconditioner) ||
        mixed$linear_preconditioner %in% c(
            'diagonal', 'triangular', 'schur'
        ),
    is.finite(mixed$relative_reduction),
    is.finite(mixed$relative_predicted_reduction),
    identical(mixed$relative_tolerance, 1e-8),
    is.integer(mixed$anderson_attempts),
    is.integer(mixed$anderson_acceptances),
    mixed$anderson_attempts >= mixed$anderson_acceptances,
    is.null(mixed$coupled_linear_solution) ||
        length(mixed$coupled_linear_solution) == length(mixed$coefficients),
    mixed$linear_iterations == 0L || (
        is.finite(mixed$linear_residual) &&
        is.finite(mixed$linear_tolerance) &&
        is.finite(mixed$linear_acceptance_tolerance) &&
        is.logical(mixed$linear_converged) &&
        length(mixed$linear_converged) == 1L &&
        mixed$linear_residual <= mixed$linear_acceptance_tolerance
    ),
    schur$linear_iterations == 0L || (
        is.finite(schur$linear_residual) &&
        is.finite(schur$linear_tolerance) &&
        is.finite(schur$linear_acceptance_tolerance) &&
        schur$linear_acceptance_tolerance >= schur$linear_tolerance &&
        schur$linear_residual <= schur$linear_acceptance_tolerance
    ),
    identical(mixed_partial$mixed_iterations, 2L),
    warm$step_mode %in% c(
        'alternating-location-scale', 'alternating-location',
        'alternating-scale', 'alternating-none',
        'alternating-location-scale-anderson',
        'alternating-location-anderson',
        'alternating-scale-anderson'
    ),
    is.finite(warm$maximum_damping),
    abs(cold$objective - warm$objective) < 1e-7,
    abs(cold$objective - joint$objective) < 1e-7,
    abs(cold$objective - mixed$objective) < 1e-7,
    abs(cold$objective - schur$objective) < 1e-7,
    abs(mixed$objective - mixed_continued$objective) < 1e-7,
    max(abs(cold$coefficients - warm$coefficients)) < 1e-4,
    max(abs(cold$linear_predictors - warm$linear_predictors)) < 1e-4,
    warm$iterations <= cold$iterations,
    cold$iterations < joint$iterations
)

cat('Sparse distributional score checks passed.\n')
