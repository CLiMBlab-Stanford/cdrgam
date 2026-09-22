library(cdrgam)

process_memory <- max(c(
    cdrgam:::.cdrgam_proc_memory('VmRSS'),
    cdrgam:::.cdrgam_r_memory()
), na.rm=TRUE)
previous_memory_limit <- getOption('cdrgam.memory_limit_bytes')
options(cdrgam.memory_limit_bytes=process_memory + 1024^2)
memory_probe <- cdrgam:::.cdrgam_memory_availability()
options(cdrgam.memory_limit_bytes=previous_memory_limit)
stopifnot(
    is.finite(process_memory),
    identical(memory_probe$source, 'option'),
    is.finite(memory_probe$available_bytes),
    memory_probe$available_bytes > 0
)

simulation <- simulate_cdr(
    list(x=function(lag) exp(-lag)),
    n_impulses=60,
    n_responses=80,
    duration=10,
    window=2,
    noise_sd=1,
    seed=9
)
formula <- response ~ irf(x, window=c(0, 2), k=5) - irf(1)
design <- prepare_cdrgam(
    formula,
    simulation$impulses,
    simulation$responses,
    quiet=TRUE
)

# A progress callback can observe structured phase/evaluation records. Raising
# an error simulates interruption after several completed objective calls.
checkpoint <- tempfile('cdrgam-sparse-checkpoint-', fileext='.rds')
events <- list()
interrupted <- tryCatch(
    cdrgam.fit(
        design,
        backend='sparse',
        checkpoint=checkpoint,
        solver_trace=function(event) {
            events[[length(events) + 1L]] <<- event
            if (event$event %in% c('new best', 'evaluation') &&
                    event$evaluation >= 12L) {
                stop('simulated interruption')
            }
        },
        sparse_control=list(gradient='exact', outer_optimizer='lbfgsb')
    ),
    error=function(error) conditionMessage(error)
)
stopifnot(
    identical(interrupted, 'simulated interruption'),
    file.exists(checkpoint),
    any(vapply(events, `[[`, character(1), 'phase') ==
        'smoothing-parameter optimization')
)
partial <- readRDS(checkpoint)
stopifnot(
    identical(partial$checkpoint_version, 1L),
    identical(partial$stage, 'optimization'),
    partial$evaluation_count >= 10L,
    is.finite(partial$best_criterion)
)

resumed_events <- list()
resumed <- cdrgam.fit(
    design,
    backend='sparse',
    checkpoint=checkpoint,
    solver_trace=function(event) {
        resumed_events[[length(resumed_events) + 1L]] <<- event
    },
    sparse_control=list(gradient='exact', outer_optimizer='lbfgsb')
)
reference <- cdrgam.fit(
    design,
    backend='sparse',
    sparse_control=list(gradient='exact', outer_optimizer='lbfgsb')
)
complete <- readRDS(checkpoint)
stopifnot(
    identical(complete$stage, 'complete'),
    isTRUE(complete$converged),
    isTRUE(resumed$sparse$convergence$resumed),
    max(abs(coef(resumed) - coef(reference))) < 1e-8,
    abs(resumed$reml - reference$reml) < 1e-8,
    any(vapply(resumed_events, `[[`, character(1), 'event') ==
        'checkpoint resumed'),
    identical(reference$sparse$gradient_requested, 'exact'),
    identical(reference$sparse$outer_optimizer_requested, 'lbfgsb'),
    identical(reference$sparse$hessian_requested, 'auto'),
    reference$sparse$hessian %in% c('analytic', 'gradient'),
    identical(
        reference$sparse$hessian,
        reference$sparse$hessian_selection$method
    ),
    is.finite(reference$sparse$hessian_selection$estimated_peak_bytes),
    length(reference$sparse$optimizer_selection$measured_exact_seconds) > 0L
)

automatic <- cdrgam.fit(
    design,
    backend='sparse',
    sparse_control=list(gradient_cores=4L, hessian='none')
)
stopifnot(
    identical(automatic$sparse$gradient_requested, 'auto'),
    identical(automatic$sparse$gradient, 'exact'),
    identical(automatic$sparse$outer_optimizer_requested, 'auto'),
    identical(automatic$sparse$outer_optimizer, 'bfgs_trust'),
    identical(automatic$sparse$hessian_requested, 'none'),
    identical(automatic$sparse$hessian, 'none'),
    automatic$sparse$control$gradient_cores == if (
        .Platform$OS.type == 'windows'
    ) 1L else 4L,
    is.list(automatic$sparse$optimizer_selection),
    automatic$sparse$optimizer_selection$policy_version == 1L
)

# Trust-region stagnation uses the resolved Hessian policy for its curvature
# assessment. An automatic analytic assessment receives the memory-derived
# allowance, while an explicit gradient strategy differentiates exact scores.
previous_element_limit <- getOption('cdrgam.max_analytic_hessian_elements')
options(cdrgam.max_analytic_hessian_elements=1)
automatic_assessment_events <- list()
automatic_assessment <- cdrgam.fit(
    design,
    backend='sparse',
    solver_trace=function(event) {
        automatic_assessment_events[[length(automatic_assessment_events) + 1L]] <<-
            event
    },
    sparse_control=list(
        gradient='exact', outer_optimizer='bfgs_trust', hessian='auto',
        optimizer_trust_radius=1e-9, optimizer_maxit=5L
    )
)
gradient_assessment_events <- list()
gradient_assessment <- cdrgam.fit(
    design,
    backend='sparse',
    solver_trace=function(event) {
        gradient_assessment_events[[length(gradient_assessment_events) + 1L]] <<-
            event
    },
    sparse_control=list(
        gradient='exact', outer_optimizer='bfgs_trust', hessian='gradient',
        optimizer_trust_radius=1e-9, optimizer_maxit=5L
    )
)
options(cdrgam.max_analytic_hessian_elements=previous_element_limit)
automatic_assessments <- Filter(
    function(event) identical(
        event$event,
        'curvature convergence assessment started'
    ),
    automatic_assessment_events
)
gradient_assessments <- Filter(
    function(event) identical(
        event$event,
        'curvature convergence assessment started'
    ),
    gradient_assessment_events
)
stopifnot(
    length(automatic_assessments) > 0L,
    all(vapply(
        automatic_assessments,
        `[[`,
        character(1),
        'method'
    ) == 'analytic'),
    !grepl(
        'current limit is 1',
        automatic_assessment$sparse$convergence$message,
        fixed=TRUE
    ),
    length(gradient_assessments) > 0L,
    all(vapply(
        gradient_assessments,
        `[[`,
        character(1),
        'method'
    ) == 'gradient')
)

select_hessian <- getFromNamespace(
    '.sparse_select_hessian_method',
    'cdrgam'
)
supports <- list(1:100, 1:50)
ample_memory <- list(
    source='test',
    limit_bytes=8 * 1024^3,
    used_bytes=0,
    available_bytes=8 * 1024^3
)
limited_memory <- list(
    source='test',
    limit_bytes=1024^2,
    used_bytes=0,
    available_bytes=1024^2
)
ample_selection <- select_hessian(
    'auto', supports, dimension=100L, chunk_size=64L,
    memory=ample_memory
)
limited_selection <- select_hessian(
    'auto', supports, dimension=100L, chunk_size=64L,
    memory=limited_memory
)
explicit_selection <- select_hessian(
    'gradient', supports, dimension=100L, chunk_size=64L,
    memory=ample_memory
)
old_element_limit <- getOption('cdrgam.max_analytic_hessian_elements')
options(cdrgam.max_analytic_hessian_elements=1)
automatic_element_limit <- cdrgam:::.sparse_analytic_hessian_element_limit(
    'auto', ample_selection
)
explicit_element_limit <- cdrgam:::.sparse_analytic_hessian_element_limit(
    'analytic', ample_selection
)
options(cdrgam.max_analytic_hessian_elements=old_element_limit)
stopifnot(
    identical(ample_selection$method, 'analytic'),
    identical(limited_selection$method, 'gradient'),
    identical(explicit_selection$method, 'gradient'),
    ample_selection$retained_elements == 15000,
    ample_selection$retained_bytes == 120000,
    ample_selection$estimated_peak_bytes > ample_selection$retained_bytes,
    automatic_element_limit == ample_selection$retained_elements,
    explicit_element_limit == 1
)

compact_call <- getFromNamespace('.cdrgam_compact_call', 'cdrgam')(
    as.call(list(cdrgam.fit, design=design)),
    'cdrgam.fit',
    'design'
)
stopifnot(
    identical(compact_call[[1L]], quote(cdrgam.fit)),
    identical(compact_call$design, quote(design))
)

# A complete checkpoint bypasses a second outer optimization while retaining
# exactly the same solution.
completed_evaluations <- complete$evaluation_count
resumed_again <- cdrgam.fit(
    design,
    backend='sparse',
    checkpoint=checkpoint,
    sparse_control=list(gradient='exact', outer_optimizer='lbfgsb')
)
stopifnot(
    readRDS(checkpoint)$evaluation_count == completed_evaluations,
    max(abs(coef(resumed_again) - coef(resumed))) < 1e-12
)

# Checkpoints are tied to the actual objective, not merely matrix dimensions.
changed_responses <- simulation$responses
changed_responses$response[[1L]] <- changed_responses$response[[1L]] + 1
changed_design <- prepare_cdrgam(
    formula,
    simulation$impulses,
    changed_responses,
    quiet=TRUE
)
mismatch <- tryCatch(
    cdrgam.fit(
        changed_design,
        backend='sparse',
        checkpoint=checkpoint,
        sparse_control=list(gradient='exact', outer_optimizer='lbfgsb')
    ),
    error=function(error) conditionMessage(error)
)
stopifnot(grepl('does not match this model', mismatch, fixed=TRUE))

# The dense reference backend obeys the same checkpoint contract.
block_checkpoint <- tempfile('cdrgam-block-checkpoint-', fileext='.rds')
block <- cdrgam.fit(
    design,
    backend='block',
    checkpoint=block_checkpoint,
    solver_trace=0
)
block_resumed <- cdrgam.fit(
    design,
    backend='block',
    checkpoint=block_checkpoint
)
stopifnot(
    identical(readRDS(block_checkpoint)$stage, 'complete'),
    max(abs(coef(block) - coef(block_resumed))) < 1e-12
)

# Central finite-difference directions are independent and give the same fit
# when evaluated by forked workers.
if (.Platform$OS.type != 'windows') {
    finite_serial <- cdrgam.fit(
        design,
        backend='sparse',
        sparse_control=list(gradient='finite', gradient_cores=1L)
    )
    finite_parallel <- cdrgam.fit(
        design,
        backend='sparse',
        sparse_control=list(gradient='finite', gradient_cores=2L)
    )
    stopifnot(
        abs(finite_serial$reml - finite_parallel$reml) < 1e-8,
        max(abs(coef(finite_serial) - coef(finite_parallel))) < 1e-8,
        finite_parallel$sparse$parallel_gradient_factorizations > 0L
    )
}

# The lower-dimensional Hessian helper recovers an analytic quadratic, and
# the Gaussian scale identity reconstructs the same outer information matrix
# as direct finite differencing of the unprofiled REML criterion.
quadratic_hessian <- matrix(c(4, 1.5, 1.5, 3), 2L, 2L)
quadratic <- cdrgam:::.central_difference_hessian(
    function(x) drop(crossprod(x, quadratic_hessian %*% x)) / 2,
    c(0.4, -0.7),
    step=1e-4
)
stopifnot(
    max(abs(quadratic$hessian - quadratic_hessian)) < 1e-6,
    identical(quadratic$evaluations, 9L)
)
quadratic_gradient_hessian <- cdrgam:::.central_difference_jacobian(
    function(x) drop(quadratic_hessian %*% x),
    c(0.4, -0.7),
    step=1e-4
)
stopifnot(
    max(abs(quadratic_gradient_hessian$hessian - quadratic_hessian)) < 1e-10,
    identical(quadratic_gradient_hessian$evaluations, 4L)
)
direct_hessian <- cdrgam.fit(
    design,
    backend='sparse',
    sparse_control=list(
        gradient='exact', outer_optimizer='lbfgsb', hessian='optimhess'
    )
)
profiled_hessian <- cdrgam.fit(
    design,
    backend='sparse',
    sparse_control=list(
        gradient='exact', outer_optimizer='lbfgsb', hessian='profiled'
    )
)
gradient_hessian <- cdrgam.fit(
    design,
    backend='sparse',
    sparse_control=list(
        gradient='exact', outer_optimizer='lbfgsb', hessian='gradient'
    )
)
analytic_hessian <- cdrgam.fit(
    design,
    backend='sparse',
    sparse_control=list(
        gradient='exact', outer_optimizer='lbfgsb', hessian='analytic'
    )
)
hessian_scale <- max(1, max(abs(direct_hessian$outer.info$hess)))
stopifnot(
    max(abs(
        direct_hessian$outer.info$hess - profiled_hessian$outer.info$hess
    )) / hessian_scale < 2e-3,
    identical(profiled_hessian$sparse$hessian, 'profiled'),
    profiled_hessian$sparse$hessian_evaluations > 0L,
    max(abs(
        direct_hessian$outer.info$hess - gradient_hessian$outer.info$hess
    )) / hessian_scale < 2e-3,
    identical(gradient_hessian$sparse$hessian, 'gradient'),
    identical(
        gradient_hessian$sparse$hessian_evaluations,
        2L * length(gradient_hessian$sp)
    ),
    max(abs(
        direct_hessian$outer.info$hess - analytic_hessian$outer.info$hess
    )) / hessian_scale < 2e-3,
    identical(analytic_hessian$sparse$hessian, 'analytic'),
    identical(analytic_hessian$sparse$hessian_evaluations, 0L),
    analytic_hessian$sparse$hessian_rhs > 0L,
    max(abs(coef(direct_hessian) - coef(profiled_hessian))) < 1e-10
)
no_hessian <- cdrgam.fit(
    design,
    backend='sparse',
    sparse_control=list(
        gradient='exact', outer_optimizer='lbfgsb', hessian='none'
    )
)
no_hessian_error <- tryCatch(
    vcov(no_hessian, unconditional=TRUE),
    error=conditionMessage
)
stopifnot(
    is.null(no_hessian$outer.info$hess),
    identical(no_hessian$sparse$hessian_evaluations, 0L),
    is.na(no_hessian$sparse$convergence$hessian_positive_definite),
    grepl('hessian="none"', no_hessian_error, fixed=TRUE),
    max(abs(coef(direct_hessian) - coef(no_hessian))) < 1e-10
)
deferred_hessian <- cdrgam.fit(
    design,
    backend='sparse',
    sparse_control=list(
        gradient='exact', outer_optimizer='lbfgsb', hessian='defer'
    )
)
stopifnot(
    is.null(deferred_hessian$outer.info$hess),
    is.environment(deferred_hessian$sparse$deferred_hessian),
    !is.null(deferred_hessian$sparse$deferred_hessian$state)
)
deferred_path <- tempfile('cdrgam-deferred-hessian-', fileext='.rds')
saveRDS(deferred_hessian, deferred_path)
deferred_hessian <- readRDS(deferred_path)
deferred_unconditional <- vcov(deferred_hessian, unconditional=TRUE)
gradient_unconditional <- vcov(gradient_hessian, unconditional=TRUE)
stopifnot(
    max(abs(deferred_unconditional - gradient_unconditional)) < 1e-8,
    is.matrix(deferred_hessian$sparse$deferred_hessian$hessian),
    is.null(deferred_hessian$sparse$deferred_hessian$state),
    identical(
        deferred_hessian$sparse$deferred_hessian$evaluations,
        2L * length(deferred_hessian$sp)
    )
)
unlink(deferred_path)

# The experimental safeguarded outer BFGS method solves a bounded analytic
# quadratic and reaches the same REML solution with fewer sparse numeric
# updates on this small model.
quadratic_center <- c(-0.4, 0.7)
quadratic_progress <- list()
quadratic_bfgs <- cdrgam:::.safeguarded_outer_bfgs(
    par=c(1.2, -1),
    fn=function(x) drop(crossprod(
        x - quadratic_center,
        quadratic_hessian %*% (x - quadratic_center)
    )) / 2,
    gr=function(x) drop(quadratic_hessian %*% (x - quadratic_center)),
    lower=c(-2, -2),
    upper=c(2, 2),
    gradient_tolerance=1e-8,
    progress=function(record) {
        quadratic_progress[[length(quadratic_progress) + 1L]] <<- record
    }
)

# Interrupted trust optimization resumes with the exact dense curvature,
# radius, counters, and iteration number. It therefore follows precisely the
# same subsequent path as the uninterrupted run.
interrupted_state <- NULL
interrupted_bfgs <- tryCatch(
    cdrgam:::.safeguarded_outer_bfgs(
        par=c(1.2, -1),
        fn=function(x) drop(crossprod(
            x - quadratic_center,
            quadratic_hessian %*% (x - quadratic_center)
        )) / 2,
        gr=function(x) drop(quadratic_hessian %*% (x - quadratic_center)),
        lower=c(-2, -2),
        upper=c(2, 2),
        gradient_tolerance=1e-8,
        progress=function(record) {
            if (record$iteration == 2L &&
                    record$event %in% c('accepted', 'rejected')) {
                interrupted_state <<- record$optimizer_state
                stop('deliberate trust-optimizer interruption')
            }
        }
    ),
    error=function(error) error
)
resumed_progress <- list()
resumed_bfgs <- cdrgam:::.safeguarded_outer_bfgs(
    par=interrupted_state$parameters,
    fn=function(x) drop(crossprod(
        x - quadratic_center,
        quadratic_hessian %*% (x - quadratic_center)
    )) / 2,
    gr=function(x) drop(quadratic_hessian %*% (x - quadratic_center)),
    lower=c(-2, -2),
    upper=c(2, 2),
    gradient_tolerance=1e-8,
    state=interrupted_state,
    progress=function(record) {
        resumed_progress[[length(resumed_progress) + 1L]] <<- record
    }
)

# At a numerical trust-radius floor the analytic Hessian either certifies a
# practically stationary solution or supplies deliberately regularized
# positive curvature for a bounded recovery attempt.
flat_assessment <- cdrgam:::.analytic_outer_convergence_assessment(
    hessian=diag(c(2, 0)),
    gradient=c(1e-8, 1e-9),
    criterion=100,
    gradient_tolerance=1e-4,
    initial_radius=2
)
unresolved_assessment <- cdrgam:::.analytic_outer_convergence_assessment(
    hessian=diag(c(2, 0)),
    gradient=c(1e-2, 1e-3),
    criterion=100,
    gradient_tolerance=1e-4,
    initial_radius=2
)
certification_calls <- 0L
floor_state <- list(
    version=1L,
    optimizer='safeguarded_outer_bfgs',
    parameters=1,
    criterion=0.5,
    gradient=1,
    hessian=matrix(1),
    trust_radius=1e-9,
    iteration=0L,
    function_evaluations=1L,
    gradient_evaluations=1L,
    accepted_steps=0L,
    rejected_steps=0L,
    consecutive_rejections=0L,
    curvature_resets=0L,
    consecutive_small_steps=0L,
    recovery_resets=0L
)
floor_fit <- cdrgam:::.safeguarded_outer_bfgs(
    par=1,
    fn=function(x) round(x^2 / 2, 8),
    gr=function(x) x,
    lower=-2,
    upper=2,
    state=floor_state,
    convergence_assessment=function(...) {
        certification_calls <<- certification_calls + 1L
        list(
            converged=TRUE,
            message='test analytic certificate',
            diagnostics=list(test=TRUE)
        )
    }
)
recovery_calls <- 0L
recovered_fit <- cdrgam:::.safeguarded_outer_bfgs(
    par=1,
    fn=function(x) round(x^2 / 2, 8),
    gr=function(x) x,
    lower=-2,
    upper=2,
    gradient_tolerance=1e-8,
    state=floor_state,
    convergence_assessment=function(...) {
        recovery_calls <<- recovery_calls + 1L
        list(
            converged=FALSE,
            message='test recovery requested',
            restart_hessian=matrix(1),
            restart_radius=0.5,
            diagnostics=list(test=TRUE)
        )
    }
)
small_step_calls <- 0L
small_step_state <- floor_state
small_step_state$trust_radius <- 1e-7
small_step_state$consecutive_small_steps <- 3L
small_step_fit <- cdrgam:::.safeguarded_outer_bfgs(
    par=1,
    fn=function(x) x^2 / 2,
    gr=function(x) x,
    lower=-2,
    upper=2,
    state=small_step_state,
    convergence_assessment=function(...) {
        small_step_calls <<- small_step_calls + 1L
        list(
            converged=TRUE,
            message='test accepted-step stagnation certificate',
            diagnostics=list(test=TRUE)
        )
    }
)
trust_fit <- cdrgam.fit(
    design,
    backend='sparse',
    sparse_control=list(
        gradient='exact',
        outer_optimizer='bfgs_trust',
        optimizer_gradient_tolerance=2e-4
    )
)
trust_checkpoint <- tempfile('cdrgam-trust-checkpoint-', fileext='.rds')
trust_backend_fit <- cdrgam.fit(
    design,
    backend='sparse',
    checkpoint=trust_checkpoint,
    sparse_control=list(
        gradient='exact',
        outer_optimizer='bfgs_trust',
        optimizer_gradient_tolerance=2e-4
    )
)
trust_checkpoint_state <- readRDS(trust_checkpoint)
interrupted_trust_checkpoint <- tempfile(
    'cdrgam-interrupted-trust-checkpoint-',
    fileext='.rds'
)
interrupted_trust_fit <- tryCatch(
    cdrgam.fit(
        design,
        backend='sparse',
        checkpoint=interrupted_trust_checkpoint,
        sparse_control=list(
            gradient='exact',
            outer_optimizer='bfgs_trust',
            optimizer_gradient_tolerance=2e-4
        ),
        solver_trace=function(record) {
            if (identical(record$event, 'outer accepted')) {
                stop('deliberate sparse-trust interruption')
            }
        }
    ),
    error=function(error) error
)
interrupted_trust_state <- readRDS(interrupted_trust_checkpoint)
resumed_trust_events <- character()
resumed_trust_fit <- cdrgam.fit(
    design,
    backend='sparse',
    checkpoint=interrupted_trust_checkpoint,
    sparse_control=list(
        gradient='exact',
        outer_optimizer='bfgs_trust',
        optimizer_gradient_tolerance=2e-4
    ),
    solver_trace=function(record) {
        resumed_trust_events <<- c(resumed_trust_events, record$event)
    }
)
stopifnot(
    identical(quadratic_bfgs$convergence, 0L),
    max(abs(quadratic_bfgs$par - quadratic_center)) < 1e-7,
    inherits(interrupted_bfgs, 'error'),
    identical(resumed_progress[[1L]]$event, 'resumed'),
    isTRUE(all.equal(resumed_bfgs$par, quadratic_bfgs$par, tolerance=0)),
    identical(resumed_bfgs$counts, quadratic_bfgs$counts),
    identical(resumed_bfgs$hessian, quadratic_bfgs$hessian),
    isTRUE(flat_assessment$converged),
    !isTRUE(unresolved_assessment$converged),
    min(eigen(
        unresolved_assessment$restart_hessian,
        symmetric=TRUE,
        only.values=TRUE
    )$values) > 0,
    identical(floor_fit$convergence, 0L),
    identical(floor_fit$message, 'test analytic certificate'),
    identical(certification_calls, 1L),
    identical(recovered_fit$convergence, 0L),
    identical(recovered_fit$recovery_resets, 1L),
    identical(recovery_calls, 1L),
    any(vapply(recovered_fit$history, function(record) {
        identical(record$event, 'curvature_recovery')
    }, logical(1))),
    identical(small_step_fit$convergence, 0L),
    identical(
        small_step_fit$message,
        'test accepted-step stagnation certificate'
    ),
    identical(small_step_calls, 1L),
    identical(small_step_fit$iterations, 0L),
    identical(trust_fit$optimizer$convergence, 0L),
    abs(trust_fit$reml - reference$reml) < 1e-6,
    trust_fit$sparse$numeric_updates < reference$sparse$numeric_updates,
    identical(trust_fit$sparse$outer_optimizer, 'bfgs_trust'),
    identical(trust_backend_fit$cdrgam$backend, 'sparse'),
    identical(trust_backend_fit$sparse$gradient, 'exact'),
    identical(trust_backend_fit$sparse$outer_optimizer, 'bfgs_trust'),
    abs(trust_backend_fit$reml - trust_fit$reml) < 1e-10,
    max(abs(coef(trust_backend_fit) - coef(trust_fit))) < 1e-10,
    inherits(interrupted_trust_fit, 'error'),
    identical(interrupted_trust_state$stage, 'optimization'),
    is.list(interrupted_trust_state$optimizer_state),
    interrupted_trust_state$optimizer_state$iteration > 0L,
    'outer resumed' %in% resumed_trust_events,
    abs(resumed_trust_fit$reml - trust_backend_fit$reml) < 1e-10,
    max(abs(coef(resumed_trust_fit) - coef(trust_backend_fit))) < 1e-10,
    identical(quadratic_progress[[1L]]$event, 'initial'),
    identical(quadratic_progress[[length(quadratic_progress)]]$event,
        'finished'),
    all(vapply(quadratic_progress, function(record) {
        all(c(
            'projected_gradient_max', 'gradient_ratio', 'step_max',
            'trust_radius', 'predicted_improvement', 'actual_improvement',
            'acceptance_ratio', 'accepted_steps', 'rejected_steps',
            'consecutive_rejections', 'curvature_resets',
            'curvature_reset', 'step_type'
        ) %in% names(record))
    }, logical(1))),
    identical(trust_checkpoint_state$stage, 'complete'),
    is.list(trust_checkpoint_state$optimizer_state),
    identical(trust_checkpoint_state$optimizer_state$version, 1L),
    is.matrix(trust_checkpoint_state$optimizer_state$hessian),
    is.null(trust_checkpoint_state$optimizer_history[[1L]]$optimizer_state),
    length(trust_checkpoint_state$optimizer_history) > 1L,
    identical(trust_checkpoint_state$optimizer_progress$event, 'finished'),
    trust_checkpoint_state$optimizer_progress[[
        'postfit_hessian_factorizations'
    ]] >= 0L,
    identical(
        trust_backend_fit$sparse$optimizer_progress,
        trust_checkpoint_state$optimizer_progress
    ),
    length(trust_backend_fit$sparse$optimizer_history) > 1L
)

invalid_trust <- tryCatch(
    cdrgam.fit(
        design,
        backend='sparse',
        sparse_control=list(
            gradient='finite',
            outer_optimizer='bfgs_trust'
        )
    ),
    error=function(error) conditionMessage(error)
)
stopifnot(grepl('requires gradient="exact"', invalid_trust, fixed=TRUE))

unlink(c(
    checkpoint,
    block_checkpoint,
    trust_checkpoint,
    interrupted_trust_checkpoint
))
