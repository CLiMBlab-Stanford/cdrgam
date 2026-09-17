library(cdrgam)

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
    fit_cdrgam(
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
        sparse_control=list(gradient='exact')
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
resumed <- fit_cdrgam(
    design,
    backend='sparse',
    checkpoint=checkpoint,
    solver_trace=function(event) {
        resumed_events[[length(resumed_events) + 1L]] <<- event
    },
    sparse_control=list(gradient='exact')
)
reference <- fit_cdrgam(
    design,
    backend='sparse',
    sparse_control=list(gradient='exact')
)
complete <- readRDS(checkpoint)
stopifnot(
    identical(complete$stage, 'complete'),
    isTRUE(complete$converged),
    isTRUE(resumed$sparse$convergence$resumed),
    max(abs(coef(resumed) - coef(reference))) < 1e-8,
    abs(resumed$reml - reference$reml) < 1e-8,
    any(vapply(resumed_events, `[[`, character(1), 'event') ==
        'checkpoint resumed')
)

# A complete checkpoint bypasses a second outer optimization while retaining
# exactly the same solution.
completed_evaluations <- complete$evaluation_count
resumed_again <- fit_cdrgam(
    design,
    backend='sparse',
    checkpoint=checkpoint,
    sparse_control=list(gradient='exact')
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
    fit_cdrgam(
        changed_design,
        backend='sparse',
        checkpoint=checkpoint,
        sparse_control=list(gradient='exact')
    ),
    error=function(error) conditionMessage(error)
)
stopifnot(grepl('does not match this model', mismatch, fixed=TRUE))

# The dense reference backend obeys the same checkpoint contract.
block_checkpoint <- tempfile('cdrgam-block-checkpoint-', fileext='.rds')
block <- fit_cdrgam(
    design,
    backend='block',
    checkpoint=block_checkpoint,
    solver_trace=0
)
block_resumed <- fit_cdrgam(
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
    finite_serial <- fit_cdrgam(
        design,
        backend='sparse',
        sparse_control=list(gradient='finite', gradient_cores=1L)
    )
    finite_parallel <- fit_cdrgam(
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
direct_hessian <- fit_cdrgam(
    design,
    backend='sparse',
    sparse_control=list(gradient='exact', hessian='optimhess')
)
profiled_hessian <- fit_cdrgam(
    design,
    backend='sparse',
    sparse_control=list(gradient='exact', hessian='profiled')
)
hessian_scale <- max(1, max(abs(direct_hessian$outer.info$hess)))
stopifnot(
    max(abs(
        direct_hessian$outer.info$hess - profiled_hessian$outer.info$hess
    )) / hessian_scale < 2e-3,
    identical(profiled_hessian$sparse$hessian, 'profiled'),
    profiled_hessian$sparse$hessian_evaluations > 0L,
    max(abs(coef(direct_hessian) - coef(profiled_hessian))) < 1e-10
)

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
trust_fit <- fit_cdrgam(
    design,
    backend='sparse',
    sparse_control=list(
        gradient='exact',
        outer_optimizer='bfgs_trust',
        optimizer_gradient_tolerance=2e-4
    )
)
trust_checkpoint <- tempfile('cdrgam-trust-checkpoint-', fileext='.rds')
trust_backend_fit <- fit_cdrgam(
    design,
    backend='sparse_trust',
    checkpoint=trust_checkpoint
)
trust_checkpoint_state <- readRDS(trust_checkpoint)
stopifnot(
    identical(quadratic_bfgs$convergence, 0L),
    max(abs(quadratic_bfgs$par - quadratic_center)) < 1e-7,
    identical(trust_fit$optimizer$convergence, 0L),
    abs(trust_fit$reml - reference$reml) < 1e-6,
    trust_fit$sparse$numeric_updates < reference$sparse$numeric_updates,
    identical(trust_fit$sparse$outer_optimizer, 'bfgs_trust'),
    identical(trust_backend_fit$cdrgam$backend, 'sparse_trust'),
    identical(trust_backend_fit$sparse$gradient, 'exact'),
    identical(trust_backend_fit$sparse$outer_optimizer, 'bfgs_trust'),
    abs(trust_backend_fit$reml - trust_fit$reml) < 1e-10,
    max(abs(coef(trust_backend_fit) - coef(trust_fit))) < 1e-10,
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
    length(trust_checkpoint_state$optimizer_history) > 1L,
    identical(trust_checkpoint_state$optimizer_progress$event, 'finished'),
    trust_checkpoint_state$optimizer_progress[[
        'postfit_hessian_factorizations'
    ]] > 0L,
    identical(
        trust_backend_fit$sparse$optimizer_progress,
        trust_checkpoint_state$optimizer_progress
    ),
    length(trust_backend_fit$sparse$optimizer_history) > 1L
)

invalid_trust <- tryCatch(
    fit_cdrgam(
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

invalid_trust_backend <- tryCatch(
    fit_cdrgam(
        design,
        backend='sparse_trust',
        sparse_control=list(gradient='finite')
    ),
    error=function(error) conditionMessage(error)
)
stopifnot(grepl('requires sparse_control$gradient="exact"',
    invalid_trust_backend, fixed=TRUE))

unlink(c(checkpoint, block_checkpoint, trust_checkpoint))
