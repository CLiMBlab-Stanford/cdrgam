.solver_trace_level <- function(trace) {
    if (is.function(trace)) {
        return(list(level=2L, callback=trace, console=FALSE))
    }
    if (is.logical(trace) && length(trace) == 1L && !is.na(trace)) {
        return(list(
            level=if (trace) 1L else 0L,
            callback=NULL,
            console=TRUE
        ))
    }
    if (is.numeric(trace) && length(trace) == 1L && is.finite(trace) &&
            trace == as.integer(trace) && trace >= 0 && trace <= 3) {
        return(list(level=as.integer(trace), callback=NULL, console=TRUE))
    }
    stop('solver_trace must be FALSE, TRUE, an integer from 0 to 3, or a function')
}

.new_solver_reporter <- function(trace, backend) {
    control <- .solver_trace_level(trace)
    started <- proc.time()[['elapsed']]
    current_phase <- 'initialization'

    emit <- function(required_level=1L, event='progress', ..., force=FALSE) {
        fields <- list(...)
        elapsed <- proc.time()[['elapsed']] - started
        record <- c(list(
            backend=backend,
            event=event,
            phase=current_phase,
            elapsed_seconds=unname(elapsed)
        ), fields)
        if (!is.null(control$callback)) {
            control$callback(record)
        }
        if (isTRUE(force) ||
                (isTRUE(control$console) && control$level >= required_level)) {
            details <- if (length(fields)) {
                paste(vapply(names(fields), function(name) {
                    value <- fields[[name]]
                    if (length(value) > 6L) {
                        value <- c(value[seq_len(6L)], '...')
                    }
                    paste0(name, '=', paste(value, collapse=','))
                }, character(1)), collapse=' ')
            } else {
                ''
            }
            message(sprintf(
                'cdrgam %s [%s; %.1fs]%s%s',
                backend,
                current_phase,
                elapsed,
                if (nzchar(event)) paste0(' ', event) else '',
                if (nzchar(details)) paste0(': ', details) else ''
            ))
        }
        invisible(record)
    }

    phase <- function(name, ..., required_level=1L) {
        current_phase <<- name
        emit(required_level, event='phase', ...)
    }

    list(
        level=control$level,
        emit=emit,
        phase=phase,
        elapsed=function() proc.time()[['elapsed']] - started
    )
}

.checkpoint_read <- function(path) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    tryCatch(
        readRDS(path),
        error=function(error) stop(
            'Could not read solver checkpoint ', sQuote(path), ': ',
            conditionMessage(error),
            call.=FALSE
        )
    )
}

.checkpoint_write <- function(state, path) {
    if (is.null(path)) return(invisible(FALSE))
    directory <- dirname(path)
    if (!dir.exists(directory)) {
        if (!dir.create(directory, recursive=TRUE, showWarnings=FALSE)) {
            stop('Could not create checkpoint directory ', sQuote(directory))
        }
    }
    temporary <- tempfile(
        pattern=paste0('.', basename(path), '-'),
        tmpdir=directory
    )
    on.exit(unlink(temporary), add=TRUE)
    saveRDS(state, temporary, version=3)
    if (!file.rename(temporary, path)) {
        if (!file.copy(temporary, path, overwrite=TRUE)) {
            stop('Could not replace solver checkpoint ', sQuote(path))
        }
        unlink(temporary)
    }
    invisible(TRUE)
}

.solver_numeric_signature <- function(value) {
    value <- as.numeric(value)
    finite <- value[is.finite(value)]
    if (!length(finite)) {
        return(c(length=length(value), finite=0, sum=0, absolute=0, square=0))
    }
    c(
        length=length(value),
        finite=length(finite),
        sum=sum(finite),
        absolute=sum(abs(finite)),
        square=sum(finite^2)
    )
}

.solver_checkpoint_signature <- function(
        backend,
        formula,
        observation_count,
        dimension,
        sp_names,
        Xty,
        system
) {
    system_values <- if (methods::is(system, 'Matrix')) {
        methods::slot(system, 'x')
    } else {
        as.numeric(system)
    }
    list(
        version=1L,
        backend=backend,
        formula=paste(deparse(formula), collapse=''),
        observations=as.integer(observation_count),
        dimension=as.integer(dimension),
        smoothing_parameters=as.character(sp_names),
        Xty=.solver_numeric_signature(Xty),
        system_dimension=dim(system),
        system_nonzeros=as.numeric(Matrix::nnzero(system)),
        system_values=.solver_numeric_signature(system_values)
    )
}

.checkpoint_validate <- function(state, signature, path) {
    if (is.null(state)) return(NULL)
    # Checkpoints from the initial prototype contained only log_sp and
    # criterion. They remain valid as unverified warm starts.
    if (is.null(state$checkpoint_version)) {
        state$legacy <- TRUE
        return(state)
    }
    if (!identical(state$checkpoint_version, 1L) ||
            !identical(state$model_signature, signature)) {
        stop(
            'Checkpoint ', sQuote(path), ' does not match this model and ',
            'cannot be resumed. Supply a different path or remove the stale ',
            'checkpoint.',
            call.=FALSE
        )
    }
    state$legacy <- FALSE
    state
}

.new_checkpoint_state <- function(signature, backend, restart_count) {
    list(
        checkpoint_version=1L,
        model_signature=signature,
        backend=backend,
        stage='optimization',
        current_restart=0L,
        restart_count=as.integer(restart_count),
        completed_runs=list(),
        optimizer_progress=NULL,
        optimizer_state=NULL,
        optimizer_history=list(),
        current_log_sp=NULL,
        current_criterion=Inf,
        best_log_sp=NULL,
        best_criterion=Inf,
        evaluation_count=0L,
        updated_at=as.character(Sys.time())
    )
}

.checkpoint_update <- function(state, ..., path=NULL) {
    changes <- list(...)
    for (name in names(changes)) state[[name]] <- changes[[name]]
    state$updated_at <- as.character(Sys.time())
    .checkpoint_write(state, path)
    state
}

.deterministic_rademacher <- function(rows, columns, seed=1729L) {
    had_seed <- exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE)
    if (had_seed) old_seed <- get('.Random.seed', envir=.GlobalEnv)
    on.exit({
        if (had_seed) {
            assign('.Random.seed', old_seed, envir=.GlobalEnv)
        } else if (exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE)) {
            rm('.Random.seed', envir=.GlobalEnv)
        }
    }, add=TRUE)
    set.seed(seed)
    matrix(
        sample(c(-1, 1), rows * columns, replace=TRUE),
        nrow=rows,
        ncol=columns
    )
}

.central_difference_hessian <- function(fn, parameters, step=1e-3) {
    count <- length(parameters)
    out <- matrix(0, nrow=count, ncol=count)
    center <- fn(parameters)
    evaluations <- 1L
    plus <- minus <- numeric(count)
    for (i in seq_len(count)) {
        forward <- backward <- parameters
        forward[[i]] <- forward[[i]] + step
        backward[[i]] <- backward[[i]] - step
        plus[[i]] <- fn(forward)
        minus[[i]] <- fn(backward)
        evaluations <- evaluations + 2L
        out[i, i] <- (plus[[i]] - 2 * center + minus[[i]]) / step^2
    }
    if (count > 1L) {
        for (i in seq_len(count - 1L)) {
            for (j in seq.int(i + 1L, count)) {
                pp <- pm <- mp <- mm <- parameters
                pp[[i]] <- pp[[i]] + step
                pp[[j]] <- pp[[j]] + step
                pm[[i]] <- pm[[i]] + step
                pm[[j]] <- pm[[j]] - step
                mp[[i]] <- mp[[i]] - step
                mp[[j]] <- mp[[j]] + step
                mm[[i]] <- mm[[i]] - step
                mm[[j]] <- mm[[j]] - step
                value <- (fn(pp) - fn(pm) - fn(mp) + fn(mm)) /
                    (4 * step^2)
                evaluations <- evaluations + 4L
                out[i, j] <- out[j, i] <- value
            }
        }
    }
    list(hessian=out, evaluations=evaluations)
}

.central_difference_jacobian <- function(fn, parameters, step=1e-3) {
    count <- length(parameters)
    out <- matrix(0, nrow=count, ncol=count)
    for (j in seq_len(count)) {
        forward <- backward <- parameters
        forward[[j]] <- forward[[j]] + step
        backward[[j]] <- backward[[j]] - step
        forward_value <- fn(forward)
        backward_value <- fn(backward)
        if (!is.numeric(forward_value) || !is.numeric(backward_value) ||
                length(forward_value) != count ||
                length(backward_value) != count ||
                any(!is.finite(forward_value)) ||
                any(!is.finite(backward_value))) {
            stop('Central-difference Jacobian received a non-finite score')
        }
        out[, j] <- (forward_value - backward_value) / (2 * step)
    }
    # Roundoff and inexact sparse solves can make a differentiated gradient
    # very slightly asymmetric. A scalar objective has a symmetric Hessian.
    out <- (out + t(out)) / 2
    list(hessian=out, evaluations=2L * count)
}

.analytic_outer_convergence_assessment <- function(
        hessian,
        gradient,
        criterion,
        gradient_tolerance,
        initial_radius,
        objective_noise=0
) {
    hessian <- (as.matrix(hessian) + t(as.matrix(hessian))) / 2
    decomposition <- eigen(hessian, symmetric=TRUE)
    values <- decomposition$values
    vectors <- decomposition$vectors
    scale <- max(1, max(abs(values)))
    curvature_tolerance <- scale * sqrt(.Machine$double.eps)
    positive <- values > curvature_tolerance
    coordinates <- drop(crossprod(vectors, gradient))
    step_coordinates <- numeric(length(values))
    if (any(positive)) {
        step_coordinates[positive] <-
            -coordinates[positive] / values[positive]
    }
    newton_step <- drop(vectors %*% step_coordinates)
    predicted_improvement <- if (any(positive)) {
        sum(coordinates[positive]^2 / values[positive]) / 2
    } else {
        0
    }
    unresolved_gradient <- if (any(!positive)) {
        max(abs(coordinates[!positive]))
    } else {
        0
    }
    step_max <- if (length(newton_step)) max(abs(newton_step)) else 0
    # REML criteria scale with the number of observations. Improvements below
    # this relative floor are not reliably distinguishable from sparse
    # factorization roundoff and cannot justify movement along a flat outer
    # direction.
    baseline_objective_floor <- max(
        1e-8,
        1e-10 * (1 + abs(criterion))
    )
    if (!is.numeric(objective_noise) || length(objective_noise) != 1L ||
            !is.finite(objective_noise) || objective_noise < 0) {
        objective_noise <- 0
    }
    objective_floor <- max(baseline_objective_floor, objective_noise)
    converged <- is.finite(predicted_improvement) &&
        predicted_improvement <= objective_floor &&
        unresolved_gradient <= gradient_tolerance
    eigenvalue_floor <- scale * 1e-6
    restart_hessian <- vectors %*% (
        pmax(values, eigenvalue_floor) * t(vectors)
    )
    restart_hessian <- (restart_hessian + t(restart_hessian)) / 2
    restart_radius <- min(
        initial_radius,
        max(0.1, 4 * sqrt(sum(newton_step^2)))
    )
    diagnostics <- list(
        predicted_improvement=predicted_improvement,
        objective_floor=objective_floor,
        baseline_objective_floor=baseline_objective_floor,
        observed_objective_noise=objective_noise,
        newton_step_max=step_max,
        unresolved_gradient=unresolved_gradient,
        curvature_tolerance=curvature_tolerance,
        minimum_eigenvalue=min(values),
        positive_directions=sum(positive),
        flat_directions=sum(!positive),
        dimension=length(values)
    )
    list(
        converged=converged,
        message=if (converged) {
            'analytic practical convergence reached in the identifiable subspace'
        } else {
            'analytic Hessian found unresolved local improvement'
        },
        restart_hessian=restart_hessian,
        restart_radius=restart_radius,
        diagnostics=diagnostics
    )
}

# Experimental bounded, safeguarded dense-BFGS outer optimizer. Unlike
# stats::optim() with a numerical score, an accepted iteration needs one
# objective factorization; the exact score is then assembled from the retained
# factor using sparse solves. The dense approximation is only q-by-q, where q
# is the number of smoothing parameters.
.safeguarded_outer_bfgs <- function(
        par,
        fn,
        gr,
        lower,
        upper,
        maxit=100L,
        gradient_tolerance=1e-4,
        initial_radius=2,
        maximum_radius=10,
        progress=NULL,
        state=NULL,
        convergence_assessment=NULL,
        maximum_recovery_resets=2L,
        minimum_radius=1e-8
) {
    par <- pmin(upper, pmax(lower, as.numeric(par)))
    dimension <- length(par)
    maximum_recovery_resets <- as.integer(maximum_recovery_resets)
    if (length(maximum_recovery_resets) != 1L ||
            is.na(maximum_recovery_resets) || maximum_recovery_resets < 0L) {
        stop('maximum_recovery_resets must be a non-negative integer')
    }
    if (!is.numeric(minimum_radius) || length(minimum_radius) != 1L ||
            !is.finite(minimum_radius) || minimum_radius <= 0) {
        stop('minimum_radius must be positive')
    }
    projected_gradient <- function(parameters, score) {
        out <- score
        at_lower <- parameters <= lower + 1e-10
        at_upper <- parameters >= upper - 1e-10
        out[at_lower & out > 0] <- 0
        out[at_upper & out < 0] <- 0
        out
    }
    scalar_finite <- function(value) {
        is.numeric(value) && length(value) == 1L && is.finite(value)
    }
    state_counters <- if (is.list(state)) c(
        state$function_evaluations,
        state$gradient_evaluations,
        state$accepted_steps,
        state$rejected_steps,
        state$consecutive_rejections,
        state$curvature_resets,
        state$consecutive_small_steps,
        state$recovery_resets,
        state$iteration
    ) else numeric()
    valid_state <- !is.null(state) && is.list(state) &&
        identical(state$version, 1L) &&
        identical(state$optimizer, 'safeguarded_outer_bfgs') &&
        length(state$parameters) == dimension &&
        length(state$gradient) == dimension &&
        identical(dim(state$hessian), c(dimension, dimension)) &&
        all(is.finite(state$parameters)) && all(is.finite(state$gradient)) &&
        all(is.finite(state$hessian)) && scalar_finite(state$criterion) &&
        scalar_finite(state$trust_radius) && state$trust_radius > 0 &&
        length(state_counters) == 9L && all(is.finite(state_counters)) &&
        all(state_counters >= 0)
    if (!is.null(state) && !valid_state) {
        warning(
            'Saved trust-optimizer state is invalid; resuming from its ',
            'parameters with fresh curvature.',
            call.=FALSE
        )
        if (is.list(state) && length(state$parameters) == dimension &&
                all(is.finite(state$parameters))) {
            par <- pmin(upper, pmax(lower, as.numeric(state$parameters)))
        }
    }
    resumed <- isTRUE(valid_state)
    if (resumed) {
        par <- pmin(upper, pmax(lower, as.numeric(state$parameters)))
        value <- as.numeric(state$criterion)
        gradient <- as.numeric(state$gradient)
        hessian <- as.matrix(state$hessian)
        radius <- as.numeric(state$trust_radius)
        function_count <- as.integer(state$function_evaluations)
        gradient_count <- as.integer(state$gradient_evaluations)
        accepted_steps <- as.integer(state$accepted_steps)
        rejected <- as.integer(state$rejected_steps)
        consecutive_rejections <- as.integer(state$consecutive_rejections)
        curvature_resets <- as.integer(state$curvature_resets)
        consecutive_small <- as.integer(state$consecutive_small_steps)
        recovery_resets <- as.integer(state$recovery_resets)
        completed_iteration <- as.integer(state$iteration)
        counters <- c(
            function_count, gradient_count, accepted_steps, rejected,
            consecutive_rejections, curvature_resets, consecutive_small,
            recovery_resets, completed_iteration
        )
        if (anyNA(counters) || any(counters < 0L)) {
            stop('Saved trust-optimizer state has invalid counters')
        }
    } else {
        function_count <- 1L
        gradient_count <- 1L
        value <- fn(par)
        gradient <- gr(par)
        scale <- max(1, max(abs(gradient)))
        hessian <- diag(scale, dimension)
        radius <- initial_radius
        accepted_steps <- 0L
        rejected <- 0L
        consecutive_rejections <- 0L
        curvature_resets <- 0L
        consecutive_small <- 0L
        recovery_resets <- 0L
        completed_iteration <- 0L
    }
    optimizer_state <- function(iteration) list(
        version=1L,
        optimizer='safeguarded_outer_bfgs',
        parameters=par,
        criterion=value,
        gradient=gradient,
        hessian=hessian,
        trust_radius=radius,
        iteration=as.integer(iteration),
        function_evaluations=function_count,
        gradient_evaluations=gradient_count,
        accepted_steps=accepted_steps,
        rejected_steps=rejected,
        consecutive_rejections=consecutive_rejections,
        curvature_resets=curvature_resets,
        consecutive_small_steps=consecutive_small,
        recovery_resets=recovery_resets
    )
    history <- list()
    emit_progress <- function(record, save=TRUE) {
        if (save) history[[length(history) + 1L]] <<- record
        emitted <- record
        emitted$optimizer_state <- optimizer_state(record$iteration)
        if (!is.null(progress)) progress(emitted)
        invisible(emitted)
    }
    observed_objective_noise <- function() {
        if (!length(history)) return(0)
        microscopic <- 1e-5 * max(1, max(abs(par)))
        parameter_tolerance <- sqrt(.Machine$double.eps) *
            max(1, max(abs(par)))
        degradation <- vapply(history, function(entry) {
            if (!identical(entry$event, 'rejected') ||
                    !is.numeric(entry$parameters) ||
                    length(entry$parameters) != length(par) ||
                    any(!is.finite(entry$parameters)) ||
                    max(abs(entry$parameters - par)) > parameter_tolerance ||
                    !is.numeric(entry$step_max) ||
                    !is.finite(entry$step_max) ||
                    entry$step_max > microscopic ||
                    !is.numeric(entry$actual_improvement) ||
                    !is.finite(entry$actual_improvement) ||
                    entry$actual_improvement >= 0) return(NA_real_)
            -entry$actual_improvement
        }, numeric(1))
        degradation <- utils::tail(degradation[is.finite(degradation)], 5L)
        if (length(degradation) < 3L) return(0)
        stats::median(degradation)
    }
    assess_stagnation <- function(record, reason) {
        assessment <- if (is.null(convergence_assessment)) NULL else {
            tryCatch(
                convergence_assessment(
                    par,
                    value,
                    gradient,
                    objective_noise=observed_objective_noise()
                ),
                error=function(error) list(
                    converged=FALSE,
                    message=paste(
                        'curvature convergence assessment failed:',
                        conditionMessage(error)
                    )
                )
            )
        }
        if (is.list(assessment) && isTRUE(assessment$converged)) {
            certified <- record
            certified$event <- 'convergence_certified'
            certified$step_type <- 'curvature_certification'
            certified$assessment <- assessment$diagnostics
            certified$stagnation_reason <- reason
            emit_progress(certified)
            return(list(
                action='converged',
                message=if (!is.null(assessment$message)) {
                    assessment$message
                } else {
                    'curvature numerical-floor convergence reached'
                }
            ))
        }
        can_recover <- recovery_resets < maximum_recovery_resets
        replacement <- if (is.list(assessment)) {
            assessment$restart_hessian
        } else {
            NULL
        }
        if (can_recover && is.matrix(replacement) &&
                identical(dim(replacement), c(dimension, dimension)) &&
                all(is.finite(replacement))) {
            hessian <<- (replacement + t(replacement)) / 2
            replacement_radius <- assessment$restart_radius
            if (!is.numeric(replacement_radius) ||
                    length(replacement_radius) != 1L ||
                    !is.finite(replacement_radius) ||
                    replacement_radius <= minimum_radius) {
                replacement_radius <- min(initial_radius, 0.1)
            }
            radius <<- min(maximum_radius, replacement_radius)
            consecutive_rejections <<- 0L
            consecutive_small <<- 0L
            recovery_resets <<- recovery_resets + 1L
            curvature_resets <<- curvature_resets + 1L
            recovery <- record
            recovery$event <- 'curvature_recovery'
            recovery$trust_radius <- radius
            recovery$consecutive_rejections <- 0L
            recovery$consecutive_small_steps <- 0L
            recovery$curvature_resets <- curvature_resets
            recovery$curvature_reset <- TRUE
            recovery$step_type <- 'curvature_hessian_reset'
            recovery$recovery_resets <- recovery_resets
            recovery$assessment <- assessment$diagnostics
            recovery$stagnation_reason <- reason
            emit_progress(recovery)
            return(list(action='recovered'))
        }
        detail <- if (is.list(assessment) && !is.null(assessment$message)) {
            assessment$message
        } else {
            'no analytic recovery was available'
        }
        list(action='failed', message=paste0(reason, '; ', detail))
    }
    projected <- projected_gradient(par, gradient)
    initial_record <- list(
        event=if (resumed) 'resumed' else 'initial',
        iteration=completed_iteration,
        parameters=par,
        criterion=value,
        candidate_criterion=value,
        gradient_max=max(abs(gradient)),
        projected_gradient_max=max(abs(projected)),
        gradient_tolerance=gradient_tolerance,
        gradient_ratio=max(abs(projected)) / gradient_tolerance,
        trust_radius=radius,
        trust_radius_before=radius,
        step_norm=0,
        step_max=0,
        predicted_improvement=NA_real_,
        actual_improvement=NA_real_,
        acceptance_ratio=NA_real_,
        accepted=TRUE,
        accepted_steps=accepted_steps,
        rejected_steps=rejected,
        consecutive_rejections=consecutive_rejections,
        curvature_resets=curvature_resets,
        curvature_reset=FALSE,
        step_type=if (resumed) 'resumed' else 'initial',
        consecutive_small_steps=consecutive_small,
        recovery_resets=recovery_resets,
        function_evaluations=function_count,
        gradient_evaluations=gradient_count
    )
    emit_progress(initial_record)
    convergence <- 1L
    message <- 'iteration limit reached'
    iterations <- completed_iteration
    continue_optimization <- TRUE
    if (resumed && (radius < minimum_radius || consecutive_small >= 3L)) {
        reason <- if (radius < minimum_radius) {
            'trust radius was below its minimum on resume'
        } else {
            'accepted steps had ceased making numerically resolvable progress'
        }
        outcome <- assess_stagnation(initial_record, reason)
        if (identical(outcome$action, 'converged')) {
            convergence <- 0L
            message <- outcome$message
            continue_optimization <- FALSE
        } else if (identical(outcome$action, 'failed')) {
            message <- outcome$message
            continue_optimization <- FALSE
        }
    }
    iteration_sequence <- if (continue_optimization &&
            completed_iteration < maxit) {
        seq.int(completed_iteration + 1L, maxit)
    } else {
        integer()
    }
    for (iteration in iteration_sequence) {
        projected <- projected_gradient(par, gradient)
        projected_norm <- max(abs(projected))
        if (projected_norm <= gradient_tolerance) {
            convergence <- 0L
            message <- 'projected gradient tolerance reached'
            break
        }
        iterations <- iteration
        step_type <- 'bfgs'
        direction <- tryCatch(
            -solve(hessian, projected),
            error=function(error) rep.int(NA_real_, dimension)
        )
        directional_derivative <- sum(projected * direction)
        descent_tolerance <- .Machine$double.eps *
            max(1, sqrt(sum(projected^2)) * sqrt(sum(direction^2)))
        if (any(!is.finite(direction)) ||
                !is.finite(directional_derivative) ||
                directional_derivative >= -descent_tolerance) {
            scale <- max(1, max(abs(projected)))
            hessian <- diag(scale, dimension)
            direction <- -projected / scale
            curvature_resets <- curvature_resets + 1L
            step_type <- 'steepest_reset'
        }
        direction_norm <- sqrt(sum(direction^2))
        if (direction_norm > radius) {
            direction <- direction * radius / direction_norm
        }
        candidate <- pmin(upper, pmax(lower, par + direction))
        step <- candidate - par
        model_reduction <- function(candidate_step) {
            reduction <- -sum(projected * candidate_step) -
                drop(crossprod(
                    candidate_step,
                    hessian %*% candidate_step
                )) / 2
            if (!is.finite(reduction)) -Inf else reduction
        }
        predicted <- model_reduction(step)
        gradient_norm <- sqrt(sum(projected^2))
        if (gradient_norm > 0) {
            gradient_curvature <- drop(crossprod(
                projected,
                hessian %*% projected
            ))
            cauchy_scale <- radius / gradient_norm
            if (is.finite(gradient_curvature) && gradient_curvature > 0) {
                cauchy_scale <- min(
                    cauchy_scale,
                    sum(projected^2) / gradient_curvature
                )
            }
            cauchy_candidate <- pmin(
                upper,
                pmax(lower, par - cauchy_scale * projected)
            )
            cauchy_step <- cauchy_candidate - par
            cauchy_predicted <- model_reduction(cauchy_step)
            if (cauchy_predicted > predicted) {
                candidate <- cauchy_candidate
                step <- cauchy_step
                predicted <- cauchy_predicted
                step_type <- 'cauchy'
            }
        }
        if (max(abs(step)) <= 1e-12) {
            convergence <- 0L
            message <- 'bounded step tolerance reached'
            break
        }
        radius_before <- radius
        if (!is.finite(predicted) || predicted <= 0) predicted <- -sum(
            gradient * step
        )
        candidate_value <- fn(candidate)
        function_count <- function_count + 1L
        actual <- value - candidate_value
        ratio <- if (is.finite(candidate_value) && predicted > 0) {
            actual / predicted
        } else {
            -Inf
        }
        accepted <- is.finite(ratio) && ratio > 1e-4 && actual > 0
        candidate_gradient <- NULL
        objective_floor <- max(1e-10, 1e-10 * (1 + abs(value)))
        noise_limited <- !accepted && is.finite(candidate_value) &&
            predicted > 0 && predicted <= objective_floor &&
            actual >= -objective_floor
        if (noise_limited) {
            candidate_gradient <- gr(candidate)
            gradient_count <- gradient_count + 1L
            candidate_projected <- projected_gradient(
                candidate,
                candidate_gradient
            )
            candidate_norm <- max(abs(candidate_projected))
            accepted <- candidate_norm <= gradient_tolerance ||
                candidate_norm < projected_norm -
                    sqrt(.Machine$double.eps) * max(1, projected_norm)
            if (accepted) step_type <- paste0(step_type, '_noise_limited')
        }
        if (!accepted) {
            rejected <- rejected + 1L
            consecutive_rejections <- consecutive_rejections + 1L
            radius <- radius / 4
            curvature_reset <- consecutive_rejections >= 3L
            if (curvature_reset) {
                scale <- max(1, max(abs(projected)))
                hessian <- diag(scale, dimension)
                curvature_resets <- curvature_resets + 1L
                consecutive_rejections <- 0L
            }
            record <- list(
                event='rejected',
                iteration=iteration,
                parameters=par,
                candidate_parameters=candidate,
                criterion=value,
                candidate_criterion=candidate_value,
                gradient_max=max(abs(gradient)),
                projected_gradient_max=projected_norm,
                gradient_tolerance=gradient_tolerance,
                gradient_ratio=projected_norm / gradient_tolerance,
                trust_radius=radius,
                trust_radius_before=radius_before,
                step_norm=sqrt(sum(step^2)),
                step_max=max(abs(step)),
                predicted_improvement=predicted,
                actual_improvement=actual,
                acceptance_ratio=ratio,
                accepted=FALSE,
                accepted_steps=accepted_steps,
                rejected_steps=rejected,
                consecutive_rejections=consecutive_rejections,
                curvature_resets=curvature_resets,
                curvature_reset=curvature_reset,
                step_type=step_type,
                consecutive_small_steps=consecutive_small,
                recovery_resets=recovery_resets,
                function_evaluations=function_count,
                gradient_evaluations=gradient_count
            )
            emit_progress(record)
            if (radius < minimum_radius) {
                outcome <- assess_stagnation(
                    record,
                    'trust radius became too small'
                )
                if (identical(outcome$action, 'converged')) {
                    convergence <- 0L
                    message <- outcome$message
                    break
                }
                if (identical(outcome$action, 'recovered')) {
                    next
                }
                message <- outcome$message
                break
            }
            next
        }
        if (is.null(candidate_gradient)) {
            candidate_gradient <- gr(candidate)
            gradient_count <- gradient_count + 1L
        }
        difference <- candidate_gradient - gradient
        hessian_step <- drop(crossprod(step, hessian %*% step))
        curvature <- sum(step * difference)
        if (is.finite(hessian_step) && hessian_step > 1e-14) {
            if (!is.finite(curvature) || curvature < 0.2 * hessian_step) {
                denominator <- hessian_step - curvature
                theta <- if (is.finite(denominator) && denominator > 0) {
                    0.8 * hessian_step / denominator
                } else {
                    0
                }
                difference <- theta * difference +
                    (1 - theta) * drop(hessian %*% step)
                curvature <- sum(step * difference)
            }
            if (is.finite(curvature) && curvature > 1e-14) {
                hessian_times_step <- drop(hessian %*% step)
                hessian <- hessian -
                    tcrossprod(hessian_times_step) / hessian_step +
                    tcrossprod(difference) / curvature
                hessian <- (hessian + t(hessian)) / 2
            }
        }
        par <- candidate
        value <- candidate_value
        gradient <- candidate_gradient
        if (ratio < 0.25) {
            radius <- radius / 2
        } else if (ratio > 0.75 &&
                sqrt(sum(step^2)) >= 0.8 * radius_before) {
            radius <- min(maximum_radius, 2 * radius)
        }
        accepted_steps <- accepted_steps + 1L
        consecutive_rejections <- 0L
        small_step <- max(abs(step)) <= 1e-3 &&
            abs(actual) <= 1e-8 * (1 + abs(value))
        consecutive_small <- if (small_step) consecutive_small + 1L else 0L
        candidate_projected <- projected_gradient(par, gradient)
        record <- list(
            event='accepted',
            iteration=iteration,
            parameters=par,
            criterion=value,
            candidate_criterion=value,
            gradient_max=max(abs(gradient)),
            projected_gradient_max=max(abs(candidate_projected)),
            previous_projected_gradient_max=projected_norm,
            gradient_tolerance=gradient_tolerance,
            gradient_ratio=max(abs(candidate_projected)) / gradient_tolerance,
            trust_radius=radius,
            trust_radius_before=radius_before,
            step_norm=sqrt(sum(step^2)),
            step_max=max(abs(step)),
            predicted_improvement=predicted,
            actual_improvement=actual,
            acceptance_ratio=ratio,
            accepted=TRUE,
            accepted_steps=accepted_steps,
            rejected_steps=rejected,
            consecutive_rejections=consecutive_rejections,
            curvature_resets=curvature_resets,
            curvature_reset=FALSE,
            step_type=step_type,
            consecutive_small_steps=consecutive_small,
            recovery_resets=recovery_resets,
            function_evaluations=function_count,
            gradient_evaluations=gradient_count
        )
        emit_progress(record)
        if (radius < minimum_radius || consecutive_small >= 3L) {
            reason <- if (radius < minimum_radius) {
                'trust radius became too small after an accepted step'
            } else {
                'accepted steps ceased making numerically resolvable progress'
            }
            outcome <- assess_stagnation(record, reason)
            if (identical(outcome$action, 'converged')) {
                convergence <- 0L
                message <- outcome$message
                break
            }
            if (identical(outcome$action, 'recovered')) next
            message <- outcome$message
            break
        }
    }
    final_projected <- projected_gradient(par, gradient)
    emit_progress(list(
        event='finished',
        iteration=iterations,
        parameters=par,
        criterion=value,
        candidate_criterion=value,
        gradient_max=max(abs(gradient)),
        projected_gradient_max=max(abs(final_projected)),
        gradient_tolerance=gradient_tolerance,
        gradient_ratio=max(abs(final_projected)) / gradient_tolerance,
        trust_radius=radius,
        trust_radius_before=radius,
        step_norm=0,
        step_max=0,
        predicted_improvement=NA_real_,
        actual_improvement=NA_real_,
        acceptance_ratio=NA_real_,
        accepted=NA,
        accepted_steps=accepted_steps,
        rejected_steps=rejected,
        consecutive_rejections=consecutive_rejections,
        curvature_resets=curvature_resets,
        curvature_reset=FALSE,
        step_type='finished',
        consecutive_small_steps=consecutive_small,
        recovery_resets=recovery_resets,
        function_evaluations=function_count,
        gradient_evaluations=gradient_count,
        convergence=convergence,
        message=message
    ), save=FALSE)
    list(
        par=par,
        value=value,
        counts=c(`function`=function_count, gradient=gradient_count),
        convergence=convergence,
        message=message,
        gradient=gradient,
        hessian=hessian,
        iterations=iterations,
        rejected_steps=rejected,
        curvature_resets=curvature_resets,
        recovery_resets=recovery_resets,
        optimizer_state=optimizer_state(iterations),
        history=history
    )
}
