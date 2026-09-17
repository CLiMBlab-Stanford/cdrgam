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
        progress=NULL
) {
    par <- pmin(upper, pmax(lower, as.numeric(par)))
    dimension <- length(par)
    function_count <- 1L
    gradient_count <- 1L
    value <- fn(par)
    gradient <- gr(par)
    scale <- max(1, max(abs(gradient)))
    hessian <- diag(scale, dimension)
    radius <- initial_radius
    projected_gradient <- function(parameters, score) {
        out <- score
        at_lower <- parameters <= lower + 1e-10
        at_upper <- parameters >= upper - 1e-10
        out[at_lower & out > 0] <- 0
        out[at_upper & out < 0] <- 0
        out
    }
    emit_progress <- function(record) {
        if (!is.null(progress)) progress(record)
        invisible(record)
    }
    accepted_steps <- 0L
    rejected <- 0L
    consecutive_rejections <- 0L
    curvature_resets <- 0L
    consecutive_small <- 0L
    projected <- projected_gradient(par, gradient)
    initial_record <- list(
        event='initial',
        iteration=0L,
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
        accepted_steps=0L,
        rejected_steps=0L,
        consecutive_rejections=0L,
        curvature_resets=0L,
        curvature_reset=FALSE,
        step_type='initial',
        consecutive_small_steps=0L,
        function_evaluations=function_count,
        gradient_evaluations=gradient_count
    )
    history <- list(initial_record)
    emit_progress(initial_record)
    convergence <- 1L
    message <- 'iteration limit reached'
    iterations <- 0L
    for (iteration in seq_len(maxit)) {
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
            reduction <- -sum(gradient * candidate_step) -
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
                function_evaluations=function_count,
                gradient_evaluations=gradient_count
            )
            history[[length(history) + 1L]] <- record
            emit_progress(record)
            if (radius < 1e-8) {
                message <- 'trust radius became too small'
                break
            }
            next
        }
        candidate_gradient <- gr(candidate)
        gradient_count <- gradient_count + 1L
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
            function_evaluations=function_count,
            gradient_evaluations=gradient_count
        )
        history[[length(history) + 1L]] <- record
        emit_progress(record)
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
        function_evaluations=function_count,
        gradient_evaluations=gradient_count,
        convergence=convergence,
        message=message
    ))
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
        history=history
    )
}
