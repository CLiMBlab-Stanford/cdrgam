.cdrgam_distributional_sparse_setup <- function(design, sparse_control, ...) {
    .fit_sparse_gaussian(
        design=design,
        family=stats::gaussian(),
        method='REML',
        trace=FALSE,
        sparse_control=sparse_control[intersect(
            names(sparse_control), c('crossprod_chunk_size', 'supernodal')
        )],
        setup_only=TRUE,
        ...
    )
}

.cdrgam_distributional_component_signature <- function(component) {
    if (inherits(component, 'cdrgam_sparse_grouped_component')) {
        return(list(
            grouped=TRUE,
            dimension=c(
                .sparse_component_nrow(component),
                .sparse_component_ncol(component)
            ),
            values=.solver_numeric_signature(component$base),
            groups=.solver_numeric_signature(component$group_index)
        ))
    }
    values <- if (methods::is(component, 'sparseMatrix')) {
        methods::slot(component, 'x')
    } else component
    list(
        grouped=FALSE,
        dimension=dim(component),
        nonzeros=as.numeric(Matrix::nnzero(component)),
        values=.solver_numeric_signature(values)
    )
}

.cdrgam_distributional_checkpoint_signature <- function(
        design, assemblies, family, control
) {
    parameter_signature <- Map(function(parameter_design, assembly) {
        list(
            formula=paste(deparse(parameter_design$effective_formula),
                collapse=''),
            observations=as.integer(assembly$observation_count),
            dimension=as.integer(assembly$dimension),
            coefficient_names=as.character(assembly$coefficient_names),
            smoothing_parameters=as.character(assembly$sp_names),
            response=.solver_numeric_signature(assembly$setup$y),
            weights=.solver_numeric_signature(assembly$setup$w),
            offset=.solver_numeric_signature(assembly$setup$offset),
            design=lapply(
                assembly$matrices,
                .cdrgam_distributional_component_signature
            ),
            penalties=lapply(assembly$penalty_components, function(penalty) {
                list(
                    dimension=dim(penalty),
                    nonzeros=as.numeric(Matrix::nnzero(penalty)),
                    values=.solver_numeric_signature(methods::slot(
                        penalty, 'x'
                    ))
                )
            })
        )
    }, design$parameters, assemblies)
    names(parameter_signature) <- names(assemblies)
    family_links <- if (is.null(family$linfo)) {
        family$link
    } else vapply(family$linfo, function(link) {
        if (!is.null(link$link)) link$link else link$name
    }, character(1))
    optimizer_signature <- list(
        gradient=control$gradient_method,
        gradient_probes=control$gradient_probes,
        optimizer_maxit=control$optimizer_maxit,
        optimizer_gradient_tolerance=
            control$optimizer_gradient_tolerance,
        inner_tolerance=control$inner_tolerance,
        inner_relative_tolerance=control$inner_relative_tolerance,
        inner_maxit=control$inner_maxit,
        lower=-25,
        upper=25
    )
    # Keep the historical default signature stable so parameter-only
    # checkpoints written before trust refinement can still be resumed.
    if (!identical(as.numeric(control$optimizer_trust_radius), 2)) {
        optimizer_signature$optimizer_trust_radius <-
            control$optimizer_trust_radius
    }
    default_adaptive_maxit <- if (identical(control$inner_maxit, 200L)) {
        600L
    } else control$inner_maxit
    if (!identical(control$inner_adaptive_maxit, default_adaptive_maxit)) {
        optimizer_signature$inner_adaptive_maxit <-
            control$inner_adaptive_maxit
    }
    if (identical(control$inner_algorithm, 'schur')) {
        optimizer_signature$inner_algorithm <- 'schur'
        optimizer_signature$inner_schur_rank <- control$inner_schur_rank
    }
    if (identical(control$criterion, 'QNCV')) {
        optimizer_signature$criterion <- 'QNCV'
        optimizer_signature$gamma <- control$gamma
        optimizer_signature$qncv_batch_size <- control$qncv_batch_size
    }
    list(
        version=1L,
        backend='distributional-sparse',
        family=family$family,
        links=family_links,
        parameters=parameter_signature,
        optimizer=optimizer_signature
    )
}

.cdrgam_sparse_penalty <- function(assembly, sp) {
    penalty <- Matrix::Matrix(
        0, assembly$dimension, assembly$dimension, sparse=TRUE
    )
    for (i in seq_along(sp)) {
        penalty <- penalty + sp[[i]] * assembly$penalty_components[[i]]
    }
    Matrix::forceSymmetric(penalty, uplo='U')
}

.cdrgam_sparse_factor <- function(system, supernodal, previous=NULL) {
    if (!is.null(previous)) {
        updated <- tryCatch(
            Matrix::update(previous, system),
            error=function(error) NULL
        )
        if (!is.null(updated)) return(updated)
    }
    Matrix::Cholesky(
        system, perm=TRUE, LDL=FALSE, super=isTRUE(supernodal)
    )
}

.cdrgam_gaulss_step_factor <- function(system, supernodal, previous=NULL) {
    attempt <- function(candidate, prior=NULL) tryCatch(
        withCallingHandlers(
            .cdrgam_sparse_factor(candidate, supernodal, prior),
            warning=function(warning) stop(conditionMessage(warning))
        ),
        error=function(error) NULL
    )
    factor <- attempt(system, previous)
    if (is.null(factor) && !is.null(previous)) factor <- attempt(system)
    if (!is.null(factor)) return(list(factor=factor, damping=0))
    diagonal_scale <- max(abs(Matrix::diag(system)), 1)
    for (relative in c(
            1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1, 10
        )) {
        damping <- relative * diagonal_scale
        adjusted <- Matrix::forceSymmetric(
            system + Matrix::Diagonal(nrow(system), x=damping),
            uplo='U'
        )
        factor <- attempt(adjusted)
        if (!is.null(factor)) {
            return(list(factor=factor, damping=damping))
        }
    }
    stop('Could not factor a damped gaulss Fisher-scoring system')
}

.cdrgam_distributional_sparse_predictor <- function(
        assembly, coefficients, offset
) {
    .cdrgam_sparse_linear_predictor(
        assembly,
        coefficients,
        offset,
        assembly$crossprod_chunk_size
    )
}

.cdrgam_distributional_sparse_predictor_batch <- function(
        assembly, coefficients
) {
    coefficients <- as.matrix(coefficients)
    output <- matrix(
        0, nrow=assembly$observation_count, ncol=ncol(coefficients)
    )
    chunk_size <- min(
        assembly$crossprod_chunk_size, assembly$observation_count
    )
    for (start in seq.int(1L, assembly$observation_count, by=chunk_size)) {
        rows <- start:min(
            assembly$observation_count, start + chunk_size - 1L
        )
        output[rows, ] <- as.matrix(
            .cdrgam_sparse_design_chunk(assembly, rows) %*% coefficients
        )
    }
    output
}

.cdrgam_distributional_cache_designs <- function(assemblies) {
    observation_count <- assemblies$location$observation_count
    sample_count <- min(
        observation_count,
        assemblies$location$crossprod_chunk_size,
        assemblies$scale$crossprod_chunk_size
    )
    sample_rows <- seq_len(sample_count)
    sample <- lapply(assemblies, function(assembly) {
        .cdrgam_sparse_design_chunk(assembly, sample_rows)
    })
    projected_bytes <- sum(vapply(sample, utils::object.size, numeric(1))) *
        observation_count / sample_count
    memory <- .cdrgam_memory_availability()
    enabled <- is.finite(memory$available_bytes) &&
        1.25 * projected_bytes <= 0.1 * memory$available_bytes
    if (enabled) {
        assemblies <- lapply(assemblies, function(assembly) {
            assembly$cached_design <- .cdrgam_sparse_design_chunk(
                assembly, seq_len(assembly$observation_count)
            )
            assembly
        })
    }
    bytes <- if (enabled) sum(vapply(assemblies, function(assembly) {
        utils::object.size(assembly$cached_design)
    }, numeric(1))) else 0
    list(
        assemblies=assemblies,
        enabled=enabled,
        bytes=bytes,
        projected_bytes=projected_bytes,
        memory_source=memory$source,
        memory_available_bytes=memory$available_bytes
    )
}

.cdrgam_gaulss_sparse_moments <- function(
        assemblies, coefficients, sigma, q, residual, weights,
        include_observed=FALSE, include_cross=include_observed,
        include_information=TRUE,
        workers=1L, blas_threads=1L,
        parameters=c('location', 'scale')
) {
    parameters <- match.arg(
        parameters, c('location', 'scale'), several.ok=TRUE
    )
    if (isTRUE(include_observed) && length(parameters) != 2L) {
        stop('Observed gaulss moments require location and scale')
    }
    if (isTRUE(include_observed) && !isTRUE(include_information)) {
        stop('Observed gaulss moments require information matrices')
    }
    use_location <- 'location' %in% parameters
    use_scale <- 'scale' %in% parameters
    location <- assemblies$location
    scale <- assemblies$scale
    chunk_sizes <- c(
        if (use_location) location$crossprod_chunk_size,
        if (use_scale) scale$crossprod_chunk_size
    )
    chunk_size <- min(chunk_sizes)
    starts <- seq.int(1L, location$observation_count, by=chunk_size)
    workers <- min(
        .cdrgam_positive_integer(workers, 'gaulss moment worker count'),
        length(starts)
    )
    evaluate <- function(local_batches) {
        location_gradient <- if (use_location) numeric(location$dimension)
        scale_gradient <- if (use_scale) numeric(scale$dimension)
        location_information <- NULL
        scale_information <- NULL
        cross_information <- NULL
        for (rows in local_batches) {
        X_location <- if (use_location) {
            .cdrgam_sparse_design_chunk(location, rows)
        }
        X_scale <- if (use_scale) {
            .cdrgam_sparse_design_chunk(scale, rows)
        }
        sigma_rows <- sigma[rows]
        q_rows <- q[rows]
        residual_rows <- residual[rows]
        weight_rows <- weights[rows]
        if (use_location) {
            location_score <- weight_rows * (-residual_rows / sigma_rows^2)
            location_gradient <- location_gradient + as.numeric(
                Matrix::crossprod(X_location, location_score)
            )
            if (isTRUE(include_information)) {
                weighted_location <- Matrix::Diagonal(
                    x=sqrt(weight_rows / sigma_rows^2)
                ) %*% X_location
                location_part <- Matrix::crossprod(weighted_location)
                location_information <- if (is.null(location_information)) {
                    location_part
                } else location_information + location_part
            }
        }
        if (use_scale) {
            scale_score <- weight_rows * q_rows * (
                1 / sigma_rows - residual_rows^2 / sigma_rows^3
            )
            scale_gradient <- scale_gradient + as.numeric(
                Matrix::crossprod(X_scale, scale_score)
            )
            if (isTRUE(include_information)) {
                scale_part <- if (isTRUE(include_observed)) {
                    scale_weight <- weight_rows * (
                        q_rows * (
                            1 / sigma_rows - residual_rows^2 / sigma_rows^3
                        ) +
                            q_rows^2 * (
                                -1 / sigma_rows^2 +
                                    3 * residual_rows^2 / sigma_rows^4
                            )
                    )
                    Matrix::crossprod(
                        X_scale,
                        Matrix::Diagonal(x=scale_weight) %*% X_scale
                    )
                } else {
                    weighted_scale <- Matrix::Diagonal(
                        x=sqrt(weight_rows * 2 * q_rows^2 / sigma_rows^2)
                    ) %*% X_scale
                    Matrix::crossprod(weighted_scale)
                }
                scale_information <- if (is.null(scale_information)) {
                    scale_part
                } else scale_information + scale_part
            }
        }
        if (isTRUE(include_observed)) {
            cross_weight <- weight_rows * 2 * q_rows * residual_rows /
                sigma_rows^3
            cross_part <- if (isTRUE(include_cross)) Matrix::crossprod(
                X_location, Matrix::Diagonal(x=cross_weight) %*% X_scale
            ) else NULL
            if (isTRUE(include_cross)) {
                cross_information <- if (is.null(cross_information)) {
                    cross_part
                } else cross_information + cross_part
            }
        }
        }
        list(
            gradients=list(
                location=location_gradient,
                scale=scale_gradient
            ),
            information=list(
                location=location_information,
                scale=scale_information,
                cross=cross_information
            )
        )
    }
    cached <- (!use_location || !is.null(location$cached_design)) &&
        (!use_scale || !is.null(scale$cached_design))
    groups <- if (cached) {
        lapply(seq_len(workers), function(worker) {
            first <- floor(
                (worker - 1L) * location$observation_count / workers
            ) + 1L
            last <- floor(worker * location$observation_count / workers)
            list(seq.int(first, last))
        })
    } else {
        lapply(
            split(starts, rep(seq_len(workers), length.out=length(starts))),
            function(local_starts) lapply(local_starts, function(start) {
                start:min(
                    location$observation_count, start + chunk_size - 1L
                )
            })
        )
    }
    parts <- if (workers > 1L && .Platform$OS.type != 'windows') {
        .cdrgam_with_blas_threads(blas_threads, suppressWarnings(
            parallel::mclapply(
                groups,
                evaluate,
                mc.cores=workers,
                mc.preschedule=TRUE,
                mc.set.seed=FALSE
            )
        ))
    } else lapply(groups, evaluate)
    valid <- length(parts) == length(groups) && all(vapply(
        parts,
        function(part) is.list(part) &&
            is.list(part$gradients) && is.list(part$information),
        logical(1)
    ))
    if (!valid) stop('A parallel gaulss moment batch failed')
    add <- function(field, parameter) {
        values <- lapply(parts, function(part) part[[field]][[parameter]])
        values <- Filter(Negate(is.null), values)
        if (!length(values)) return(NULL)
        Reduce(`+`, values)
    }
    list(
        gradients=list(
            location=add('gradients', 'location'),
            scale=add('gradients', 'scale')
        ),
        information=list(
            location=add('information', 'location'),
            scale=add('information', 'scale'),
            cross=if (isTRUE(include_observed) && isTRUE(include_cross)) {
                add('information', 'cross')
            } else NULL
        ),
        chunks=length(starts),
        batches=sum(lengths(groups)),
        workers=workers
    )
}

.cdrgam_gaulss_numerically_converged <- function(
        directional_derivative, best_change, objective,
        minimum_step=2^-20
) {
    values <- c(
        directional_derivative, best_change, objective, minimum_step
    )
    if (any(!is.finite(values)) || directional_derivative >= 0 ||
            minimum_step <= 0) return(FALSE)
    numerical_floor <- 1000 * .Machine$double.eps * (1 + abs(objective))
    minimum_step * -directional_derivative <= numerical_floor &&
        best_change <= numerical_floor
}

.cdrgam_gaulss_relative_convergence <- function(
        previous, current, directional_derivative, tolerance=1e-8
) {
    values <- c(previous, current, directional_derivative, tolerance)
    if (any(!is.finite(values)) || tolerance <= 0 || current > previous ||
            directional_derivative >= 0) {
        return(list(
            converged=FALSE,
            actual=Inf,
            predicted=Inf
        ))
    }
    scale <- 1 + abs(previous)
    actual <- (previous - current) / scale
    predicted <- -directional_derivative / scale
    list(
        converged=actual <= tolerance && predicted <= tolerance,
        actual=actual,
        predicted=predicted
    )
}

.cdrgam_anderson_candidate <- function(
        states, residuals, block_sizes, depth=2L, ridge=1e-10
) {
    count <- min(length(states), length(residuals))
    depth <- min(as.integer(depth), count - 1L)
    if (depth < 1L || length(block_sizes) < 1L ||
            sum(block_sizes) != length(states[[count]])) return(NULL)
    indices <- seq.int(count - depth, count)
    local_states <- states[indices]
    local_residuals <- residuals[indices]
    if (any(vapply(
            c(local_states, local_residuals),
            function(value) length(value) != sum(block_sizes) ||
                any(!is.finite(value)),
            logical(1)
        ))) return(NULL)
    state_differences <- vapply(seq_len(depth), function(index) {
        local_states[[index + 1L]] - local_states[[index]]
    }, numeric(sum(block_sizes)))
    residual_differences <- vapply(seq_len(depth), function(index) {
        local_residuals[[index + 1L]] - local_residuals[[index]]
    }, numeric(sum(block_sizes)))
    current_state <- local_states[[depth + 1L]]
    current_residual <- local_residuals[[depth + 1L]]
    block <- rep(seq_along(block_sizes), block_sizes)
    scales <- vapply(seq_along(block_sizes), function(index) {
        rows <- block == index
        residual_scale <- sqrt(mean(current_residual[rows]^2))
        state_scale <- sqrt(mean(current_state[rows]^2))
        max(residual_scale, sqrt(.Machine$double.eps) * (1 + state_scale))
    }, numeric(1))
    scaled_differences <- residual_differences / scales[block]
    scaled_residual <- current_residual / scales[block]
    gram <- crossprod(scaled_differences)
    stabilization <- max(ridge * max(1, sum(diag(gram))), .Machine$double.eps)
    coefficients <- tryCatch(
        solve(
            gram + diag(stabilization, nrow=depth),
            as.numeric(crossprod(scaled_differences, scaled_residual))
        ),
        error=function(error) NULL
    )
    if (is.null(coefficients) || any(!is.finite(coefficients)) ||
            max(abs(coefficients)) > 10) return(NULL)
    candidate <- current_state + current_residual -
        as.numeric((state_differences + residual_differences) %*% coefficients)
    if (any(!is.finite(candidate))) NULL else candidate
}

.cdrgam_gaulss_block_update <- function(
        parameter, assembly, coefficients, predictors, gradient, information,
        penalty, factor, objective, value
) {
    system <- Matrix::forceSymmetric(information + penalty, uplo='U')
    factorization <- .cdrgam_gaulss_step_factor(
        system, assembly$supernodal, factor
    )
    step <- -as.numeric(.cdr_factor_solve(factorization$factor, gradient))
    predictor_step <- .cdrgam_distributional_sparse_predictor(
        assembly, step, numeric(assembly$observation_count)
    )
    directional_derivative <- sum(gradient * step)
    step_size <- 1
    minimum_step <- 2^-20
    best_change <- Inf
    while (step_size >= minimum_step) {
        candidate_coefficients <- coefficients
        candidate_predictors <- predictors
        candidate_coefficients[[parameter]] <-
            coefficients[[parameter]] + step_size * step
        candidate_predictors[[parameter]] <-
            predictors[[parameter]] + step_size * predictor_step
        candidate_value <- objective(
            candidate_coefficients, candidate_predictors
        )
        best_change <- min(best_change, candidate_value - value, na.rm=TRUE)
        if (is.finite(candidate_value) && candidate_value < value) {
            return(list(
                coefficients=candidate_coefficients,
                predictors=candidate_predictors,
                value=candidate_value,
                factor=factorization$factor,
                damping=factorization$damping,
                accepted=TRUE,
                step=step_size * step,
                directional_derivative=directional_derivative,
                best_change=best_change,
                numerically_converged=FALSE
            ))
        }
        step_size <- step_size / 2
    }
    list(
        coefficients=coefficients,
        predictors=predictors,
        value=value,
        factor=factorization$factor,
        damping=factorization$damping,
        accepted=FALSE,
        step=numeric(length(step)),
        directional_derivative=directional_derivative,
        best_change=best_change,
        numerically_converged=.cdrgam_gaulss_numerically_converged(
            directional_derivative, best_change, value, minimum_step
        )
    )
}

.cdrgam_gaulss_cross_product <- function(
        assemblies, cross_weight, vectors, workers=1L, blas_threads=1L,
        outputs=c('location', 'scale')
) {
    outputs <- match.arg(
        outputs, c('location', 'scale'), several.ok=TRUE
    )
    observation_count <- assemblies$location$observation_count
    matrix_input <- list(
        location=is.matrix(vectors$location),
        scale=is.matrix(vectors$scale)
    )
    predictions <- list(
        location=if ('scale' %in% outputs) {
            .cdrgam_distributional_sparse_predictor_batch(
                assemblies$location, vectors$location
            )
        } else NULL,
        scale=if ('location' %in% outputs) {
            .cdrgam_distributional_sparse_predictor_batch(
                assemblies$scale, vectors$scale
            )
        } else NULL
    )
    if (length(outputs) > 1L &&
            ncol(predictions$location) != ncol(predictions$scale)) {
        stop('Joint gaulss cross-products require equal batch widths')
    }
    chunk_size <- min(
        assemblies$location$crossprod_chunk_size,
        assemblies$scale$crossprod_chunk_size
    )
    starts <- seq.int(1L, observation_count, by=chunk_size)
    workers <- min(
        .cdrgam_positive_integer(workers, 'gaulss cross-product worker count'),
        length(starts)
    )
    groups <- split(starts, rep(seq_len(workers), length.out=length(starts)))
    evaluate <- function(local_starts) {
        location <- if ('location' %in% outputs) {
            matrix(
                0, nrow=assemblies$location$dimension,
                ncol=ncol(predictions$scale)
            )
        } else NULL
        scale <- if ('scale' %in% outputs) {
            matrix(
                0, nrow=assemblies$scale$dimension,
                ncol=ncol(predictions$location)
            )
        } else NULL
        for (start in local_starts) {
            rows <- start:min(observation_count, start + chunk_size - 1L)
            X_location <- if ('location' %in% outputs) {
                .cdrgam_sparse_design_chunk(assemblies$location, rows)
            } else NULL
            X_scale <- if ('scale' %in% outputs) {
                .cdrgam_sparse_design_chunk(assemblies$scale, rows)
            } else NULL
            local_weight <- cross_weight[rows]
            if ('location' %in% outputs) {
                location <- location + as.matrix(Matrix::crossprod(
                    X_location,
                    local_weight * predictions$scale[rows, , drop=FALSE]
                ))
            }
            if ('scale' %in% outputs) {
                scale <- scale + as.matrix(Matrix::crossprod(
                    X_scale,
                    local_weight * predictions$location[rows, , drop=FALSE]
                ))
            }
        }
        list(location=location, scale=scale)
    }
    parts <- if (workers > 1L && .Platform$OS.type != 'windows') {
        .cdrgam_with_blas_threads(blas_threads, suppressWarnings(
            parallel::mclapply(
                groups, evaluate, mc.cores=workers,
                mc.preschedule=TRUE, mc.set.seed=FALSE
            )
        ))
    } else lapply(groups, evaluate)
    list(
        location=if ('location' %in% outputs) {
            value <- Reduce(`+`, lapply(parts, `[[`, 'location'))
            if (matrix_input$scale) value else as.numeric(value)
        } else NULL,
        scale=if ('scale' %in% outputs) {
            value <- Reduce(`+`, lapply(parts, `[[`, 'scale'))
            if (matrix_input$location) value else as.numeric(value)
        } else NULL
    )
}

.cdrgam_bicgstab <- function(
        right_hand_side, product, preconditioner, tolerance=1e-6,
        maxit=50L, residual_measure=NULL, initial_solution=NULL,
        progress=NULL
) {
    notify <- function(event, ...) {
        if (is.function(progress)) progress(event, ...)
        invisible(NULL)
    }
    timed_application <- function(kind, iteration, application, operation) {
        notify(
            paste0('Krylov ', kind, ' started'),
            iteration=iteration, application=application
        )
        started <- proc.time()[['elapsed']]
        value <- operation()
        seconds <- proc.time()[['elapsed']] - started
        if (identical(kind, 'operator')) {
            operator_seconds <<- operator_seconds + seconds
            operator_applications <<- operator_applications + 1L
        } else {
            preconditioner_seconds <<- preconditioner_seconds + seconds
            preconditioner_applications <<- preconditioner_applications + 1L
        }
        notify(
            paste0('Krylov ', kind, ' complete'),
            iteration=iteration,
            application=application,
            seconds=format(seconds, digits=5),
            cumulative_seconds=format(
                if (identical(kind, 'operator')) {
                    operator_seconds
                } else preconditioner_seconds,
                digits=5
            )
        )
        value
    }
    right_hand_side <- as.numeric(right_hand_side)
    if (is.null(residual_measure)) {
        residual_measure <- function(value) sqrt(sum(value^2))
    }
    operator_seconds <- 0
    preconditioner_seconds <- 0
    operator_applications <- 0L
    preconditioner_applications <- 0L
    residual_scale <- residual_measure(right_hand_side)
    if (!is.finite(residual_scale) || residual_scale < 0) {
        stop('Krylov residual measure was not finite and nonnegative')
    }
    if (residual_scale <= .Machine$double.eps) residual_scale <- 1
    warm_started <- length(initial_solution) == length(right_hand_side) &&
        all(is.finite(initial_solution))
    solution <- if (warm_started) {
        as.numeric(initial_solution)
    } else numeric(length(right_hand_side))
    residual <- if (warm_started) {
        right_hand_side - timed_application(
            'operator', 0L, 1L, function() product(solution)
        )
    } else right_hand_side
    residual_norm <- residual_measure(residual) / residual_scale
    initial_residual <- residual_norm
    if (warm_started && (!is.finite(residual_norm) || residual_norm >= 1)) {
        notify(
            'Krylov warm start rejected',
            residual=format(residual_norm, digits=5)
        )
        solution <- numeric(length(right_hand_side))
        residual <- right_hand_side
        residual_norm <- 1
        warm_started <- FALSE
    }
    shadow <- residual
    notify(
        'Krylov solve started',
        dimension=length(right_hand_side),
        tolerance=format(tolerance, digits=5),
        maximum_iterations=maxit,
        residual=format(residual_norm, digits=5),
        warm_start=warm_started
    )
    if (residual_norm <= tolerance) {
        notify(
            'Krylov solve complete',
            iterations=0L,
            residual=format(residual_norm, digits=5),
            termination='initial residual',
            operator_seconds='0',
            preconditioner_seconds='0'
        )
        return(list(
            solution=solution, converged=TRUE, iterations=0L,
            residual=residual_norm, termination='initial residual',
            warm_started=warm_started,
            initial_residual=initial_residual,
            operator_applications=operator_applications,
            preconditioner_applications=preconditioner_applications,
            operator_seconds=operator_seconds,
            preconditioner_seconds=preconditioner_seconds
        ))
    }
    rho <- alpha <- omega <- 1
    direction <- image <- numeric(length(residual))
    termination <- 'iteration limit'
    iteration <- 0L
    for (iteration in seq_len(maxit)) {
        next_rho <- sum(shadow * residual)
        if (!is.finite(next_rho) || abs(next_rho) <= .Machine$double.eps) {
            termination <- 'rho breakdown'
            break
        }
        beta <- (next_rho / rho) * (alpha / omega)
        direction <- residual + beta * (direction - omega * image)
        preconditioned_direction <- timed_application(
            'preconditioner', iteration, 1L,
            function() preconditioner(direction)
        )
        image <- timed_application(
            'operator', iteration, 1L,
            function() product(preconditioned_direction)
        )
        denominator <- sum(shadow * image)
        if (!is.finite(denominator) ||
                abs(denominator) <= .Machine$double.eps) {
            termination <- 'alpha breakdown'
            break
        }
        alpha <- next_rho / denominator
        intermediate <- residual - alpha * image
        intermediate_norm <- residual_measure(intermediate) / residual_scale
        if (is.finite(intermediate_norm) && intermediate_norm <= tolerance) {
            solution <- solution + alpha * preconditioned_direction
            residual_norm <- intermediate_norm
            termination <- 'tolerance'
            output <- list(
                solution=solution, converged=TRUE, iterations=iteration,
                residual=residual_norm, termination=termination,
                warm_started=warm_started,
                initial_residual=initial_residual,
                operator_applications=operator_applications,
                preconditioner_applications=preconditioner_applications,
                operator_seconds=operator_seconds,
                preconditioner_seconds=preconditioner_seconds
            )
            notify(
                'Krylov solve complete',
                iterations=iteration,
                residual=format(residual_norm, digits=5),
                termination=termination,
                operator_seconds=format(operator_seconds, digits=5),
                preconditioner_seconds=format(
                    preconditioner_seconds, digits=5
                )
            )
            return(output)
        }
        preconditioned_intermediate <- timed_application(
            'preconditioner', iteration, 2L,
            function() preconditioner(intermediate)
        )
        image_intermediate <- timed_application(
            'operator', iteration, 2L,
            function() product(preconditioned_intermediate)
        )
        denominator <- sum(image_intermediate^2)
        if (!is.finite(denominator) ||
                denominator <= .Machine$double.eps) {
            termination <- 'omega breakdown'
            break
        }
        omega <- sum(image_intermediate * intermediate) / denominator
        if (!is.finite(omega) || abs(omega) <= .Machine$double.eps) {
            termination <- 'omega breakdown'
            break
        }
        solution <- solution + alpha * preconditioned_direction +
            omega * preconditioned_intermediate
        residual <- intermediate - omega * image_intermediate
        residual_norm <- residual_measure(residual) / residual_scale
        if (is.finite(residual_norm) && residual_norm <= tolerance) {
            termination <- 'tolerance'
            output <- list(
                solution=solution, converged=TRUE, iterations=iteration,
                residual=residual_norm, termination=termination,
                warm_started=warm_started,
                initial_residual=initial_residual,
                operator_applications=operator_applications,
                preconditioner_applications=preconditioner_applications,
                operator_seconds=operator_seconds,
                preconditioner_seconds=preconditioner_seconds
            )
            notify(
                'Krylov solve complete',
                iterations=iteration,
                residual=format(residual_norm, digits=5),
                termination=termination,
                operator_seconds=format(operator_seconds, digits=5),
                preconditioner_seconds=format(
                    preconditioner_seconds, digits=5
                )
            )
            return(output)
        }
        notify(
            'Krylov iteration complete',
            iteration=iteration,
            residual=format(residual_norm, digits=5),
            operator_seconds=format(operator_seconds, digits=5),
            preconditioner_seconds=format(
                preconditioner_seconds, digits=5
            )
        )
        rho <- next_rho
    }
    output <- list(
        solution=solution, converged=FALSE,
        iterations=iteration,
        residual=residual_norm, termination=termination,
        warm_started=warm_started,
        initial_residual=initial_residual,
        operator_applications=operator_applications,
        preconditioner_applications=preconditioner_applications,
        operator_seconds=operator_seconds,
        preconditioner_seconds=preconditioner_seconds
    )
    notify(
        'Krylov solve complete',
        iterations=iteration,
        residual=format(residual_norm, digits=5),
        termination=termination,
        operator_seconds=format(operator_seconds, digits=5),
        preconditioner_seconds=format(preconditioner_seconds, digits=5)
    )
    output
}

.cdrgam_fgmres <- function(
        right_hand_side, product, preconditioner, tolerance=1e-6,
        maxit=20L, residual_measure=NULL, initial_solution=NULL,
        progress=NULL
) {
    notify <- function(event, ...) {
        if (is.function(progress)) progress(event, ...)
        invisible(NULL)
    }
    timed_application <- function(kind, iteration, application, operation) {
        notify(
            paste0('Krylov ', kind, ' started'),
            iteration=iteration, application=application
        )
        started <- proc.time()[['elapsed']]
        value <- operation()
        seconds <- proc.time()[['elapsed']] - started
        if (identical(kind, 'operator')) {
            operator_seconds <<- operator_seconds + seconds
            operator_applications <<- operator_applications + 1L
        } else {
            preconditioner_seconds <<- preconditioner_seconds + seconds
            preconditioner_applications <<- preconditioner_applications + 1L
        }
        notify(
            paste0('Krylov ', kind, ' complete'),
            iteration=iteration,
            application=application,
            seconds=format(seconds, digits=5),
            cumulative_seconds=format(
                if (identical(kind, 'operator')) {
                    operator_seconds
                } else preconditioner_seconds,
                digits=5
            )
        )
        value
    }
    least_squares <- function(matrix, right_hand_side) {
        decomposition <- svd(matrix, nu=nrow(matrix), nv=ncol(matrix))
        if (!length(decomposition$d)) return(numeric(ncol(matrix)))
        threshold <- max(dim(matrix)) * max(decomposition$d) *
            .Machine$double.eps
        keep <- decomposition$d > threshold
        if (!any(keep)) return(numeric(ncol(matrix)))
        indices <- which(keep)
        projected <- as.numeric(crossprod(
            decomposition$u[, indices, drop=FALSE], right_hand_side
        )) / decomposition$d[indices]
        as.numeric(
            decomposition$v[, indices, drop=FALSE] %*% projected
        )
    }
    right_hand_side <- as.numeric(right_hand_side)
    dimension <- length(right_hand_side)
    maxit <- min(as.integer(maxit), dimension)
    if (is.null(residual_measure)) {
        residual_measure <- function(value) sqrt(sum(value^2))
    }
    operator_seconds <- 0
    preconditioner_seconds <- 0
    operator_applications <- 0L
    preconditioner_applications <- 0L
    residual_scale <- residual_measure(right_hand_side)
    if (!is.finite(residual_scale) || residual_scale < 0) {
        stop('Krylov residual measure was not finite and nonnegative')
    }
    if (residual_scale <= .Machine$double.eps) residual_scale <- 1
    warm_started <- length(initial_solution) == dimension &&
        all(is.finite(initial_solution))
    solution <- if (warm_started) {
        as.numeric(initial_solution)
    } else numeric(dimension)
    residual <- if (warm_started) {
        right_hand_side - timed_application(
            'operator', 0L, 1L, function() product(solution)
        )
    } else right_hand_side
    residual_norm <- residual_measure(residual) / residual_scale
    initial_residual <- residual_norm
    if (warm_started && (!is.finite(residual_norm) || residual_norm >= 1)) {
        notify(
            'Krylov warm start rejected',
            residual=format(residual_norm, digits=5)
        )
        solution <- numeric(dimension)
        residual <- right_hand_side
        residual_norm <- 1
        warm_started <- FALSE
    }
    notify(
        'Krylov solve started',
        method='fgmres',
        dimension=dimension,
        tolerance=format(tolerance, digits=5),
        maximum_iterations=maxit,
        residual=format(residual_norm, digits=5),
        warm_start=warm_started
    )
    finish <- function(converged, iterations, termination) {
        output <- list(
            solution=solution,
            converged=converged,
            iterations=as.integer(iterations),
            residual=residual_norm,
            termination=termination,
            warm_started=warm_started,
            initial_residual=initial_residual,
            operator_applications=operator_applications,
            preconditioner_applications=preconditioner_applications,
            operator_seconds=operator_seconds,
            preconditioner_seconds=preconditioner_seconds
        )
        notify(
            'Krylov solve complete',
            method='fgmres',
            iterations=iterations,
            residual=format(residual_norm, digits=5),
            termination=termination,
            operator_seconds=format(operator_seconds, digits=5),
            preconditioner_seconds=format(preconditioner_seconds, digits=5)
        )
        output
    }
    if (residual_norm <= tolerance) {
        return(finish(TRUE, 0L, 'initial residual'))
    }
    euclidean_norm <- sqrt(sum(residual^2))
    if (!is.finite(euclidean_norm) ||
            euclidean_norm <= .Machine$double.eps) {
        return(finish(FALSE, 0L, 'initial residual breakdown'))
    }
    initial_solution_value <- solution
    basis <- matrix(0, nrow=dimension, ncol=maxit + 1L)
    images <- matrix(0, nrow=dimension, ncol=maxit)
    hessenberg <- matrix(0, nrow=maxit + 1L, ncol=maxit)
    basis[, 1L] <- residual / euclidean_norm
    termination <- 'iteration limit'
    iteration <- 0L
    for (iteration in seq_len(maxit)) {
        images[, iteration] <- timed_application(
            'preconditioner', iteration, 1L,
            function() preconditioner(basis[, iteration])
        )
        candidate_image <- timed_application(
            'operator', iteration, 1L,
            function() product(images[, iteration])
        )
        for (index in seq_len(iteration)) {
            coefficient <- sum(basis[, index] * candidate_image)
            hessenberg[index, iteration] <- coefficient
            candidate_image <- candidate_image -
                coefficient * basis[, index]
        }
        for (index in seq_len(iteration)) {
            correction <- sum(basis[, index] * candidate_image)
            hessenberg[index, iteration] <-
                hessenberg[index, iteration] + correction
            candidate_image <- candidate_image - correction * basis[, index]
        }
        next_norm <- sqrt(sum(candidate_image^2))
        hessenberg[iteration + 1L, iteration] <- next_norm
        breakdown <- !is.finite(next_norm) ||
            next_norm <= sqrt(.Machine$double.eps)
        if (!breakdown) {
            basis[, iteration + 1L] <- candidate_image / next_norm
        }
        local_hessenberg <- hessenberg[
            seq_len(iteration + 1L), seq_len(iteration), drop=FALSE
        ]
        local_right_hand_side <- c(
            euclidean_norm, numeric(iteration)
        )
        coefficients <- least_squares(
            local_hessenberg, local_right_hand_side
        )
        solution <- initial_solution_value + as.numeric(
            images[, seq_len(iteration), drop=FALSE] %*% coefficients
        )
        residual_coordinates <- local_right_hand_side -
            as.numeric(local_hessenberg %*% coefficients)
        candidate_residual <- as.numeric(
            basis[, seq_len(iteration + 1L), drop=FALSE] %*%
                residual_coordinates
        )
        residual_norm <- residual_measure(candidate_residual) / residual_scale
        notify(
            'Krylov iteration complete',
            method='fgmres',
            iteration=iteration,
            residual=format(residual_norm, digits=5),
            operator_seconds=format(operator_seconds, digits=5),
            preconditioner_seconds=format(
                preconditioner_seconds, digits=5
            )
        )
        if (is.finite(residual_norm) && residual_norm <= tolerance) {
            return(finish(TRUE, iteration, 'tolerance'))
        }
        if (breakdown) {
            termination <- 'Arnoldi breakdown'
            break
        }
    }
    finish(FALSE, iteration, termination)
}

.cdrgam_gaulss_block_precondition <- function(
        vector, factors, split_vector,
        mode=c('diagonal', 'triangular'), cross_scale=NULL
) {
    mode <- match.arg(mode)
    vectors <- split_vector(vector)
    location <- as.numeric(.cdr_factor_solve(
        factors$location, vectors$location
    ))
    scale_rhs <- vectors$scale
    if (identical(mode, 'triangular')) {
        if (!is.function(cross_scale)) {
            stop('A cross-scale action is required for triangular solves')
        }
        scale_rhs <- scale_rhs - cross_scale(location)
    }
    scale <- as.numeric(.cdr_factor_solve(factors$scale, scale_rhs))
    c(location, scale)
}

.cdrgam_gaulss_schur_preconditioner <- function(
        diagonal_solve, correction, right_hand_side, rank=6L, progress=NULL
) {
    notify <- function(event, ...) {
        if (is.function(progress)) progress(event, ...)
        invisible(NULL)
    }
    baseline <- list(
        apply=diagonal_solve,
        rank=0L,
        middle_minimum=NA_real_,
        middle_maximum=NA_real_
    )
    dimension <- length(right_hand_side)
    rank <- min(as.integer(rank), dimension)
    if (!is.finite(rank) || rank < 1L || dimension < 1L) return(baseline)
    notify(
        'Schur preconditioner setup started',
        dimension=dimension,
        requested_rank=rank
    )
    started <- proc.time()[['elapsed']]
    sketch <- matrix(0, nrow=dimension, ncol=rank)
    leading <- diagonal_solve(right_hand_side)
    leading_norm <- sqrt(sum(leading^2))
    if (!is.finite(leading_norm) || leading_norm <= .Machine$double.eps) {
        leading <- right_hand_side
        leading_norm <- sqrt(sum(leading^2))
    }
    if (!is.finite(leading_norm) || leading_norm <= .Machine$double.eps) {
        return(baseline)
    }
    sketch[, 1L] <- leading / leading_norm
    if (rank > 1L) {
        indices <- seq_len(dimension)
        for (column in 2:rank) {
            phase <- indices * (sqrt(2) * column + sqrt(3)) + column^2
            sketch[, column] <- ifelse(sin(phase) >= 0, 1, -1)
        }
    }
    sketch <- qr.Q(qr(sketch), complete=FALSE)
    apply_batch <- function(operation, values) {
        batched <- tryCatch(operation(values), error=function(error) NULL)
        if (is.matrix(batched) && identical(dim(batched), dim(values))) {
            return(batched)
        }
        vapply(seq_len(ncol(values)), function(column) {
            operation(values[, column])
        }, numeric(nrow(values)))
    }
    image <- apply_batch(correction, sketch)
    gram <- (crossprod(sketch, image) + crossprod(image, sketch)) / 2
    decomposition <- eigen(gram, symmetric=TRUE)
    scale <- max(abs(decomposition$values), 1)
    keep <- decomposition$values > 1e-8 * scale
    if (!any(keep)) {
        notify(
            'Schur preconditioner setup complete',
            seconds=format(proc.time()[['elapsed']] - started, digits=5),
            retained_rank=0L,
            reason='no positive sampled correction directions'
        )
        return(baseline)
    }
    vectors <- decomposition$vectors[, keep, drop=FALSE]
    values <- decomposition$values[keep]
    low_rank <- image %*% sweep(vectors, 2L, sqrt(values), `/`)
    diagonal_low_rank <- apply_batch(diagonal_solve, low_rank)
    middle <- diag(ncol(low_rank)) - crossprod(
        low_rank, diagonal_low_rank
    )
    middle <- (middle + t(middle)) / 2
    middle_decomposition <- eigen(middle, symmetric=TRUE)
    middle_scale <- max(abs(middle_decomposition$values), 1)
    floor <- 1e-6 * middle_scale
    inverse_values <- ifelse(
        middle_decomposition$values >= 0, 1, -1
    ) / pmax(abs(middle_decomposition$values), floor)
    apply <- function(vector) {
        base <- diagonal_solve(vector)
        coordinates <- crossprod(low_rank, base)
        correction_coordinates <- middle_decomposition$vectors %*% (
            inverse_values * crossprod(
                middle_decomposition$vectors, coordinates
            )
        )
        as.numeric(base + diagonal_low_rank %*% correction_coordinates)
    }
    notify(
        'Schur preconditioner setup complete',
        seconds=format(proc.time()[['elapsed']] - started, digits=5),
        retained_rank=ncol(low_rank),
        middle_minimum=format(
            min(abs(middle_decomposition$values)), digits=5
        ),
        middle_maximum=format(
            max(abs(middle_decomposition$values)), digits=5
        )
    )
    list(
        apply=apply,
        rank=as.integer(ncol(low_rank)),
        middle_minimum=min(abs(middle_decomposition$values)),
        middle_maximum=max(abs(middle_decomposition$values))
    )
}

.cdrgam_gaulss_schur_solve <- function(
        systems, factors, gradients, cross_location, cross_scale,
        tolerance, maxit, initial_solution=NULL,
        eliminate=c('auto', 'location', 'scale'), preconditioner_rank=6L,
        progress=NULL
) {
    eliminate <- match.arg(eliminate)
    dimensions <- vapply(systems, nrow, integer(1))
    if (identical(eliminate, 'auto')) {
        eliminate <- if (dimensions[['location']] >=
                dimensions[['scale']]) 'location' else 'scale'
    }
    remaining <- if (identical(eliminate, 'location')) 'scale' else 'location'
    solve_block <- function(parameter, right_hand_side) {
        solved <- .cdr_factor_solve(factors[[parameter]], right_hand_side)
        if (is.matrix(right_hand_side)) as.matrix(solved) else as.numeric(solved)
    }
    cross_from <- function(parameter, vector) {
        if (identical(parameter, 'location')) {
            cross_scale(vector)
        } else cross_location(vector)
    }
    eliminated_gradient <- gradients[[eliminate]]
    remaining_gradient <- gradients[[remaining]]
    right_hand_side <- -remaining_gradient + cross_from(
        eliminate, solve_block(eliminate, eliminated_gradient)
    )
    product <- function(vector) {
        crossed <- cross_from(remaining, vector)
        as.numeric(systems[[remaining]] %*% vector) -
            cross_from(eliminate, solve_block(eliminate, crossed))
    }
    correction <- function(vector) {
        cross_from(
            eliminate,
            solve_block(eliminate, cross_from(remaining, vector))
        )
    }
    diagonal_preconditioner <- function(vector) {
        solve_block(remaining, vector)
    }
    approximation <- .cdrgam_gaulss_schur_preconditioner(
        diagonal_preconditioner,
        correction,
        right_hand_side,
        rank=preconditioner_rank,
        progress=progress
    )
    preconditioner <- approximation$apply
    residual_measure <- function(vector) {
        sqrt(sum(diagonal_preconditioner(vector)^2))
    }
    initial_remaining <- NULL
    if (length(initial_solution) == sum(dimensions) &&
            all(is.finite(initial_solution))) {
        ranges <- split(
            seq_len(sum(dimensions)),
            rep(names(dimensions), dimensions)
        )
        initial_remaining <- initial_solution[ranges[[remaining]]]
    }
    solved <- .cdrgam_fgmres(
        right_hand_side, product, preconditioner,
        tolerance=tolerance,
        maxit=min(as.integer(maxit), dimensions[[remaining]]),
        residual_measure=residual_measure,
        initial_solution=initial_remaining,
        progress=progress
    )
    if (any(!is.finite(solved$solution))) {
        solved$eliminate <- eliminate
        solved$remaining <- remaining
        solved$schur_preconditioner_rank <- approximation$rank
        return(solved)
    }
    steps <- list(location=NULL, scale=NULL)
    steps[[remaining]] <- solved$solution
    steps[[eliminate]] <- solve_block(
        eliminate,
        -eliminated_gradient - cross_from(remaining, solved$solution)
    )
    solved$solution <- unlist(steps, use.names=FALSE)
    solved$eliminate <- eliminate
    solved$remaining <- remaining
    solved$schur_preconditioner_rank <- approximation$rank
    solved
}

.cdrgam_gaulss_coupled_update <- function(
        assemblies, coefficients, predictors, gradients, moments, penalties,
        factors, cross_weight, objective, value, workers=1L, blas_threads=1L,
        previous_nonlinear_norm=NA_real_, initial_solution=NULL,
        preconditioner_mode=c('diagonal', 'triangular', 'schur'),
        schur_rank=6L, progress=NULL
) {
    preconditioner_mode <- match.arg(preconditioner_mode)
    if (is.function(progress)) progress(
        'coupled factorization started',
        preconditioner=preconditioner_mode
    )
    factor_started <- proc.time()[['elapsed']]
    systems <- list(
        location=Matrix::forceSymmetric(
            moments$information$location + penalties$location, uplo='U'
        ),
        scale=Matrix::forceSymmetric(
            moments$information$scale + penalties$scale, uplo='U'
        )
    )
    factorizations <- Map(function(system, assembly, previous) {
        .cdrgam_gaulss_step_factor(system, assembly$supernodal, previous)
    }, systems, assemblies, factors)
    factors <- lapply(factorizations, `[[`, 'factor')
    if (is.function(progress)) progress(
        'coupled factorization complete',
        seconds=format(proc.time()[['elapsed']] - factor_started, digits=5)
    )
    dimensions <- vapply(assemblies, `[[`, integer(1), 'dimension')
    ranges <- split(seq_len(sum(dimensions)), rep(names(dimensions), dimensions))
    split_vector <- function(vector) list(
        location=vector[ranges$location], scale=vector[ranges$scale]
    )
    gradient <- unlist(gradients, use.names=FALSE)
    product <- function(vector) {
        vectors <- split_vector(vector)
        cross <- .cdrgam_gaulss_cross_product(
            assemblies, cross_weight, vectors,
            workers=workers, blas_threads=blas_threads
        )
        unlist(list(
            location=as.numeric(systems$location %*% vectors$location) +
                cross$location,
            scale=as.numeric(systems$scale %*% vectors$scale) + cross$scale
        ), use.names=FALSE)
    }
    cross_scale <- function(location) {
        .cdrgam_gaulss_cross_product(
            assemblies, cross_weight,
            list(location=location, scale=numeric(dimensions[['scale']])),
            workers=workers, blas_threads=blas_threads,
            outputs='scale'
        )$scale
    }
    cross_location <- function(scale) {
        .cdrgam_gaulss_cross_product(
            assemblies, cross_weight,
            list(location=numeric(dimensions[['location']]), scale=scale),
            workers=workers, blas_threads=blas_threads,
            outputs='location'
        )$location
    }
    diagonal_preconditioner <- function(vector) {
        .cdrgam_gaulss_block_precondition(
            vector, factors, split_vector, mode='diagonal'
        )
    }
    preconditioner <- if (identical(preconditioner_mode, 'diagonal')) {
        diagonal_preconditioner
    } else if (identical(preconditioner_mode, 'triangular')) function(vector) {
        .cdrgam_gaulss_block_precondition(
            vector, factors, split_vector, mode='triangular',
            cross_scale=cross_scale
        )
    } else diagonal_preconditioner
    residual_measure <- function(vector) {
        sqrt(sum(diagonal_preconditioner(vector)^2))
    }
    nonlinear_norm <- residual_measure(gradient)
    linear_tolerance <- .cdrgam_gaulss_forcing_tolerance(
        nonlinear_norm, previous_nonlinear_norm
    )
    linear_acceptance_tolerance <- if (identical(
        preconditioner_mode, 'schur'
    )) {
        .cdrgam_gaulss_relaxed_forcing_tolerance(linear_tolerance)
    } else linear_tolerance
    solved <- if (identical(preconditioner_mode, 'schur')) {
        .cdrgam_gaulss_schur_solve(
            systems, factors, gradients,
            cross_location, cross_scale,
            tolerance=linear_acceptance_tolerance,
            maxit=12L,
            initial_solution=initial_solution,
            preconditioner_rank=schur_rank,
            progress=progress
        )
    } else .cdrgam_bicgstab(
        -gradient, product, preconditioner,
        tolerance=linear_tolerance, maxit=min(12L, length(gradient)),
        residual_measure=residual_measure,
        initial_solution=initial_solution,
        progress=progress
    )
    if (any(!is.finite(solved$solution))) {
        stop(
            'Matrix-free coupled solve failed: ', solved$termination,
            '; relative residual ', format(solved$residual, digits=5)
        )
    }
    linear_converged <- FALSE
    if (identical(preconditioner_mode, 'schur')) {
        full_residual <- -gradient - product(solved$solution)
        solved$residual <- residual_measure(full_residual) /
            max(residual_measure(-gradient), .Machine$double.eps)
        linear_converged <- is.finite(solved$residual) &&
            solved$residual <= linear_tolerance
        if (!is.finite(solved$residual) ||
                solved$residual > linear_acceptance_tolerance) {
            stop(
                'Matrix-free Schur solve failed full-system certification; ',
                'relative residual ', format(solved$residual, digits=5),
                '; acceptance tolerance ',
                format(linear_acceptance_tolerance, digits=5)
            )
        }
    } else if (!isTRUE(solved$converged)) {
        stop(
            'Matrix-free coupled solve failed: ', solved$termination,
            '; relative residual ', format(solved$residual, digits=5)
        )
    } else linear_converged <- TRUE
    step <- solved$solution
    steps <- split_vector(step)
    predictor_steps <- Map(function(assembly, parameter_step) {
        .cdrgam_distributional_sparse_predictor(
            assembly, parameter_step,
            numeric(assembly$observation_count)
        )
    }, assemblies, steps)
    directional_derivative <- sum(gradient * step)
    if (!is.finite(directional_derivative) || directional_derivative >= 0) {
        stop(
            'Matrix-free coupled solve did not produce a descent direction; ',
            'directional derivative ',
            format(directional_derivative, digits=5)
        )
    }
    if (!linear_converged && is.function(progress)) progress(
        'coupled inexact direction eligible',
        residual=format(solved$residual, digits=5),
        forcing_tolerance=format(linear_tolerance, digits=5),
        acceptance_tolerance=format(
            linear_acceptance_tolerance, digits=5
        )
    )
    step_size <- 1
    minimum_step <- 2^-20
    best_change <- Inf
    while (step_size >= minimum_step) {
        candidate_coefficients <- Map(function(current, candidate_step) {
            current + step_size * candidate_step
        }, coefficients, steps)
        candidate_predictors <- Map(function(current, candidate_step) {
            current + step_size * candidate_step
        }, predictors, predictor_steps)
        candidate_value <- objective(
            candidate_coefficients, candidate_predictors
        )
        best_change <- min(best_change, candidate_value - value, na.rm=TRUE)
        if (is.finite(candidate_value) && candidate_value < value) {
            return(list(
                coefficients=candidate_coefficients,
                predictors=candidate_predictors,
                value=candidate_value,
                factors=factors,
                damping=max(vapply(
                    factorizations, `[[`, numeric(1), 'damping'
                )),
                accepted=TRUE,
                step=step_size * step,
                linear_iterations=solved$iterations,
                linear_residual=solved$residual,
                linear_tolerance=linear_tolerance,
                linear_acceptance_tolerance=linear_acceptance_tolerance,
                linear_converged=linear_converged,
                linear_warm_started=solved$warm_started,
                linear_initial_residual=solved$initial_residual,
                linear_solution=solved$solution,
                preconditioner_mode=preconditioner_mode,
                nonlinear_norm=nonlinear_norm,
                directional_derivative=step_size * directional_derivative,
                best_change=best_change,
                numerically_converged=FALSE
            ))
        }
        step_size <- step_size / 2
    }
    list(
        coefficients=coefficients,
        predictors=predictors,
        value=value,
        factors=factors,
        damping=max(vapply(
            factorizations, `[[`, numeric(1), 'damping'
        )),
        accepted=FALSE,
        step=numeric(length(step)),
        linear_iterations=solved$iterations,
        linear_residual=solved$residual,
        linear_tolerance=linear_tolerance,
        linear_acceptance_tolerance=linear_acceptance_tolerance,
        linear_converged=linear_converged,
        linear_warm_started=solved$warm_started,
        linear_initial_residual=solved$initial_residual,
        linear_solution=solved$solution,
        preconditioner_mode=preconditioner_mode,
        nonlinear_norm=nonlinear_norm,
        directional_derivative=directional_derivative,
        best_change=best_change,
        numerically_converged=.cdrgam_gaulss_numerically_converged(
            directional_derivative, best_change, value, minimum_step
        )
    )
}

.cdrgam_gaulss_inner_algorithm <- function(algorithm='auto', dimension) {
    algorithm <- match.arg(
        algorithm, c('auto', 'mixed', 'schur', 'alternating', 'joint')
    )
    if (!identical(algorithm, 'auto')) return(algorithm)
    if (dimension >= 256L) 'mixed' else 'joint'
}

.cdrgam_gaulss_forcing_tolerance <- function(
        nonlinear_norm, previous_norm=NA_real_, minimum=1e-6,
        maximum=0.1, gamma=0.9, exponent=1.5
) {
    if (!is.finite(nonlinear_norm) || nonlinear_norm < 0 ||
            !is.finite(previous_norm) || previous_norm <= 0) {
        return(maximum)
    }
    ratio <- nonlinear_norm / previous_norm
    min(maximum, max(minimum, gamma * ratio^exponent))
}

.cdrgam_gaulss_relaxed_forcing_tolerance <- function(
        forcing_tolerance, maximum=0.5
) {
    min(maximum, sqrt(max(0, forcing_tolerance)))
}

.cdrgam_gaulss_retry_delay <- function(consecutive_failures, maximum=8L) {
    consecutive_failures <- max(0L, as.integer(consecutive_failures))
    maximum <- max(0L, as.integer(maximum))
    if (consecutive_failures == 0L || maximum == 0L) return(0L)
    as.integer(min(maximum, 2^min(consecutive_failures - 1L, 30L)))
}

.cdrgam_gaulss_resume_state <- function(
        coefficients, predictors, algorithm, mixed_iterations=0L,
        moment_workers=1L, anderson_states=list(), anderson_residuals=list(),
        anderson_attempts=0L, anderson_acceptances=0L,
        coupled_disabled=FALSE, coupled_failures=0L,
        coupled_consecutive_failures=0L, coupled_retry_delay=0L,
        coupled_linear_solution=NULL, nonlinear_norm=NA_real_
) {
    list(
        parameter_coefficients=coefficients,
        linear_predictors=do.call(cbind, predictors),
        algorithm=algorithm,
        mixed_iterations=as.integer(mixed_iterations),
        moment_workers=as.integer(moment_workers),
        anderson_states=anderson_states,
        anderson_residuals=anderson_residuals,
        anderson_attempts=as.integer(anderson_attempts),
        anderson_acceptances=as.integer(anderson_acceptances),
        coupled_disabled=isTRUE(coupled_disabled),
        coupled_failures=as.integer(coupled_failures),
        coupled_consecutive_failures=as.integer(
            coupled_consecutive_failures
        ),
        coupled_retry_delay=as.integer(coupled_retry_delay),
        coupled_linear_solution=coupled_linear_solution,
        nonlinear_norm=nonlinear_norm
    )
}

.cdrgam_gaulss_solution_resume_state <- function(solution) {
    fields <- c(
        'parameter_coefficients', 'linear_predictors', 'algorithm',
        'mixed_iterations', 'moment_workers', 'anderson_states',
        'anderson_residuals', 'anderson_attempts', 'anderson_acceptances',
        'coupled_disabled', 'coupled_failures',
        'coupled_consecutive_failures', 'coupled_retry_delay',
        'coupled_linear_solution', 'nonlinear_norm'
    )
    solution[intersect(fields, names(solution))]
}

.cdrgam_gaulss_sparse_fixed <- function(
        assemblies, smoothing_parameters, b, tolerance=1e-5,
        relative_tolerance=1e-8, maxit=200L, initial=NULL, cores=1L,
        moment_workers=NULL, schur_rank=6L,
        algorithm=c('auto', 'mixed', 'schur', 'alternating', 'joint'),
        progress=NULL, checkpoint=NULL
) {
    notify <- function(event, ...) {
        if (is.function(progress)) progress(event, ...)
        invisible(NULL)
    }
    dimension <- sum(vapply(assemblies, `[[`, integer(1), 'dimension'))
    dimensions <- vapply(assemblies, `[[`, integer(1), 'dimension')
    algorithm <- .cdrgam_gaulss_inner_algorithm(
        match.arg(algorithm), dimension
    )
    counts <- vapply(
        assemblies, function(assembly) length(assembly$penalty_components),
        integer(1)
    )
    location_indices <- seq_len(counts[['location']])
    scale_indices <- counts[['location']] + seq_len(counts[['scale']])
    penalties <- list(
        location=.cdrgam_sparse_penalty(
            assemblies$location, smoothing_parameters[location_indices]
        ),
        scale=.cdrgam_sparse_penalty(
            assemblies$scale, smoothing_parameters[scale_indices]
        )
    )
    y <- assemblies$location$setup$y
    weights <- assemblies$location$setup$w
    if (is.null(weights)) weights <- rep.int(1, length(y))
    cores <- .cdrgam_available_cores(cores)
    requested_moment_workers <- .cdrgam_parallel_plan(
        cores,
        max(1L, ceiling(length(y) / min(
            assemblies$location$crossprod_chunk_size,
            assemblies$scale$crossprod_chunk_size
        ))),
        moment_workers
    )$workers
    resolved_moment_workers <- 1L
    offsets <- lapply(assemblies, function(assembly) {
        value <- assembly$setup$offset
        if (is.null(value)) numeric(length(y)) else value
    })
    objective <- function(coefficients, predictors=NULL) {
        if (is.null(predictors)) {
            predictors <- Map(
                .cdrgam_distributional_sparse_predictor,
                assemblies,
                coefficients,
                offsets
            )
        }
        sigma <- b + exp(predictors$scale)
        if (any(!is.finite(sigma))) return(Inf)
        residual <- y - predictors$location
        sum(weights * (log(sigma) + 0.5 * (residual / sigma)^2)) +
            0.5 * as.numeric(Matrix::crossprod(
                coefficients$location, penalties$location %*% coefficients$location
            )) +
            0.5 * as.numeric(Matrix::crossprod(
                coefficients$scale, penalties$scale %*% coefficients$scale
            ))
    }
    coefficients <- if (!is.null(initial) &&
            is.list(initial$parameter_coefficients)) {
        initial$parameter_coefficients
    } else if (is.list(initial) &&
            all(c('location', 'scale') %in% names(initial))) {
        initial[c('location', 'scale')]
    } else NULL
    valid_initial <- !is.null(coefficients) &&
        length(coefficients$location) == assemblies$location$dimension &&
        length(coefficients$scale) == assemblies$scale$dimension &&
        all(is.finite(unlist(coefficients, use.names=FALSE)))
    notify(
        'solve started',
        algorithm=algorithm,
        dimension=dimension,
        observations=length(y),
        warm_start=valid_initial,
        maximum_iterations=maxit
    )
    resolved_moment_workers <- if (valid_initial &&
            length(initial$moment_workers) == 1L &&
            is.finite(initial$moment_workers)) {
        min(requested_moment_workers, as.integer(initial$moment_workers))
    } else 1L
    calibrate_moment_workers <- function(moments, output_multiplier=1) {
        if (resolved_moment_workers > 1L || requested_moment_workers <= 1L) {
            return(invisible(resolved_moment_workers))
        }
        output_bytes <- sum(vapply(
            c(moments$gradients, moments$information),
            function(value) if (is.null(value)) 0 else utils::object.size(value),
            numeric(1)
        )) * output_multiplier
        plan <- .cdrgam_memory_parallel_plan(
            cores=cores,
            tasks=max(1L, moments$chunks),
            workers=requested_moment_workers,
            per_worker_bytes=max(1, 2 * output_bytes),
            memory_fraction=0.35,
            reserve_bytes=1024^3
        )
        resolved_moment_workers <<- plan$workers
        invisible(resolved_moment_workers)
    }
    retained_predictors <- if (valid_initial &&
            is.matrix(initial$linear_predictors) &&
            nrow(initial$linear_predictors) == length(y) &&
            all(c('location', 'scale') %in%
                colnames(initial$linear_predictors)) &&
            all(is.finite(initial$linear_predictors))) {
        list(
            location=initial$linear_predictors[, 'location'],
            scale=initial$linear_predictors[, 'scale']
        )
    } else NULL
    retained_factors <- if (valid_initial &&
            is.list(initial$inner_factors) &&
            all(c('location', 'scale') %in% names(initial$inner_factors))) {
        initial$inner_factors
    } else list(location=NULL, scale=NULL)
    if (!valid_initial) {
        notify('initialization started')
        initial_sigma <- sqrt(stats::weighted.mean(
            (y - stats::weighted.mean(y, weights))^2, weights
        ))
        initial_eta_scale <- rep.int(
            log(max(initial_sigma - b, 1e-4)), length(y)
        )
        coefficients <- list(
            location=numeric(assemblies$location$dimension),
            scale=numeric(assemblies$scale$dimension)
        )
        predictors <- list(
            location=offsets$location,
            scale=initial_eta_scale
        )
        q <- exp(predictors$scale)
        sigma <- b + q
        initial_moments_started <- proc.time()[['elapsed']]
        moments <- .cdrgam_gaulss_sparse_moments(
            assemblies, coefficients, sigma, q, y - predictors$location,
            weights, parameters='location'
        )
        notify(
            'initial moments complete',
            seconds=format(
                proc.time()[['elapsed']] - initial_moments_started,
                digits=5
            ),
            chunks=moments$chunks
        )
        calibrate_moment_workers(moments, output_multiplier=2)
        location_system <- Matrix::forceSymmetric(
            moments$information$location + penalties$location, uplo='U'
        )
        location_rhs <- -moments$gradients$location
        factor_location <- .cdrgam_gaulss_step_factor(
            location_system, assemblies$location$supernodal
        )$factor
        coefficients$location <- as.numeric(
            .cdr_factor_solve(factor_location, location_rhs)
        )
        scale_system <- NULL
        scale_rhs <- numeric(assemblies$scale$dimension)
        scale_starts <- seq.int(
            1L, length(y), by=assemblies$scale$crossprod_chunk_size
        )
        notify('initial scale system started', chunks=length(scale_starts))
        scale_started <- proc.time()[['elapsed']]
        for (chunk_index in seq_along(scale_starts)) {
            start <- scale_starts[[chunk_index]]
            rows <- start:min(
                length(y),
                start + assemblies$scale$crossprod_chunk_size - 1L
            )
            X <- .cdrgam_sparse_design_chunk(assemblies$scale, rows)
            part <- Matrix::crossprod(
                Matrix::Diagonal(x=sqrt(weights[rows])) %*% X
            )
            scale_system <- if (is.null(scale_system)) {
                part
            } else scale_system + part
            scale_rhs <- scale_rhs + as.numeric(Matrix::crossprod(
                X,
                weights[rows] *
                    (initial_eta_scale[rows] - offsets$scale[rows])
            ))
            notify(
                'initial scale chunk complete',
                chunk=chunk_index,
                chunks=length(scale_starts),
                seconds=format(
                    proc.time()[['elapsed']] - scale_started,
                    digits=5
                )
            )
        }
        scale_system <- Matrix::forceSymmetric(
            scale_system + penalties$scale, uplo='U'
        )
        factor_scale <- .cdrgam_gaulss_step_factor(
            scale_system, assemblies$scale$supernodal
        )$factor
        coefficients$scale <- as.numeric(
            .cdr_factor_solve(factor_scale, scale_rhs)
        )
        notify(
            'initialization complete',
            seconds=format(
                proc.time()[['elapsed']] - initial_moments_started,
                digits=5
            )
        )
    }
    predictors <- if (!is.null(retained_predictors)) {
        retained_predictors
    } else Map(
        .cdrgam_distributional_sparse_predictor,
        assemblies,
        coefficients,
        offsets
    )
    value <- objective(coefficients, predictors)
    if (!is.finite(value) && valid_initial) {
        return(.cdrgam_gaulss_sparse_fixed(
            assemblies, smoothing_parameters, b,
            tolerance=tolerance, relative_tolerance=relative_tolerance,
            maxit=maxit, initial=NULL,
            cores=cores, moment_workers=moment_workers,
            schur_rank=schur_rank, algorithm=algorithm, progress=progress,
            checkpoint=checkpoint
        ))
    }
    converged <- FALSE
    gradient_norm <- Inf
    termination <- 'iteration limit'
    step_norm <- NA_real_
    directional_derivative <- NA_real_
    best_line_search_change <- NA_real_
    step_mode <- NA_character_
    linear_iterations <- 0L
    linear_residual <- NA_real_
    linear_tolerance <- NA_real_
    linear_acceptance_tolerance <- NA_real_
    linear_converged <- NA
    linear_warm_started <- FALSE
    linear_initial_residual <- NA_real_
    linear_preconditioner <- NA_character_
    relative_reduction <- Inf
    relative_predicted_reduction <- Inf
    anderson_states <- if (valid_initial &&
            is.list(initial$anderson_states)) {
        initial$anderson_states
    } else list()
    anderson_residuals <- if (valid_initial &&
            is.list(initial$anderson_residuals)) {
        initial$anderson_residuals
    } else list()
    anderson_attempts <- if (valid_initial &&
            length(initial$anderson_attempts) == 1L &&
            is.finite(initial$anderson_attempts)) {
        as.integer(initial$anderson_attempts)
    } else 0L
    anderson_acceptances <- if (valid_initial &&
            length(initial$anderson_acceptances) == 1L &&
            is.finite(initial$anderson_acceptances)) {
        as.integer(initial$anderson_acceptances)
    } else 0L
    maximum_damping <- 0
    factor_location <- retained_factors$location
    factor_scale <- retained_factors$scale
    factor_joint <- NULL
    coupled_disabled <- !identical(algorithm, 'schur') && valid_initial &&
        isTRUE(initial$coupled_disabled)
    coupled_failures <- if (valid_initial &&
            length(initial$coupled_failures) == 1L &&
            is.finite(initial$coupled_failures)) {
        as.integer(initial$coupled_failures)
    } else 0L
    coupled_consecutive_failures <- if (valid_initial &&
            length(initial$coupled_consecutive_failures) == 1L &&
            is.finite(initial$coupled_consecutive_failures)) {
        as.integer(initial$coupled_consecutive_failures)
    } else 0L
    coupled_retry_delay <- if (identical(algorithm, 'schur') &&
            valid_initial &&
            length(initial$coupled_retry_delay) == 1L &&
            is.finite(initial$coupled_retry_delay) &&
            initial$coupled_retry_delay >= 0) {
        as.integer(initial$coupled_retry_delay)
    } else 0L
    nonlinear_norm <- if (valid_initial &&
            length(initial$nonlinear_norm) == 1L &&
            is.finite(initial$nonlinear_norm) &&
            initial$nonlinear_norm >= 0) {
        initial$nonlinear_norm
    } else NA_real_
    coupled_linear_solution <- if (valid_initial &&
            length(initial$coupled_linear_solution) == dimension &&
            all(is.finite(initial$coupled_linear_solution))) {
        as.numeric(initial$coupled_linear_solution)
    } else NULL
    mixed_iteration_offset <- if (valid_initial &&
            initial$algorithm %in% c('mixed', 'schur') &&
            length(initial$mixed_iterations) == 1L &&
            is.finite(initial$mixed_iterations)) {
        as.integer(initial$mixed_iterations)
    } else 0L
    for (iteration in seq_len(maxit)) {
        global_iteration <- mixed_iteration_offset + iteration
        if (is.function(checkpoint)) checkpoint(
            .cdrgam_gaulss_resume_state(
                coefficients, predictors, algorithm,
                mixed_iterations=global_iteration - 1L,
                moment_workers=resolved_moment_workers,
                anderson_states=anderson_states,
                anderson_residuals=anderson_residuals,
                anderson_attempts=anderson_attempts,
                anderson_acceptances=anderson_acceptances,
                coupled_disabled=coupled_disabled,
                coupled_failures=coupled_failures,
                coupled_consecutive_failures=
                    coupled_consecutive_failures,
                coupled_retry_delay=coupled_retry_delay,
                coupled_linear_solution=coupled_linear_solution,
                nonlinear_norm=nonlinear_norm
            )
        )
        coupled_eligible <- algorithm %in% c('mixed', 'schur') &&
            global_iteration > 2L && !coupled_disabled
        coupled_step <- coupled_eligible && coupled_retry_delay <= 0L
        if (coupled_eligible && !coupled_step) {
            notify(
                'coupled update skipped',
                iteration=global_iteration,
                retry_delay=coupled_retry_delay
            )
            coupled_retry_delay <- max(0L, coupled_retry_delay - 1L)
        }
        notify(
            'iteration started',
            iteration=global_iteration,
            objective=format(value, digits=10),
            coupled=coupled_step,
            moment_workers=resolved_moment_workers
        )
        q <- exp(predictors$scale)
        sigma <- b + q
        residual <- y - predictors$location
        moments_started <- proc.time()[['elapsed']]
        notify(
            'moments started',
            iteration=mixed_iteration_offset + iteration,
            coupled=coupled_step,
            workers=resolved_moment_workers
        )
        moments <- .cdrgam_gaulss_sparse_moments(
            assemblies, coefficients, sigma, q, residual, weights,
            include_observed=coupled_step, include_cross=FALSE,
            workers=resolved_moment_workers,
            blas_threads=max(1L, cores %/% resolved_moment_workers)
        )
        notify(
            'moments complete',
            iteration=mixed_iteration_offset + iteration,
            seconds=format(
                proc.time()[['elapsed']] - moments_started, digits=5
            ),
            chunks=moments$chunks,
            workers=resolved_moment_workers
        )
        calibrate_moment_workers(moments)
        gradients <- list(
            location=moments$gradients$location +
                as.numeric(penalties$location %*% coefficients$location),
            scale=moments$gradients$scale +
                as.numeric(penalties$scale %*% coefficients$scale)
        )
        gradient_norm <- max(abs(unlist(gradients, use.names=FALSE)))
        if (gradient_norm < tolerance) {
            converged <- TRUE
            termination <- 'score tolerance'
            break
        }
        coupled_failed <- FALSE
        if (coupled_step) {
            coupled_value <- value
            requested_preconditioner <- if (identical(
                algorithm, 'schur'
            )) 'schur' else 'diagonal'
            coupled_progress <- function(event, ...) notify(
                event,
                inner_iteration=mixed_iteration_offset + iteration,
                ...
            )
            coupled_update <- tryCatch(
                .cdrgam_gaulss_coupled_update(
                    assemblies, coefficients, predictors, gradients,
                    moments, penalties,
                    list(location=factor_location, scale=factor_scale),
                    weights * 2 * q * residual / sigma^3,
                    objective, value,
                    workers=requested_moment_workers,
                    blas_threads=max(1L, cores %/% requested_moment_workers),
                    previous_nonlinear_norm=nonlinear_norm,
                    initial_solution=coupled_linear_solution,
                    preconditioner_mode=requested_preconditioner,
                    schur_rank=schur_rank,
                    progress=coupled_progress
                ),
                error=function(error) {
                    notify(
                        'coupled update failed',
                        iteration=mixed_iteration_offset + iteration,
                        reason=conditionMessage(error)
                    )
                    NULL
                }
            )
            if (!is.null(coupled_update)) {
                factor_location <- coupled_update$factors$location
                factor_scale <- coupled_update$factors$scale
                maximum_damping <- max(
                    maximum_damping, coupled_update$damping
                )
                step_norm <- max(abs(coupled_update$step))
                directional_derivative <-
                    coupled_update$directional_derivative
                best_line_search_change <- coupled_update$best_change
                linear_iterations <- coupled_update$linear_iterations
                linear_residual <- coupled_update$linear_residual
                linear_tolerance <- coupled_update$linear_tolerance
                linear_acceptance_tolerance <-
                    coupled_update$linear_acceptance_tolerance
                linear_converged <- coupled_update$linear_converged
                linear_warm_started <- coupled_update$linear_warm_started
                linear_initial_residual <-
                    coupled_update$linear_initial_residual
                linear_preconditioner <-
                    coupled_update$preconditioner_mode
                nonlinear_norm <- coupled_update$nonlinear_norm
                step_mode <- 'coupled'
                if (coupled_update$accepted) {
                    coupled_consecutive_failures <- 0L
                    coupled_retry_delay <- 0L
                    coupled_linear_solution <-
                        coupled_update$linear_solution
                    coefficients <- coupled_update$coefficients
                    predictors <- coupled_update$predictors
                    value <- coupled_update$value
                    relative <- .cdrgam_gaulss_relative_convergence(
                        coupled_value,
                        value,
                        directional_derivative,
                        relative_tolerance
                    )
                    relative_reduction <- relative$actual
                    relative_predicted_reduction <- relative$predicted
                    notify(
                        'iteration complete',
                        iteration=mixed_iteration_offset + iteration,
                        mode='coupled',
                        objective=format(value, digits=10),
                        gradient=format(gradient_norm, digits=5),
                        linear_iterations=linear_iterations,
                        linear_residual=format(linear_residual, digits=5),
                        linear_tolerance=format(linear_tolerance, digits=5),
                        linear_acceptance_tolerance=format(
                            linear_acceptance_tolerance, digits=5
                        ),
                        linear_converged=linear_converged,
                        linear_warm_start=linear_warm_started,
                        linear_preconditioner=linear_preconditioner,
                        linear_initial_residual=format(
                            linear_initial_residual, digits=5
                        ),
                        relative_reduction=format(
                            relative_reduction, digits=5
                        ),
                        relative_predicted_reduction=format(
                            relative_predicted_reduction, digits=5
                        )
                    )
                    if (isTRUE(relative$converged)) {
                        converged <- TRUE
                        termination <- 'relative objective tolerance'
                        break
                    }
                    next
                }
                if (coupled_update$numerically_converged) {
                    converged <- TRUE
                    termination <- 'numerical precision'
                    break
                }
            }
            coupled_failures <- coupled_failures + 1L
            coupled_consecutive_failures <-
                coupled_consecutive_failures + 1L
            coupled_linear_solution <- NULL
            if (identical(algorithm, 'schur')) {
                coupled_retry_delay <- .cdrgam_gaulss_retry_delay(
                    coupled_consecutive_failures
                )
            } else {
                coupled_disabled <- coupled_consecutive_failures >= 3L
            }
            notify(
                if (coupled_disabled) {
                    'coupled updates disabled'
                } else 'coupled update deferred',
                iteration=mixed_iteration_offset + iteration,
                failures=coupled_failures,
                consecutive_failures=coupled_consecutive_failures,
                retry_delay=coupled_retry_delay,
                reason=if (is.null(coupled_update)) {
                    'bounded Krylov probe failed'
                } else 'coupled line search rejected the probe direction'
            )
            coupled_failed <- TRUE
        }
        if (algorithm %in% c('alternating', 'mixed', 'schur')) {
            block_started <- proc.time()[['elapsed']]
            iteration_value <- value
            iteration_coefficients <- unlist(coefficients, use.names=FALSE)
            notify(
                'alternating update started',
                iteration=mixed_iteration_offset + iteration
            )
            location_update <- .cdrgam_gaulss_block_update(
                'location', assemblies$location, coefficients, predictors,
                gradients$location, moments$information$location,
                penalties$location, factor_location, objective, value
            )
            coefficients <- location_update$coefficients
            predictors <- location_update$predictors
            value <- location_update$value
            factor_location <- location_update$factor

            q <- exp(predictors$scale)
            sigma <- b + q
            residual <- y - predictors$location
            scale_moments <- .cdrgam_gaulss_sparse_moments(
                assemblies, coefficients, sigma, q, residual, weights,
                include_information=coupled_step,
                workers=resolved_moment_workers,
                blas_threads=max(1L, cores %/% resolved_moment_workers),
                parameters='scale'
            )
            scale_gradient <- scale_moments$gradients$scale +
                as.numeric(penalties$scale %*% coefficients$scale)
            scale_update <- .cdrgam_gaulss_block_update(
                'scale', assemblies$scale, coefficients, predictors,
                scale_gradient,
                if (coupled_step) {
                    scale_moments$information$scale
                } else moments$information$scale,
                penalties$scale, factor_scale, objective, value
            )
            coefficients <- scale_update$coefficients
            predictors <- scale_update$predictors
            value <- scale_update$value
            factor_scale <- scale_update$factor

            accepted <- c(
                location=location_update$accepted,
                scale=scale_update$accepted
            )
            maximum_damping <- max(
                maximum_damping,
                location_update$damping,
                scale_update$damping
            )
            accepted_steps <- c(location_update$step, scale_update$step)
            step_norm <- max(abs(accepted_steps))
            directional_derivative <-
                sum(gradients$location * location_update$step) +
                sum(scale_gradient * scale_update$step)
            best_line_search_change <- min(
                location_update$best_change,
                scale_update$best_change,
                na.rm=TRUE
            )
            accepted_names <- names(accepted)[accepted]
            step_mode <- if (length(accepted_names)) {
                paste0(
                    if (coupled_failed) 'coupled-fallback-' else {
                        'alternating-'
                    },
                    paste(accepted_names, collapse='-')
                )
            } else 'alternating-none'
            if (!any(accepted)) {
                numerical <- c(
                    location_update$numerically_converged,
                    scale_update$numerically_converged
                )
                if (all(numerical)) {
                    converged <- TRUE
                    termination <- 'numerical precision'
                } else {
                    termination <- 'alternating line search failed'
                }
                break
            }
            mapped_coefficients <- unlist(coefficients, use.names=FALSE)
            acceleration_eligible <- identical(algorithm, 'alternating') ||
                coupled_disabled
            if (acceleration_eligible) {
                anderson_states <- c(
                    anderson_states, list(iteration_coefficients)
                )
                anderson_residuals <- c(
                    anderson_residuals,
                    list(mapped_coefficients - iteration_coefficients)
                )
                retained <- tail(seq_along(anderson_states), 3L)
                anderson_states <- anderson_states[retained]
                anderson_residuals <- anderson_residuals[retained]
            } else {
                anderson_states <- list()
                anderson_residuals <- list()
            }
            global_iteration <- mixed_iteration_offset + iteration
            attempt_acceleration <- acceleration_eligible &&
                length(anderson_states) >= 3L &&
                global_iteration %% 3L == 0L
            if (attempt_acceleration) {
                anderson_attempts <- anderson_attempts + 1L
                candidate_vector <- .cdrgam_anderson_candidate(
                    anderson_states,
                    anderson_residuals,
                    dimensions,
                    depth=2L
                )
                if (!is.null(candidate_vector)) {
                    candidate_coefficients <- list(
                        location=candidate_vector[
                            seq_len(dimensions[['location']])
                        ],
                        scale=tail(candidate_vector, dimensions[['scale']])
                    )
                    candidate_predictors <- Map(
                        .cdrgam_distributional_sparse_predictor,
                        assemblies,
                        candidate_coefficients,
                        offsets
                    )
                    candidate_value <- objective(
                        candidate_coefficients, candidate_predictors
                    )
                    if (is.finite(candidate_value) &&
                            candidate_value < value) {
                        coefficients <- candidate_coefficients
                        predictors <- candidate_predictors
                        value <- candidate_value
                        best_line_search_change <- min(
                            best_line_search_change,
                            candidate_value - iteration_value
                        )
                        accepted_step <-
                            candidate_vector - iteration_coefficients
                        directional_derivative <- sum(
                            unlist(gradients, use.names=FALSE) * accepted_step
                        )
                        step_norm <- max(abs(accepted_step))
                        step_mode <- paste0(step_mode, '-anderson')
                        anderson_acceptances <- anderson_acceptances + 1L
                    }
                }
            }
            relative <- .cdrgam_gaulss_relative_convergence(
                iteration_value,
                value,
                directional_derivative,
                relative_tolerance
            )
            relative_reduction <- relative$actual
            relative_predicted_reduction <- relative$predicted
            notify(
                'iteration complete',
                iteration=mixed_iteration_offset + iteration,
                mode=step_mode,
                seconds=format(
                    proc.time()[['elapsed']] - block_started, digits=5
                ),
                objective=format(value, digits=10),
                gradient=format(gradient_norm, digits=5),
                relative_reduction=format(relative_reduction, digits=5),
                relative_predicted_reduction=format(
                    relative_predicted_reduction, digits=5
                ),
                anderson_attempts=anderson_attempts,
                anderson_acceptances=anderson_acceptances
            )
            if (isTRUE(relative$converged)) {
                converged <- TRUE
                termination <- 'relative objective tolerance'
                break
            }
            next
        }
        systems <- list(
            location=Matrix::forceSymmetric(
                moments$information$location + penalties$location, uplo='U'
            ),
            scale=Matrix::forceSymmetric(
                moments$information$scale + penalties$scale, uplo='U'
            )
        )
        location_factorization <- .cdrgam_gaulss_step_factor(
            systems$location,
            assemblies$location$supernodal,
            factor_location
        )
        scale_factorization <- .cdrgam_gaulss_step_factor(
            systems$scale,
            assemblies$scale$supernodal,
            factor_scale
        )
        factor_location <- location_factorization$factor
        factor_scale <- scale_factorization$factor
        maximum_damping <- max(
            maximum_damping,
            location_factorization$damping,
            scale_factorization$damping
        )
        steps <- list(
            location=-as.numeric(.cdr_factor_solve(
                factor_location, gradients$location
            )),
            scale=-as.numeric(.cdr_factor_solve(factor_scale, gradients$scale))
        )
        gradient_vector <- unlist(gradients, use.names=FALSE)
        step_vector <- unlist(steps, use.names=FALSE)
        step_norm <- max(abs(step_vector))
        directional_derivative <- sum(gradient_vector * step_vector)
        predictor_steps <- Map(function(assembly, step) {
            .cdrgam_distributional_sparse_predictor(
                assembly,
                step,
                numeric(assembly$observation_count)
            )
        }, assemblies, steps)
        step_size <- 1
        minimum_step <- 2^-20
        accepted <- FALSE
        step_mode <- 'joint'
        accepted_step <- NULL
        best_line_search_change <- Inf
        while (step_size >= minimum_step) {
            candidate <- Map(function(current, step) {
                current + step_size * step
            }, coefficients, steps)
            candidate_predictors <- Map(function(current, step) {
                current + step_size * step
            }, predictors, predictor_steps)
            candidate_value <- objective(candidate, candidate_predictors)
            best_line_search_change <- min(
                best_line_search_change, candidate_value - value,
                na.rm=TRUE
            )
            if (is.finite(candidate_value) && candidate_value < value) {
                coefficients <- candidate
                predictors <- candidate_predictors
                value <- candidate_value
                accepted <- TRUE
                accepted_step <- step_size * step_vector
                break
            }
            step_size <- step_size / 2
        }
        if (!accepted) {
            block_derivatives <- vapply(names(steps), function(parameter) {
                sum(gradients[[parameter]] * steps[[parameter]])
            }, numeric(1))
            for (parameter in names(sort(block_derivatives))) {
                block_step_size <- 1
                while (block_step_size >= minimum_step) {
                    candidate <- coefficients
                    candidate_predictors <- predictors
                    candidate[[parameter]] <- coefficients[[parameter]] +
                        block_step_size * steps[[parameter]]
                    candidate_predictors[[parameter]] <-
                        predictors[[parameter]] +
                        block_step_size * predictor_steps[[parameter]]
                    candidate_value <- objective(
                        candidate, candidate_predictors
                    )
                    best_line_search_change <- min(
                        best_line_search_change,
                        candidate_value - value,
                        na.rm=TRUE
                    )
                    if (is.finite(candidate_value) && candidate_value < value) {
                        coefficients <- candidate
                        predictors <- candidate_predictors
                        value <- candidate_value
                        accepted <- TRUE
                        step_mode <- paste0(parameter, '-only')
                        accepted_step <- numeric(length(step_vector))
                        parameter_offset <- if (identical(
                                parameter, 'location'
                            )) 0L else length(steps$location)
                        indices <- parameter_offset +
                            seq_along(steps[[parameter]])
                        accepted_step[indices] <-
                            block_step_size * steps[[parameter]]
                        break
                    }
                    block_step_size <- block_step_size / 2
                }
                if (accepted) break
            }
        }
        if (accepted && max(abs(accepted_step)) < tolerance / 100) {
            converged <- TRUE
            termination <- 'step tolerance'
            break
        }
        if (!accepted) {
            if (.cdrgam_gaulss_numerically_converged(
                    directional_derivative,
                    best_line_search_change,
                    value,
                    minimum_step
                )) {
                converged <- TRUE
                termination <- 'numerical precision'
                break
            }
            termination <- 'line search failed'
            break
        }
    }
    q <- exp(predictors$scale)
    sigma <- b + q
    residual <- y - predictors$location
    notify(
        'final Hessian moments started',
        iterations=iteration,
        termination=termination
    )
    final_started <- proc.time()[['elapsed']]
    observed <- .cdrgam_gaulss_sparse_moments(
        assemblies, coefficients, sigma, q, residual, weights,
        include_observed=TRUE,
        workers=resolved_moment_workers,
        blas_threads=max(1L, cores %/% resolved_moment_workers)
    )
    notify(
        'final Hessian moments complete',
        seconds=format(proc.time()[['elapsed']] - final_started, digits=5),
        chunks=observed$chunks
    )
    likelihood_hessian <- rbind(
        cbind(observed$information$location, observed$information$cross),
        cbind(
            Matrix::t(observed$information$cross),
            observed$information$scale
        )
    )
    penalty <- Matrix::bdiag(penalties$location, penalties$scale)
    hessian <- Matrix::forceSymmetric(
        likelihood_hessian + penalty, uplo='U'
    )
    joint_factorization <- .cdrgam_gaulss_step_factor(
        hessian,
        assemblies$location$supernodal || assemblies$scale$supernodal,
        factor_joint
    )
    joint_factor <- joint_factorization$factor
    maximum_damping <- max(maximum_damping, joint_factorization$damping)
    notify(
        'solve complete',
        iterations=iteration,
        converged=converged,
        termination=termination,
        objective=format(value, digits=10),
        gradient=format(gradient_norm, digits=5)
    )
    list(
        coefficients=unlist(coefficients, use.names=FALSE),
        parameter_coefficients=coefficients,
        linear_predictors=do.call(cbind, predictors),
        fitted_values=cbind(location=predictors$location, scale=1 / sigma),
        sigma=sigma,
        residuals=residual,
        objective=value,
        likelihood_hessian=likelihood_hessian,
        hessian=hessian,
        factor=joint_factor,
        penalty=penalty,
        parameter_penalties=penalties,
        inner_factors=list(
            location=factor_location,
            scale=factor_scale
        ),
        converged=converged,
        iterations=iteration,
        gradient_norm=gradient_norm,
        termination=termination,
        step_norm=step_norm,
        directional_derivative=directional_derivative,
        best_line_search_change=best_line_search_change,
        step_mode=step_mode,
        linear_iterations=linear_iterations,
        linear_residual=linear_residual,
        linear_tolerance=linear_tolerance,
        linear_acceptance_tolerance=linear_acceptance_tolerance,
        linear_converged=linear_converged,
        linear_warm_started=linear_warm_started,
        linear_initial_residual=linear_initial_residual,
        linear_preconditioner=linear_preconditioner,
        coupled_linear_solution=coupled_linear_solution,
        relative_reduction=relative_reduction,
        relative_predicted_reduction=relative_predicted_reduction,
        relative_tolerance=relative_tolerance,
        anderson_attempts=anderson_attempts,
        anderson_acceptances=anderson_acceptances,
        anderson_states=anderson_states,
        anderson_residuals=anderson_residuals,
        nonlinear_norm=nonlinear_norm,
        maximum_damping=maximum_damping,
        prior_weights=weights,
        chunks=observed$chunks,
        moment_workers=resolved_moment_workers,
        algorithm=algorithm,
        coupled_disabled=coupled_disabled,
        coupled_failures=coupled_failures,
        coupled_consecutive_failures=coupled_consecutive_failures,
        coupled_retry_delay=coupled_retry_delay,
        mixed_iterations=if (algorithm %in% c('mixed', 'schur')) {
            mixed_iteration_offset + iteration
        } else 0L,
        warm_started=valid_initial
    )
}

.cdrgam_sparse_row_inner <- function(design, coefficient_by_row) {
    coefficient_by_row <- as.matrix(coefficient_by_row)
    if (ncol(design) != nrow(coefficient_by_row) ||
            nrow(design) != ncol(coefficient_by_row)) {
        stop('Row-inner-product dimensions do not match')
    }
    if (!methods::is(design, 'sparseMatrix')) {
        return(rowSums(design * t(coefficient_by_row)))
    }
    entries <- Matrix::summary(design)
    if (!nrow(entries)) return(numeric(nrow(design)))
    values <- entries$x * coefficient_by_row[cbind(
        entries$j, entries$i
    )]
    summed <- rowsum(values, entries$i, reorder=FALSE)
    output <- numeric(nrow(design))
    output[as.integer(rownames(summed))] <- summed[, 1L]
    output
}

.cdrgam_gaulss_sparse_qncv <- function(
        assemblies, solution, gamma=1, batch_size=NULL,
        coefficient_derivatives=NULL, hessian_derivatives=NULL,
        workers=1L, blas_threads=1L, resume=NULL, checkpoint=NULL
) {
    if (length(gamma) != 1L || !is.numeric(gamma) ||
            !is.finite(gamma) || gamma <= 0) {
        stop('gamma must be one positive finite number')
    }
    dimensions <- vapply(assemblies, `[[`, integer(1), 'dimension')
    ranges <- split(seq_len(sum(dimensions)), rep(names(dimensions), dimensions))
    derivative_count <- if (is.null(coefficient_derivatives)) 0L else {
        coefficient_derivatives <- as.matrix(coefficient_derivatives)
        ncol(coefficient_derivatives)
    }
    if (derivative_count && (
            nrow(coefficient_derivatives) != sum(dimensions) ||
            length(hessian_derivatives) != derivative_count)) {
        stop('QNCV derivatives do not match the joint sparse system')
    }
    observation_count <- assemblies$location$observation_count
    if (is.null(batch_size)) batch_size <- min(256L, observation_count)
    batch_size <- min(as.integer(batch_size), observation_count)
    if (!is.finite(batch_size) || batch_size < 1L) {
        stop('QNCV batch size must be positive')
    }
    q <- exp(solution$linear_predictors[, 'scale'])
    sigma <- solution$sigma
    residual <- solution$residuals
    weights <- solution$prior_weights
    location_score <- weights * residual / sigma^2
    scale_score <- weights * q * (
        -1 / sigma + residual^2 / sigma^3
    )
    location_curvature <- weights / sigma^2
    cross_curvature <- weights * 2 * q * residual / sigma^3
    scale_curvature <- weights * (
        q * (1 / sigma - residual^2 / sigma^3) +
            q^2 * (-1 / sigma^2 + 3 * residual^2 / sigma^4)
    )
    starts <- seq.int(1L, observation_count, by=batch_size)
    workers <- min(
        .cdrgam_positive_integer(workers, 'QNCV worker count'),
        length(starts)
    )
    resume_valid <- is.list(resume) &&
        identical(resume$batch_size, batch_size) &&
        identical(resume$observation_count, observation_count) &&
        identical(resume$derivative_count, derivative_count) &&
        is.list(resume$results)
    results <- if (resume_valid) resume$results else list()
    completed_starts <- unlist(lapply(
        results, function(result) result$starts
    ), use.names=FALSE)
    pending_starts <- setdiff(starts, completed_starts)
    evaluate_batches <- function(local_starts) {
        local_rows <- integer()
        local_delta_location <- numeric()
        local_delta_scale <- numeric()
        local_score <- numeric(derivative_count)
        local_minimum_determinant <- Inf
        for (start in local_starts) {
        rows <- start:min(observation_count, start + batch_size - 1L)
        count <- length(rows)
        X_location <- .cdrgam_sparse_design_chunk(
            assemblies$location, rows
        )
        X_scale <- .cdrgam_sparse_design_chunk(assemblies$scale, rows)
        zero_location <- Matrix::Matrix(
            0, dimensions[['location']], count, sparse=TRUE
        )
        zero_scale <- Matrix::Matrix(
            0, dimensions[['scale']], count, sparse=TRUE
        )
        right_hand_sides <- rbind(
            cbind(Matrix::t(X_location), zero_location),
            cbind(zero_scale, Matrix::t(X_scale))
        )
        inverse_rows <- as.matrix(.cdr_factor_solve(
            solution$factor, right_hand_sides
        ))
        location_columns <- seq_len(count)
        scale_columns <- count + seq_len(count)
        a11 <- .cdrgam_sparse_row_inner(
            X_location,
            inverse_rows[ranges$location, location_columns, drop=FALSE]
        )
        a12_location <- .cdrgam_sparse_row_inner(
            X_location,
            inverse_rows[ranges$location, scale_columns, drop=FALSE]
        )
        a12_scale <- .cdrgam_sparse_row_inner(
            X_scale,
            inverse_rows[ranges$scale, location_columns, drop=FALSE]
        )
        a12 <- (a12_location + a12_scale) / 2
        a22 <- .cdrgam_sparse_row_inner(
            X_scale,
            inverse_rows[ranges$scale, scale_columns, drop=FALSE]
        )
        w11 <- location_curvature[rows]
        w12 <- cross_curvature[rows]
        w22 <- scale_curvature[rows]
        c11 <- 1 - a11 * w11 - a12 * w12
        c12 <- -a11 * w12 - a12 * w22
        c21 <- -a12 * w11 - a22 * w12
        c22 <- 1 - a12 * w12 - a22 * w22
        determinant <- c11 * c22 - c12 * c21
        local_minimum_determinant <- min(
            local_minimum_determinant, determinant
        )
        if (any(!is.finite(determinant)) || any(abs(determinant) < 1e-10)) {
            stop('QNCV leave-out system is numerically singular')
        }
        v1 <- a11 * location_score[rows] + a12 * scale_score[rows]
        v2 <- a12 * location_score[rows] + a22 * scale_score[rows]
        delta1 <- -(c22 * v1 - c12 * v2) / determinant
        delta2 <- -(-c21 * v1 + c11 * v2) / determinant
        local_rows <- c(local_rows, rows)
        local_delta_location <- c(local_delta_location, delta1)
        local_delta_scale <- c(local_delta_scale, delta2)
        if (derivative_count) {
            predictor_derivatives <- list(
                location=as.matrix(X_location %*%
                    coefficient_derivatives[
                        ranges$location, , drop=FALSE
                    ]),
                scale=as.matrix(X_scale %*%
                    coefficient_derivatives[ranges$scale, , drop=FALSE])
            )
            z_location <- inverse_rows[, location_columns, drop=FALSE]
            z_scale <- inverse_rows[, scale_columns, drop=FALSE]
            l1 <- location_score[rows]
            l2 <- scale_score[rows]
            for (index in seq_len(derivative_count)) {
                hessian_product <- as.matrix(
                    hessian_derivatives[[index]] %*% inverse_rows
                )
                product_location <- hessian_product[
                    , location_columns, drop=FALSE
                ]
                product_scale <- hessian_product[
                    , scale_columns, drop=FALSE
                ]
                da11 <- -colSums(z_location * product_location)
                da12 <- -0.5 * (
                    colSums(z_location * product_scale) +
                        colSums(z_scale * product_location)
                )
                da22 <- -colSums(z_scale * product_scale)
                deta1 <- predictor_derivatives$location[, index]
                deta2 <- predictor_derivatives$scale[, index]
                dw11 <- -2 * weights[rows] * q[rows] / sigma[rows]^3 *
                    deta2
                dw12 <- -2 * weights[rows] * q[rows] / sigma[rows]^3 *
                    deta1 + 2 * weights[rows] * (
                        q[rows] * residual[rows] / sigma[rows]^3 -
                            3 * q[rows]^2 * residual[rows] /
                                sigma[rows]^4
                    ) * deta2
                scale_residual_weight <- weights[rows] * (
                    -2 * q[rows] * residual[rows] / sigma[rows]^3 +
                        6 * q[rows]^2 * residual[rows] / sigma[rows]^4
                )
                scale_q_weight <- weights[rows] * (
                    1 / sigma[rows] - 3 * q[rows] / sigma[rows]^2 -
                        residual[rows]^2 / sigma[rows]^3 +
                        2 * q[rows]^2 / sigma[rows]^3 +
                        9 * q[rows] * residual[rows]^2 / sigma[rows]^4 -
                        12 * q[rows]^2 * residual[rows]^2 /
                            sigma[rows]^5
                )
                dw22 <- -scale_residual_weight * deta1 +
                    scale_q_weight * q[rows] * deta2
                dl1 <- -w11 * deta1 - w12 * deta2
                dl2 <- -w12 * deta1 - w22 * deta2
                dc11 <- -(da11 * w11 + da12 * w12 +
                    a11 * dw11 + a12 * dw12)
                dc12 <- -(da11 * w12 + da12 * w22 +
                    a11 * dw12 + a12 * dw22)
                dc21 <- -(da12 * w11 + da22 * w12 +
                    a12 * dw11 + a22 * dw12)
                dc22 <- -(da12 * w12 + da22 * w22 +
                    a12 * dw12 + a22 * dw22)
                dv1 <- da11 * l1 + da12 * l2 + a11 * dl1 + a12 * dl2
                dv2 <- da12 * l1 + da22 * l2 + a12 * dl1 + a22 * dl2
                right1 <- dv1 + dc11 * delta1 + dc12 * delta2
                right2 <- dv2 + dc21 * delta1 + dc22 * delta2
                ddelta1 <- -(c22 * right1 - c12 * right2) / determinant
                ddelta2 <- -(-c21 * right1 + c11 * right2) / determinant
                full_derivative <- -(l1 * deta1 + l2 * deta2)
                correction_derivative <-
                    dl1 * delta1 + dl2 * delta2 +
                    l1 * ddelta1 + l2 * ddelta2 -
                    ddelta1 * (w11 * delta1 + w12 * delta2) -
                    ddelta2 * (w12 * delta1 + w22 * delta2) -
                    0.5 * (
                        dw11 * delta1^2 +
                            2 * dw12 * delta1 * delta2 +
                            dw22 * delta2^2
                    )
                local_score[[index]] <- local_score[[index]] + sum(
                    full_derivative - gamma * correction_derivative
                )
            }
        }
        }
        list(
            starts=local_starts,
            rows=local_rows,
            delta_location=local_delta_location,
            delta_scale=local_delta_scale,
            score=local_score,
            minimum_determinant=local_minimum_determinant
        )
    }
    wave_size <- max(workers, 8L * workers)
    waves <- split(
        pending_starts,
        ceiling(seq_along(pending_starts) / wave_size)
    )
    for (wave in waves) {
        wave_workers <- min(workers, length(wave))
        worker_starts <- split(
            wave,
            rep(seq_len(wave_workers), length.out=length(wave))
        )
        current <- if (wave_workers > 1L && .Platform$OS.type != 'windows') {
            .cdrgam_with_blas_threads(blas_threads, suppressWarnings(
                parallel::mclapply(
                    worker_starts,
                    evaluate_batches,
                    mc.cores=wave_workers,
                    mc.preschedule=TRUE,
                    mc.set.seed=FALSE
                )
            ))
        } else {
            lapply(worker_starts, evaluate_batches)
        }
        valid <- length(current) == length(worker_starts) && all(vapply(
            current,
            function(result) is.list(result) &&
                is.numeric(result$minimum_determinant),
            logical(1)
        ))
        if (!valid) stop('A parallel QNCV response batch failed')
        results <- c(results, current)
        if (is.function(checkpoint)) checkpoint(list(
            batch_size=batch_size,
            observation_count=observation_count,
            derivative_count=derivative_count,
            results=results
        ))
    }
    valid <- length(results) > 0L && all(vapply(
        results,
        function(result) is.list(result) &&
            is.numeric(result$minimum_determinant),
        logical(1)
    ))
    if (!valid) stop('A parallel QNCV response batch failed')
    delta_location <- numeric(observation_count)
    delta_scale <- numeric(observation_count)
    score <- numeric(derivative_count)
    for (result in results) {
        delta_location[result$rows] <- result$delta_location
        delta_scale[result$rows] <- result$delta_scale
        score <- score + result$score
    }
    minimum_determinant <- min(vapply(
        results, `[[`, numeric(1), 'minimum_determinant'
    ))
    log_likelihood <- -weights * (
        log(sigma) + 0.5 * (residual / sigma)^2 + 0.5 * log(2 * pi)
    )
    linear_change <- location_score * delta_location +
        scale_score * delta_scale
    quadratic_change <- -(
        location_curvature * delta_location^2 +
            2 * cross_curvature * delta_location * delta_scale +
            scale_curvature * delta_scale^2
    )
    criterion <- -sum(log_likelihood) - gamma * sum(
        linear_change + 0.5 * quadratic_change
    )
    list(
        criterion=criterion,
        linear_predictors=solution$linear_predictors + cbind(
            location=delta_location,
            scale=delta_scale
        ),
        predictor_change=cbind(
            location=delta_location,
            scale=delta_scale
        ),
        score=if (derivative_count) score else NULL,
        minimum_determinant=minimum_determinant,
        batch_size=batch_size,
        workers=workers
    )
}

.cdrgam_gaulss_qncv_direction_plan <- function(
        assemblies, direction_count, cores, workers,
        response_batch_size=256L,
        derivative_workspace_bytes=0,
        memory=.cdrgam_memory_availability()
) {
    direction_count <- .cdrgam_positive_integer(
        direction_count, 'QNCV direction count'
    )
    plan <- .cdrgam_parallel_plan(cores, direction_count, workers)
    requested_workers <- plan$workers
    observation_count <- assemblies$location$observation_count
    sample_count <- min(
        observation_count,
        assemblies$location$crossprod_chunk_size,
        assemblies$scale$crossprod_chunk_size
    )
    sample_rows <- seq_len(sample_count)
    sample <- lapply(assemblies, function(assembly) {
        .cdrgam_sparse_design_chunk(assembly, sample_rows)
    })
    projected_design_bytes <- sum(vapply(
        sample, utils::object.size, numeric(1)
    )) * observation_count / sample_count
    direction_bytes <- 16 * as.double(observation_count) * direction_count
    dimension <- sum(vapply(assemblies, `[[`, integer(1), 'dimension'))
    qncv_worker_bytes <- 80 * as.double(dimension) * response_batch_size
    derivative_workspace_bytes <- max(
        0, as.double(derivative_workspace_bytes)
    )
    rm(sample)
    memory_budget_bytes <- if (is.finite(memory$available_bytes)) {
        0.25 * memory$available_bytes
    } else NA_real_
    memory_workers <- if (is.finite(memory_budget_bytes) &&
            projected_design_bytes > 0) {
        max(0L, floor(
            (memory_budget_bytes - direction_bytes) /
                projected_design_bytes - 1
        ))
    } else 0L
    resolved_workers <- min(requested_workers, memory_workers)
    cache_parallel <- resolved_workers > 1L &&
        .Platform$OS.type != 'windows'
    if (!cache_parallel) resolved_workers <- 1L
    plan$workers <- as.integer(resolved_workers)
    plan$blas_threads <- max(1L, plan$cores %/% plan$workers)
    response_budget_bytes <- if (is.finite(memory$available_bytes)) {
        0.15 * memory$available_bytes
    } else NA_real_
    response_memory_workers <- if (is.finite(response_budget_bytes) &&
            qncv_worker_bytes > 0) {
        max(0L, floor(
            (response_budget_bytes - derivative_workspace_bytes) /
                qncv_worker_bytes
        ))
    } else 0L
    plan$response_workers <- as.integer(max(
        1L, min(requested_workers, response_memory_workers)
    ))
    if (.Platform$OS.type == 'windows') plan$response_workers <- 1L
    plan$response_blas_threads <- max(
        1L, plan$cores %/% plan$response_workers
    )
    plan$requested_workers <- as.integer(requested_workers)
    plan$memory_workers <- as.integer(memory_workers)
    plan$memory_limited <- plan$workers < requested_workers
    plan$memory_source <- memory$source
    plan$memory_available_bytes <- memory$available_bytes
    plan$memory_budget_bytes <- memory_budget_bytes
    plan$projected_design_bytes <- as.numeric(projected_design_bytes)
    plan$direction_bytes <- as.numeric(direction_bytes)
    plan$qncv_worker_bytes <- as.numeric(qncv_worker_bytes)
    plan$derivative_workspace_bytes <- derivative_workspace_bytes
    plan$response_memory_budget_bytes <- response_budget_bytes
    plan$response_memory_workers <- as.integer(response_memory_workers)
    plan$cache_parallel <- cache_parallel
    plan
}

.cdrgam_gaulss_sparse_qncv_score <- function(
        assemblies, solution, sp, family, gamma=1, batch_size=NULL,
        cores=1L, workers=1L, score_batch_size=NULL,
        resume=NULL, checkpoint=NULL
) {
    dimensions <- vapply(assemblies, `[[`, integer(1), 'dimension')
    ranges <- split(seq_len(sum(dimensions)), rep(names(dimensions), dimensions))
    counts <- vapply(
        assemblies,
        function(assembly) length(assembly$penalty_components),
        integer(1)
    )
    tasks <- list()
    offset <- 0L
    for (parameter in names(assemblies)) {
        for (local_index in seq_len(counts[[parameter]])) {
            global_index <- offset + local_index
            tasks[[global_index]] <- list(
                parameter=parameter,
                local_index=local_index,
                global_index=global_index
            )
        }
        offset <- offset + counts[[parameter]]
    }
    if (!length(tasks)) return(numeric())
    derivative_bytes <- max(
        1,
        32 * as.double(Matrix::nnzero(solution$hessian)) +
            16 * as.double(assemblies$location$observation_count)
    )
    plan <- .cdrgam_sparse_score_batch_plan(
        length(tasks), workers, derivative_bytes, score_batch_size,
        memory_fraction=0.25
    )
    direction_plan <- .cdrgam_gaulss_qncv_direction_plan(
        assemblies,
        direction_count=plan$batch_size,
        cores=cores,
        workers=plan$workers,
        response_batch_size=if (is.null(batch_size)) 256L else batch_size,
        derivative_workspace_bytes=derivative_bytes * plan$batch_size
    )
    plan$workers <- direction_plan$workers
    resume_valid <- is.list(resume) &&
        identical(resume$task_count, length(tasks)) &&
        identical(resume$groups, plan$groups) &&
        is.list(resume$results)
    results <- if (resume_valid) resume$results else list()
    active <- if (resume_valid) resume$active else NULL
    evaluate <- function(indices, qncv_resume=NULL, qncv_checkpoint=NULL) tryCatch({
        started <- proc.time()[['elapsed']]
        selected <- tasks[indices]
        penalties <- lapply(selected, function(task) {
            .cdrgam_distributional_global_penalty(
                assemblies, ranges, task$parameter, task$local_index,
                sp[[task$global_index]]
            )
        })
        right_hand_sides <- do.call(cbind, lapply(
            penalties,
            function(penalty) as.numeric(
                penalty %*% solution$coefficients
            )
        ))
        coefficient_derivatives <- -as.matrix(.cdr_factor_solve(
            solution$factor, right_hand_sides
        ))
        likelihood_derivatives <-
            .cdrgam_gaulss_sparse_hessian_directions(
                assemblies,
                solution,
                coefficient_derivatives,
                .cdrgam_gaulss_b(family),
                workers=plan$workers,
                cache_parallel=direction_plan$cache_parallel
            )
        hessian_derivatives <- Map(function(penalty, likelihood) {
            Matrix::forceSymmetric(penalty + likelihood, uplo='U')
        }, penalties, likelihood_derivatives)
        evaluated <- .cdrgam_gaulss_sparse_qncv(
            assemblies,
            solution,
            gamma=gamma,
            batch_size=batch_size,
            coefficient_derivatives=coefficient_derivatives,
            hessian_derivatives=hessian_derivatives,
            workers=direction_plan$response_workers,
            blas_threads=direction_plan$response_blas_threads,
            resume=qncv_resume,
            checkpoint=qncv_checkpoint
        )
        list(
            indices=indices,
            score=evaluated$score,
            criterion=evaluated$criterion,
            minimum_determinant=evaluated$minimum_determinant,
            seconds=proc.time()[['elapsed']] - started,
            error=NULL
        )
    }, error=function(error) list(
        indices=indices,
        score=NULL,
        criterion=NA_real_,
        minimum_determinant=NA_real_,
        seconds=NA_real_,
        error=conditionMessage(error)
    ))
    completed <- unlist(lapply(
        results, function(result) result$indices
    ), use.names=FALSE)
    for (indices in plan$groups) {
        if (all(indices %in% completed)) next
        same_active <- is.list(active) && identical(active$indices, indices)
        qncv_checkpoint <- if (is.function(checkpoint)) {
            function(state) {
                active <<- list(indices=indices, qncv=state)
                checkpoint(list(
                    task_count=length(tasks), groups=plan$groups,
                    results=results, active=active
                ))
            }
        } else NULL
        result <- .cdrgam_with_blas_threads(
            direction_plan$blas_threads,
            evaluate(
                indices,
                qncv_resume=if (same_active) active$qncv else NULL,
                qncv_checkpoint=qncv_checkpoint
            )
        )
        if (!is.null(result$error)) stop(
            'Exact QNCV score batch ',
            paste(result$indices, collapse=', '),
            ' failed: ', result$error
        )
        results[[length(results) + 1L]] <- result
        completed <- c(completed, result$indices)
        active <- NULL
        if (is.function(checkpoint)) checkpoint(list(
            task_count=length(tasks), groups=plan$groups,
            results=results, active=NULL
        ))
    }
    failed <- which(vapply(
        results, function(result) !is.null(result$error), logical(1)
    ))
    if (length(failed)) {
        result <- results[[failed[[1L]]]]
        stop(
            'Exact QNCV score batch ',
            paste(result$indices, collapse=', '),
            ' failed: ', result$error
        )
    }
    output <- numeric(length(tasks))
    for (result in results) output[result$indices] <- result$score
    attr(output, 'criterion') <- results[[1L]]$criterion
    attr(output, 'minimum_determinant') <- min(vapply(
        results, `[[`, numeric(1), 'minimum_determinant'
    ))
    attr(output, 'score_plan') <- list(
        workers=plan$workers,
        response_workers=direction_plan$response_workers,
        requested_workers=direction_plan$requested_workers,
        blas_threads=direction_plan$blas_threads,
        response_blas_threads=direction_plan$response_blas_threads,
        memory_limited=direction_plan$memory_limited,
        design_cache=direction_plan$cache_parallel,
        parallel_axis='QNCV Hessian directions and response batches',
        batches=length(plan$groups),
        batch_size=plan$batch_size,
        response_batch_size=batch_size,
        derivative_bytes=plan$derivative_bytes,
        memory_source=plan$memory$source,
        memory_available_bytes=plan$memory$available_bytes,
        memory_budget_bytes=direction_plan$memory_budget_bytes,
        projected_design_bytes=direction_plan$projected_design_bytes,
        direction_bytes=direction_plan$direction_bytes,
        qncv_worker_bytes=direction_plan$qncv_worker_bytes,
        derivative_workspace_bytes=
            direction_plan$derivative_workspace_bytes,
        response_memory_budget_bytes=
            direction_plan$response_memory_budget_bytes,
        batch_seconds=vapply(results, `[[`, numeric(1), 'seconds')
    )
    output
}

.cdrgam_distributional_global_penalty <- function(
        assemblies, ranges, parameter, index, sp
) {
    dimension <- sum(vapply(assemblies, `[[`, integer(1), 'dimension'))
    derivative <- Matrix::Matrix(
        0, dimension, dimension, sparse=TRUE
    )
    range <- ranges[[parameter]]
    derivative[range, range] <- sp *
        assemblies[[parameter]]$penalty_components[[index]]
    Matrix::forceSymmetric(derivative, uplo='U')
}

.cdrgam_gaulss_sparse_hessian_direction <- function(
        assemblies, solution, coefficient_derivative, b
) {
    dimensions <- vapply(assemblies, `[[`, integer(1), 'dimension')
    ranges <- split(seq_len(sum(dimensions)), rep(names(dimensions), dimensions))
    direction <- list(
        location=.cdrgam_distributional_sparse_predictor(
            assemblies$location,
            coefficient_derivative[ranges$location],
            numeric(assemblies$location$observation_count)
        ),
        scale=.cdrgam_distributional_sparse_predictor(
            assemblies$scale,
            coefficient_derivative[ranges$scale],
            numeric(assemblies$scale$observation_count)
        )
    )
    q <- exp(solution$linear_predictors[, 'scale'])
    sigma <- b + q
    residual <- solution$residuals
    weights <- solution$prior_weights
    location <- NULL
    cross <- NULL
    scale <- NULL
    chunk_size <- min(
        assemblies$location$crossprod_chunk_size,
        assemblies$scale$crossprod_chunk_size
    )
    for (start in seq.int(1L, length(residual), by=chunk_size)) {
        rows <- start:min(length(residual), start + chunk_size - 1L)
        X_location <- .cdrgam_sparse_design_chunk(assemblies$location, rows)
        X_scale <- .cdrgam_sparse_design_chunk(assemblies$scale, rows)
        q_rows <- q[rows]
        sigma_rows <- sigma[rows]
        residual_rows <- residual[rows]
        weight_rows <- weights[rows]
        d_location <- direction$location[rows]
        d_scale <- direction$scale[rows]
        location_scale_weight <-
            -2 * weight_rows * q_rows / sigma_rows^3
        cross_location_weight <-
            -2 * weight_rows * q_rows / sigma_rows^3
        cross_scale_weight <- 2 * weight_rows * (
            q_rows * residual_rows / sigma_rows^3 -
                3 * q_rows^2 * residual_rows / sigma_rows^4
        )
        scale_residual_weight <- weight_rows * (
            -2 * q_rows * residual_rows / sigma_rows^3 +
                6 * q_rows^2 * residual_rows / sigma_rows^4
        )
        scale_q_weight <- weight_rows * (
            1 / sigma_rows - 3 * q_rows / sigma_rows^2 -
                residual_rows^2 / sigma_rows^3 +
                    2 * q_rows^2 / sigma_rows^3 +
                    9 * q_rows * residual_rows^2 / sigma_rows^4 -
                    12 * q_rows^2 * residual_rows^2 / sigma_rows^5
        )
        d_location_weight <- location_scale_weight * d_scale
        d_cross_weight <- cross_location_weight * d_location +
            cross_scale_weight * d_scale
        d_scale_weight <- -scale_residual_weight * d_location +
            scale_q_weight * q_rows * d_scale
        location_part <- Matrix::crossprod(
            X_location,
            Matrix::Diagonal(x=d_location_weight) %*% X_location
        )
        cross_part <- Matrix::crossprod(
            X_location,
            Matrix::Diagonal(x=d_cross_weight) %*% X_scale
        )
        scale_part <- Matrix::crossprod(
            X_scale,
            Matrix::Diagonal(x=d_scale_weight) %*% X_scale
        )
        location <- if (is.null(location)) {
            location_part
        } else location + location_part
        cross <- if (is.null(cross)) cross_part else cross + cross_part
        scale <- if (is.null(scale)) scale_part else scale + scale_part
    }
    Matrix::forceSymmetric(rbind(
        cbind(location, cross),
        cbind(Matrix::t(cross), scale)
    ), uplo='U')
}

.cdrgam_gaulss_sparse_hessian_directions <- function(
        assemblies, solution, coefficient_derivatives, b, workers=1L,
        cache_parallel=NULL
) {
    coefficient_derivatives <- as.matrix(coefficient_derivatives)
    dimensions <- vapply(assemblies, `[[`, integer(1), 'dimension')
    if (nrow(coefficient_derivatives) != sum(dimensions)) {
        stop('Coefficient derivatives do not match the joint sparse system')
    }
    ranges <- split(seq_len(sum(dimensions)), rep(names(dimensions), dimensions))
    coefficient_directions <- list(
        location=coefficient_derivatives[
            ranges$location, , drop=FALSE
        ],
        scale=coefficient_derivatives[ranges$scale, , drop=FALSE]
    )
    count <- ncol(coefficient_derivatives)
    workers <- min(as.integer(workers), max(1L, count))
    parallel <- workers > 1L && .Platform$OS.type != 'windows'
    cached_design <- NULL
    cache_bytes <- 0
    if (parallel) {
        observation_count <- assemblies$location$observation_count
        sample_count <- min(
            observation_count,
            assemblies$location$crossprod_chunk_size,
            assemblies$scale$crossprod_chunk_size
        )
        sample_rows <- seq_len(sample_count)
        sample <- lapply(assemblies, function(assembly) {
            .cdrgam_sparse_design_chunk(assembly, sample_rows)
        })
        projected_bytes <- sum(vapply(sample, utils::object.size, numeric(1))) *
            observation_count / sample_count
        direction_bytes <- 16 * as.double(observation_count) * count
        cache_safe <- if (is.null(cache_parallel)) {
            memory <- .cdrgam_memory_availability()
            projected_working_bytes <-
                projected_bytes * (workers + 1) + direction_bytes
            is.finite(memory$available_bytes) &&
                projected_working_bytes <= 0.25 * memory$available_bytes
        } else isTRUE(cache_parallel)
        if (cache_safe) {
            cached_design <- lapply(assemblies, function(assembly) {
                .cdrgam_sparse_design_chunk(
                    assembly, seq_len(assembly$observation_count)
                )
            })
            cache_bytes <- sum(vapply(
                cached_design, utils::object.size, numeric(1)
            ))
        }
        rm(sample)
    }
    location <- cross <- scale <- vector('list', count)
    q <- exp(solution$linear_predictors[, 'scale'])
    sigma <- b + q
    residual <- solution$residuals
    weights <- solution$prior_weights
    if (!is.null(cached_design)) {
        direction <- list(
            location=as.matrix(
                cached_design$location %*% coefficient_directions$location
            ),
            scale=as.matrix(
                cached_design$scale %*% coefficient_directions$scale
            )
        )
        location_scale_weight <- -2 * weights * q / sigma^3
        cross_location_weight <- location_scale_weight
        cross_scale_weight <- 2 * weights * (
            q * residual / sigma^3 - 3 * q^2 * residual / sigma^4
        )
        scale_residual_weight <- weights * (
            -2 * q * residual / sigma^3 +
                6 * q^2 * residual / sigma^4
        )
        scale_q_weight <- weights * (
            1 / sigma - 3 * q / sigma^2 - residual^2 / sigma^3 +
                2 * q^2 / sigma^3 + 9 * q * residual^2 / sigma^4 -
                12 * q^2 * residual^2 / sigma^5
        )
        evaluate <- function(index) {
            d_location <- direction$location[, index]
            d_scale <- direction$scale[, index]
            location_part <- Matrix::crossprod(
                cached_design$location,
                Matrix::Diagonal(
                    x=location_scale_weight * d_scale
                ) %*% cached_design$location
            )
            cross_part <- Matrix::crossprod(
                cached_design$location,
                Matrix::Diagonal(x=
                    cross_location_weight * d_location +
                        cross_scale_weight * d_scale
                ) %*% cached_design$scale
            )
            scale_part <- Matrix::crossprod(
                cached_design$scale,
                Matrix::Diagonal(x=
                    -scale_residual_weight * d_location +
                        scale_q_weight * q * d_scale
                ) %*% cached_design$scale
            )
            Matrix::forceSymmetric(rbind(
                cbind(location_part, cross_part),
                cbind(Matrix::t(cross_part), scale_part)
            ), uplo='U')
        }
        output <- suppressWarnings(parallel::mclapply(
            seq_len(count), evaluate, mc.cores=workers,
            mc.preschedule=TRUE, mc.set.seed=FALSE
        ))
        valid <- length(output) == count && all(vapply(
            output,
            function(value) methods::is(value, 'sparseMatrix'),
            logical(1)
        ))
        if (valid) {
            attr(output, 'direction_workers') <- workers
            attr(output, 'design_cache_bytes') <- cache_bytes
            return(output)
        }
        warning(
            'Parallel likelihood derivatives failed; retrying with streamed ',
            'serial derivatives',
            call.=FALSE
        )
    }
    chunk_size <- min(
        assemblies$location$crossprod_chunk_size,
        assemblies$scale$crossprod_chunk_size
    )
    for (start in seq.int(1L, length(residual), by=chunk_size)) {
        rows <- start:min(length(residual), start + chunk_size - 1L)
        X_location <- .cdrgam_sparse_design_chunk(assemblies$location, rows)
        X_scale <- .cdrgam_sparse_design_chunk(assemblies$scale, rows)
        direction <- list(
            location=as.matrix(X_location %*%
                coefficient_directions$location),
            scale=as.matrix(X_scale %*% coefficient_directions$scale)
        )
        q_rows <- q[rows]
        sigma_rows <- sigma[rows]
        residual_rows <- residual[rows]
        weight_rows <- weights[rows]
        location_scale_weight <-
            -2 * weight_rows * q_rows / sigma_rows^3
        cross_location_weight <-
            -2 * weight_rows * q_rows / sigma_rows^3
        cross_scale_weight <- 2 * weight_rows * (
            q_rows * residual_rows / sigma_rows^3 -
                3 * q_rows^2 * residual_rows / sigma_rows^4
        )
        scale_residual_weight <- weight_rows * (
            -2 * q_rows * residual_rows / sigma_rows^3 +
                6 * q_rows^2 * residual_rows / sigma_rows^4
        )
        scale_q_weight <- weight_rows * (
            1 / sigma_rows - 3 * q_rows / sigma_rows^2 -
                residual_rows^2 / sigma_rows^3 +
                    2 * q_rows^2 / sigma_rows^3 +
                    9 * q_rows * residual_rows^2 / sigma_rows^4 -
                    12 * q_rows^2 * residual_rows^2 / sigma_rows^5
        )
        for (index in seq_len(count)) {
            d_location <- direction$location[, index]
            d_scale <- direction$scale[, index]
            d_location_weight <- location_scale_weight * d_scale
            d_cross_weight <- cross_location_weight * d_location +
                cross_scale_weight * d_scale
            d_scale_weight <- -scale_residual_weight * d_location +
                scale_q_weight * q_rows * d_scale
            location_part <- Matrix::crossprod(
                X_location,
                Matrix::Diagonal(x=d_location_weight) %*% X_location
            )
            cross_part <- Matrix::crossprod(
                X_location,
                Matrix::Diagonal(x=d_cross_weight) %*% X_scale
            )
            scale_part <- Matrix::crossprod(
                X_scale,
                Matrix::Diagonal(x=d_scale_weight) %*% X_scale
            )
            location[[index]] <- if (is.null(location[[index]])) {
                location_part
            } else location[[index]] + location_part
            cross[[index]] <- if (is.null(cross[[index]])) {
                cross_part
            } else cross[[index]] + cross_part
            scale[[index]] <- if (is.null(scale[[index]])) {
                scale_part
            } else scale[[index]] + scale_part
        }
    }
    output <- lapply(seq_len(count), function(index) {
        Matrix::forceSymmetric(rbind(
            cbind(location[[index]], cross[[index]]),
            cbind(Matrix::t(cross[[index]]), scale[[index]])
        ), uplo='U')
    })
    attr(output, 'direction_workers') <- 1L
    attr(output, 'design_cache_bytes') <- 0
    output
}

.cdrgam_gaulss_sparse_score_reference <- function(
        assemblies, solution, sp, family, workers=1L
) {
    dimensions <- vapply(assemblies, `[[`, integer(1), 'dimension')
    ranges <- split(seq_len(sum(dimensions)), rep(names(dimensions), dimensions))
    counts <- vapply(
        assemblies, function(assembly) length(assembly$penalty_components),
        integer(1)
    )
    penalty_scores <- list()
    offset <- 0L
    for (parameter in names(assemblies)) {
        count <- counts[[parameter]]
        indices <- offset + seq_len(count)
        penalty_scores[[parameter]] <- .sparse_penalty_logdet_score(
            assemblies[[parameter]]$blocks,
            sp[indices],
            count
        )
        offset <- offset + count
    }
    beta <- solution$coefficients
    b <- .cdrgam_gaulss_b(family)
    global_index <- 0L
    tasks <- list()
    for (parameter in names(assemblies)) {
        for (local_index in seq_len(counts[[parameter]])) {
            global_index <- global_index + 1L
            tasks[[global_index]] <- list(
                parameter=parameter,
                local_index=local_index,
                global_index=global_index
            )
        }
    }
    evaluate <- function(task) {
        tryCatch({
            penalty_derivative <- .cdrgam_distributional_global_penalty(
                assemblies,
                ranges,
                task$parameter,
                task$local_index,
                sp[[task$global_index]]
            )
            coefficient_derivative <- -as.numeric(.cdr_factor_solve(
                solution$factor,
                penalty_derivative %*% beta
            ))
            likelihood_derivative <-
                .cdrgam_gaulss_sparse_hessian_direction(
                    assemblies, solution, coefficient_derivative, b
                )
            hessian_derivative <- Matrix::forceSymmetric(
                penalty_derivative + likelihood_derivative,
                uplo='U'
            )
            value <- as.numeric(Matrix::crossprod(
                beta, penalty_derivative %*% beta
            )) + .sparse_logdet_score(
                solution$factor, hessian_derivative, chunk_size=256L
            ) - penalty_scores[[task$parameter]][[task$local_index]]
            list(value=value, error=NULL)
        }, error=function(error) {
            list(value=NA_real_, error=conditionMessage(error))
        })
    }
    workers <- min(as.integer(workers), length(tasks))
    parallel <- workers > 1L && .Platform$OS.type != 'windows'
    results <- if (parallel) {
        parallel::mclapply(
            tasks, evaluate, mc.cores=workers, mc.preschedule=FALSE
        )
    } else lapply(tasks, evaluate)
    failed <- which(vapply(
        results, function(result) !is.null(result$error), logical(1)
    ))
    if (length(failed)) {
        first <- failed[[1L]]
        stop(
            'Exact distributional score component ', first,
            ' failed: ', results[[first]]$error
        )
    }
    vapply(results, `[[`, numeric(1), 'value')
}

.cdrgam_gaulss_sparse_score_batch <- function(
        assemblies, solution, sp, family, tasks, penalty_scores,
        trace_workers=1L
) {
    dimensions <- vapply(assemblies, `[[`, integer(1), 'dimension')
    ranges <- split(seq_len(sum(dimensions)), rep(names(dimensions), dimensions))
    beta <- solution$coefficients
    penalty_derivatives <- lapply(tasks, function(task) {
        .cdrgam_distributional_global_penalty(
            assemblies,
            ranges,
            task$parameter,
            task$local_index,
            sp[[task$global_index]]
        )
    })
    .cdrgam_sparse_exact_score_batch(
        solution$factor,
        beta,
        penalty_derivatives,
        vapply(tasks, function(task) {
            penalty_scores[[task$parameter]][[task$local_index]]
        }, numeric(1)),
        likelihood_derivatives=function(coefficient_derivatives) {
            .cdrgam_gaulss_sparse_hessian_directions(
                assemblies,
                solution,
                coefficient_derivatives,
                .cdrgam_gaulss_b(family),
                workers=trace_workers
            )
        },
        trace_workers=trace_workers
    )
}

.cdrgam_gaulss_sparse_score <- function(
        assemblies, solution, sp, family, workers=1L, batch_size=NULL
) {
    counts <- vapply(
        assemblies, function(assembly) length(assembly$penalty_components),
        integer(1)
    )
    penalty_scores <- list()
    offset <- 0L
    tasks <- list()
    for (parameter in names(assemblies)) {
        count <- counts[[parameter]]
        indices <- offset + seq_len(count)
        penalty_scores[[parameter]] <- .sparse_penalty_logdet_score(
            assemblies[[parameter]]$blocks,
            sp[indices],
            count
        )
        for (local_index in seq_len(count)) {
            global_index <- offset + local_index
            tasks[[global_index]] <- list(
                parameter=parameter,
                local_index=local_index,
                global_index=global_index
            )
        }
        offset <- offset + count
    }
    if (!length(tasks)) {
        output <- numeric()
        attr(output, 'score_plan') <- list(
            workers=0L,
            parallel_axis='not required',
            batches=0L,
            batch_size=0L,
            derivative_bytes=0,
            memory_source='not required',
            memory_available_bytes=NA_real_
        )
        return(output)
    }
    derivative_bytes <- max(
        1,
        32 * as.double(Matrix::nnzero(solution$hessian)) +
            16 * as.double(assemblies$location$observation_count)
    )
    plan <- .cdrgam_sparse_score_batch_plan(
        length(tasks), workers, derivative_bytes, batch_size
    )
    groups <- plan$groups
    workers <- plan$workers
    evaluate <- function(indices) tryCatch(
        list(
            indices=indices,
            values=.cdrgam_gaulss_sparse_score_batch(
                assemblies,
                solution,
                sp,
                family,
                tasks[indices],
                penalty_scores,
                trace_workers=workers
            ),
            error=NULL
        ),
        error=function(error) {
            list(indices=indices, values=NULL, error=conditionMessage(error))
        }
    )
    results <- lapply(groups, evaluate)
    missing <- which(!vapply(results, function(result) {
        is.list(result) &&
            all(c('indices', 'values', 'error') %in% names(result))
    }, logical(1)))
    if (length(missing)) {
        stop(
            'Exact distributional score worker ', missing[[1L]],
            ' did not return a result; it may have exceeded its memory limit'
        )
    }
    failed <- which(vapply(
        results, function(result) !is.null(result$error), logical(1)
    ))
    if (length(failed)) {
        first <- results[[failed[[1L]]]]
        stop(
            'Exact distributional score batch ',
            paste(first$indices, collapse=', '),
            ' failed: ', first$error
        )
    }
    output <- numeric(length(tasks))
    for (result in results) output[result$indices] <- result$values
    attr(output, 'score_plan') <- list(
        workers=workers,
        parallel_axis='likelihood directions and inverse columns',
        batches=length(groups),
        batch_size=plan$batch_size,
        derivative_bytes=plan$derivative_bytes,
        memory_source=plan$memory$source,
        memory_available_bytes=plan$memory$available_bytes,
        batch_timings=lapply(results, function(result) {
            attr(result$values, 'timing')
        })
    )
    output
}

.cdrgam_gaulss_sparse_stochastic_score <- function(
        assemblies, solution, sp, family, probes
) {
    probes <- as.matrix(probes)
    dimensions <- vapply(assemblies, `[[`, integer(1), 'dimension')
    dimension <- sum(dimensions)
    if (nrow(probes) != dimension || !ncol(probes)) {
        stop('Stochastic score probes do not match the joint sparse system')
    }
    ranges <- split(seq_len(dimension), rep(names(dimensions), dimensions))
    counts <- vapply(
        assemblies, function(assembly) length(assembly$penalty_components),
        integer(1)
    )
    tasks <- list()
    penalty_scores <- numeric()
    offset <- 0L
    for (parameter in names(assemblies)) {
        count <- counts[[parameter]]
        indices <- offset + seq_len(count)
        local_scores <- .sparse_penalty_logdet_score(
            assemblies[[parameter]]$blocks, sp[indices], count
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
    if (!length(tasks)) return(numeric())

    started <- proc.time()[['elapsed']]
    penalties <- lapply(tasks, function(task) {
        .cdrgam_distributional_global_penalty(
            assemblies, ranges, task$parameter, task$local_index,
            sp[[task$global_index]]
        )
    })
    right_hand_sides <- do.call(cbind, lapply(penalties, function(value) {
        as.numeric(value %*% solution$coefficients)
    }))
    coefficient_derivatives <- -as.matrix(.cdr_factor_solve(
        solution$factor, right_hand_sides
    ))
    coefficient_seconds <- proc.time()[['elapsed']] - started

    started <- proc.time()[['elapsed']]
    inverse_probes <- as.matrix(.cdr_factor_solve(solution$factor, probes))
    probe_solve_seconds <- proc.time()[['elapsed']] - started
    started <- proc.time()[['elapsed']]
    penalty_traces <- vapply(penalties, function(penalty) {
        mean(colSums(
            inverse_probes * as.matrix(penalty %*% probes)
        ))
    }, numeric(1))
    penalty_probe_seconds <- proc.time()[['elapsed']] - started

    location_directions <- coefficient_derivatives[
        ranges$location, , drop=FALSE
    ]
    scale_directions <- coefficient_derivatives[
        ranges$scale, , drop=FALSE
    ]
    likelihood_traces <- numeric(length(tasks))
    q <- exp(solution$linear_predictors[, 'scale'])
    b <- .cdrgam_gaulss_b(family)
    sigma <- b + q
    residual <- solution$residuals
    weights <- solution$prior_weights
    chunk_size <- min(
        assemblies$location$crossprod_chunk_size,
        assemblies$scale$crossprod_chunk_size
    )
    started <- proc.time()[['elapsed']]
    for (start in seq.int(1L, length(residual), by=chunk_size)) {
        rows <- start:min(length(residual), start + chunk_size - 1L)
        X_location <- .cdrgam_sparse_design_chunk(
            assemblies$location, rows
        )
        X_scale <- .cdrgam_sparse_design_chunk(assemblies$scale, rows)
        direction_location <- as.matrix(
            X_location %*% location_directions
        )
        direction_scale <- as.matrix(X_scale %*% scale_directions)
        probe_location <- as.matrix(
            X_location %*% probes[ranges$location, , drop=FALSE]
        )
        probe_scale <- as.matrix(
            X_scale %*% probes[ranges$scale, , drop=FALSE]
        )
        inverse_location <- as.matrix(
            X_location %*% inverse_probes[ranges$location, , drop=FALSE]
        )
        inverse_scale <- as.matrix(
            X_scale %*% inverse_probes[ranges$scale, , drop=FALSE]
        )
        location_leverage <- rowMeans(
            inverse_location * probe_location
        )
        cross_leverage <- rowMeans(
            inverse_location * probe_scale + inverse_scale * probe_location
        )
        scale_leverage <- rowMeans(inverse_scale * probe_scale)

        q_rows <- q[rows]
        sigma_rows <- sigma[rows]
        residual_rows <- residual[rows]
        weight_rows <- weights[rows]
        location_scale_weight <-
            -2 * weight_rows * q_rows / sigma_rows^3
        cross_location_weight <- location_scale_weight
        cross_scale_weight <- 2 * weight_rows * (
            q_rows * residual_rows / sigma_rows^3 -
                3 * q_rows^2 * residual_rows / sigma_rows^4
        )
        scale_residual_weight <- weight_rows * (
            -2 * q_rows * residual_rows / sigma_rows^3 +
                6 * q_rows^2 * residual_rows / sigma_rows^4
        )
        scale_q_weight <- weight_rows * (
            1 / sigma_rows - 3 * q_rows / sigma_rows^2 -
                residual_rows^2 / sigma_rows^3 +
                2 * q_rows^2 / sigma_rows^3 +
                9 * q_rows * residual_rows^2 / sigma_rows^4 -
                12 * q_rows^2 * residual_rows^2 / sigma_rows^5
        )
        for (index in seq_along(tasks)) {
            d_location <- direction_location[, index]
            d_scale <- direction_scale[, index]
            d_location_weight <- location_scale_weight * d_scale
            d_cross_weight <- cross_location_weight * d_location +
                cross_scale_weight * d_scale
            d_scale_weight <- -scale_residual_weight * d_location +
                scale_q_weight * q_rows * d_scale
            likelihood_traces[[index]] <- likelihood_traces[[index]] + sum(
                d_location_weight * location_leverage +
                    d_cross_weight * cross_leverage +
                    d_scale_weight * scale_leverage
            )
        }
    }
    likelihood_seconds <- proc.time()[['elapsed']] - started
    started <- proc.time()[['elapsed']]
    quadratic <- vapply(seq_along(penalties), function(index) {
        as.numeric(Matrix::crossprod(
            solution$coefficients,
            penalties[[index]] %*% solution$coefficients
        ))
    }, numeric(1))
    output <- quadratic + penalty_traces + likelihood_traces - penalty_scores
    attr(output, 'score_plan') <- list(
        method='stochastic matrix-free',
        probes=ncol(probes),
        timing=c(
            coefficient_derivatives=coefficient_seconds,
            probe_solves=probe_solve_seconds,
            penalty_probe_products=penalty_probe_seconds,
            likelihood_quadratic_forms=likelihood_seconds,
            quadratic=proc.time()[['elapsed']] - started
        )
    )
    output
}

.cdrgam_gaulss_cache_candidate <- function(
        parameters, evaluation, best_parameters=NULL, best_evaluation=NULL
) {
    same_as_best <- !is.null(best_parameters) &&
        identical(as.numeric(parameters), as.numeric(best_parameters))
    if (is.null(evaluation$solution) && same_as_best &&
            !is.null(best_evaluation$solution)) {
        return(best_evaluation)
    }
    evaluation
}

.cdrgam_gaulss_sparse_laml <- function(
        assemblies, family, control, reporter=NULL, checkpoint=NULL,
        checkpoint_signature=NULL
) {
    qncv <- identical(control$criterion, 'QNCV')
    counts <- vapply(
        assemblies, function(assembly) length(assembly$penalty_components),
        integer(1)
    )
    initial <- unlist(lapply(
        assemblies, .cdrgam_sparse_initial_sp
    ), use.names=FALSE)
    initial_log_sp <- log(initial)
    raw_checkpoint_state <- .checkpoint_read(checkpoint)
    failed_checkpoint <- is.list(raw_checkpoint_state) &&
        identical(raw_checkpoint_state$backend, 'distributional-sparse') &&
        identical(raw_checkpoint_state$stage, 'optimization') &&
        (length(raw_checkpoint_state$best_criterion) != 1L ||
            !is.finite(raw_checkpoint_state$best_criterion)) &&
        length(raw_checkpoint_state$current_criterion) == 1L &&
        is.finite(raw_checkpoint_state$current_criterion) &&
        raw_checkpoint_state$current_criterion >= 1e49
    checkpoint_state <- if (failed_checkpoint) NULL else {
        .checkpoint_validate(
            raw_checkpoint_state,
            checkpoint_signature,
            checkpoint
        )
    }
    if (!is.null(checkpoint_state) && isTRUE(checkpoint_state$legacy)) {
        stop(
            'Legacy checkpoints cannot resume a distributional fit. Supply ',
            'a different path or remove the stale checkpoint.',
            call.=FALSE
        )
    }
    resumed <- !is.null(checkpoint_state)
    if (is.null(checkpoint_state)) {
        checkpoint_state <- .new_checkpoint_state(
            checkpoint_signature, 'distributional-sparse', 0L
        )
        checkpoint_state$phase <- if (
            identical(control$gradient_method, 'hybrid')
        ) 'stochastic' else 'exact'
        checkpoint_state$stochastic_complete <- FALSE
        checkpoint_state$stochastic_warmup <- NULL
        checkpoint_state$optimizer_counts <- c(`function`=0L, gradient=0L)
        checkpoint_state$exact_restart <- 0L
    }
    b <- .cdrgam_gaulss_b(family)
    checkpoint_count <- function(name) {
        value <- checkpoint_state[[name]]
        if (is.null(value) || !is.finite(value)) 0L else as.integer(value)
    }
    evaluations <- checkpoint_count('evaluation_count')
    score_evaluations <- checkpoint_count('score_evaluations')
    stochastic_score_evaluations <- checkpoint_count(
        'stochastic_score_evaluations'
    )
    cached_parameters <- NULL
    cached_evaluation <- NULL
    best_parameters <- NULL
    best_evaluation <- NULL
    warm_solution <- NULL
    accepted_parameters <- NULL
    accepted_evaluation <- NULL
    warm_starts <- 0L
    cold_fallbacks <- 0L
    hybrid_warmup_active <- FALSE
    hybrid_rejected_solves <- 0L
    hybrid_rejected_seconds <- 0
    hybrid_reference_seconds <- numeric()
    hybrid_transition <- NULL
    score_plan <- NULL
    stochastic_score_plan <- NULL
    checkpoint_phase <- if (is.null(checkpoint_state$phase)) {
        if (identical(control$gradient_method, 'hybrid')) {
            'stochastic'
        } else 'exact'
    } else checkpoint_state$phase
    write_checkpoint <- function(
            parameters=NULL, evaluation=NULL, phase=checkpoint_phase,
            stage='optimization', optimization=NULL
    ) {
        if (is.null(checkpoint)) return(invisible(FALSE))
        checkpoint_state$stage <<- stage
        checkpoint_state$phase <<- phase
        checkpoint_state$evaluation_count <<- evaluations
        checkpoint_state$score_evaluations <<- score_evaluations
        checkpoint_state$stochastic_score_evaluations <<-
            stochastic_score_evaluations
        checkpoint_state$warm_starts <<- warm_starts
        checkpoint_state$cold_fallbacks <<- cold_fallbacks
        if (!is.null(parameters)) {
            checkpoint_state$current_log_sp <<- as.numeric(parameters)
        }
        if (!is.null(evaluation) && is.finite(evaluation$criterion)) {
            checkpoint_state$current_criterion <<- evaluation$criterion
        }
        if (!is.null(best_parameters) && !is.null(best_evaluation) &&
                is.finite(best_evaluation$criterion)) {
            checkpoint_state$best_log_sp <<- best_parameters
            checkpoint_state$best_criterion <<- best_evaluation$criterion
        }
        if (!is.null(optimization)) {
            checkpoint_state$optimization <<- optimization
        }
        checkpoint_state$updated_at <<- as.character(Sys.time())
        .checkpoint_write(checkpoint_state, checkpoint)
    }
    inflight_matches <- function(parameters) {
        state <- checkpoint_state$inflight
        is.list(state) && identical(
            as.numeric(state$parameters), as.numeric(parameters)
        )
    }
    write_inflight <- function(
            parameters, stage, inner=NULL, criterion_qncv=NULL,
            criterion_result=NULL, score=NULL
    ) {
        if (is.null(checkpoint)) return(invisible(FALSE))
        checkpoint_state$inflight <<- list(
            parameters=as.numeric(parameters),
            stage=stage,
            inner=inner,
            criterion_qncv=criterion_qncv,
            criterion_result=criterion_result,
            score=score,
            updated_at=as.character(Sys.time())
        )
        checkpoint_state$updated_at <<- as.character(Sys.time())
        .checkpoint_write(checkpoint_state, checkpoint)
    }
    if (resumed && !is.null(reporter)) reporter$emit(
        1L,
        'checkpoint resumed',
        path=checkpoint,
        stage=checkpoint_state$stage,
        phase=checkpoint_phase,
        evaluations=evaluations,
        best=checkpoint_state$best_criterion
    )
    if (failed_checkpoint && !is.null(reporter)) reporter$emit(
        1L,
        'failed checkpoint discarded',
        path=checkpoint,
        criterion=raw_checkpoint_state$current_criterion
    )
    dimension <- sum(vapply(assemblies, `[[`, integer(1), 'dimension'))
    gradient_probes <- if (identical(control$gradient_method, 'hybrid')) {
        .deterministic_rademacher(dimension, control$gradient_probes)
    } else NULL
    remember <- function(log_sp, value) {
        cached_parameters <<- as.numeric(log_sp)
        cached_evaluation <<- value
        value
    }
    invalid_evaluation <- function(log_sp, reason) {
        score <- numeric(length(log_sp))
        finite <- is.finite(log_sp)
        excess <- pmax(abs(log_sp[finite]) - 25, 0)
        score[finite] <- 2 * sign(log_sp[finite]) * excess
        if (!is.null(best_parameters) && !any(excess > 0) && all(finite)) {
            score <- 2 * (log_sp - best_parameters)
        }
        value <- list(
            criterion=1e50 + sum(excess^2),
            score=score,
            solution=NULL,
            invalid_reason=reason
        )
        cached <- .cdrgam_gaulss_cache_candidate(
            log_sp, value, best_parameters, best_evaluation
        )
        remember(log_sp, cached)
        value
    }
    score_evaluation <- function(log_sp, value) {
        if (!is.null(value$score) || is.null(value$solution)) return(value)
        score_evaluations <<- score_evaluations + 1L
        score_started <- proc.time()[['elapsed']]
        inflight <- if (inflight_matches(log_sp)) {
            checkpoint_state$inflight
        } else NULL
        inner_state <- .cdrgam_gaulss_solution_resume_state(value$solution)
        score_checkpoint <- function(state) write_inflight(
            log_sp, 'exact_score', inner=inner_state,
            criterion_result=value$solution$qncv, score=state
        )
        scored <- tryCatch(
            list(
                value=if (qncv) {
                    .cdrgam_gaulss_sparse_qncv_score(
                        assemblies,
                        value$solution,
                        value$solution$sp,
                        family,
                        gamma=control$gamma,
                        batch_size=control$qncv_batch_size,
                        cores=control$cores,
                        workers=control$gradient_workers,
                        score_batch_size=control$score_batch_size,
                        resume=if (is.list(inflight)) inflight$score else NULL,
                        checkpoint=score_checkpoint
                    )
                } else {
                    .cdrgam_with_blas_threads(
                        control$gradient_blas_threads,
                        .cdrgam_gaulss_sparse_score(
                            assemblies,
                            value$solution,
                            value$solution$sp,
                            family,
                            workers=control$gradient_workers,
                            batch_size=control$score_batch_size
                        )
                    )
                },
                error=NULL
            ),
            error=function(error) list(
                value=NULL, error=conditionMessage(error)
            )
        )
        score <- scored$value
        if (is.null(score) || any(!is.finite(score))) {
            reason <- if (is.null(score)) {
                paste('exact outer score failed:', scored$error)
            } else 'exact outer score was not finite'
            if (!is.null(reporter)) reporter$emit(
                1L,
                'outer exact score failed',
                evaluation=evaluations,
                reason=reason
            )
            stop(reason, call.=FALSE)
        }
        score_plan <<- attr(score, 'score_plan')
        value$score <- as.numeric(score)
        remember(log_sp, value)
        if (!is.null(best_parameters) && identical(
                as.numeric(log_sp), best_parameters)) {
            best_evaluation <<- value
        }
        checkpoint_state$inflight <<- NULL
        write_checkpoint(log_sp, value, phase='exact')
        if (!is.null(reporter)) reporter$emit(
            1L,
            'outer exact score complete',
            evaluation=evaluations,
            seconds=format(
                proc.time()[['elapsed']] - score_started,
                digits=5
            ),
            batches=score_plan$batches,
            batch_size=score_plan$batch_size,
            workers=score_plan$workers,
            criterion=format(value$criterion, digits=10),
            maximum_score=if (length(value$score)) {
                format(max(abs(value$score)), digits=5)
            } else 'none'
        )
        value
    }
    evaluate <- function(log_sp, need_score=FALSE) {
        if (!is.null(cached_parameters) && identical(
                as.numeric(log_sp), cached_parameters)) {
            if (need_score) {
                return(score_evaluation(log_sp, cached_evaluation))
            }
            return(cached_evaluation)
        }
        evaluations <<- evaluations + 1L
        evaluation_started <- proc.time()[['elapsed']]
        inflight <- if (inflight_matches(log_sp)) {
            checkpoint_state$inflight
        } else NULL
        if (any(!is.finite(log_sp)) || any(abs(log_sp) > 25)) {
            return(invalid_evaluation(
                log_sp, 'smoothing parameters left the finite search region'
            ))
        }
        last_inner_checkpoint <- -Inf
        inner_checkpoint <- function(state) {
            now <- proc.time()[['elapsed']]
            if (!is.finite(last_inner_checkpoint) ||
                    now - last_inner_checkpoint >= 60) {
                write_inflight(log_sp, 'inner', inner=state)
                last_inner_checkpoint <<- now
            }
            invisible(NULL)
        }
        solve_fixed <- function(initial) {
            if (hybrid_warmup_active) {
                policy <- .hybrid_rejection_policy(
                    hybrid_rejected_solves,
                    hybrid_rejected_seconds,
                    hybrid_reference_seconds
                )
                if (policy$transition) {
                    hybrid_transition <<- policy
                    if (!is.null(reporter)) reporter$emit(
                        1L,
                        'stochastic warm-up transition requested',
                        reason=policy$reason,
                        rejected_solves=policy$rejected_solves,
                        rejected_seconds=format(
                            policy$rejected_seconds, digits=5
                        ),
                        reference_seconds=format(
                            policy$reference_seconds, digits=5
                        ),
                        exact_score_budget_seconds=format(
                            policy$seconds_limit, digits=5
                        )
                    )
                    stop(.hybrid_transition_condition(policy))
                }
            }
            solve_started <- proc.time()[['elapsed']]
            inner_progress <- if (is.null(reporter)) NULL else {
                function(event, ...) reporter$emit(
                    1L,
                    paste('inner', event),
                    evaluation=evaluations,
                    ...
                )
            }
            result <- tryCatch({
                solution <- .cdrgam_gaulss_sparse_fixed(
                    assemblies, exp(log_sp), b,
                    tolerance=control$inner_tolerance,
                    relative_tolerance=control$inner_relative_tolerance,
                    maxit=control$inner_maxit,
                    initial=initial,
                    cores=control$cores,
                    schur_rank=control$inner_schur_rank,
                    algorithm=control$inner_algorithm,
                    progress=inner_progress,
                    checkpoint=inner_checkpoint
                )
                total_iterations <- max(
                    solution$iterations,
                    if (length(solution$mixed_iterations) == 1L) {
                        solution$mixed_iterations
                    } else 0L
                )
                maximum_damping <- solution$maximum_damping
                continuations <- 0L
                while (!isTRUE(solution$converged) &&
                        identical(solution$termination, 'iteration limit') &&
                        total_iterations < control$inner_adaptive_maxit &&
                        is.finite(solution$directional_derivative) &&
                        solution$directional_derivative < 0 &&
                        is.finite(solution$best_line_search_change) &&
                        solution$best_line_search_change < 0) {
                    remaining <- control$inner_adaptive_maxit -
                        total_iterations
                    continuations <- continuations + 1L
                    if (!is.null(reporter)) reporter$emit(
                        1L,
                        'inner solve continuation',
                        evaluation=evaluations,
                        completed_iterations=total_iterations,
                        remaining_iterations=remaining,
                        gradient=format(solution$gradient_norm, digits=5),
                        best_change=format(
                            solution$best_line_search_change,
                            digits=5
                        )
                    )
                    continued <- .cdrgam_gaulss_sparse_fixed(
                        assemblies, exp(log_sp), b,
                        tolerance=control$inner_tolerance,
                        relative_tolerance=control$inner_relative_tolerance,
                        maxit=min(control$inner_maxit, remaining),
                        initial=solution,
                        cores=control$cores,
                        schur_rank=control$inner_schur_rank,
                        algorithm=control$inner_algorithm,
                        progress=inner_progress,
                        checkpoint=inner_checkpoint
                    )
                    total_iterations <- total_iterations +
                        continued$iterations
                    maximum_damping <- max(
                        maximum_damping,
                        continued$maximum_damping
                    )
                    solution <- continued
                }
                solution$iterations <- total_iterations
                solution$continuations <- continuations
                solution$maximum_damping <- maximum_damping
                list(solution=solution, error=NULL)
            }, error=function(error) {
                list(solution=NULL, error=conditionMessage(error))
            })
            if (hybrid_warmup_active) {
                hybrid_rejected_solves <<- hybrid_rejected_solves + 1L
                hybrid_rejected_seconds <<- hybrid_rejected_seconds +
                    proc.time()[['elapsed']] - solve_started
            }
            result
        }
        fixed_failure <- function(result) {
            if (is.null(result$solution)) {
                return(paste('failed:', result$error))
            }
            solution <- result$solution
            paste0(
                solution$termination,
                ' after ', solution$iterations, ' iterations',
                '; maximum penalized score ',
                format(solution$gradient_norm, digits=5),
                '; proposed-step norm ',
                format(solution$step_norm, digits=5),
                '; step mode ', solution$step_mode,
                '; maximum damping ',
                format(solution$maximum_damping, digits=5)
            )
        }
        resume_solution <- if (is.list(inflight) &&
                is.list(inflight$inner)) inflight$inner else NULL
        initial_solution <- if (!is.null(resume_solution)) {
            resume_solution
        } else warm_solution
        used_warm_start <- !is.null(initial_solution)
        if (used_warm_start) warm_starts <<- warm_starts + 1L
        fixed <- solve_fixed(initial_solution)
        warm_failure <- if (is.null(fixed$solution) ||
                !isTRUE(fixed$solution$converged)) {
            fixed_failure(fixed)
        } else NULL
        if (used_warm_start && (is.null(fixed$solution) ||
                !isTRUE(fixed$solution$converged))) {
            cold_fallbacks <<- cold_fallbacks + 1L
            fixed <- solve_fixed(NULL)
        }
        solution <- fixed$solution
        if (is.null(solution)) {
            reason <- paste('inner solve failed:', fixed$error)
            if (!is.null(warm_failure)) {
                reason <- paste(reason, 'warm-start result:', warm_failure)
            }
            if (!is.null(reporter)) reporter$emit(
                1L,
                'outer inner solve rejected',
                evaluation=evaluations,
                warm_start=used_warm_start,
                reason=reason
            )
            return(invalid_evaluation(log_sp, reason))
        }
        if (!isTRUE(solution$converged)) {
            reason <- paste0(
                'inner solve did not converge after ', solution$iterations,
                ' iterations; maximum penalized score ',
                format(solution$gradient_norm, digits=5),
                '; termination ', solution$termination,
                '; proposed-step norm ',
                format(solution$step_norm, digits=5),
                '; directional derivative ',
                format(solution$directional_derivative, digits=5),
                '; best line-search change ',
                format(solution$best_line_search_change, digits=5),
                '; objective ', format(solution$objective, digits=10),
                '; step mode ', solution$step_mode,
                '; maximum damping ',
                format(solution$maximum_damping, digits=5)
            )
            if (!is.null(warm_failure)) {
                reason <- paste(reason, 'warm-start result:', warm_failure)
            }
            if (!is.null(reporter)) reporter$emit(
                1L,
                'outer inner solve rejected',
                evaluation=evaluations,
                warm_start=used_warm_start,
                reason=reason
            )
            return(invalid_evaluation(log_sp, reason))
        }
        warm_solution <<- solution
        inner_state <- .cdrgam_gaulss_solution_resume_state(solution)
        write_inflight(log_sp, 'criterion_ready', inner=inner_state)
        inner_seconds <- proc.time()[['elapsed']] - evaluation_started
        if (!is.null(reporter)) reporter$emit(
            1L,
            'outer inner solve complete',
            evaluation=evaluations,
            seconds=format(inner_seconds, digits=5),
            iterations=solution$iterations,
            moment_workers=solution$moment_workers,
            linear_iterations=solution$linear_iterations,
            linear_residual=format(solution$linear_residual, digits=5),
            warm_start=used_warm_start,
            cold_fallback=used_warm_start && !isTRUE(solution$warm_started)
        )
        criterion_evaluation <- if (qncv) {
            if (is.list(inflight) &&
                    is.list(inflight$criterion_result)) {
                inflight$criterion_result
            } else {
                criterion_checkpoint <- function(state) write_inflight(
                    log_sp, 'criterion', inner=inner_state,
                    criterion_qncv=state
                )
                result <- .cdrgam_gaulss_sparse_qncv(
                    assemblies,
                    solution,
                    gamma=control$gamma,
                    batch_size=control$qncv_batch_size,
                    workers=control$gradient_workers,
                    blas_threads=control$gradient_blas_threads,
                    resume=if (is.list(inflight)) {
                        inflight$criterion_qncv
                    } else NULL,
                    checkpoint=criterion_checkpoint
                )
                write_inflight(
                    log_sp, 'score_ready', inner=inner_state,
                    criterion_result=result,
                    score=if (is.list(inflight)) inflight$score else NULL
                )
                result
            }
        } else NULL
        if (qncv) {
            criterion <- criterion_evaluation$criterion
        } else {
            penalty_determinant <- 0
            offset <- 0L
            for (parameter in names(assemblies)) {
                count <- counts[[parameter]]
                indices <- offset + seq_len(count)
                penalty_determinant <- penalty_determinant +
                    .sparse_penalty_logdet(
                        assemblies[[parameter]]$blocks,
                        exp(log_sp[indices])
                    )$value
                offset <- offset + count
            }
            criterion <- 2 * solution$objective +
                .cdr_factor_logdet(solution$factor) - penalty_determinant
        }
        if (!is.finite(criterion)) {
            return(invalid_evaluation(
                log_sp,
                paste(control$criterion, 'criterion was not finite')
            ))
        }
        solution$criterion <- criterion
        solution$sp <- exp(log_sp)
        if (qncv) solution$qncv <- criterion_evaluation
        value <- list(
            criterion=criterion,
            score=NULL,
            solution=solution,
            inner_seconds=inner_seconds
        )
        remember(log_sp, value)
        if (is.null(best_evaluation) || is.null(best_evaluation$solution) ||
                criterion < best_evaluation$criterion) {
            best_parameters <<- as.numeric(log_sp)
            best_evaluation <<- value
        }
        if (need_score) score_evaluation(log_sp, value) else value
    }
    if (!length(initial_log_sp)) {
        evaluated <- evaluate(numeric())
        solution <- evaluated$solution
        optimization <- list(
            par=numeric(), value=solution$criterion,
            convergence=0L,
            counts=c(`function`=1L, gradient=NA_integer_)
        )
        write_checkpoint(
            numeric(), evaluated, phase='exact', stage='complete',
            optimization=optimization
        )
        return(list(
            optimization=optimization,
            solution=solution,
            evaluations=evaluations,
            score_evaluations=score_evaluations,
            stochastic_score_evaluations=stochastic_score_evaluations,
            warm_starts=warm_starts,
            cold_fallbacks=cold_fallbacks,
            score_plan=score_plan,
            stochastic_score_plan=stochastic_score_plan
        ))
    }
    objective <- function(parameters) evaluate(parameters)$criterion
    objective_validity <- function(parameters, criterion) {
        value <- evaluate(parameters)
        is.finite(criterion) && !is.null(value$solution)
    }
    score <- function(parameters) evaluate(parameters, need_score=TRUE)$score
    stochastic_score <- function(parameters) {
        value <- evaluate(parameters)
        if (is.null(value$solution)) return(value$score)
        stochastic_score_evaluations <<- stochastic_score_evaluations + 1L
        started <- proc.time()[['elapsed']]
        stochastic <- tryCatch(
            list(
                value=.cdrgam_gaulss_sparse_stochastic_score(
                    assemblies,
                    value$solution,
                    value$solution$sp,
                    family,
                    gradient_probes
                ),
                error=NULL
            ),
            error=function(error) list(
                value=NULL, error=conditionMessage(error)
            )
        )
        if (is.null(stochastic$value) || any(!is.finite(stochastic$value))) {
            reason <- if (is.null(stochastic$value)) {
                paste('stochastic outer score failed:', stochastic$error)
            } else 'stochastic outer score was not finite'
            return(invalid_evaluation(parameters, reason)$score)
        }
        stochastic_score_plan <<- attr(stochastic$value, 'score_plan')
        hybrid_reference_seconds <<- c(
            utils::tail(hybrid_reference_seconds, 8L),
            value$inner_seconds
        )
        hybrid_rejected_solves <<- 0L
        hybrid_rejected_seconds <<- 0
        write_checkpoint(parameters, value, phase='stochastic')
        if (!is.null(reporter)) reporter$emit(
            1L,
            'outer stochastic score complete',
            evaluation=evaluations,
            probes=stochastic_score_plan$probes,
            seconds=format(proc.time()[['elapsed']] - started, digits=5),
            criterion=format(value$criterion, digits=10),
            maximum_score=format(max(abs(stochastic$value)), digits=5)
        )
        as.numeric(stochastic$value)
    }
    parameters <- if (!is.null(checkpoint_state$best_log_sp)) {
        as.numeric(checkpoint_state$best_log_sp)
    } else if (!is.null(checkpoint_state$current_log_sp)) {
        as.numeric(checkpoint_state$current_log_sp)
    } else initial_log_sp
    total_counts <- checkpoint_state$optimizer_counts
    if (is.null(total_counts) || length(total_counts) != 2L) {
        total_counts <- c(`function`=0L, gradient=0L)
    }
    stochastic_warmup <- checkpoint_state$stochastic_warmup
    completed_optimization <- if (
        identical(checkpoint_state$stage, 'complete')
    ) checkpoint_state$optimization else NULL
    run_stochastic <- !qncv &&
        identical(control$gradient_method, 'hybrid') &&
        !identical(checkpoint_phase, 'exact') &&
        !identical(checkpoint_state$stage, 'complete')
    if (run_stochastic) {
        warmup_limit <- min(25L, control$optimizer_maxit)
        prior_warmup_scores <- stochastic_score_evaluations
        remaining_warmup <- max(0L, warmup_limit - prior_warmup_scores)
        if (!is.null(reporter)) reporter$phase(
            'stochastic smoothing-parameter warm-up',
            probes=control$gradient_probes,
            maximum_iterations=remaining_warmup
        )
        warmup_evaluations <- evaluations
        warmup_scores <- stochastic_score_evaluations
        hybrid_warmup_active <- remaining_warmup > 0L
        stochastic_warmup <- if (!remaining_warmup) {
            list(
                par=parameters,
                value=evaluate(parameters)$criterion,
                counts=c(`function`=0L, gradient=0L),
                convergence=1L,
                message='stochastic warm-up iteration budget was exhausted'
            )
        } else tryCatch(
            stats::optim(
                par=parameters,
                fn=objective,
                gr=stochastic_score,
                method='BFGS',
                control=list(
                    maxit=remaining_warmup,
                    reltol=1e-7
                )
            ),
            cdrgam_hybrid_transition=function(condition) list(
                par=best_parameters,
                value=best_evaluation$criterion,
                counts=c(
                    `function`=evaluations - warmup_evaluations,
                    gradient=stochastic_score_evaluations - warmup_scores
                ),
                convergence=1L,
                message=condition$message,
                transition=hybrid_transition
            ),
            finally={
                hybrid_warmup_active <- FALSE
            }
        )
        total_counts <- total_counts + stochastic_warmup$counts
        if (!is.null(best_parameters)) {
            parameters <- best_parameters
            warm_solution <- best_evaluation$solution
            cached_parameters <- best_parameters
            cached_evaluation <- best_evaluation
            best_evaluation <- score_evaluation(
                best_parameters, best_evaluation
            )
        }
        checkpoint_phase <- 'exact'
        checkpoint_state$stochastic_complete <- TRUE
        checkpoint_state$stochastic_warmup <- stochastic_warmup
        checkpoint_state$optimizer_counts <- total_counts
        write_checkpoint(parameters, best_evaluation, phase='exact')
        if (!is.null(reporter)) reporter$phase(
            'exact smoothing-parameter refinement',
            criterion=format(evaluate(parameters)$criterion, digits=10),
            transition=if (is.null(stochastic_warmup$transition)) {
                'stochastic warm-up completed'
            } else stochastic_warmup$transition$reason
        )
    }
    lower_bound <- rep.int(-25, length(parameters))
    upper_bound <- rep.int(25, length(parameters))
    exact_base_counts <- checkpoint_state$exact_base_counts
    if (is.null(exact_base_counts) || length(exact_base_counts) != 2L) {
        exact_base_counts <- total_counts
    }
    checkpoint_state$exact_base_counts <- exact_base_counts
    directional_diagnostic <- NULL
    optimizer_progress_callback <- function(record) {
        if (isTRUE(record$accepted)) {
            matching_cache <- !is.null(cached_parameters) && identical(
                as.numeric(record$parameters), cached_parameters
            ) && !is.null(cached_evaluation$solution)
            if (matching_cache) {
                accepted_parameters <<- cached_parameters
                accepted_evaluation <<- cached_evaluation
            } else if (!is.null(best_parameters) && identical(
                    as.numeric(record$parameters), best_parameters
                ) && !is.null(best_evaluation$solution)) {
                accepted_parameters <<- best_parameters
                accepted_evaluation <<- best_evaluation
            }
        } else if (identical(record$event, 'rejected') &&
                !is.null(accepted_evaluation$solution)) {
            warm_solution <<- accepted_evaluation$solution
        }
        optimizer_state <- record$optimizer_state
        record$optimizer_state <- NULL
        checkpoint_state$optimizer_state <<- optimizer_state
        checkpoint_state$optimizer_progress <<- record
        state_counts <- c(
            `function`=optimizer_state$function_evaluations,
            gradient=optimizer_state$gradient_evaluations
        )
        checkpoint_state$optimizer_counts <<- exact_base_counts + state_counts
        optimizer_history <- checkpoint_state$optimizer_history
        if (is.null(optimizer_history)) optimizer_history <- list()
        optimizer_history[[length(optimizer_history) + 1L]] <- record
        if (length(optimizer_history) > 1000L) {
            optimizer_history <- utils::tail(optimizer_history, 1000L)
        }
        checkpoint_state$optimizer_history <<- optimizer_history
        write_checkpoint(
            record$parameters,
            list(criterion=record$criterion),
            phase='exact'
        )
        if (!is.null(reporter)) reporter$emit(
            1L,
            paste('outer', record$event),
            iteration=record$iteration,
            criterion=format(record$criterion, digits=10),
            projected_gradient=format(
                record$projected_gradient_max,
                digits=5
            ),
            gradient_ratio=format(record$gradient_ratio, digits=5),
            step_max=format(record$step_max, digits=5),
            trust_radius=format(record$trust_radius, digits=5),
            actual_improvement=format(
                record$actual_improvement,
                digits=5
            ),
            acceptance_ratio=format(record$acceptance_ratio, digits=5),
            accepted=record$accepted,
            rejected=record$rejected_steps,
            consecutive_rejections=record$consecutive_rejections,
            step_type=record$step_type,
            curvature_resets=record$curvature_resets,
            curvature_reset=record$curvature_reset,
            recovery_resets=record$recovery_resets
        )
        invisible(record)
    }
    curvature_recovery_assessment <- function(
            log_sp,
            criterion,
            gradient,
            objective_noise=0
    ) {
        scale <- max(1, max(abs(gradient)))
        list(
            converged=FALSE,
            message='diagonal outer-curvature recovery requested',
            restart_hessian=diag(scale, length(gradient)),
            restart_radius=min(control$optimizer_trust_radius, 0.1),
            diagnostics=list(
                criterion=criterion,
                maximum_score=max(abs(gradient)),
                objective_noise=objective_noise
            )
        )
    }
    stagnation_diagnostic <- function(
            log_sp,
            criterion,
            gradient,
            objective_noise=0
    ) {
        diagnostic <- .outer_directional_derivative_check(
            log_sp,
            criterion,
            gradient,
            objective,
            lower_bound,
            upper_bound,
            objective_noise=objective_noise,
            objective_validity=objective_validity
        )
        directional_diagnostic <<- diagnostic
        if (!is.null(reporter)) reporter$emit(
            1L,
            'outer directional derivative check',
            criterion=format(criterion, digits=10),
            analytic=format(diagnostic$analytic, digits=5),
            finite_difference=format(
                diagnostic$finite_difference,
                digits=5
            ),
            discrepancy=format(diagnostic$discrepancy, digits=5),
            tolerance=format(diagnostic$tolerance, digits=5),
            consistent=diagnostic$consistent,
            step=diagnostic$step
        )
        diagnostic$message <- paste0(
            diagnostic$message,
            ' (analytic ', format(diagnostic$analytic, digits=5),
            ', finite difference ',
            format(diagnostic$finite_difference, digits=5), ')'
        )
        diagnostic
    }
    optimization <- completed_optimization
    if (!is.null(completed_optimization)) {
        parameters <- completed_optimization$par
        if (!is.null(reporter)) reporter$phase(
            'completed smoothing-parameter optimization recovered',
            criterion=format(completed_optimization$value, digits=10)
        )
    } else {
        saved_optimizer_state <- checkpoint_state$optimizer_state
        use_optimizer_state <- is.list(saved_optimizer_state) &&
            length(saved_optimizer_state$parameters) == length(parameters) &&
            all(is.finite(saved_optimizer_state$parameters)) &&
            length(checkpoint_state$best_criterion) == 1L &&
            is.finite(checkpoint_state$best_criterion)
        if (use_optimizer_state) {
            parameters <- as.numeric(saved_optimizer_state$parameters)
        }
        optimization <- .safeguarded_outer_bfgs(
            par=parameters,
            fn=objective,
            gr=score,
            lower=lower_bound,
            upper=upper_bound,
            maxit=control$optimizer_maxit,
            gradient_tolerance=control$optimizer_gradient_tolerance,
            initial_radius=control$optimizer_trust_radius,
            progress=optimizer_progress_callback,
            state=if (use_optimizer_state) saved_optimizer_state else NULL,
            convergence_assessment=curvature_recovery_assessment,
            stagnation_diagnostic=stagnation_diagnostic,
            objective_validity=objective_validity
        )
        total_counts <- exact_base_counts + optimization$counts
        optimization$counts <- total_counts
        optimization$directional_diagnostic <- directional_diagnostic
        checkpoint_state$optimizer_counts <- total_counts
    }
    optimization$stochastic_warmup <- stochastic_warmup
    retained <- evaluate(optimization$par, need_score=TRUE)
    recovered <- is.null(retained$solution) ||
        (!is.null(best_evaluation) &&
            best_evaluation$criterion < retained$criterion)
    if (recovered && !is.null(best_evaluation)) {
        best_evaluation <- score_evaluation(best_parameters, best_evaluation)
        retained <- best_evaluation
        optimization$par <- best_parameters
        optimization$value <- best_evaluation$criterion
        optimization$gradient <- best_evaluation$score
        optimization$convergence <- 1L
        optimization$message <- paste0(
            optimization$message,
            '; optimizer endpoint was replaced by its best valid evaluation'
        )
    }
    projected_gradient <- optimization$gradient
    at_lower <- optimization$par <= lower_bound + 1e-10
    at_upper <- optimization$par >= upper_bound - 1e-10
    projected_gradient[at_lower & projected_gradient > 0] <- 0
    projected_gradient[at_upper & projected_gradient < 0] <- 0
    optimization$projected_gradient <- projected_gradient
    if (max(abs(projected_gradient)) >
            control$optimizer_gradient_tolerance) {
        optimization$convergence <- 1L
        score_message <- paste0(
            'outer score ',
            format(max(abs(projected_gradient)), digits=5),
            ' exceeds tolerance ',
            format(control$optimizer_gradient_tolerance, digits=5)
        )
        if (is.null(optimization$message) || !nzchar(optimization$message)) {
            optimization$message <- score_message
        } else if (!grepl(score_message, optimization$message, fixed=TRUE)) {
            optimization$message <- paste(
                optimization$message,
                score_message,
                sep='; '
            )
        }
    } else {
        optimization$convergence <- 0L
    }
    solution <- retained$solution
    if (is.null(solution)) {
        stop(
            'Sparse gaulss ', control$criterion,
            ' optimization ended at an invalid fit: ',
            retained$invalid_reason
        )
    }
    checkpoint_complete <- identical(optimization$convergence, 0L)
    if (!checkpoint_complete) {
        checkpoint_state$optimizer_state <- NULL
        checkpoint_state$optimization <- NULL
        checkpoint_state$optimizer_history <- list()
        checkpoint_state$optimizer_progress <- NULL
        checkpoint_state$exact_restart <-
            checkpoint_count('exact_restart') + 1L
    }
    write_checkpoint(
        optimization$par,
        retained,
        phase='exact',
        stage=if (checkpoint_complete) 'complete' else 'optimization',
        optimization=if (checkpoint_complete) optimization else NULL
    )
    list(
        optimization=optimization,
        solution=solution,
        evaluations=evaluations,
        score_evaluations=score_evaluations,
        stochastic_score_evaluations=stochastic_score_evaluations,
        warm_starts=warm_starts,
        cold_fallbacks=cold_fallbacks,
        score_plan=score_plan,
        stochastic_score_plan=stochastic_score_plan
    )
}

.cdrgam_distributional_sparse_terms <- function(design, assemblies) {
    metadata <- list()
    labels <- character()
    parameter_terms <- vector('list', length(assemblies))
    names(parameter_terms) <- names(assemblies)
    offset <- 0L
    for (parameter in names(assemblies)) {
        terms <- design$parameters[[parameter]]$terms
        parameter_terms[[parameter]] <- integer(length(terms))
        local <- assemblies[[parameter]]$term_metadata
        for (term_index in seq_along(terms)) {
            info <- local[[term_index]]
            info$coefficient_index <- info$coefficient_index + offset
            metadata[[length(metadata) + 1L]] <- info
            labels[[length(labels) + 1L]] <- paste0(
                parameter, ':', terms[[term_index]]$name
            )
            parameter_terms[[parameter]][[term_index]] <- length(metadata)
        }
        offset <- offset + assemblies[[parameter]]$dimension
    }
    list(metadata=metadata, labels=labels, parameter_terms=parameter_terms)
}

.cdrgam_distributional_sparse_smooths <- function(assemblies) {
    output <- list()
    offset <- 0L
    for (parameter in names(assemblies)) {
        for (smooth in assemblies[[parameter]]$smooth) {
            smooth$first.para <- smooth$first.para + offset
            smooth$last.para <- smooth$last.para + offset
            smooth$label <- paste0(parameter, ':', smooth$label)
            output[[length(output) + 1L]] <- smooth
        }
        offset <- offset + assemblies[[parameter]]$dimension
    }
    output
}

.fit_distributional_sparse <- function(
        design, family, method=NULL, checkpoint=NULL, trace=FALSE,
        sparse_control=list(), gamma=1, ...
) {
    dots <- list(...)
    if (!is.null(dots$weights) && any(dots$weights != 1)) {
        stop(
            'gaulss does not support non-unit prior weights; mgcv accepts ',
            'but ignores them'
        )
    }
    if (!is.null(method) && !(method %in% c('REML', 'fREML', 'QNCV'))) {
        stop('The distributional sparse backend supports REML and QNCV')
    }
    criterion_method <- if (identical(method, 'QNCV')) 'QNCV' else 'REML'
    qncv <- identical(criterion_method, 'QNCV')
    if (length(gamma) != 1L || !is.numeric(gamma) ||
            !is.finite(gamma) || gamma <= 0) {
        stop('gamma must be one positive finite number')
    }
    if (qncv && !is.null(dots$nei)) {
        stop(
            'Sparse QNCV currently supports its default ',
            'leave-one-response-out neighborhoods only'
        )
    }
    allowed <- c(
        'crossprod_chunk_size', 'supernodal', 'optimizer_maxit',
        'optimizer_gradient_tolerance', 'optimizer_trust_radius',
        'inner_tolerance', 'inner_relative_tolerance', 'inner_maxit',
        'inner_adaptive_maxit',
        'inner_algorithm', 'inner_schur_rank',
        'cores', 'gradient',
        'gradient_probes', 'gradient_workers', 'score_batch_size',
        'qncv_batch_size'
    )
    unknown <- setdiff(names(sparse_control), allowed)
    if (length(unknown)) {
        stop(
            'Unsupported distributional sparse_control entries: ',
            paste(unknown, collapse=', ')
        )
    }
    control <- list(
        optimizer_maxit=if (is.null(sparse_control$optimizer_maxit)) {
            500L
        } else as.integer(sparse_control$optimizer_maxit),
        optimizer_gradient_tolerance=if (is.null(
            sparse_control$optimizer_gradient_tolerance
        )) 2e-3 else sparse_control$optimizer_gradient_tolerance,
        optimizer_trust_radius=if (is.null(
            sparse_control$optimizer_trust_radius
        )) 2 else sparse_control$optimizer_trust_radius,
        inner_tolerance=if (is.null(sparse_control$inner_tolerance)) {
            1e-5
        } else sparse_control$inner_tolerance,
        inner_relative_tolerance=if (is.null(
            sparse_control$inner_relative_tolerance
        )) 1e-8 else sparse_control$inner_relative_tolerance,
        inner_maxit=if (is.null(sparse_control$inner_maxit)) {
            200L
        } else as.integer(sparse_control$inner_maxit),
        inner_adaptive_maxit=if (is.null(
            sparse_control$inner_adaptive_maxit
        )) {
            if (is.null(sparse_control$inner_maxit)) 600L else {
                as.integer(sparse_control$inner_maxit)
            }
        } else as.integer(sparse_control$inner_adaptive_maxit)
    )
    numeric_control <- unlist(control, use.names=FALSE)
    if (any(!is.finite(numeric_control)) || any(numeric_control <= 0)) {
        stop('Distributional sparse optimizer controls must be positive')
    }
    if (control$inner_adaptive_maxit < control$inner_maxit) {
        stop('inner_adaptive_maxit must be at least inner_maxit')
    }
    requested_inner_algorithm <- if (is.null(
        sparse_control$inner_algorithm
    )) {
        'auto'
    } else match.arg(
        sparse_control$inner_algorithm,
        c('auto', 'mixed', 'schur', 'alternating', 'joint')
    )
    control$inner_schur_rank <- if (is.null(
        sparse_control$inner_schur_rank
    )) 6L else .cdrgam_positive_integer(
        sparse_control$inner_schur_rank,
        'sparse_control$inner_schur_rank'
    )
    cores <- .cdrgam_available_cores(sparse_control$cores)
    gradient_requested <- if (is.null(sparse_control$gradient)) {
        'auto'
    } else {
        match.arg(sparse_control$gradient, c('auto', 'exact', 'hybrid'))
    }
    if (qncv && identical(gradient_requested, 'hybrid')) {
        stop('Sparse QNCV supports gradient="auto" or "exact"')
    }
    gradient_probes <- if (is.null(sparse_control$gradient_probes)) {
        64L
    } else {
        .cdrgam_positive_integer(
            sparse_control$gradient_probes,
            'sparse_control$gradient_probes'
        )
    }
    previous_blas_threads <- .cdrgam_blas_threads()
    on.exit(
        RhpcBLASctl::blas_set_num_threads(previous_blas_threads),
        add=TRUE
    )
    RhpcBLASctl::blas_set_num_threads(cores)
    reporter <- .new_solver_reporter(trace, 'distributional sparse')
    reporter$phase('parameter model and sparse-pattern setup')
    assemblies <- lapply(design$parameters, function(parameter_design) {
        do.call(
            .cdrgam_distributional_sparse_setup,
            c(
                list(design=parameter_design, sparse_control=sparse_control),
                dots
            )
        )
    })
    design_cache <- .cdrgam_distributional_cache_designs(assemblies)
    assemblies <- design_cache$assemblies
    reporter$emit(
        1L,
        'distributional design cache',
        enabled=design_cache$enabled,
        megabytes=format(design_cache$bytes / 1024^2, digits=5),
        projected_megabytes=format(
            design_cache$projected_bytes / 1024^2, digits=5
        ),
        memory_source=design_cache$memory_source,
        available_megabytes=format(
            design_cache$memory_available_bytes / 1024^2, digits=5
        )
    )
    y <- assemblies$location$setup$y
    if (!identical(y, assemblies$scale$setup$y)) {
        stop('Distributional parameter setups retained different responses')
    }
    score_tasks <- sum(vapply(
        assemblies,
        function(assembly) length(assembly$penalty_components),
        integer(1)
    ))
    dimension <- sum(vapply(assemblies, `[[`, integer(1), 'dimension'))
    control$inner_algorithm <- .cdrgam_gaulss_inner_algorithm(
        requested_inner_algorithm, dimension
    )
    gradient_method <- if (qncv) {
        'exact'
    } else if (identical(gradient_requested, 'auto')) {
        if (dimension > max(512L, 4L * gradient_probes)) 'hybrid' else 'exact'
    } else gradient_requested
    gradient_plan <- .cdrgam_parallel_plan(
        cores,
        max(1L, score_tasks),
        sparse_control$gradient_workers
    )
    if (.Platform$OS.type == 'windows' &&
            !is.null(sparse_control$gradient_workers) &&
            sparse_control$gradient_workers > 1L) {
        warning(
            'Parallel exact distributional scores require fork support; ',
            'using one gradient worker on Windows',
            call.=FALSE
        )
    }
    control$cores <- cores
    control$inner_blas_threads <- cores
    control$gradient_method <- gradient_method
    control$gradient_probes <- gradient_probes
    control$gradient_workers <- gradient_plan$workers
    control$gradient_blas_threads <- gradient_plan$blas_threads
    control$gradient_worker_source <- gradient_plan$worker_source
    control$score_batch_size <- if (is.null(sparse_control$score_batch_size)) {
        NULL
    } else {
        .cdrgam_positive_integer(
            sparse_control$score_batch_size,
            'sparse_control$score_batch_size'
        )
    }
    control$qncv_batch_size <- if (qncv) {
        memory <- .cdrgam_memory_availability()
        target_bytes <- if (is.finite(memory$available_bytes)) {
            min(256 * 1024^2, 0.05 * memory$available_bytes)
        } else 128 * 1024^2
        automatic <- max(1L, min(
            length(y),
            assemblies$location$crossprod_chunk_size,
            assemblies$scale$crossprod_chunk_size,
            floor(target_bytes / max(1, 32 * as.double(dimension)))
        ))
        if (is.null(sparse_control$qncv_batch_size)) {
            as.integer(automatic)
        } else min(
            length(y),
            .cdrgam_positive_integer(
                sparse_control$qncv_batch_size,
                'sparse_control$qncv_batch_size'
            )
        )
    } else NULL
    control$criterion <- criterion_method
    control$gamma <- gamma
    reporter$phase(
        'joint smoothing-parameter optimization',
        criterion=criterion_method,
        smoothing_parameters=score_tasks,
        requested_gradient=gradient_requested,
        gradient=control$gradient_method,
        gradient_probes=if (identical(
            control$gradient_method, 'hybrid'
        )) control$gradient_probes else 0L,
        cores=control$cores,
        inner_blas_threads=control$inner_blas_threads,
        gradient_workers=control$gradient_workers,
        gradient_blas_threads=control$gradient_blas_threads,
        exact_score_batching=if (qncv) {
            'memory-bounded Hessian directions and leave-out solves'
        } else 'memory-bounded shared traces',
        qncv_response_batch_size=if (qncv) {
            control$qncv_batch_size
        } else 0L,
        inner_schur_rank=if (identical(
            control$inner_algorithm, 'schur'
        )) control$inner_schur_rank else 0L,
        inner_solver=paste(
            'streamed sparse', control$inner_algorithm, 'Fisher scoring'
        )
    )
    checkpoint_signature <- .cdrgam_distributional_checkpoint_signature(
        design, assemblies, family, control
    )
    result <- .cdrgam_with_blas_threads(
        control$inner_blas_threads,
        .cdrgam_gaulss_sparse_laml(
            assemblies, family, control, reporter=reporter,
            checkpoint=checkpoint,
            checkpoint_signature=checkpoint_signature
        )
    )
    if (qncv && !is.null(result$score_plan)) {
        control$gradient_workers_requested <-
            result$score_plan$requested_workers
        control$gradient_workers <- result$score_plan$workers
        control$gradient_blas_threads <- result$score_plan$blas_threads
        control$gradient_memory_limited <-
            result$score_plan$memory_limited
        control$gradient_memory_source <- result$score_plan$memory_source
        control$gradient_memory_available_bytes <-
            result$score_plan$memory_available_bytes
        control$gradient_memory_budget_bytes <-
            result$score_plan$memory_budget_bytes
        control$gradient_design_cache <- result$score_plan$design_cache
    }
    solution <- result$solution
    dimensions <- vapply(assemblies, `[[`, integer(1), 'dimension')
    ranges <- split(seq_len(sum(dimensions)), rep(names(dimensions), dimensions))
    coefficient_names <- unlist(Map(function(parameter, assembly) {
        paste0(parameter, ':', assembly$coefficient_names)
    }, names(assemblies), assemblies), use.names=FALSE)
    coefficients <- stats::setNames(solution$coefficients, coefficient_names)
    likelihood_trace <- .sparse_logdet_score(
        solution$factor, solution$likelihood_hessian, chunk_size=256L
    )
    effective_df <- likelihood_trace
    term_info <- .cdrgam_distributional_sparse_terms(design, assemblies)
    smooths <- .cdrgam_distributional_sparse_smooths(assemblies)
    counts <- vapply(
        assemblies, function(assembly) length(assembly$penalty_components),
        integer(1)
    )
    sp_names <- unlist(Map(function(parameter, assembly) {
        paste0(parameter, ':', assembly$sp_names)
    }, names(assemblies), assemblies), use.names=FALSE)
    smoothing_parameters <- stats::setNames(solution$sp, sp_names)
    penalty_trace <- vapply(seq_along(solution$sp), function(i) {
        parameter <- if (i <= counts[['location']]) 'location' else 'scale'
        local_index <- if (identical(parameter, 'location')) {
            i
        } else i - counts[['location']]
        derivative <- Matrix::Matrix(
            0, sum(dimensions), sum(dimensions), sparse=TRUE
        )
        range <- ranges[[parameter]]
        derivative[range, range] <- solution$sp[[i]] *
            assemblies[[parameter]]$penalty_components[[local_index]]
        .sparse_logdet_score(solution$factor, derivative, chunk_size=256L)
    }, numeric(1))
    parameter_preparation <- lapply(
        design$parameters, .cdrgam_distributional_preparation
    )
    identifiability <- lapply(design$parameters, `[[`, 'identifiability')
    output <- list(
        coefficients=coefficients,
        fitted.values=solution$fitted_values,
        residuals=solution$residuals,
        linear.predictors=solution$linear_predictors,
        family=family,
        sp=smoothing_parameters,
        scale=1,
        sig2=1,
        method=criterion_method,
        smooth=smooths,
        converged=isTRUE(solution$converged) &&
            identical(result$optimization$convergence, 0L),
        edf=effective_df,
        df.residual=length(y) - effective_df,
        y=y,
        prior.weights=solution$prior_weights,
        deviance=sum(
            solution$prior_weights *
                (solution$residuals / solution$sigma)^2
        ),
        loglik=-sum(solution$prior_weights * (
            log(solution$sigma) +
                0.5 * (solution$residuals / solution$sigma)^2 +
                0.5 * log(2 * pi)
        )),
        reml=if (qncv) NULL else solution$criterion,
        gcv.ubre=if (qncv) stats::setNames(
            solution$criterion, 'QNCV'
        ) else NULL,
        qncv=if (qncv) solution$criterion else NULL,
        optimizer=result$optimization,
        distributional=list(
            sigma=solution$sigma,
            coefficient_ranges=ranges,
            inner_iterations=solution$iterations,
            inner_continuations=solution$continuations,
            inner_algorithm=solution$algorithm,
            inner_moment_workers=solution$moment_workers,
            gradient_norm=solution$gradient_norm,
            likelihood_hessian=solution$likelihood_hessian,
            evaluations=result$evaluations,
            score_evaluations=result$score_evaluations,
            stochastic_score_evaluations=
                result$stochastic_score_evaluations,
            warm_starts=result$warm_starts,
            cold_fallbacks=result$cold_fallbacks,
            score_plan=result$score_plan,
            stochastic_score_plan=result$stochastic_score_plan,
            qncv=if (qncv) list(
                minimum_determinant=solution$qncv$minimum_determinant,
                response_batch_size=solution$qncv$batch_size
            ) else NULL
        ),
        sparse=list(
            factor=solution$factor,
            dimension=sum(dimensions),
            effective_df=effective_df,
            penalty_trace=penalty_trace,
            control=c(
                control,
                list(
                    gradient=control$gradient_method,
                    crossprod_chunk_size=vapply(
                        assemblies, `[[`, integer(1),
                        'crossprod_chunk_size'
                    )
                )
            ),
            factor_class=class(solution$factor)[[1L]],
            factor_nonzeros=.cdr_factor_nonzeros(solution$factor),
            condition_indicator=.cdr_factor_condition_indicator(
                solution$factor,
                system=solution$system
            ),
            design_cache=design_cache[setdiff(
                names(design_cache), 'assemblies'
            )],
            crossprod_chunks=solution$chunks,
            observation_count=length(y)
        ),
        cdrgam=list(
            schema_version=1L,
            engine='sparse',
            backend='sparse',
            distributional=TRUE,
            parameter_names=names(assemblies),
            parameter_terms=term_info$parameter_terms,
            parametric_indices=unlist(Map(function(range, assembly) {
                range[seq_len(assembly$setup$nsdf)]
            }, ranges, assemblies), use.names=FALSE),
            formula=list(
                user=design$formula,
                normalized=design$normalized_formula,
                effective=design$effective_formula,
                mgcv=lapply(assemblies, function(assembly) {
                    assembly$setup$formula
                })
            ),
            preparation=list(parameters=parameter_preparation),
            scaling=NULL,
            identifiability=identifiability,
            term_labels=term_info$labels,
            terms=term_info$metadata,
            prediction=list(
                setups=lapply(assemblies, function(assembly) {
                    .cdr_prediction_setup(assembly$setup)
                }),
                random_effects=lapply(assemblies, `[[`, 'random_prediction'),
                coefficient_ranges=ranges
            ),
            rank=list(
                action='error', resolution='identified', regularization=0,
                tolerance=.rank_tolerance(NULL)
            ),
            solver=paste(
                if (qncv) {
                    'sparse Gaussian location-scale QNCV solver'
                } else 'sparse Gaussian location-scale LAML solver',
                if (qncv) {
                    paste(
                        '(streamed', control$inner_algorithm,
                        'Fisher scoring, low-rank leave-out',
                        'updates, exact-score BFGS)'
                    )
                } else if (identical(control$gradient_method, 'hybrid')) {
                    paste(
                        '(streamed', control$inner_algorithm,
                        'Fisher scoring, matrix-free stochastic',
                        'warm-up, exact-score BFGS refinement)'
                    )
                } else {
                    paste(
                        '(streamed', control$inner_algorithm,
                        'Fisher scoring, exact-score BFGS)'
                    )
                }
            )
        )
    )
    class(output) <- c(
        'cdrgam_distributional_sparse', 'cdrgam_sparse', 'cdrgam'
    )
    reporter$emit(
        1L,
        'fit complete',
        criterion=format(solution$criterion, digits=10),
        evaluations=result$evaluations,
        score_evaluations=result$score_evaluations,
        stochastic_score_evaluations=result$stochastic_score_evaluations,
        warm_starts=result$warm_starts,
        cold_fallbacks=result$cold_fallbacks,
        converged=output$converged
    )
    output
}

.cdrgam_distributional_sparse_lpmatrices <- function(
        object, impulses, responses, chunk_size
) {
    output <- list()
    for (parameter in object$cdrgam$parameter_names) {
        term_indices <- object$cdrgam$parameter_terms[[parameter]]
        range <- object$cdrgam$prediction$coefficient_ranges[[parameter]]
        shim <- object
        shim$coefficients <- object$coefficients[range]
        shim$cdrgam$preparation <-
            object$cdrgam$preparation$parameters[[parameter]]
        shim$cdrgam$scaling <- shim$cdrgam$preparation$scaling
        shim$cdrgam$terms <- object$cdrgam$terms[term_indices]
        shim$cdrgam$prediction <- list(
            setup=object$cdrgam$prediction$setups[[parameter]],
            random_effects=object$cdrgam$prediction$random_effects[[parameter]]
        )
        class(shim) <- c('cdrgam_sparse', 'cdrgam')
        X <- .predict_cdrgam_streams(
            shim,
            impulses=impulses,
            responses=responses,
            type='lpmatrix',
            chunk_size=chunk_size
        )
        ordinary <- .cdr_setup_lpmatrix(
            object$cdrgam$prediction$setups[[parameter]], responses
        )
        output[[parameter]] <- list(X=X, offset=ordinary$offset)
    }
    output
}

#' @rdname predict.cdrgam
#' @export
predict.cdrgam_distributional_sparse <- function(
        object,
        newdata=NULL,
        type=c('response', 'link', 'lpmatrix'),
        se.fit=FALSE,
        chunk_size=10000,
        unconditional=FALSE,
        ...
) {
    type <- match.arg(type)
    if (isTRUE(unconditional)) {
        stop('Smoothing-parameter covariance is unavailable')
    }
    if (is.null(newdata)) {
        if (identical(type, 'link')) return(object$linear.predictors)
        if (identical(type, 'response')) return(object$fitted.values)
        stop('Training lpmatrices are not retained by the sparse backend')
    }
    if (!is.list(newdata) ||
            !all(c('impulses', 'responses') %in% names(newdata))) {
        stop('newdata must contain impulse and response data frames')
    }
    matrices <- .cdrgam_distributional_sparse_lpmatrices(
        object, newdata$impulses, newdata$responses, chunk_size
    )
    if (identical(type, 'lpmatrix')) {
        return(lapply(matrices, `[[`, 'X'))
    }
    links <- matrix(
        NA_real_, nrow(newdata$responses), length(matrices),
        dimnames=list(NULL, names(matrices))
    )
    standard_errors <- links
    for (parameter in names(matrices)) {
        range <- object$cdrgam$prediction$coefficient_ranges[[parameter]]
        assembled <- matrices[[parameter]]
        X <- assembled$X
        links[, parameter] <- as.numeric(
            X %*% object$coefficients[range] + assembled$offset
        )
        if (isTRUE(se.fit)) {
            selector <- Matrix::Matrix(
                0, nrow=length(object$coefficients), ncol=nrow(X), sparse=TRUE
            )
            selector[range, ] <- Matrix::t(X)
            solved <- .cdr_factor_solve(object$sparse$factor, selector)
            standard_errors[, parameter] <- sqrt(pmax(
                0, Matrix::rowSums(X * Matrix::t(solved[range, , drop=FALSE]))
            ))
        }
    }
    fit <- links
    if (identical(type, 'response')) {
        fit[, 'scale'] <- object$family$linfo[[2L]]$linkinv(links[, 'scale'])
        if (isTRUE(se.fit)) {
            standard_errors[, 'scale'] <- standard_errors[, 'scale'] *
                abs(object$family$linfo[[2L]]$mu.eta(links[, 'scale']))
        }
    }
    if (!isTRUE(se.fit)) return(fit)
    list(fit=fit, se.fit=standard_errors)
}

#' @export
logLik.cdrgam_distributional_sparse <- function(object, ...) {
    value <- object$loglik
    attr(value, 'df') <- object$edf
    attr(value, 'nobs') <- length(object$y)
    class(value) <- 'logLik'
    value
}

#' @export
summary.cdrgam_distributional_sparse <- function(
        object, dispersion=NULL, freq=FALSE, re.test=TRUE,
        all.coefficients=FALSE, ...
) {
    shim <- object
    dimension <- length(object$coefficients)
    if (dimension^2 > getOption('cdrgam.max_vcov_elements', 2e8)) {
        stop(
            'The summary covariance exceeds cdrgam.max_vcov_elements; ',
            'increase the option or request selected inference with estimate_irf()'
        )
    }
    shim$Vp <- as.matrix(.cdr_factor_solve(
        object$sparse$factor, Matrix::Diagonal(dimension)
    ))
    shim$coefficient.edf <- diag(
        shim$Vp %*% as.matrix(object$distributional$likelihood_hessian)
    )
    class(shim) <- c('cdrgam_distributional_block', 'cdrgam_block', 'cdrgam')
    summary.cdrgam_distributional_block(
        shim,
        dispersion=dispersion,
        freq=freq,
        re.test=re.test,
        all.coefficients=all.coefficients,
        ...
    )
}
