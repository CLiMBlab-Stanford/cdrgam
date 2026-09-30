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
    for (relative in c(1e-10, 1e-8, 1e-6, 1e-4)) {
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
        include_observed=FALSE
) {
    location <- assemblies$location
    scale <- assemblies$scale
    location_gradient <- numeric(location$dimension)
    scale_gradient <- numeric(scale$dimension)
    location_information <- NULL
    scale_information <- NULL
    cross_information <- NULL
    chunk_size <- min(
        location$crossprod_chunk_size,
        scale$crossprod_chunk_size
    )
    starts <- seq.int(1L, location$observation_count, by=chunk_size)
    for (start in starts) {
        rows <- start:min(location$observation_count, start + chunk_size - 1L)
        X_location <- .cdrgam_sparse_design_chunk(location, rows)
        X_scale <- .cdrgam_sparse_design_chunk(scale, rows)
        sigma_rows <- sigma[rows]
        q_rows <- q[rows]
        residual_rows <- residual[rows]
        weight_rows <- weights[rows]
        location_score <- weight_rows * (-residual_rows / sigma_rows^2)
        scale_score <- weight_rows * q_rows * (
            1 / sigma_rows - residual_rows^2 / sigma_rows^3
        )
        location_gradient <- location_gradient + as.numeric(
            Matrix::crossprod(X_location, location_score)
        )
        scale_gradient <- scale_gradient + as.numeric(
            Matrix::crossprod(X_scale, scale_score)
        )
        weighted_location <- Matrix::Diagonal(
            x=sqrt(weight_rows / sigma_rows^2)
        ) %*% X_location
        weighted_scale <- Matrix::Diagonal(
            x=sqrt(weight_rows * 2 * q_rows^2 / sigma_rows^2)
        ) %*% X_scale
        location_part <- Matrix::crossprod(weighted_location)
        scale_part <- Matrix::crossprod(weighted_scale)
        location_information <- if (is.null(location_information)) {
            location_part
        } else location_information + location_part
        scale_information <- if (is.null(scale_information)) {
            scale_part
        } else scale_information + scale_part
        if (isTRUE(include_observed)) {
            cross_weight <- weight_rows * 2 * q_rows * residual_rows /
                sigma_rows^3
            scale_weight <- weight_rows * (
                q_rows * (1 / sigma_rows - residual_rows^2 / sigma_rows^3) +
                    q_rows^2 * (
                        -1 / sigma_rows^2 +
                            3 * residual_rows^2 / sigma_rows^4
                    )
            )
            cross_part <- Matrix::crossprod(
                X_location,
                Matrix::Diagonal(x=cross_weight) %*% X_scale
            )
            observed_scale <- Matrix::crossprod(
                X_scale,
                Matrix::Diagonal(x=scale_weight) %*% X_scale
            )
            cross_information <- if (is.null(cross_information)) {
                cross_part
            } else cross_information + cross_part
            scale_information <- scale_information - scale_part + observed_scale
        }
    }
    list(
        gradients=list(location=location_gradient, scale=scale_gradient),
        information=list(
            location=location_information,
            scale=scale_information,
            cross=cross_information
        ),
        chunks=length(starts)
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

.cdrgam_gaulss_sparse_fixed <- function(
        assemblies, smoothing_parameters, b, tolerance=1e-5, maxit=200L,
        initial=NULL
) {
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
        moments <- .cdrgam_gaulss_sparse_moments(
            assemblies, coefficients, sigma, q, y - predictors$location,
            weights
        )
        location_system <- Matrix::forceSymmetric(
            moments$information$location + penalties$location, uplo='U'
        )
        location_rhs <- numeric(assemblies$location$dimension)
        for (start in seq.int(
                1L, length(y),
                by=assemblies$location$crossprod_chunk_size)) {
            rows <- start:min(
                length(y),
                start + assemblies$location$crossprod_chunk_size - 1L
            )
            X <- .cdrgam_sparse_design_chunk(assemblies$location, rows)
            location_rhs <- location_rhs + as.numeric(Matrix::crossprod(
                X,
                weights[rows] * (y[rows] - offsets$location[rows]) /
                    sigma[rows]^2
            ))
        }
        factor_location <- .cdrgam_gaulss_step_factor(
            location_system, assemblies$location$supernodal
        )$factor
        coefficients$location <- as.numeric(
            .cdr_factor_solve(factor_location, location_rhs)
        )
        scale_system <- NULL
        scale_rhs <- numeric(assemblies$scale$dimension)
        for (start in seq.int(
                1L, length(y), by=assemblies$scale$crossprod_chunk_size)) {
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
            tolerance=tolerance, maxit=maxit, initial=NULL
        ))
    }
    converged <- FALSE
    gradient_norm <- Inf
    termination <- 'iteration limit'
    step_norm <- NA_real_
    directional_derivative <- NA_real_
    best_line_search_change <- NA_real_
    step_mode <- NA_character_
    maximum_damping <- 0
    factor_location <- retained_factors$location
    factor_scale <- retained_factors$scale
    for (iteration in seq_len(maxit)) {
        q <- exp(predictors$scale)
        sigma <- b + q
        residual <- y - predictors$location
        moments <- .cdrgam_gaulss_sparse_moments(
            assemblies, coefficients, sigma, q, residual, weights
        )
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
    observed <- .cdrgam_gaulss_sparse_moments(
        assemblies, coefficients, sigma, q, residual, weights,
        include_observed=TRUE
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
    joint_factor <- .cdrgam_sparse_factor(
        hessian,
        assemblies$location$supernodal || assemblies$scale$supernodal
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
        maximum_damping=maximum_damping,
        prior_weights=weights,
        chunks=observed$chunks,
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
        coefficient_derivatives=NULL, hessian_derivatives=NULL
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
    delta_location <- numeric(observation_count)
    delta_scale <- numeric(observation_count)
    score <- numeric(derivative_count)
    minimum_determinant <- Inf
    for (start in seq.int(1L, observation_count, by=batch_size)) {
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
        minimum_determinant <- min(minimum_determinant, determinant)
        if (any(!is.finite(determinant)) || any(abs(determinant) < 1e-10)) {
            stop('QNCV leave-out system is numerically singular')
        }
        v1 <- a11 * location_score[rows] + a12 * scale_score[rows]
        v2 <- a12 * location_score[rows] + a22 * scale_score[rows]
        delta_location[rows] <- -(c22 * v1 - c12 * v2) / determinant
        delta_scale[rows] <- -(-c21 * v1 + c11 * v2) / determinant
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
            delta1 <- delta_location[rows]
            delta2 <- delta_scale[rows]
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
                score[[index]] <- score[[index]] + sum(
                    full_derivative - gamma * correction_derivative
                )
            }
        }
    }
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
        batch_size=batch_size
    )
}

.cdrgam_gaulss_qncv_direction_plan <- function(
        assemblies, direction_count, cores, workers,
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
    plan$requested_workers <- as.integer(requested_workers)
    plan$memory_workers <- as.integer(memory_workers)
    plan$memory_limited <- plan$workers < requested_workers
    plan$memory_source <- memory$source
    plan$memory_available_bytes <- memory$available_bytes
    plan$memory_budget_bytes <- memory_budget_bytes
    plan$projected_design_bytes <- as.numeric(projected_design_bytes)
    plan$direction_bytes <- as.numeric(direction_bytes)
    plan$cache_parallel <- cache_parallel
    plan
}

.cdrgam_gaulss_sparse_qncv_score <- function(
        assemblies, solution, sp, family, gamma=1, batch_size=NULL,
        cores=1L, workers=1L, score_batch_size=NULL
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
        length(tasks), workers, derivative_bytes, score_batch_size
    )
    direction_plan <- .cdrgam_gaulss_qncv_direction_plan(
        assemblies,
        direction_count=plan$batch_size,
        cores=cores,
        workers=plan$workers
    )
    plan$workers <- direction_plan$workers
    evaluate <- function(indices) tryCatch({
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
            hessian_derivatives=hessian_derivatives
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
    results <- .cdrgam_with_blas_threads(
        direction_plan$blas_threads,
        lapply(plan$groups, evaluate)
    )
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
        requested_workers=direction_plan$requested_workers,
        blas_threads=direction_plan$blas_threads,
        memory_limited=direction_plan$memory_limited,
        design_cache=direction_plan$cache_parallel,
        parallel_axis='QNCV Hessian directions and leave-out solves',
        batches=length(plan$groups),
        batch_size=plan$batch_size,
        response_batch_size=batch_size,
        derivative_bytes=plan$derivative_bytes,
        memory_source=plan$memory$source,
        memory_available_bytes=plan$memory$available_bytes,
        memory_budget_bytes=direction_plan$memory_budget_bytes,
        projected_design_bytes=direction_plan$projected_design_bytes,
        direction_bytes=direction_plan$direction_bytes,
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
    checkpoint_state <- .checkpoint_validate(
        .checkpoint_read(checkpoint),
        checkpoint_signature,
        checkpoint
    )
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
    if (resumed && !is.null(reporter)) reporter$emit(
        1L,
        'checkpoint resumed',
        path=checkpoint,
        stage=checkpoint_state$stage,
        phase=checkpoint_phase,
        evaluations=evaluations,
        best=checkpoint_state$best_criterion
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
                        score_batch_size=control$score_batch_size
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
            return(invalid_evaluation(log_sp, reason))
        }
        score_plan <<- attr(score, 'score_plan')
        value$score <- as.numeric(score)
        remember(log_sp, value)
        if (!is.null(best_parameters) && identical(
                as.numeric(log_sp), best_parameters)) {
            best_evaluation <<- value
        }
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
        if (any(!is.finite(log_sp)) || any(abs(log_sp) > 25)) {
            return(invalid_evaluation(
                log_sp, 'smoothing parameters left the finite search region'
            ))
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
            result <- tryCatch(
                list(solution=.cdrgam_gaulss_sparse_fixed(
                    assemblies, exp(log_sp), b,
                    tolerance=control$inner_tolerance,
                    maxit=control$inner_maxit,
                    initial=initial
                ), error=NULL),
                error=function(error) {
                    list(solution=NULL, error=conditionMessage(error))
                }
            )
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
        used_warm_start <- !is.null(warm_solution)
        if (used_warm_start) warm_starts <<- warm_starts + 1L
        fixed <- solve_fixed(warm_solution)
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
        inner_seconds <- proc.time()[['elapsed']] - evaluation_started
        if (!is.null(reporter)) reporter$emit(
            1L,
            'outer inner solve complete',
            evaluation=evaluations,
            seconds=format(inner_seconds, digits=5),
            iterations=solution$iterations,
            warm_start=used_warm_start,
            cold_fallback=used_warm_start && !isTRUE(solution$warm_started)
        )
        criterion_evaluation <- if (qncv) {
            .cdrgam_gaulss_sparse_qncv(
                assemblies,
                solution,
                gamma=control$gamma,
                batch_size=control$qncv_batch_size
            )
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
            all(is.finite(saved_optimizer_state$parameters))
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
        'inner_tolerance', 'inner_maxit', 'cores', 'gradient',
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
        inner_maxit=if (is.null(sparse_control$inner_maxit)) {
            200L
        } else as.integer(sparse_control$inner_maxit)
    )
    numeric_control <- unlist(control, use.names=FALSE)
    if (any(!is.finite(numeric_control)) || any(numeric_control <= 0)) {
        stop('Distributional sparse optimizer controls must be positive')
    }
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
        inner_solver='streamed sparse Fisher scoring'
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
                solution$factor
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
                        '(streamed Fisher scoring, low-rank leave-out',
                        'updates, exact-score BFGS)'
                    )
                } else if (identical(control$gradient_method, 'hybrid')) {
                    paste(
                        '(streamed Fisher scoring, matrix-free stochastic',
                        'warm-up, exact-score BFGS refinement)'
                    )
                } else {
                    '(streamed Fisher scoring, exact-score BFGS)'
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
