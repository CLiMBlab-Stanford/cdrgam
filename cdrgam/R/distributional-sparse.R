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

.cdrgam_gaulss_sparse_fixed <- function(
        assemblies, smoothing_parameters, b, tolerance=1e-5, maxit=200L
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
    initial_sigma <- sqrt(stats::weighted.mean(
        (y - stats::weighted.mean(y, weights))^2, weights
    ))
    initial_eta_scale <- rep.int(log(max(initial_sigma - b, 1e-4)), length(y))
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
            1L, length(y), by=assemblies$location$crossprod_chunk_size)) {
        rows <- start:min(
            length(y), start + assemblies$location$crossprod_chunk_size - 1L
        )
        X <- .cdrgam_sparse_design_chunk(assemblies$location, rows)
        location_rhs <- location_rhs + as.numeric(Matrix::crossprod(
            X, weights[rows] * (y[rows] - offsets$location[rows]) / sigma[rows]^2
        ))
    }
    factor_location <- .cdrgam_sparse_factor(
        location_system, assemblies$location$supernodal
    )
    coefficients$location <- as.numeric(
        .cdr_factor_solve(factor_location, location_rhs)
    )
    scale_system <- NULL
    scale_rhs <- numeric(assemblies$scale$dimension)
    for (start in seq.int(
            1L, length(y), by=assemblies$scale$crossprod_chunk_size)) {
        rows <- start:min(
            length(y), start + assemblies$scale$crossprod_chunk_size - 1L
        )
        X <- .cdrgam_sparse_design_chunk(assemblies$scale, rows)
        part <- Matrix::crossprod(
            Matrix::Diagonal(x=sqrt(weights[rows])) %*% X
        )
        scale_system <- if (is.null(scale_system)) part else scale_system + part
        scale_rhs <- scale_rhs + as.numeric(Matrix::crossprod(
            X, weights[rows] * (initial_eta_scale[rows] - offsets$scale[rows])
        ))
    }
    scale_system <- Matrix::forceSymmetric(
        scale_system + penalties$scale, uplo='U'
    )
    factor_scale <- .cdrgam_sparse_factor(
        scale_system, assemblies$scale$supernodal
    )
    coefficients$scale <- as.numeric(
        .cdr_factor_solve(factor_scale, scale_rhs)
    )
    objective <- function(coefficients) {
        eta_location <- .cdrgam_distributional_sparse_predictor(
            assemblies$location, coefficients$location, offsets$location
        )
        eta_scale <- .cdrgam_distributional_sparse_predictor(
            assemblies$scale, coefficients$scale, offsets$scale
        )
        sigma <- b + exp(eta_scale)
        if (any(!is.finite(sigma))) return(Inf)
        residual <- y - eta_location
        sum(weights * (log(sigma) + 0.5 * (residual / sigma)^2)) +
            0.5 * as.numeric(Matrix::crossprod(
                coefficients$location, penalties$location %*% coefficients$location
            )) +
            0.5 * as.numeric(Matrix::crossprod(
                coefficients$scale, penalties$scale %*% coefficients$scale
            ))
    }
    value <- objective(coefficients)
    converged <- FALSE
    gradient_norm <- Inf
    factor_location <- NULL
    factor_scale <- NULL
    for (iteration in seq_len(maxit)) {
        predictors <- Map(
            .cdrgam_distributional_sparse_predictor,
            assemblies,
            coefficients,
            offsets
        )
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
        systems <- list(
            location=Matrix::forceSymmetric(
                moments$information$location + penalties$location, uplo='U'
            ),
            scale=Matrix::forceSymmetric(
                moments$information$scale + penalties$scale, uplo='U'
            )
        )
        factor_location <- .cdrgam_sparse_factor(
            systems$location, assemblies$location$supernodal, factor_location
        )
        factor_scale <- .cdrgam_sparse_factor(
            systems$scale, assemblies$scale$supernodal, factor_scale
        )
        steps <- list(
            location=-as.numeric(.cdr_factor_solve(
                factor_location, gradients$location
            )),
            scale=-as.numeric(.cdr_factor_solve(factor_scale, gradients$scale))
        )
        step_size <- 1
        accepted <- FALSE
        while (step_size >= 2^-20) {
            candidate <- Map(function(current, step) {
                current + step_size * step
            }, coefficients, steps)
            candidate_value <- objective(candidate)
            if (is.finite(candidate_value) && candidate_value < value) {
                coefficients <- candidate
                value <- candidate_value
                accepted <- TRUE
                break
            }
            step_size <- step_size / 2
        }
        if (gradient_norm < tolerance ||
                (accepted && max(abs(step_size * unlist(
                    steps, use.names=FALSE
                ))) < tolerance / 100)) {
            converged <- TRUE
            break
        }
        if (!accepted) break
    }
    predictors <- Map(
        .cdrgam_distributional_sparse_predictor,
        assemblies,
        coefficients,
        offsets
    )
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
        converged=converged,
        iterations=iteration,
        gradient_norm=gradient_norm,
        prior_weights=weights,
        chunks=observed$chunks
    )
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
        d_q <- q_rows * d_scale
        d_sigma <- d_q
        d_residual <- -d_location
        d_location_weight <- -2 * weight_rows * d_sigma / sigma_rows^3
        d_cross_weight <- 2 * weight_rows * (
            (d_q * residual_rows + q_rows * d_residual) / sigma_rows^3 -
                3 * q_rows * residual_rows * d_sigma / sigma_rows^4
        )
        d_scale_dr <- weight_rows * (
            -2 * q_rows * residual_rows / sigma_rows^3 +
                6 * q_rows^2 * residual_rows / sigma_rows^4
        )
        d_scale_dq <- weight_rows * (
            1 / sigma_rows - q_rows / sigma_rows^2 -
                residual_rows^2 / sigma_rows^3 +
                3 * q_rows * residual_rows^2 / sigma_rows^4 -
                2 * q_rows / sigma_rows^2 +
                2 * q_rows^2 / sigma_rows^3 +
                6 * q_rows * residual_rows^2 / sigma_rows^4 -
                12 * q_rows^2 * residual_rows^2 / sigma_rows^5
        )
        d_scale_weight <- d_scale_dr * d_residual + d_scale_dq * d_q
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

.cdrgam_gaulss_sparse_score <- function(assemblies, solution, sp, family) {
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
    scores <- numeric(sum(counts))
    beta <- solution$coefficients
    b <- .cdrgam_gaulss_b(family)
    global_index <- 0L
    for (parameter in names(assemblies)) {
        for (local_index in seq_len(counts[[parameter]])) {
            global_index <- global_index + 1L
            penalty_derivative <- .cdrgam_distributional_global_penalty(
                assemblies,
                ranges,
                parameter,
                local_index,
                sp[[global_index]]
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
            scores[[global_index]] <- as.numeric(Matrix::crossprod(
                beta, penalty_derivative %*% beta
            )) + .sparse_logdet_score(
                solution$factor, hessian_derivative, chunk_size=256L
            ) - penalty_scores[[parameter]][[local_index]]
        }
    }
    scores
}

.cdrgam_gaulss_sparse_laml <- function(assemblies, family, control) {
    counts <- vapply(
        assemblies, function(assembly) length(assembly$penalty_components),
        integer(1)
    )
    initial <- unlist(lapply(
        assemblies, .cdrgam_sparse_initial_sp
    ), use.names=FALSE)
    initial_log_sp <- log(initial)
    b <- .cdrgam_gaulss_b(family)
    evaluations <- 0L
    cached_parameters <- NULL
    cached_evaluation <- NULL
    evaluate <- function(log_sp) {
        if (!is.null(cached_parameters) && identical(
                as.numeric(log_sp), cached_parameters)) {
            return(cached_evaluation)
        }
        evaluations <<- evaluations + 1L
        if (any(!is.finite(log_sp)) || any(abs(log_sp) > 25)) {
            return(list(criterion=1e100, score=rep.int(0, length(log_sp))))
        }
        solution <- tryCatch(
            .cdrgam_gaulss_sparse_fixed(
                assemblies, exp(log_sp), b,
                tolerance=control$inner_tolerance,
                maxit=control$inner_maxit
            ),
            error=function(error) NULL
        )
        if (is.null(solution) || !isTRUE(solution$converged)) {
            return(list(criterion=1e100, score=rep.int(0, length(log_sp))))
        }
        penalty_determinant <- 0
        offset <- 0L
        for (parameter in names(assemblies)) {
            count <- counts[[parameter]]
            indices <- offset + seq_len(count)
            penalty_determinant <- penalty_determinant +
                .sparse_penalty_logdet(
                    assemblies[[parameter]]$blocks, exp(log_sp[indices])
                )$value
            offset <- offset + count
        }
        criterion <- 2 * solution$objective +
            .cdr_factor_logdet(solution$factor) - penalty_determinant
        if (!is.finite(criterion)) {
            return(list(criterion=1e100, score=rep.int(0, length(log_sp))))
        }
        solution$criterion <- criterion
        solution$sp <- exp(log_sp)
        value <- list(
            criterion=criterion,
            score=.cdrgam_gaulss_sparse_score(
                assemblies, solution, solution$sp, family
            ),
            solution=solution
        )
        cached_parameters <<- as.numeric(log_sp)
        cached_evaluation <<- value
        value
    }
    if (!length(initial_log_sp)) {
        evaluated <- evaluate(numeric())
        solution <- evaluated$solution
        return(list(
            optimization=list(
                par=numeric(), value=solution$criterion,
                convergence=0L,
                counts=c(`function`=1L, gradient=NA_integer_)
            ),
            solution=solution,
            evaluations=evaluations
        ))
    }
    optimization <- stats::optim(
        par=initial_log_sp,
        fn=function(parameters) evaluate(parameters)$criterion,
        gr=function(parameters) evaluate(parameters)$score,
        method='L-BFGS-B',
        lower=rep.int(-25, length(initial_log_sp)),
        upper=rep.int(25, length(initial_log_sp)),
        control=list(
            maxit=control$optimizer_maxit,
            pgtol=control$optimizer_gradient_tolerance,
            factr=1e7
        )
    )
    optimization$gradient <- evaluate(optimization$par)$score
    retained <- evaluate(optimization$par)
    solution <- retained$solution
    if (is.null(solution)) {
        stop('Sparse gaulss LAML optimization ended at an invalid fit')
    }
    list(
        optimization=optimization,
        solution=solution,
        evaluations=evaluations
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
        design, family, method=NULL, trace=FALSE, sparse_control=list(), ...
) {
    dots <- list(...)
    if (!is.null(dots$weights) && any(dots$weights != 1)) {
        stop(
            'gaulss does not support non-unit prior weights; mgcv accepts ',
            'but ignores them'
        )
    }
    if (!is.null(method) && !(method %in% c('REML', 'fREML'))) {
        stop('The distributional sparse backend currently supports only REML')
    }
    allowed <- c(
        'crossprod_chunk_size', 'supernodal', 'optimizer_maxit',
        'optimizer_gradient_tolerance', 'inner_tolerance', 'inner_maxit'
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
        )) 1e-4 else sparse_control$optimizer_gradient_tolerance,
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
    y <- assemblies$location$setup$y
    if (!identical(y, assemblies$scale$setup$y)) {
        stop('Distributional parameter setups retained different responses')
    }
    reporter$phase(
        'joint smoothing-parameter optimization',
        smoothing_parameters=sum(vapply(
            assemblies,
            function(assembly) length(assembly$penalty_components),
            integer(1)
        )),
        gradient='exact',
        inner_solver='streamed sparse Fisher scoring'
    )
    result <- .cdrgam_gaulss_sparse_laml(assemblies, family, control)
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
        method='REML',
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
        reml=solution$criterion,
        optimizer=result$optimization,
        distributional=list(
            sigma=solution$sigma,
            coefficient_ranges=ranges,
            inner_iterations=solution$iterations,
            gradient_norm=solution$gradient_norm,
            likelihood_hessian=solution$likelihood_hessian
        ),
        sparse=list(
            factor=solution$factor,
            dimension=sum(dimensions),
            effective_df=effective_df,
            penalty_trace=penalty_trace,
            control=c(
                control,
                list(
                    gradient='exact',
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
                'sparse Gaussian location-scale LAML reference solver',
                '(streamed Fisher scoring, exact outer score)'
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
