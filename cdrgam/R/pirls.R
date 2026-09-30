.cdrgam_valid_family_state <- function(family, eta, mu) {
    all(is.finite(eta)) && all(is.finite(mu)) &&
        (is.null(family$valideta) || family$valideta(eta)) &&
        (is.null(family$validmu) || family$validmu(mu))
}

.cdrgam_initial_mu <- function(family, y, weights) {
    nobs <- length(y)
    n <- rep.int(1, nobs)
    mustart <- NULL
    etastart <- start <- NULL
    eval(family$initialize)
    if (is.null(mustart) || length(mustart) != nobs ||
            any(!is.finite(mustart))) {
        stop('The family did not produce valid initial fitted means')
    }
    mustart
}

.cdrgam_gamma_saturated_loglik <- function(y, weights, dispersion) {
    positive <- weights > 0
    y <- y[positive]
    weights <- weights[positive]
    local_dispersion <- dispersion / weights
    inverse <- 1 / local_dispersion
    value <- sum(
        -lgamma(inverse) - log(local_dispersion) * inverse - inverse -
            log(y)
    )
    derivative <- sum(
        (digamma(inverse) + log(local_dispersion)) /
            local_dispersion^2 / weights
    )
    list(value=value, derivative=derivative)
}

.cdrgam_gamma_loglik <- function(y, mu, weights, dispersion) {
    positive <- weights > 0
    sum(weights[positive] * stats::dgamma(
        y[positive],
        shape=1 / dispersion,
        scale=mu[positive] * dispersion,
        log=TRUE
    ))
}

.cdrgam_gamma_reported_scale <- function(
        y, mu, weights, effective_df
) {
    pearson <- sum(weights * (y - mu)^2 / mu^2) /
        (length(y) - effective_df)
    fletcher_adjustment <- max(-0.9, mean(2 * (y - mu) / mu))
    pearson / (1 + fletcher_adjustment)
}

# Dense penalized IRLS for a fixed collection of smoothing parameters. This
# deliberately contains no smoothing-parameter logic: it is the numerical
# reference inner solve shared by generalized marginal-likelihood tests.
.cdrgam_dense_pirls <- function(
        setup,
        family,
        smoothing_parameters,
        fixed_ridge=0,
        initial_coefficients=NULL,
        tolerance=1e-8,
        max_iterations=100L,
        maximum_halvings=25L
) {
    family <- .as_family(family)
    X <- setup$X
    y <- setup$y
    offset <- setup$offset
    if (is.null(offset)) offset <- numeric(length(y))
    prior_weights <- setup$w
    if (is.null(prior_weights)) prior_weights <- rep.int(1, length(y))
    if (any(!is.finite(prior_weights)) || any(prior_weights < 0)) {
        stop('PIRLS requires finite nonnegative prior weights')
    }
    penalty <- .embed_penalties(setup, smoothing_parameters)
    dimension <- ncol(X)
    coefficients <- if (is.null(initial_coefficients)) {
        numeric(dimension)
    } else {
        if (length(initial_coefficients) != dimension ||
                any(!is.finite(initial_coefficients))) {
            stop('initial_coefficients has the wrong dimension or values')
        }
        as.numeric(initial_coefficients)
    }
    initial_mu <- .cdrgam_initial_mu(family, y, prior_weights)
    eta <- family$linkfun(initial_mu)
    objective <- Inf
    converged <- FALSE
    factor <- NULL
    working_weights <- numeric(length(y))
    for (iteration in seq_len(max_iterations)) {
        mu <- family$linkinv(eta)
        derivative <- family$mu.eta(eta)
        variance <- family$variance(mu)
        valid <- .cdrgam_valid_family_state(family, eta, mu) &&
            all(is.finite(derivative)) && all(is.finite(variance)) &&
            all(variance > 0) && all(derivative != 0)
        if (!valid) stop('PIRLS encountered an invalid family state')
        working_weights <- prior_weights * derivative^2 / variance
        working_response <- eta - offset + (y - mu) / derivative
        root_weights <- sqrt(working_weights)
        weighted_X <- X * root_weights
        system <- crossprod(weighted_X) + penalty +
            diag(fixed_ridge, dimension)
        rhs <- crossprod(weighted_X, working_response * root_weights)
        factor <- tryCatch(chol(system), error=function(error) NULL)
        if (is.null(factor)) {
            stop('The penalized PIRLS system is not positive definite')
        }
        proposal <- drop(backsolve(factor, forwardsolve(t(factor), rhs)))
        accepted <- FALSE
        halvings <- 0L
        repeat {
            candidate_eta <- drop(offset + X %*% proposal)
            candidate_mu <- family$linkinv(candidate_eta)
            candidate_objective <- if (.cdrgam_valid_family_state(
                    family, candidate_eta, candidate_mu)) {
                sum(family$dev.resids(y, candidate_mu, prior_weights)) +
                    drop(crossprod(proposal, penalty %*% proposal)) +
                    fixed_ridge * sum(proposal^2)
            } else Inf
            if (is.finite(candidate_objective) &&
                    (!is.finite(objective) ||
                        candidate_objective <= objective +
                            tolerance * (1 + abs(objective)))) {
                accepted <- TRUE
                break
            }
            if (halvings >= maximum_halvings) break
            proposal <- (proposal + coefficients) / 2
            halvings <- halvings + 1L
        }
        if (!accepted) {
            stop('PIRLS step halving failed to improve the penalized deviance')
        }
        change <- if (is.finite(objective)) {
            abs(objective - candidate_objective) / (1 + abs(objective))
        } else Inf
        coefficients <- proposal
        eta <- candidate_eta
        objective <- candidate_objective
        if (is.finite(change) && change < tolerance) {
            converged <- TRUE
            break
        }
    }
    mu <- family$linkinv(eta)
    derivative <- family$mu.eta(eta)
    variance <- family$variance(mu)
    working_weights <- prior_weights * derivative^2 / variance
    weighted_X <- X * sqrt(working_weights)
    system <- crossprod(weighted_X) + penalty +
        diag(fixed_ridge, dimension)
    factor <- chol(system)
    list(
        coefficients=coefficients,
        linear_predictors=eta,
        fitted_values=mu,
        residuals=y - mu,
        deviance=sum(family$dev.resids(y, mu, prior_weights)),
        penalized_deviance=objective,
        penalty=penalty,
        system=system,
        factor=factor,
        working_weights=working_weights,
        prior_weights=prior_weights,
        iterations=iteration,
        converged=converged
    )
}

# Sparse fixed-smoothing PIRLS reference. The current caller supplies an mgcv
# setup matrix; the production integration will replace that input with the
# streamed component assembler already used by the Gaussian sparse backend.
.cdrgam_sparse_pirls <- function(
        setup,
        family,
        smoothing_parameters,
        fixed_ridge=0,
        initial_coefficients=NULL,
        tolerance=1e-8,
        max_iterations=100L,
        maximum_halvings=25L,
        supernodal=FALSE
) {
    family <- .as_family(family)
    X <- Matrix::Matrix(setup$X, sparse=TRUE)
    y <- setup$y
    offset <- setup$offset
    if (is.null(offset)) offset <- numeric(length(y))
    prior_weights <- setup$w
    if (is.null(prior_weights)) prior_weights <- rep.int(1, length(y))
    if (any(!is.finite(prior_weights)) || any(prior_weights < 0)) {
        stop('Sparse PIRLS requires finite nonnegative prior weights')
    }
    penalty <- Matrix::Matrix(
        .embed_penalties(setup, smoothing_parameters),
        sparse=TRUE
    )
    dimension <- ncol(X)
    coefficients <- if (is.null(initial_coefficients)) {
        numeric(dimension)
    } else {
        if (length(initial_coefficients) != dimension ||
                any(!is.finite(initial_coefficients))) {
            stop('initial_coefficients has the wrong dimension or values')
        }
        as.numeric(initial_coefficients)
    }
    eta <- family$linkfun(.cdrgam_initial_mu(family, y, prior_weights))
    objective <- Inf
    converged <- FALSE
    factor <- NULL
    numeric_updates <- 0L
    for (iteration in seq_len(max_iterations)) {
        mu <- family$linkinv(eta)
        derivative <- family$mu.eta(eta)
        variance <- family$variance(mu)
        valid <- .cdrgam_valid_family_state(family, eta, mu) &&
            all(is.finite(derivative)) && all(is.finite(variance)) &&
            all(variance > 0) && all(derivative != 0)
        if (!valid) stop('Sparse PIRLS encountered an invalid family state')
        working_weights <- prior_weights * derivative^2 / variance
        working_response <- eta - offset + (y - mu) / derivative
        weighted_X <- Matrix::Diagonal(x=sqrt(working_weights)) %*% X
        system <- Matrix::forceSymmetric(
            Matrix::crossprod(weighted_X) + penalty +
                Matrix::Diagonal(dimension, x=fixed_ridge),
            uplo='U'
        )
        rhs <- as.numeric(Matrix::crossprod(
            weighted_X,
            working_response * sqrt(working_weights)
        ))
        factor <- if (is.null(factor)) {
            Matrix::Cholesky(
                system,
                perm=TRUE,
                LDL=FALSE,
                super=supernodal
            )
        } else {
            numeric_updates <- numeric_updates + 1L
            tryCatch(
                Matrix::update(factor, system),
                error=function(error) Matrix::Cholesky(
                    system,
                    perm=TRUE,
                    LDL=FALSE,
                    super=supernodal
                )
            )
        }
        proposal <- as.numeric(.cdr_factor_solve(factor, rhs))
        accepted <- FALSE
        halvings <- 0L
        repeat {
            candidate_eta <- as.numeric(offset + X %*% proposal)
            candidate_mu <- family$linkinv(candidate_eta)
            candidate_objective <- if (.cdrgam_valid_family_state(
                    family, candidate_eta, candidate_mu)) {
                sum(family$dev.resids(y, candidate_mu, prior_weights)) +
                    as.numeric(Matrix::crossprod(
                        proposal,
                        penalty %*% proposal
                    )) + fixed_ridge * sum(proposal^2)
            } else Inf
            if (is.finite(candidate_objective) &&
                    (!is.finite(objective) ||
                        candidate_objective <= objective +
                            tolerance * (1 + abs(objective)))) {
                accepted <- TRUE
                break
            }
            if (halvings >= maximum_halvings) break
            proposal <- (proposal + coefficients) / 2
            halvings <- halvings + 1L
        }
        if (!accepted) {
            stop(
                'Sparse PIRLS step halving failed to improve the ',
                'penalized deviance'
            )
        }
        change <- if (is.finite(objective)) {
            abs(objective - candidate_objective) / (1 + abs(objective))
        } else Inf
        coefficients <- proposal
        eta <- candidate_eta
        objective <- candidate_objective
        if (is.finite(change) && change < tolerance) {
            converged <- TRUE
            break
        }
    }
    mu <- family$linkinv(eta)
    list(
        coefficients=coefficients,
        linear_predictors=eta,
        fitted_values=mu,
        residuals=y - mu,
        deviance=sum(family$dev.resids(y, mu, prior_weights)),
        penalized_deviance=objective,
        penalty=penalty,
        system=system,
        factor=factor,
        working_weights=working_weights,
        prior_weights=prior_weights,
        iterations=iteration,
        numeric_updates=numeric_updates,
        converged=converged
    )
}

.cdrgam_streamed_sparse_laml <- function(
        assembly,
        family,
        log_sp,
        initial_coefficients=NULL,
        initial_solution=NULL,
        fixed_ridge=0,
        tolerance=1e-9,
        chunk_size=assembly$crossprod_chunk_size,
        score=FALSE,
        score_workers=1L,
        score_batch_size=NULL,
        score_blas_threads=1L
) {
    family <- .as_family(family)
    canonical <- (identical(family$family, 'binomial') &&
            identical(family$link, 'logit')) ||
        (identical(family$family, 'poisson') &&
            identical(family$link, 'log'))
    estimated_gamma <- identical(family$family, 'Gamma') &&
        identical(family$link, 'log')
    if (!canonical && !estimated_gamma) {
        stop(
            'Sparse generalized LAML currently requires binomial(logit) ',
            ', poisson(log), or Gamma(log)'
        )
    }
    penalty_count <- length(assembly$penalty_components)
    parameter_count <- penalty_count + as.integer(estimated_gamma)
    if (length(log_sp) != parameter_count) {
        stop('log_sp has the wrong dimension for this family')
    }
    smoothing_log_parameters <- log_sp[seq_len(penalty_count)]
    sp <- exp(smoothing_log_parameters)
    dispersion <- if (estimated_gamma) {
        exp(log_sp[[parameter_count]])
    } else 1
    solution <- .cdrgam_streamed_sparse_pirls(
        assembly,
        family,
        sp,
        fixed_ridge=fixed_ridge,
        initial_coefficients=if (!is.null(initial_solution)) {
            initial_solution$coefficients
        } else initial_coefficients,
        initial_factor=if (!is.null(initial_solution)) {
            initial_solution$factor
        } else NULL,
        initial_linear_predictors=if (!is.null(initial_solution)) {
            initial_solution$linear_predictors
        } else NULL,
        tolerance=tolerance,
        chunk_size=chunk_size
    )
    if (!isTRUE(solution$converged)) {
        stop('Sparse PIRLS did not converge')
    }
    penalty_determinant <- .sparse_penalty_logdet(assembly$blocks, sp)
    criterion <- solution$penalized_deviance / dispersion +
        .cdr_factor_logdet(solution$factor) - penalty_determinant$value
    saturated <- NULL
    if (estimated_gamma) {
        saturated <- .cdrgam_gamma_saturated_loglik(
            assembly$setup$y,
            solution$prior_weights,
            dispersion
        )
        nullity <- assembly$dimension - penalty_determinant$rank
        criterion <- criterion - 2 * saturated$value -
            nullity * log(2 * pi * dispersion)
    }
    output <- list(
        criterion=criterion,
        sp=sp,
        scale=dispersion,
        solution=solution,
        saturated=saturated,
        penalty_logdet=penalty_determinant$value,
        penalty_rank=penalty_determinant$rank
    )
    if (!isTRUE(score)) return(output)
    .cdrgam_streamed_sparse_laml_score(
        assembly,
        family,
        output,
        workers=score_workers,
        batch_size=score_batch_size,
        blas_threads=score_blas_threads
    )
}

.cdrgam_streamed_sparse_laml_score <- function(
        assembly, family, evaluation, workers=1L, batch_size=NULL,
        blas_threads=1L
) {
    scores <- .cdrgam_with_blas_threads(
        blas_threads,
        .cdrgam_streamed_sparse_score(
            assembly,
            evaluation$solution,
            evaluation$sp,
            family,
            evaluation$scale,
            workers=workers,
            batch_size=batch_size
        )
    )
    estimated_gamma <- identical(family$family, 'Gamma') &&
        identical(family$link, 'log')
    if (estimated_gamma) {
        parameter_count <- length(assembly$penalty_components) + 1L
        nullity <- assembly$dimension - evaluation$penalty_rank
        scores[[parameter_count]] <-
            -evaluation$solution$penalized_deviance / evaluation$scale -
            2 * evaluation$saturated$derivative * evaluation$scale - nullity
    }
    evaluation$score_plan <- attr(scores, 'score_plan')
    evaluation$score <- as.numeric(scores)
    evaluation
}

.cdrgam_streamed_sparse_score <- function(
        assembly, solution, sp, family, dispersion=1,
        workers=1L, batch_size=NULL
) {
    count <- length(assembly$penalty_components)
    if (!count) return(numeric())
    penalty_scores <- .sparse_penalty_logdet_score(
        assembly$blocks, sp, count
    )
    mu <- solution$fitted_values
    prior_weights <- solution$prior_weights
    weight_derivative <- if (identical(family$family, 'poisson')) {
        prior_weights * mu
    } else if (identical(family$family, 'binomial')) {
        prior_weights * mu * (1 - mu) * (1 - 2 * mu)
    } else numeric(length(mu))
    derivative_bytes <- max(
        1,
        16 * as.double(Matrix::nnzero(solution$system)) +
            8 * as.double(assembly$observation_count)
    )
    plan <- .cdrgam_sparse_score_batch_plan(
        count, workers, derivative_bytes, batch_size
    )
    groups <- plan$groups
    workers <- plan$workers
    evaluate <- function(indices) tryCatch({
        derivatives <- lapply(indices, function(index) {
            sp[[index]] * assembly$penalty_components[[index]]
        })
        values <- .cdrgam_sparse_exact_score_batch(
            solution$factor,
            solution$coefficients,
            derivatives,
            penalty_scores[indices],
            likelihood_derivatives=function(coefficient_derivatives) {
                output <- vector('list', ncol(coefficient_derivatives))
                starts <- seq.int(
                    1L, assembly$observation_count,
                    by=min(
                        assembly$crossprod_chunk_size,
                        assembly$observation_count
                    )
                )
                for (start in starts) {
                    rows <- start:min(
                        assembly$observation_count,
                        start + assembly$crossprod_chunk_size - 1L
                    )
                    X <- .cdrgam_sparse_design_chunk(assembly, rows)
                    eta_derivatives <- as.matrix(
                        X %*% coefficient_derivatives
                    )
                    for (direction in seq_len(ncol(coefficient_derivatives))) {
                        contribution <- Matrix::crossprod(
                            X,
                            Matrix::Diagonal(x=
                                weight_derivative[rows] *
                                eta_derivatives[, direction]
                            ) %*% X
                        )
                        output[[direction]] <- if (is.null(
                                output[[direction]])) {
                            contribution
                        } else output[[direction]] + contribution
                    }
                }
                output
            },
            quadratic_scale=1 / dispersion,
            trace_workers=workers
        )
        list(indices=indices, values=values, error=NULL)
    }, error=function(error) {
        list(indices=indices, values=NULL, error=conditionMessage(error))
    })
    results <- lapply(groups, evaluate)
    missing <- which(!vapply(results, function(result) {
        is.list(result) &&
            all(c('indices', 'values', 'error') %in% names(result))
    }, logical(1)))
    if (length(missing)) {
        stop(
            'Exact generalized score worker ', missing[[1L]],
            ' did not return a result; it may have exceeded its memory limit'
        )
    }
    failed <- which(vapply(
        results, function(result) !is.null(result$error), logical(1)
    ))
    if (length(failed)) {
        first <- results[[failed[[1L]]]]
        stop(
            'Exact generalized score batch ',
            paste(first$indices, collapse=', '),
            ' failed: ', first$error
        )
    }
    output <- numeric(count)
    for (result in results) output[result$indices] <- result$values
    attr(output, 'score_plan') <- list(
        workers=workers,
        parallel_axis='inverse columns',
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

.cdrgam_optimize_streamed_sparse_laml <- function(
        assembly,
        family,
        initial_log_sp=NULL,
        fixed_ridge=0,
        tolerance=1e-9,
        chunk_size=assembly$crossprod_chunk_size,
        max_iterations=100L,
        gradient_tolerance=1e-4,
        trust_radius=2,
        score_workers=1L,
        score_batch_size=NULL,
        score_blas_threads=1L,
        reporter=NULL
) {
    family <- .as_family(family)
    estimated_gamma <- identical(family$family, 'Gamma') &&
        identical(family$link, 'log')
    parameter_count <- length(assembly$penalty_components) +
        as.integer(estimated_gamma)
    if (is.null(initial_log_sp)) {
        initial_log_sp <- numeric(parameter_count)
        if (estimated_gamma && parameter_count > 1L) {
            initial_log_sp[seq_len(parameter_count - 1L)] <-
                log(.cdrgam_sparse_initial_sp(assembly))
        }
    }
    if (length(initial_log_sp) != parameter_count ||
            any(!is.finite(initial_log_sp))) {
        stop('initial_log_sp has the wrong dimension or values')
    }
    cached_parameters <- NULL
    cached_evaluation <- NULL
    warm_solution <- NULL
    warm_starts <- 0L
    cold_fallbacks <- 0L
    evaluations <- 0L
    score_evaluations <- 0L
    score_plan <- NULL
    score_evaluation <- function(parameters, value) {
        if (!is.null(value$score)) return(value)
        score_evaluations <<- score_evaluations + 1L
        started <- proc.time()[['elapsed']]
        value <- .cdrgam_streamed_sparse_laml_score(
            assembly,
            family,
            value,
            workers=score_workers,
            batch_size=score_batch_size,
            blas_threads=score_blas_threads
        )
        score_plan <<- value$score_plan
        cached_parameters <<- as.numeric(parameters)
        cached_evaluation <<- value
        if (!is.null(reporter)) reporter$emit(
            1L,
            'outer exact score complete',
            evaluation=evaluations,
            score_evaluation=score_evaluations,
            seconds=format(proc.time()[['elapsed']] - started, digits=5),
            maximum_score=if (length(value$score)) {
                format(max(abs(value$score)), digits=5)
            } else 'none'
        )
        value
    }
    evaluate <- function(parameters, need_score=FALSE) {
        if (!is.null(cached_parameters) && identical(
                as.numeric(parameters), cached_parameters)) {
            if (need_score) {
                return(score_evaluation(parameters, cached_evaluation))
            }
            return(cached_evaluation)
        }
        evaluations <<- evaluations + 1L
        evaluation_started <- proc.time()[['elapsed']]
        solve <- function(initial) tryCatch(
            list(value=.cdrgam_streamed_sparse_laml(
                assembly,
                family,
                parameters,
                initial_solution=initial,
                fixed_ridge=fixed_ridge,
                tolerance=tolerance,
                chunk_size=chunk_size,
                score=FALSE
            ), error=NULL),
            error=function(error) {
                list(value=NULL, error=conditionMessage(error))
            }
        )
        used_warm_start <- !is.null(warm_solution)
        if (used_warm_start) warm_starts <<- warm_starts + 1L
        solved <- solve(warm_solution)
        if (used_warm_start && is.null(solved$value)) {
            cold_fallbacks <<- cold_fallbacks + 1L
            solved <- solve(NULL)
        }
        if (is.null(solved$value)) stop(solved$error)
        value <- solved$value
        warm_solution <<- value$solution
        if (!is.null(reporter)) reporter$emit(
            1L,
            'outer inner solve complete',
            evaluation=evaluations,
            seconds=format(
                proc.time()[['elapsed']] - evaluation_started,
                digits=5
            ),
            iterations=value$solution$iterations,
            warm_start=used_warm_start,
            cold_fallback=used_warm_start &&
                !isTRUE(value$solution$warm_started),
            criterion=format(value$criterion, digits=10)
        )
        cached_parameters <<- as.numeric(parameters)
        cached_evaluation <<- value
        if (need_score) score_evaluation(parameters, value) else value
    }
    objective <- function(parameters) evaluate(parameters)$criterion
    gradient <- function(parameters) {
        evaluate(parameters, need_score=TRUE)$score
    }
    optimization <- if (estimated_gamma) {
        value <- stats::optim(
            par=initial_log_sp,
            fn=objective,
            gr=gradient,
            method='L-BFGS-B',
            lower=rep.int(-25, parameter_count),
            upper=rep.int(25, parameter_count),
            control=list(maxit=max_iterations, pgtol=gradient_tolerance)
        )
        value$gradient <- gradient(value$par)
        value
    } else {
        .safeguarded_outer_bfgs(
            par=initial_log_sp,
            fn=objective,
            gr=gradient,
            lower=rep.int(-25, parameter_count),
            upper=rep.int(25, parameter_count),
            maxit=max_iterations,
            gradient_tolerance=gradient_tolerance,
            initial_radius=trust_radius
        )
    }
    retained <- evaluate(optimization$par, need_score=TRUE)
    list(
        optimization=optimization,
        retained=retained,
        evaluations=evaluations,
        score_evaluations=score_evaluations,
        warm_starts=warm_starts,
        cold_fallbacks=cold_fallbacks,
        score_plan=score_plan
    )
}

.cdrgam_streamed_sparse_prediction_error <- function(
        assembly,
        family,
        log_sp,
        initial_solution=NULL,
        fixed_ridge=0,
        gamma=1,
        tolerance=1e-9,
        chunk_size=assembly$crossprod_chunk_size
) {
    family <- .as_family(family)
    canonical <- (identical(family$family, 'binomial') &&
            identical(family$link, 'logit')) ||
        (identical(family$family, 'poisson') &&
            identical(family$link, 'log'))
    estimated_gamma <- identical(family$family, 'Gamma') &&
        identical(family$link, 'log')
    if (!canonical && !estimated_gamma) {
        stop(
            'Sparse generalized GCV.Cp currently requires ',
            'binomial(logit), poisson(log), or Gamma(log)'
        )
    }
    if (length(log_sp) != length(assembly$penalty_components)) {
        stop('log_sp has the wrong smoothing-parameter dimension')
    }
    if (length(gamma) != 1L || !is.numeric(gamma) ||
            !is.finite(gamma) || gamma <= 0) {
        stop('gamma must be one positive finite number')
    }
    sp <- exp(log_sp)
    solution <- .cdrgam_streamed_sparse_pirls(
        assembly,
        family,
        sp,
        fixed_ridge=fixed_ridge,
        initial_coefficients=if (!is.null(initial_solution)) {
            initial_solution$coefficients
        } else NULL,
        initial_factor=if (!is.null(initial_solution)) {
            initial_solution$factor
        } else NULL,
        initial_linear_predictors=if (!is.null(initial_solution)) {
            initial_solution$linear_predictors
        } else NULL,
        tolerance=tolerance,
        chunk_size=chunk_size
    )
    if (!isTRUE(solution$converged)) {
        stop('Sparse PIRLS did not converge')
    }
    penalty_trace <- vapply(seq_along(sp), function(i) {
        .sparse_logdet_score(
            solution$factor,
            sp[[i]] * assembly$penalty_components[[i]],
            chunk_size=256L
        )
    }, numeric(1))
    effective_df <- assembly$dimension - sum(penalty_trace)
    observation_count <- assembly$observation_count
    denominator <- observation_count - gamma * effective_df
    if (!is.finite(effective_df) || denominator <= 0) {
        stop('The generalized prediction-error denominator is not positive')
    }
    criterion_name <- if (estimated_gamma) 'GCV' else 'UBRE'
    criterion <- if (estimated_gamma) {
        observation_count * solution$deviance / denominator^2
    } else {
        solution$deviance / observation_count +
            2 * gamma * effective_df / observation_count - 1
    }
    if (!is.finite(criterion)) {
        stop('The generalized prediction-error criterion is not finite')
    }
    scale <- if (estimated_gamma) {
        .cdrgam_gamma_reported_scale(
            assembly$setup$y,
            solution$fitted_values,
            solution$prior_weights,
            effective_df
        )
    } else 1
    list(
        criterion=criterion,
        criterion_name=criterion_name,
        sp=sp,
        scale=scale,
        solution=solution,
        penalty_trace=penalty_trace,
        effective_df=effective_df
    )
}

.cdrgam_generalized_gcv_gradient_selection <- function(
        requested, assembly, family=NULL, solution=NULL,
        memory=.cdrgam_memory_availability()
) {
    supports <- .sparse_penalty_supports(assembly$penalty_components)
    union_count <- length(unique(unlist(supports, use.names=FALSE)))
    retained_elements <- as.double(union_count)^2
    inverse_chunk_size <- .cdrgam_sparse_trace_chunk_size(
        assembly$dimension, 1L, memory
    )
    retained_bytes <- 8 * retained_elements
    transient_bytes <- 16 * as.double(assembly$dimension) *
        inverse_chunk_size +
        24 * as.double(assembly$observation_count)
    estimated_peak_bytes <- retained_bytes + transient_bytes
    available <- memory$available_bytes
    budget <- if (is.finite(available)) {
        max(0, min(available / 2, available - 1024^3))
    } else NA_real_
    memory_safe <- is.finite(budget) && estimated_peak_bytes <= budget
    cost_resolved <- !is.null(solution) && !is.null(family)
    exact_rhs <- finite_factorizations <- NA_real_
    exact_design_passes <- finite_design_passes <- NA_real_
    cost_safe <- NA
    if (cost_resolved) {
        penalty_count <- length(assembly$penalty_components)
        iterations <- max(1L, solution$iterations)
        gamma_family <- identical(family$family, 'Gamma')
        exact_rhs <- (if (gamma_family) 1 else 2) * union_count +
            penalty_count
        finite_factorizations <- 2 * penalty_count * iterations
        exact_design_passes <- if (gamma_family) {
            1 + penalty_count
        } else ceiling(union_count / inverse_chunk_size) + penalty_count
        finite_design_passes <- 2 * penalty_count * iterations
        cost_safe <- exact_rhs <= 32 * finite_factorizations &&
            exact_design_passes <= finite_design_passes
    }
    if (!identical(requested, 'auto')) {
        method <- requested
        reason <- 'explicit generalized GCV gradient method'
    } else if (!is.finite(budget)) {
        method <- 'finite'
        reason <- 'available memory could not be determined'
    } else if (!memory_safe) {
        method <- 'finite'
        reason <- 'exact generalized GCV gradient exceeds the memory budget'
    } else if (cost_resolved && !cost_safe) {
        method <- 'finite'
        reason <- 'finite differences have lower predicted structural cost'
    } else {
        method <- 'exact'
        reason <- if (cost_resolved) {
            'exact generalized GCV gradient has lower predicted cost'
        } else 'exact generalized GCV gradient fits the memory budget'
    }
    list(
        policy_version=1L,
        requested=requested,
        method=method,
        reason=reason,
        retained_elements=retained_elements,
        retained_bytes=retained_bytes,
        transient_bytes=transient_bytes,
        estimated_peak_bytes=estimated_peak_bytes,
        inverse_chunk_size=inverse_chunk_size,
        memory_source=memory$source,
        memory_available_bytes=available,
        memory_budget_bytes=budget,
        cost_resolved=cost_resolved,
        exact_rhs=exact_rhs,
        finite_factorizations=finite_factorizations,
        exact_design_passes=exact_design_passes,
        finite_design_passes=finite_design_passes
    )
}

.cdrgam_generalized_observed_factor <- function(
        assembly, family, solution, fixed_ridge, chunk_size
) {
    if (!identical(family$family, 'Gamma')) return(solution$factor)
    observed_information <- NULL
    y <- assembly$setup$y
    weights <- solution$prior_weights * y / solution$fitted_values
    for (start in seq.int(
            1L, assembly$observation_count,
            by=min(chunk_size, assembly$observation_count))) {
        rows <- start:min(
            assembly$observation_count, start + chunk_size - 1L
        )
        X <- .cdrgam_sparse_design_chunk(assembly, rows)
        weighted <- Matrix::Diagonal(x=sqrt(weights[rows])) %*% X
        contribution <- Matrix::crossprod(weighted)
        observed_information <- if (is.null(observed_information)) {
            contribution
        } else observed_information + contribution
    }
    system <- Matrix::forceSymmetric(
        observed_information + solution$penalty +
            Matrix::Diagonal(assembly$dimension, x=fixed_ridge),
        uplo='U'
    )
    tryCatch(
        Matrix::update(solution$factor, system),
        error=function(error) Matrix::Cholesky(
            system, perm=TRUE, LDL=FALSE, super=assembly$supernodal
        )
    )
}

.cdrgam_generalized_exact_prediction_error_gradient <- function(
        assembly, family, evaluation, fixed_ridge=0, gamma=1,
        chunk_size=assembly$crossprod_chunk_size, batch_size=NULL,
        inverse_chunk_size=256L,
        element_limit=getOption('cdrgam.max_analytic_hessian_elements', 2e8)
) {
    solution <- evaluation$solution
    sp <- evaluation$sp
    components <- assembly$penalty_components
    count <- length(components)
    if (!count) return(numeric())
    supports <- .sparse_penalty_supports(components)
    union_support <- sort(unique(unlist(supports, use.names=FALSE)))
    support_positions <- lapply(supports, match, table=union_support)
    retained_elements <- as.double(length(union_support))^2
    if (!is.numeric(element_limit) || length(element_limit) != 1L ||
            !is.finite(element_limit) || element_limit < 1) {
        stop('The exact generalized GCV-gradient element limit must be positive')
    }
    if (retained_elements > element_limit) {
        stop(
            'The exact generalized GCV gradient would retain ',
            format(retained_elements, scientific=FALSE, big.mark=','),
            ' selected-inverse elements; the current limit is ',
            format(element_limit, scientific=FALSE, big.mark=','),
            '. Use gradient="finite" or raise option ',
            'cdrgam.max_analytic_hessian_elements explicitly.'
        )
    }
    inverse_chunk_size <- .cdrgam_positive_integer(
        inverse_chunk_size, 'exact generalized GCV inverse chunk size'
    )
    derivatives <- lapply(seq_len(count), function(index) {
        sp[[index]] * components[[index]]
    })
    penalty <- Reduce(`+`, derivatives)
    selected_inverse <- matrix(
        0, nrow=length(union_support), ncol=length(union_support)
    )
    weight_changes <- !identical(family$family, 'Gamma')
    penalty_leverage <- if (weight_changes) {
        numeric(assembly$observation_count)
    } else NULL
    inverse_rhs <- 0L
    starts <- seq.int(
        1L, assembly$observation_count,
        by=min(chunk_size, assembly$observation_count)
    )
    for (start in seq.int(
            1L, length(union_support), by=inverse_chunk_size)) {
        selected <- start:min(
            length(union_support), start + inverse_chunk_size - 1L
        )
        columns <- union_support[selected]
        selector <- Matrix::sparseMatrix(
            i=columns, j=seq_along(columns), x=1,
            dims=c(assembly$dimension, length(columns))
        )
        right_hand_sides <- if (weight_changes) {
            cbind(selector, penalty[, columns, drop=FALSE])
        } else selector
        solved <- as.matrix(.cdr_factor_solve(
            solution$factor, right_hand_sides
        ))
        inverse_rhs <- inverse_rhs + ncol(solved)
        inverse_columns <- solved[, seq_along(columns), drop=FALSE]
        selected_inverse[, selected] <- inverse_columns[
            union_support, , drop=FALSE
        ]
        if (weight_changes) {
            penalty_columns <- solved[
                , length(columns) + seq_along(columns), drop=FALSE
            ]
            for (row_start in starts) {
                rows <- row_start:min(
                    assembly$observation_count,
                    row_start + chunk_size - 1L
                )
                X <- .cdrgam_sparse_design_chunk(assembly, rows)
                left <- as.matrix(X %*% inverse_columns)
                right <- as.matrix(X %*% penalty_columns)
                penalty_leverage[rows] <- penalty_leverage[rows] +
                    rowSums(left * right)
            }
        }
    }
    inverse_derivatives <- vector('list', count)
    system_scores <- numeric(count)
    for (index in seq_len(count)) {
        support <- supports[[index]]
        positions <- support_positions[[index]]
        local_penalty <- as.matrix(derivatives[[index]][
            support, support, drop=FALSE
        ])
        value <- selected_inverse[, positions, drop=FALSE] %*%
            local_penalty
        inverse_derivatives[[index]] <- value
        system_scores[[index]] <- sum(value[cbind(
            positions, seq_along(support)
        )])
    }
    trace_cross_sums <- numeric(count)
    for (i in seq_len(count)) {
        positions_i <- support_positions[[i]]
        for (j in seq_len(i)) {
            support_j <- supports[[j]]
            positions_j <- support_positions[[j]]
            left <- inverse_derivatives[[j]][
                positions_i, , drop=FALSE
            ]
            right <- inverse_derivatives[[i]][
                positions_j, , drop=FALSE
            ]
            trace_cross <- sum(left * t(right))
            trace_cross_sums[[i]] <- trace_cross_sums[[i]] + trace_cross
            if (i != j) {
                trace_cross_sums[[j]] <-
                    trace_cross_sums[[j]] + trace_cross
            }
        }
    }
    edf_gradient <- -(system_scores - trace_cross_sums)
    observed_factor <- .cdrgam_generalized_observed_factor(
        assembly, family, solution, fixed_ridge, chunk_size
    )
    derivative_bytes <- max(
        1, 8 * as.double(assembly$dimension) +
            16 * as.double(assembly$observation_count)
    )
    plan <- .cdrgam_sparse_score_batch_plan(
        count, 1L, derivative_bytes, batch_size
    )
    deviance_gradient <- numeric(count)
    weight_gradient <- numeric(count)
    mu <- solution$fitted_values
    y <- assembly$setup$y
    prior_weights <- solution$prior_weights
    score_eta <- if (identical(family$family, 'Gamma')) {
        prior_weights * (1 - y / mu)
    } else prior_weights * (mu - y)
    working_weight_derivative <- if (identical(
            family$family, 'poisson'
        )) {
        prior_weights * mu
    } else if (identical(family$family, 'binomial')) {
        prior_weights * mu * (1 - mu) * (1 - 2 * mu)
    } else numeric(length(mu))
    for (indices in plan$groups) {
        right_hand_sides <- do.call(cbind, lapply(indices, function(index) {
            as.numeric(derivatives[[index]] %*% solution$coefficients)
        }))
        coefficient_derivatives <- -as.matrix(.cdr_factor_solve(
            observed_factor, right_hand_sides
        ))
        for (row_start in starts) {
            rows <- row_start:min(
                assembly$observation_count,
                row_start + chunk_size - 1L
            )
            X <- .cdrgam_sparse_design_chunk(assembly, rows)
            eta_derivatives <- as.matrix(X %*% coefficient_derivatives)
            deviance_gradient[indices] <-
                deviance_gradient[indices] +
                2 * colSums(eta_derivatives * score_eta[rows])
            if (weight_changes) {
                weight_gradient[indices] <-
                    weight_gradient[indices] + colSums(
                        eta_derivatives * working_weight_derivative[rows] *
                            penalty_leverage[rows]
                    )
            }
        }
    }
    edf_gradient <- edf_gradient + weight_gradient
    observation_count <- assembly$observation_count
    if (identical(family$family, 'Gamma')) {
        denominator <- observation_count - gamma * evaluation$effective_df
        gradient <- observation_count * (
            deviance_gradient / denominator^2 +
                2 * gamma * solution$deviance * edf_gradient /
                    denominator^3
        )
    } else {
        gradient <- deviance_gradient / observation_count +
            2 * gamma * edf_gradient / observation_count
    }
    attr(gradient, 'score_plan') <- list(
        workers=1L,
        parallel_axis='selected inverse and smoothing directions',
        batches=length(plan$groups),
        batch_size=plan$batch_size,
        inverse_chunk_size=inverse_chunk_size,
        inverse_rhs=inverse_rhs,
        retained_elements=retained_elements,
        memory_source=plan$memory$source,
        memory_available_bytes=plan$memory$available_bytes
    )
    gradient
}

.cdrgam_optimize_streamed_sparse_prediction_error <- function(
        assembly,
        family,
        initial_log_sp=NULL,
        fixed_ridge=0,
        gamma=1,
        tolerance=1e-9,
        chunk_size=assembly$crossprod_chunk_size,
        max_iterations=100L,
        gradient_tolerance=1e-4,
        finite_difference_step=1e-3,
        cores=1L,
        gradient_workers_requested=NULL,
        gradient_requested='auto',
        score_batch_size=NULL,
        reporter=NULL
) {
    family <- .as_family(family)
    parameter_count <- length(assembly$penalty_components)
    if (is.null(initial_log_sp)) {
        initial_log_sp <- log(.cdrgam_sparse_initial_sp(assembly))
    }
    if (length(initial_log_sp) != parameter_count ||
            any(!is.finite(initial_log_sp))) {
        stop('initial_log_sp has the wrong dimension or values')
    }
    if (!is.numeric(finite_difference_step) ||
            length(finite_difference_step) != 1L ||
            !is.finite(finite_difference_step) ||
            finite_difference_step <= 0) {
        stop('finite_difference_step must be positive')
    }
    lower <- rep.int(-25, parameter_count)
    upper <- rep.int(35, parameter_count)
    cached_parameters <- NULL
    cached_evaluation <- NULL
    warm_solution <- NULL
    evaluations <- 0L
    gradient_evaluations <- 0L
    warm_starts <- 0L
    cold_fallbacks <- 0L
    gradient_selection <- .cdrgam_generalized_gcv_gradient_selection(
        gradient_requested, assembly
    )
    gradient_method <- gradient_selection$method
    gradient_plan <- if (identical(gradient_method, 'finite')) {
        plan <- .cdrgam_parallel_plan(
            cores,
            max(1L, 2L * parameter_count),
            gradient_workers_requested
        )
        plan$requested_workers <- plan$workers
        plan$memory_limited <- FALSE
        plan$per_worker_bytes <- NA_real_
        plan$memory_source <- 'not-yet-measured'
        plan$memory_available_bytes <- NA_real_
        plan$memory_budget_bytes <- NA_real_
        plan
    } else list(
        cores=cores,
        workers=1L,
        blas_threads=cores,
        forked=FALSE,
        worker_source='exact selected-inverse gradient',
        requested_workers=1L,
        memory_limited=FALSE,
        per_worker_bytes=gradient_selection$estimated_peak_bytes,
        memory_source=gradient_selection$memory_source,
        memory_available_bytes=gradient_selection$memory_available_bytes,
        memory_budget_bytes=gradient_selection$memory_budget_bytes
    )
    resolve_gradient_plan <- function(solution) {
        sample_count <- min(
            assembly$observation_count,
            assembly$crossprod_chunk_size
        )
        sample_design <- .cdrgam_sparse_design_chunk(
            assembly, seq_len(sample_count)
        )
        transient_bytes <- 4 * as.double(utils::object.size(sample_design))
        rm(sample_design)
        gradient_plan <<- .cdrgam_memory_parallel_plan(
            cores,
            max(1L, 2L * parameter_count),
            gradient_workers_requested,
            .cdrgam_sparse_worker_bytes(
                solution$system,
                solution$factor,
                transient_bytes=transient_bytes
            )
        )
    }
    solve <- function(parameters, initial) {
        attempt <- tryCatch(
            .cdrgam_streamed_sparse_prediction_error(
                assembly,
                family,
                parameters,
                initial_solution=initial,
                fixed_ridge=fixed_ridge,
                gamma=gamma,
                tolerance=tolerance,
                chunk_size=chunk_size
            ),
            error=function(error) structure(
                list(message=conditionMessage(error)),
                class='cdrgam_prediction_error_failure'
            )
        )
        if (inherits(attempt, 'cdrgam_prediction_error_failure') &&
                !is.null(initial)) {
            attempt <- tryCatch(
                .cdrgam_streamed_sparse_prediction_error(
                    assembly,
                    family,
                    parameters,
                    initial_solution=NULL,
                    fixed_ridge=fixed_ridge,
                    gamma=gamma,
                    tolerance=tolerance,
                    chunk_size=chunk_size
                ),
                error=function(error) structure(
                    list(message=conditionMessage(error)),
                    class='cdrgam_prediction_error_failure'
                )
            )
        }
        attempt
    }
    evaluate <- function(parameters) {
        if (!is.null(cached_parameters) && identical(
                as.numeric(parameters), cached_parameters)) {
            return(cached_evaluation)
        }
        evaluations <<- evaluations + 1L
        used_warm_start <- !is.null(warm_solution)
        if (used_warm_start) warm_starts <<- warm_starts + 1L
        started <- proc.time()[['elapsed']]
        value <- solve(parameters, warm_solution)
        if (inherits(value, 'cdrgam_prediction_error_failure')) {
            if (used_warm_start) cold_fallbacks <<- cold_fallbacks + 1L
            stop(value$message)
        }
        warm_solution <<- value$solution
        cached_parameters <<- as.numeric(parameters)
        cached_evaluation <<- value
        if (!is.null(reporter)) reporter$emit(
            1L,
            'outer inner solve complete',
            evaluation=evaluations,
            seconds=format(proc.time()[['elapsed']] - started, digits=5),
            iterations=value$solution$iterations,
            warm_start=used_warm_start,
            criterion=format(value$criterion, digits=10)
        )
        value
    }
    objective <- function(parameters) evaluate(parameters)$criterion
    finite_gradient <- function(parameters) {
        started <- proc.time()[['elapsed']]
        base <- evaluate(parameters)
        resolve_gradient_plan(base$solution)
        points <- vector('list', 2L * parameter_count)
        widths <- numeric(parameter_count)
        for (i in seq_len(parameter_count)) {
            lower_point <- upper_point <- parameters
            lower_point[[i]] <- max(
                lower[[i]], parameters[[i]] - finite_difference_step
            )
            upper_point[[i]] <- min(
                upper[[i]], parameters[[i]] + finite_difference_step
            )
            points[[2L * i - 1L]] <- lower_point
            points[[2L * i]] <- upper_point
            widths[[i]] <- upper_point[[i]] - lower_point[[i]]
        }
        worker <- function(point) {
            value <- solve(point, base$solution)
            if (inherits(value, 'cdrgam_prediction_error_failure')) {
                stop(value$message)
            }
            value$criterion
        }
        values <- .cdrgam_with_blas_threads(gradient_plan$blas_threads, {
            results <- if (gradient_plan$workers > 1L &&
                    .Platform$OS.type != 'windows') {
                parallel::mclapply(
                    points,
                    worker,
                    mc.cores=gradient_plan$workers,
                    mc.preschedule=TRUE,
                    mc.set.seed=FALSE
                )
            } else lapply(points, worker)
            valid <- vapply(results, function(value) {
                is.numeric(value) && length(value) == 1L &&
                    is.finite(value)
            }, logical(1))
            if (!all(valid)) {
                stop('Parallel finite-difference gradient evaluation failed')
            }
            vapply(results, as.numeric, numeric(1))
        })
        evaluations <<- evaluations + length(points)
        gradient_evaluations <<- gradient_evaluations + 1L
        score <- (values[seq.int(2L, length(values), by=2L)] -
            values[seq.int(1L, length(values), by=2L)]) / widths
        if (!is.null(reporter)) reporter$emit(
            1L,
            'outer finite gradient complete',
            gradient_evaluation=gradient_evaluations,
            evaluations=length(points),
            workers=gradient_plan$workers,
            blas_threads=gradient_plan$blas_threads,
            memory_limited=gradient_plan$memory_limited,
            seconds=format(proc.time()[['elapsed']] - started, digits=5),
            maximum=format(max(abs(score)), digits=5)
        )
        score
    }
    exact_gradient <- function(parameters) {
        started <- proc.time()[['elapsed']]
        base <- evaluate(parameters)
        score <- .cdrgam_with_blas_threads(
            cores,
            .cdrgam_generalized_exact_prediction_error_gradient(
                assembly,
                family,
                base,
                fixed_ridge=fixed_ridge,
                gamma=gamma,
                chunk_size=chunk_size,
                batch_size=score_batch_size,
                inverse_chunk_size=gradient_selection$inverse_chunk_size,
                element_limit=if (identical(
                    gradient_requested, 'auto'
                )) {
                    max(1, gradient_selection$retained_elements)
                } else getOption(
                    'cdrgam.max_analytic_hessian_elements', 2e8
                )
            )
        )
        gradient_evaluations <<- gradient_evaluations + 1L
        score_plan <- attr(score, 'score_plan')
        gradient_plan <<- c(
            gradient_plan[setdiff(names(gradient_plan), names(score_plan))],
            score_plan
        )
        if (!is.null(reporter)) reporter$emit(
            1L,
            'outer exact gradient complete',
            gradient_evaluation=gradient_evaluations,
            inverse_rhs=score_plan$inverse_rhs,
            batch_size=score_plan$batch_size,
            blas_threads=cores,
            seconds=format(proc.time()[['elapsed']] - started, digits=5),
            maximum=format(max(abs(score)), digits=5)
        )
        as.numeric(score)
    }
    gradient <- function(parameters) {
        if (identical(gradient_method, 'exact')) {
            exact_gradient(parameters)
        } else finite_gradient(parameters)
    }
    if (!parameter_count) {
        retained <- evaluate(numeric())
        optimization <- list(
            par=numeric(),
            value=retained$criterion,
            gradient=numeric(),
            convergence=0L,
            counts=c(`function`=1L, gradient=0L),
            message='no smoothing parameters'
        )
    } else {
        # Prediction-error criteria can put a null-space solution beyond a
        # shallow local basin. Common shifts cheaply identify its attraction
        # region before each finite-difference optimization becomes costly.
        probe_shifts <- c(0, -4, 4, -8, 8)
        probe_points <- lapply(probe_shifts, function(shift) {
            pmin(upper, pmax(lower, initial_log_sp + shift))
        })
        probe_objective <- function(parameters) tryCatch(
            objective(parameters),
            error=function(error) Inf
        )
        probe_values <- vapply(probe_points, probe_objective, numeric(1))
        if (!any(is.finite(probe_values))) {
            stop('Could not evaluate the generalized prediction-error fit')
        }
        optimization_start <- probe_points[[which.min(probe_values)]]
        if (identical(gradient_requested, 'auto')) {
            gradient_selection <-
                .cdrgam_generalized_gcv_gradient_selection(
                    gradient_requested,
                    assembly,
                    family=family,
                    solution=evaluate(optimization_start)$solution
                )
            gradient_method <- gradient_selection$method
        }
        if (!is.null(reporter)) reporter$emit(
            1L,
            'outer initialization probe complete',
            evaluations=length(probe_points),
            selected_shift=probe_shifts[[which.min(probe_values)]],
            criterion=format(min(probe_values), digits=10),
            gradient=gradient_method,
            gradient_reason=gradient_selection$reason
        )
        optimization <- stats::optim(
            par=optimization_start,
            fn=objective,
            gr=gradient,
            method='L-BFGS-B',
            lower=lower,
            upper=upper,
            control=list(
                maxit=max_iterations,
                pgtol=gradient_tolerance,
                factr=1e7
            )
        )
        # Audit each coordinate beyond the converged basin. One restart is
        # enough to cross a local barrier without turning this into an
        # open-ended multistart search.
        escape_points <- list(optimization$par)
        for (i in seq_len(parameter_count)) {
            for (shift in c(-8, 8)) {
                point <- optimization$par
                point[[i]] <- min(upper[[i]], max(
                    lower[[i]], point[[i]] + shift
                ))
                escape_points[[length(escape_points) + 1L]] <- point
            }
        }
        escape_values <- vapply(escape_points, probe_objective, numeric(1))
        escape <- which.min(escape_values)
        improvement_tolerance <- sqrt(.Machine$double.eps) *
            (1 + abs(optimization$value))
        basin_restart <- escape_values[[escape]] <
            optimization$value - improvement_tolerance
        if (basin_restart) {
            if (!is.null(reporter)) reporter$emit(
                1L,
                'outer basin escape detected',
                evaluations=length(escape_points),
                previous_criterion=format(optimization$value, digits=10),
                candidate_criterion=format(
                    escape_values[[escape]], digits=10
                )
            )
            previous_counts <- optimization$counts
            optimization <- stats::optim(
                par=escape_points[[escape]],
                fn=objective,
                gr=gradient,
                method='L-BFGS-B',
                lower=lower,
                upper=upper,
                control=list(
                    maxit=max_iterations,
                    pgtol=gradient_tolerance,
                    factr=1e7
                )
            )
            optimization$counts <- optimization$counts + previous_counts
        }
        optimization$basin_restart <- basin_restart
        retained <- evaluate(optimization$par)
        optimization$gradient <- gradient(optimization$par)
        projected_gradient <- optimization$gradient
        at_lower <- optimization$par - lower <= sqrt(.Machine$double.eps)
        at_upper <- upper - optimization$par <= sqrt(.Machine$double.eps)
        projected_gradient[at_lower & projected_gradient > 0] <- 0
        projected_gradient[at_upper & projected_gradient < 0] <- 0
        optimization$projected_gradient <- projected_gradient
        if (max(abs(projected_gradient)) > gradient_tolerance) {
            optimization$convergence <- 1L
            optimization$message <- paste0(
                optimization$message,
                if (nzchar(optimization$message)) '; ' else '',
                'outer projected gradient exceeds tolerance'
            )
        }
    }
    list(
        optimization=optimization,
        retained=retained,
        evaluations=evaluations,
        gradient_evaluations=gradient_evaluations,
        warm_starts=warm_starts,
        cold_fallbacks=cold_fallbacks,
        gradient_plan=c(
            gradient_plan,
            if (is.null(gradient_plan$parallel_axis)) {
                list(parallel_axis='finite-difference smoothing parameters')
            } else list()
        ),
        gradient_method=gradient_method,
        gradient_selection=gradient_selection
    )
}

.cdrgam_sparse_design_chunk <- function(assembly, rows) {
    if (!is.null(assembly$cached_design)) {
        if (length(rows) == assembly$observation_count &&
                identical(rows, seq_len(assembly$observation_count))) {
            return(assembly$cached_design)
        }
        return(assembly$cached_design[rows, , drop=FALSE])
    }
    pieces <- lapply(
        assembly$matrices,
        .sparse_component_rows,
        rows=rows
    )
    if (length(pieces) == 1L) pieces[[1L]] else do.call(cbind, pieces)
}

.cdrgam_sparse_linear_predictor <- function(
        assembly, coefficients, offset, chunk_size
) {
    output <- numeric(assembly$observation_count)
    for (start in seq.int(
            1L, assembly$observation_count,
            by=min(chunk_size, assembly$observation_count))) {
        rows <- start:min(
            assembly$observation_count,
            start + chunk_size - 1L
        )
        X <- .cdrgam_sparse_design_chunk(assembly, rows)
        output[rows] <- as.numeric(offset[rows] + X %*% coefficients)
    }
    output
}

.cdrgam_initial_sp_from_diagonal <- function(
        information_diagonal,
        penalty_components
) {
    count <- length(penalty_components)
    if (!count) return(numeric())
    penalty_diagonals <- lapply(
        penalty_components,
        function(penalty) abs(Matrix::diag(penalty))
    )
    smoothing_parameters <- vapply(penalty_diagonals, function(diagonal) {
        threshold <- max(diagonal) * .Machine$double.eps^0.8
        active <- diagonal > threshold
        if (!any(active)) return(1)
        mean(information_diagonal[active]) / mean(diagonal[active])
    }, numeric(1))
    combined <- numeric(length(information_diagonal))
    for (i in seq_len(count)) {
        combined <- combined + smoothing_parameters[[i]] *
            penalty_diagonals[[i]]
    }
    active <- information_diagonal > 0 & combined > 0
    if (any(active)) {
        ratio <- function() mean(
            information_diagonal[active] /
                (information_diagonal[active] + combined[active])
        )
        while (ratio() > 0.4) {
            smoothing_parameters <- smoothing_parameters * 10
            combined <- combined * 10
        }
        while (ratio() < 0.4) {
            smoothing_parameters <- smoothing_parameters / 10
            combined <- combined / 10
        }
    }
    pmax(smoothing_parameters, .Machine$double.eps)
}

.cdrgam_sparse_initial_sp <- function(assembly) {
    information_diagonal <- numeric(assembly$dimension)
    chunk_size <- assembly$crossprod_chunk_size
    for (start in seq.int(
            1L, assembly$observation_count,
            by=min(chunk_size, assembly$observation_count))) {
        rows <- start:min(
            assembly$observation_count,
            start + chunk_size - 1L
        )
        X <- .cdrgam_sparse_design_chunk(assembly, rows)
        information_diagonal <- information_diagonal +
            Matrix::colSums(X^2)
    }
    .cdrgam_initial_sp_from_diagonal(
        information_diagonal,
        assembly$penalty_components
    )
}

# Streamed sparse PIRLS over the component representation used by the sparse
# Gaussian backend. No global response-by-coefficient matrix is materialized.
.cdrgam_streamed_sparse_pirls <- function(
        assembly,
        family,
        smoothing_parameters,
        fixed_ridge=0,
        initial_coefficients=NULL,
        initial_factor=NULL,
        initial_linear_predictors=NULL,
        tolerance=1e-8,
        max_iterations=100L,
        maximum_halvings=25L,
        chunk_size=assembly$crossprod_chunk_size
) {
    family <- .as_family(family)
    y <- assembly$setup$y
    offset <- assembly$setup$offset
    if (is.null(offset)) offset <- numeric(length(y))
    prior_weights <- assembly$setup$w
    if (is.null(prior_weights)) prior_weights <- rep.int(1, length(y))
    if (any(!is.finite(prior_weights)) || any(prior_weights < 0)) {
        stop('Streamed sparse PIRLS requires finite nonnegative prior weights')
    }
    if (length(smoothing_parameters) !=
            length(assembly$penalty_components)) {
        stop('smoothing_parameters has the wrong dimension')
    }
    penalty <- Matrix::Matrix(
        0,
        nrow=assembly$dimension,
        ncol=assembly$dimension,
        sparse=TRUE
    )
    for (i in seq_along(smoothing_parameters)) {
        penalty <- penalty + smoothing_parameters[[i]] *
            assembly$penalty_components[[i]]
    }
    warm_started <- !is.null(initial_coefficients)
    coefficients <- if (!warm_started) {
        numeric(assembly$dimension)
    } else {
        if (length(initial_coefficients) != assembly$dimension ||
                any(!is.finite(initial_coefficients))) {
            stop('initial_coefficients has the wrong dimension or values')
        }
        as.numeric(initial_coefficients)
    }
    retained_predictor <- warm_started &&
        length(initial_linear_predictors) == assembly$observation_count &&
        all(is.finite(initial_linear_predictors))
    eta <- if (retained_predictor) {
        as.numeric(initial_linear_predictors)
    } else if (warm_started) {
        .cdrgam_sparse_linear_predictor(
            assembly, coefficients, offset, chunk_size
        )
    } else family$linkfun(.cdrgam_initial_mu(family, y, prior_weights))
    initial_mu <- family$linkinv(eta)
    objective <- if (warm_started &&
            .cdrgam_valid_family_state(family, eta, initial_mu)) {
        sum(family$dev.resids(y, initial_mu, prior_weights)) +
            as.numeric(Matrix::crossprod(
                coefficients, penalty %*% coefficients
            )) + fixed_ridge * sum(coefficients^2)
    } else Inf
    converged <- FALSE
    factor <- initial_factor
    numeric_updates <- 0L
    starts <- seq.int(
        1L,
        assembly$observation_count,
        by=min(chunk_size, assembly$observation_count)
    )
    for (iteration in seq_len(max_iterations)) {
        mu <- family$linkinv(eta)
        derivative <- family$mu.eta(eta)
        variance <- family$variance(mu)
        valid <- .cdrgam_valid_family_state(family, eta, mu) &&
            all(is.finite(derivative)) && all(is.finite(variance)) &&
            all(variance > 0) && all(derivative != 0)
        if (!valid) {
            stop('Streamed sparse PIRLS encountered an invalid family state')
        }
        working_weights <- prior_weights * derivative^2 / variance
        working_response <- eta - offset + (y - mu) / derivative
        information <- NULL
        rhs <- numeric(assembly$dimension)
        for (start in starts) {
            rows <- start:min(
                assembly$observation_count,
                start + chunk_size - 1L
            )
            X <- .cdrgam_sparse_design_chunk(assembly, rows)
            root_weights <- sqrt(working_weights[rows])
            weighted_X <- Matrix::Diagonal(x=root_weights) %*% X
            contribution <- Matrix::crossprod(weighted_X)
            information <- if (is.null(information)) {
                contribution
            } else information + contribution
            rhs <- rhs + as.numeric(Matrix::crossprod(
                weighted_X,
                working_response[rows] * root_weights
            ))
        }
        system <- Matrix::forceSymmetric(
            information + penalty +
                Matrix::Diagonal(assembly$dimension, x=fixed_ridge),
            uplo='U'
        )
        factor <- if (is.null(factor)) {
            Matrix::Cholesky(
                system,
                perm=TRUE,
                LDL=FALSE,
                super=assembly$supernodal
            )
        } else {
            numeric_updates <- numeric_updates + 1L
            tryCatch(
                Matrix::update(factor, system),
                error=function(error) Matrix::Cholesky(
                    system,
                    perm=TRUE,
                    LDL=FALSE,
                    super=assembly$supernodal
                )
            )
        }
        proposal <- as.numeric(.cdr_factor_solve(factor, rhs))
        accepted <- FALSE
        halvings <- 0L
        repeat {
            candidate_eta <- .cdrgam_sparse_linear_predictor(
                assembly,
                proposal,
                offset,
                chunk_size
            )
            candidate_mu <- family$linkinv(candidate_eta)
            candidate_objective <- if (.cdrgam_valid_family_state(
                    family, candidate_eta, candidate_mu)) {
                sum(family$dev.resids(y, candidate_mu, prior_weights)) +
                    as.numeric(Matrix::crossprod(
                        proposal,
                        penalty %*% proposal
                    )) + fixed_ridge * sum(proposal^2)
            } else Inf
            if (is.finite(candidate_objective) &&
                    (!is.finite(objective) ||
                        candidate_objective <= objective +
                            tolerance * (1 + abs(objective)))) {
                accepted <- TRUE
                break
            }
            if (halvings >= maximum_halvings) break
            proposal <- (proposal + coefficients) / 2
            halvings <- halvings + 1L
        }
        if (!accepted) {
            stop(
                'Streamed sparse PIRLS step halving failed to improve the ',
                'penalized deviance'
            )
        }
        change <- if (is.finite(objective)) {
            abs(objective - candidate_objective) / (1 + abs(objective))
        } else Inf
        coefficients <- proposal
        eta <- candidate_eta
        objective <- candidate_objective
        if (is.finite(change) && change < tolerance) {
            converged <- TRUE
            break
        }
    }
    mu <- family$linkinv(eta)
    derivative <- family$mu.eta(eta)
    variance <- family$variance(mu)
    working_weights <- prior_weights * derivative^2 / variance
    information <- NULL
    for (start in starts) {
        rows <- start:min(
            assembly$observation_count,
            start + chunk_size - 1L
        )
        X <- .cdrgam_sparse_design_chunk(assembly, rows)
        weighted_X <- Matrix::Diagonal(
            x=sqrt(working_weights[rows])
        ) %*% X
        contribution <- Matrix::crossprod(weighted_X)
        information <- if (is.null(information)) {
            contribution
        } else information + contribution
    }
    system <- Matrix::forceSymmetric(
        information + penalty +
            Matrix::Diagonal(assembly$dimension, x=fixed_ridge),
        uplo='U'
    )
    factor <- tryCatch(
        Matrix::update(factor, system),
        error=function(error) Matrix::Cholesky(
            system,
            perm=TRUE,
            LDL=FALSE,
            super=assembly$supernodal
        )
    )
    numeric_updates <- numeric_updates + 1L
    list(
        coefficients=coefficients,
        linear_predictors=eta,
        fitted_values=mu,
        residuals=y - mu,
        deviance=sum(family$dev.resids(y, mu, prior_weights)),
        penalized_deviance=objective,
        penalty=penalty,
        system=system,
        factor=factor,
        working_weights=working_weights,
        prior_weights=prior_weights,
        iterations=iteration,
        numeric_updates=numeric_updates,
        chunks=length(starts),
        converged=converged,
        warm_started=warm_started
    )
}

# Dense Laplace-approximated marginal-likelihood reference. The scalar
# criterion is intentionally finite-differenced; this path validates the
# statistical target before an analytic sparse score is introduced.
.cdrgam_dense_initial_sp <- function(setup) {
    count <- length(setup$S)
    if (!count) return(numeric())
    information_diagonal <- colSums(setup$X^2)
    combined <- numeric(ncol(setup$X))
    smoothing_parameters <- numeric(count)
    for (i in seq_len(count)) {
        diagonal <- abs(diag(setup$S[[i]]))
        threshold <- max(diagonal) * .Machine$double.eps^0.8
        active <- diagonal > threshold
        indices <- setup$off[[i]] + seq_along(diagonal) - 1L
        smoothing_parameters[[i]] <- if (any(active)) {
            mean(information_diagonal[indices[active]]) /
                mean(diagonal[active])
        } else 1
        combined[indices] <- combined[indices] +
            smoothing_parameters[[i]] * diagonal
    }
    active <- information_diagonal > 0 & combined > 0
    if (any(active)) {
        ratio <- function() mean(
            information_diagonal[active] /
                (information_diagonal[active] + combined[active])
        )
        while (ratio() > 0.4) {
            smoothing_parameters <- smoothing_parameters * 10
            combined <- combined * 10
        }
        while (ratio() < 0.4) {
            smoothing_parameters <- smoothing_parameters / 10
            combined <- combined / 10
        }
    }
    pmax(smoothing_parameters, .Machine$double.eps)
}

.cdrgam_dense_laml <- function(
        setup,
        family,
        initial_log_sp=NULL,
        fixed_ridge=0,
        lower=-25,
        upper=25,
        pirls_tolerance=1e-9,
        optimizer_control=list(factr=1e7)
) {
    family <- .as_family(family)
    canonical <- (identical(family$family, 'binomial') &&
            identical(family$link, 'logit')) ||
        (identical(family$family, 'poisson') &&
            identical(family$link, 'log'))
    estimated_gamma <- identical(family$family, 'Gamma') &&
        identical(family$link, 'log')
    if (!canonical && !estimated_gamma) {
        stop(
            'The dense LAML reference currently requires binomial(logit) ',
            ', poisson(log), or Gamma(log)'
        )
    }
    penalty_count <- length(setup$S)
    parameter_count <- penalty_count + as.integer(estimated_gamma)
    if (is.null(initial_log_sp)) {
        initial_log_sp <- numeric(parameter_count)
        if (estimated_gamma && penalty_count) {
            initial_log_sp[seq_len(penalty_count)] <-
                log(.cdrgam_dense_initial_sp(setup))
        }
    }
    if (length(initial_log_sp) != parameter_count ||
            any(!is.finite(initial_log_sp))) {
        stop('initial_log_sp has the wrong dimension or values for this family')
    }
    evaluations <- 0L
    evaluate <- function(parameters, retain=FALSE) {
        evaluations <<- evaluations + 1L
        log_sp <- parameters[seq_len(penalty_count)]
        dispersion <- if (estimated_gamma) {
            exp(parameters[[parameter_count]])
        } else 1
        solution <- tryCatch(
            .cdrgam_dense_pirls(
                setup,
                family,
                exp(log_sp),
                fixed_ridge=fixed_ridge,
                tolerance=pirls_tolerance
            ),
            error=function(error) NULL
        )
        if (is.null(solution) || !isTRUE(solution$converged)) {
            return(if (retain) NULL else 1e100)
        }
        penalty_determinant <- .positive_log_determinant(solution$penalty)
        log_system_determinant <- 2 * sum(log(diag(solution$factor)))
        criterion <- solution$penalized_deviance / dispersion +
            log_system_determinant - penalty_determinant$value
        if (estimated_gamma) {
            saturated <- .cdrgam_gamma_saturated_loglik(
                setup$y,
                solution$prior_weights,
                dispersion
            )
            nullity <- ncol(setup$X) - penalty_determinant$rank
            criterion <- criterion - 2 * saturated$value -
                nullity * log(2 * pi * dispersion)
        }
        if (!is.finite(criterion)) {
            return(if (retain) NULL else 1e100)
        }
        if (!retain) return(criterion)
        solution$criterion <- criterion
        solution$sp <- exp(log_sp)
        solution$scale <- dispersion
        solution$penalty_rank <- penalty_determinant$rank
        solution
    }
    if (!parameter_count) {
        solution <- evaluate(numeric(), retain=TRUE)
        return(list(
            optimization=list(
                par=numeric(), value=solution$criterion,
                convergence=0L, counts=c(`function`=1L, gradient=NA_integer_)
            ),
            solution=solution,
            evaluations=evaluations
        ))
    }
    optimization <- stats::optim(
        initial_log_sp,
        evaluate,
        method='L-BFGS-B',
        lower=rep.int(lower, parameter_count),
        upper=rep.int(upper, parameter_count),
        control=optimizer_control
    )
    solution <- evaluate(optimization$par, retain=TRUE)
    if (is.null(solution)) stop('Dense LAML optimization ended at an invalid fit')
    list(
        optimization=optimization,
        solution=solution,
        evaluations=evaluations
    )
}
