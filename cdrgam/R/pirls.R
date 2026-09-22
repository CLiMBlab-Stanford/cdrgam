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
        fixed_ridge=0,
        tolerance=1e-9,
        chunk_size=assembly$crossprod_chunk_size,
        score=FALSE
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
        initial_coefficients=initial_coefficients,
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
        penalty_logdet=penalty_determinant$value,
        penalty_rank=penalty_determinant$rank
    )
    if (!isTRUE(score)) return(output)
    penalty_score <- .sparse_penalty_logdet_score(
        assembly$blocks,
        sp,
        penalty_count
    )
    beta <- solution$coefficients
    mu <- solution$fitted_values
    prior_weights <- solution$prior_weights
    weight_derivative <- if (identical(family$family, 'poisson')) {
        prior_weights * mu
    } else if (identical(family$family, 'binomial')) {
        prior_weights * mu * (1 - mu) * (1 - 2 * mu)
    } else numeric(length(mu))
    scores <- numeric(parameter_count)
    zero_offset <- numeric(assembly$observation_count)
    starts <- seq.int(
        1L,
        assembly$observation_count,
        by=min(chunk_size, assembly$observation_count)
    )
    for (j in seq_len(penalty_count)) {
        penalty_derivative <- sp[[j]] * assembly$penalty_components[[j]]
        coefficient_derivative <- -as.numeric(.cdr_factor_solve(
            solution$factor,
            penalty_derivative %*% beta
        ))
        eta_derivative <- .cdrgam_sparse_linear_predictor(
            assembly,
            coefficient_derivative,
            zero_offset,
            chunk_size
        )
        curvature_derivative <- NULL
        for (start in starts) {
            rows <- start:min(
                assembly$observation_count,
                start + chunk_size - 1L
            )
            X <- .cdrgam_sparse_design_chunk(assembly, rows)
            contribution <- Matrix::crossprod(
                X,
                Matrix::Diagonal(
                    x=weight_derivative[rows] * eta_derivative[rows]
                ) %*% X
            )
            curvature_derivative <- if (is.null(curvature_derivative)) {
                contribution
            } else curvature_derivative + contribution
        }
        system_derivative <- Matrix::forceSymmetric(
            penalty_derivative + curvature_derivative,
            uplo='U'
        )
        scores[[j]] <- as.numeric(Matrix::crossprod(
            beta,
            penalty_derivative %*% beta
        )) / dispersion + .sparse_logdet_score(
            solution$factor,
            system_derivative
        ) - penalty_score[[j]]
    }
    if (estimated_gamma) {
        nullity <- assembly$dimension - penalty_determinant$rank
        scores[[parameter_count]] <-
            -solution$penalized_deviance / dispersion -
            2 * saturated$derivative * dispersion - nullity
    }
    output$score <- scores
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
        trust_radius=2
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
    evaluate <- function(parameters) {
        if (!is.null(cached_parameters) && identical(
                as.numeric(parameters), cached_parameters)) {
            return(cached_evaluation)
        }
        value <- .cdrgam_streamed_sparse_laml(
            assembly,
            family,
            parameters,
            fixed_ridge=fixed_ridge,
            tolerance=tolerance,
            chunk_size=chunk_size,
            score=TRUE
        )
        cached_parameters <<- as.numeric(parameters)
        cached_evaluation <<- value
        value
    }
    objective <- function(parameters) evaluate(parameters)$criterion
    gradient <- function(parameters) evaluate(parameters)$score
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
    retained <- evaluate(optimization$par)
    list(optimization=optimization, retained=retained)
}

.cdrgam_sparse_design_chunk <- function(assembly, rows) {
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

.cdrgam_sparse_initial_sp <- function(assembly) {
    count <- length(assembly$penalty_components)
    if (!count) return(numeric())
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
    penalty_diagonals <- lapply(
        assembly$penalty_components,
        function(penalty) abs(Matrix::diag(penalty))
    )
    smoothing_parameters <- vapply(penalty_diagonals, function(diagonal) {
        threshold <- max(diagonal) * .Machine$double.eps^0.8
        active <- diagonal > threshold
        if (!any(active)) return(1)
        mean(information_diagonal[active]) / mean(diagonal[active])
    }, numeric(1))
    combined <- numeric(assembly$dimension)
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

# Streamed sparse PIRLS over the component representation used by the sparse
# Gaussian backend. No global response-by-coefficient matrix is materialized.
.cdrgam_streamed_sparse_pirls <- function(
        assembly,
        family,
        smoothing_parameters,
        fixed_ridge=0,
        initial_coefficients=NULL,
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
    coefficients <- if (is.null(initial_coefficients)) {
        numeric(assembly$dimension)
    } else {
        if (length(initial_coefficients) != assembly$dimension ||
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
        converged=converged
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
            return(if (retain) NULL else .Machine$double.xmax / 100)
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
            return(if (retain) NULL else .Machine$double.xmax / 100)
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
