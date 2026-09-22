.cdrgam_distributional_block_setup <- function(design, ...) {
    .fit_compressed_mgcv(
        y=design$responses[[design$response_name]],
        terms=design$terms,
        family=stats::gaussian(),
        method='REML',
        engine='gam',
        base_formula=design$ordinary_formula,
        response_data=design$responses,
        user_formula=design$formula,
        preparation=.cdrgam_distributional_preparation(design),
        setup_only=TRUE,
        ...
    )
}

.cdrgam_gaulss_b <- function(family) {
    value <- environment(family$ll)$b
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 0) {
        stop('Could not recover the gaulss minimum standard deviation')
    }
    value
}

.cdrgam_distributional_solve <- function(system, rhs) {
    scale <- max(abs(diag(system)), 1)
    stabilized <- system + diag(1e-10 * scale, nrow(system))
    factor <- chol(stabilized)
    backsolve(factor, forwardsolve(t(factor), rhs))
}

.cdrgam_gaulss_fixed <- function(
        setups,
        smoothing_parameters,
        b,
        tolerance=1e-5,
        maxit=200L
) {
    location <- setups$location
    scale <- setups$scale
    X_location <- location$X
    X_scale <- scale$X
    y <- location$y
    weights <- location$w
    if (is.null(weights)) weights <- rep.int(1, length(y))
    counts <- vapply(setups, function(setup) length(setup$S), integer(1))
    location_indices <- seq_len(counts[['location']])
    scale_indices <- counts[['location']] + seq_len(counts[['scale']])
    location_penalty <- .embed_penalties(
        location, smoothing_parameters[location_indices]
    )
    scale_penalty <- .embed_penalties(
        scale, smoothing_parameters[scale_indices]
    )
    location_system <- crossprod(X_location * sqrt(weights)) +
        location_penalty
    location_coefficients <- .cdrgam_distributional_solve(
        location_system,
        crossprod(X_location, weights * (y - location$offset))
    )
    weighted_mean <- stats::weighted.mean(y, weights)
    initial_sigma <- sqrt(stats::weighted.mean(
        (y - weighted_mean)^2, weights
    ))
    scale_target <- rep.int(
        log(max(initial_sigma - b, 1e-4)), length(y)
    ) - scale$offset
    scale_system <- crossprod(X_scale * sqrt(weights)) + scale_penalty
    scale_coefficients <- .cdrgam_distributional_solve(
        scale_system,
        crossprod(X_scale, weights * scale_target)
    )
    objective <- function(location_coefficients, scale_coefficients) {
        location_predictor <- drop(
            X_location %*% location_coefficients + location$offset
        )
        scale_predictor <- drop(
            X_scale %*% scale_coefficients + scale$offset
        )
        sigma <- b + exp(scale_predictor)
        if (any(!is.finite(sigma))) return(Inf)
        residual <- y - location_predictor
        sum(weights * (log(sigma) + 0.5 * (residual / sigma)^2)) +
            0.5 * drop(crossprod(
                location_coefficients,
                location_penalty %*% location_coefficients
            )) +
            0.5 * drop(crossprod(
                scale_coefficients,
                scale_penalty %*% scale_coefficients
            ))
    }
    value <- objective(location_coefficients, scale_coefficients)
    converged <- FALSE
    gradient_norm <- Inf
    for (iteration in seq_len(maxit)) {
        location_predictor <- drop(
            X_location %*% location_coefficients + location$offset
        )
        scale_predictor <- drop(
            X_scale %*% scale_coefficients + scale$offset
        )
        q <- exp(scale_predictor)
        sigma <- b + q
        residual <- y - location_predictor
        if (any(!is.finite(sigma)) || any(sigma <= 0)) break
        location_gradient <- crossprod(
            X_location, weights * (-residual / sigma^2)
        ) + location_penalty %*% location_coefficients
        scale_gradient <- crossprod(
            X_scale,
            weights * q * (1 / sigma - residual^2 / sigma^3)
        ) + scale_penalty %*% scale_coefficients
        gradient_norm <- max(abs(c(location_gradient, scale_gradient)))
        if (!is.finite(gradient_norm)) break
        location_information <- crossprod(
            X_location * sqrt(weights / sigma^2)
        ) + location_penalty
        scale_information <- crossprod(
            X_scale * sqrt(weights * 2 * q^2 / sigma^2)
        ) + scale_penalty
        location_step <- -.cdrgam_distributional_solve(
            location_information, location_gradient
        )
        scale_step <- -.cdrgam_distributional_solve(
            scale_information, scale_gradient
        )
        step_size <- 1
        accepted <- FALSE
        while (step_size >= 2^-20) {
            candidate_location <- location_coefficients +
                step_size * location_step
            candidate_scale <- scale_coefficients + step_size * scale_step
            candidate_value <- objective(candidate_location, candidate_scale)
            if (is.finite(candidate_value) && candidate_value < value) {
                location_coefficients <- candidate_location
                scale_coefficients <- candidate_scale
                value <- candidate_value
                accepted <- TRUE
                break
            }
            step_size <- step_size / 2
        }
        if (gradient_norm < tolerance ||
                (accepted && max(abs(step_size * c(
                    location_step, scale_step
                ))) < tolerance / 100)) {
            converged <- TRUE
            break
        }
        if (!accepted) break
    }
    location_predictor <- drop(
        X_location %*% location_coefficients + location$offset
    )
    scale_predictor <- drop(X_scale %*% scale_coefficients + scale$offset)
    q <- exp(scale_predictor)
    sigma <- b + q
    residual <- y - location_predictor
    cross_weight <- weights * 2 * q * residual / sigma^3
    scale_weight <- weights * (
        q * (1 / sigma - residual^2 / sigma^3) +
            q^2 * (-1 / sigma^2 + 3 * residual^2 / sigma^4)
    )
    location_hessian <- crossprod(
        X_location * sqrt(weights / sigma^2)
    )
    cross_hessian <- crossprod(X_location, X_scale * cross_weight)
    scale_hessian <- crossprod(X_scale, X_scale * scale_weight)
    likelihood_hessian <- rbind(
        cbind(location_hessian, cross_hessian),
        cbind(t(cross_hessian), scale_hessian)
    )
    penalty <- as.matrix(Matrix::bdiag(
        location_penalty, scale_penalty
    ))
    hessian <- likelihood_hessian + penalty
    list(
        coefficients=c(location_coefficients, scale_coefficients),
        parameter_coefficients=list(
            location=location_coefficients,
            scale=scale_coefficients
        ),
        linear_predictors=cbind(
            location=location_predictor,
            scale=scale_predictor
        ),
        fitted_values=cbind(
            location=location_predictor,
            scale=1 / sigma
        ),
        sigma=sigma,
        residuals=residual,
        objective=value,
        likelihood_hessian=likelihood_hessian,
        hessian=hessian,
        penalty=penalty,
        parameter_penalties=list(
            location=location_penalty,
            scale=scale_penalty
        ),
        converged=converged,
        iterations=iteration,
        gradient_norm=gradient_norm,
        prior_weights=weights
    )
}

.cdrgam_gaulss_dense_laml <- function(setups, family) {
    counts <- vapply(setups, function(setup) length(setup$S), integer(1))
    parameter_count <- sum(counts)
    initial <- unlist(lapply(setups, .cdrgam_dense_initial_sp), use.names=FALSE)
    initial_log_sp <- log(initial)
    b <- .cdrgam_gaulss_b(family)
    evaluations <- 0L
    evaluate <- function(log_sp, retain=FALSE) {
        evaluations <<- evaluations + 1L
        if (any(!is.finite(log_sp)) || any(abs(log_sp) > 25)) {
            return(if (retain) NULL else 1e100)
        }
        solution <- tryCatch(
            .cdrgam_gaulss_fixed(setups, exp(log_sp), b),
            error=function(error) NULL
        )
        if (is.null(solution) || !isTRUE(solution$converged)) {
            return(if (retain) NULL else 1e100)
        }
        factor <- tryCatch(chol(solution$hessian), error=function(error) NULL)
        if (is.null(factor)) return(if (retain) NULL else 1e100)
        penalty_determinant <- sum(vapply(
            solution$parameter_penalties,
            function(penalty) .positive_log_determinant(penalty)$value,
            numeric(1)
        ))
        criterion <- 2 * solution$objective +
            2 * sum(log(diag(factor))) - penalty_determinant
        if (!is.finite(criterion)) {
            return(if (retain) NULL else 1e100)
        }
        if (!retain) return(criterion)
        solution$criterion <- criterion
        solution$sp <- exp(log_sp)
        solution$factor <- factor
        solution
    }
    if (!parameter_count) {
        solution <- .cdrgam_gaulss_fixed(setups, numeric(), b)
        solution$criterion <- 2 * solution$objective +
            as.numeric(determinant(solution$hessian, logarithm=TRUE)$modulus)
        solution$sp <- numeric()
        solution$factor <- chol(solution$hessian)
        return(list(
            optimization=list(
                par=numeric(), value=solution$criterion,
                convergence=0L,
                counts=c(`function`=1L, gradient=NA_integer_)
            ),
            solution=solution,
            evaluations=1L
        ))
    }
    optimization <- stats::optim(
        initial_log_sp,
        evaluate,
        method='BFGS',
        control=list(
            maxit=500L,
            reltol=1e-9,
            ndeps=rep.int(1e-3, parameter_count)
        )
    )
    solution <- evaluate(optimization$par, retain=TRUE)
    if (is.null(solution)) {
        stop('Dense gaulss LAML optimization ended at an invalid fit')
    }
    list(
        optimization=optimization,
        solution=solution,
        evaluations=evaluations
    )
}

.cdrgam_distributional_smooths <- function(setups) {
    output <- list()
    offset <- 0L
    for (parameter in names(setups)) {
        for (smooth in setups[[parameter]]$smooth) {
            smooth$first.para <- smooth$first.para + offset
            smooth$last.para <- smooth$last.para + offset
            smooth$label <- paste0(parameter, ':', smooth$label)
            output[[length(output) + 1L]] <- smooth
        }
        offset <- offset + ncol(setups[[parameter]]$X)
    }
    output
}

.cdrgam_distributional_block_terms <- function(design, setups) {
    metadata <- list()
    labels <- character()
    parameter_terms <- vector('list', length(setups))
    names(parameter_terms) <- names(setups)
    offset <- 0L
    for (parameter in names(setups)) {
        terms <- design$parameters[[parameter]]$terms
        parameter_terms[[parameter]] <- integer(length(terms))
        for (term_index in seq_along(terms)) {
            info <- .cdrgam_distributional_term_metadata(
                terms[[term_index]],
                list(smooth=setups[[parameter]]$smooth),
                paste0('cdr_term_', term_index),
                parameter
            )
            info$coefficient_index <- info$coefficient_index + offset
            metadata[[length(metadata) + 1L]] <- info
            labels[[length(labels) + 1L]] <- paste0(
                parameter, ':', terms[[term_index]]$name
            )
            parameter_terms[[parameter]][[term_index]] <- length(metadata)
        }
        offset <- offset + ncol(setups[[parameter]]$X)
    }
    list(metadata=metadata, labels=labels, parameter_terms=parameter_terms)
}

.fit_distributional_block <- function(
        design,
        family,
        method=NULL,
        trace=FALSE,
        ...
) {
    dots <- list(...)
    if (!is.null(dots$weights) && any(dots$weights != 1)) {
        stop(
            'gaulss does not support non-unit prior weights; mgcv accepts ',
            'but ignores them'
        )
    }
    if (!is.null(method) && !(method %in% c('REML', 'fREML'))) {
        stop('The distributional block backend currently supports only REML')
    }
    reporter <- .new_solver_reporter(trace, 'distributional block')
    reporter$phase('parameter model setup')
    design$parameters <- lapply(
        design$parameters,
        .materialize_cdr_design,
        sparse=FALSE
    )
    setups <- lapply(design$parameters, function(parameter_design) {
        do.call(
            .cdrgam_distributional_block_setup,
            c(list(design=parameter_design), dots)
        )
    })
    setups <- lapply(setups, function(setup) {
        .audit_mgcv_setup(setup, .rank_tolerance(NULL))$setup
    })
    y <- setups$location$y
    if (!identical(y, setups$scale$y)) {
        stop('Distributional parameter setups retained different responses')
    }
    weights <- lapply(setups, function(setup) {
        if (is.null(setup$w)) rep.int(1, length(y)) else setup$w
    })
    if (!isTRUE(all.equal(weights$location, weights$scale))) {
        stop('Distributional parameter setups retained different weights')
    }
    reporter$phase(
        'joint smoothing-parameter optimization',
        smoothing_parameters=sum(vapply(
            setups, function(setup) length(setup$S), integer(1)
        )),
        gradient='finite',
        inner_solver='Fisher scoring'
    )
    result <- .cdrgam_gaulss_dense_laml(setups, family)
    solution <- result$solution
    dimensions <- vapply(setups, function(setup) ncol(setup$X), integer(1))
    ranges <- split(seq_len(sum(dimensions)), rep(names(dimensions), dimensions))
    coefficient_names <- unlist(Map(function(parameter, setup) {
        paste0(parameter, ':', colnames(setup$X))
    }, names(setups), setups), use.names=FALSE)
    coefficients <- stats::setNames(solution$coefficients, coefficient_names)
    covariance <- chol2inv(solution$factor)
    dimnames(covariance) <- list(coefficient_names, coefficient_names)
    influence <- diag(covariance %*% solution$likelihood_hessian)
    smooths <- .cdrgam_distributional_smooths(setups)
    term_info <- .cdrgam_distributional_block_terms(design, setups)
    sp_names <- unlist(Map(function(parameter, setup) {
        names <- names(setup$sp)
        if (is.null(names)) names <- paste0('sp', seq_along(setup$S))
        paste0(parameter, ':', names)
    }, names(setups), setups), use.names=FALSE)
    smoothing_parameters <- stats::setNames(solution$sp, sp_names)
    parametric_indices <- unlist(Map(function(range, setup) {
        range[seq_len(setup$nsdf)]
    }, ranges, setups), use.names=FALSE)
    preparation <- lapply(
        design$parameters,
        .cdrgam_distributional_preparation
    )
    prediction_setups <- lapply(setups, .cdr_prediction_setup)
    identifiability <- lapply(
        design$parameters, `[[`, 'identifiability'
    )
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
        Vp=covariance,
        Vc=NULL,
        edf=sum(influence),
        coefficient.edf=influence,
        df.residual=length(y) - sum(influence),
        y=y,
        prior.weights=solution$prior_weights,
        X=do.call(cbind, lapply(setups, `[[`, 'X')),
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
        converged=isTRUE(solution$converged) &&
            identical(result$optimization$convergence, 0L),
        distributional=list(
            hessian=solution$hessian,
            likelihood_hessian=solution$likelihood_hessian,
            sigma=solution$sigma,
            coefficient_ranges=ranges,
            inner_iterations=solution$iterations,
            gradient_norm=solution$gradient_norm
        ),
        cdrgam=list(
            schema_version=1L,
            engine='block',
            backend='block',
            distributional=TRUE,
            parameter_names=names(setups),
            parameter_terms=term_info$parameter_terms,
            parametric_indices=parametric_indices,
            formula=list(
                user=design$formula,
                normalized=design$normalized_formula,
                effective=design$effective_formula,
                mgcv=lapply(setups, `[[`, 'formula')
            ),
            preparation=list(parameters=preparation),
            scaling=NULL,
            identifiability=identifiability,
            term_labels=term_info$labels,
            terms=term_info$metadata,
            prediction=list(
                setups=prediction_setups,
                coefficient_ranges=ranges
            ),
            rank=list(
                action='error',
                resolution='identified',
                regularization=0,
                tolerance=.rank_tolerance(NULL)
            ),
            solver='dense Gaussian location-scale LAML reference solver'
        )
    )
    class(output) <- c(
        'cdrgam_distributional_block', 'cdrgam_block', 'cdrgam'
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

.cdrgam_distributional_block_lpmatrices <- function(
        object,
        impulses,
        responses,
        chunk_size
) {
    output <- list()
    for (parameter in object$cdrgam$parameter_names) {
        term_indices <- object$cdrgam$parameter_terms[[parameter]]
        shim <- object
        shim$cdrgam$preparation <-
            object$cdrgam$preparation$parameters[[parameter]]
        shim$cdrgam$scaling <- shim$cdrgam$preparation$scaling
        shim$cdrgam$terms <- object$cdrgam$terms[term_indices]
        matrices <- .cdr_predict_irf_matrices(
            shim,
            impulses,
            responses,
            chunk_size,
            source_impulses=impulses,
            source_responses=responses
        )
        names(matrices) <- paste0('cdr_term_', seq_along(matrices))
        output[[parameter]] <- .cdr_setup_lpmatrix(
            object$cdrgam$prediction$setups[[parameter]],
            responses,
            matrices
        )
    }
    output
}

#' @rdname predict.cdrgam
#' @export
predict.cdrgam_distributional_block <- function(
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
        matrices <- lapply(object$cdrgam$parameter_names, function(parameter) {
            range <- object$cdrgam$prediction$coefficient_ranges[[parameter]]
            object$X[, range, drop=FALSE]
        })
        names(matrices) <- object$cdrgam$parameter_names
        return(matrices)
    }
    if (!is.list(newdata) ||
            !all(c('impulses', 'responses') %in% names(newdata))) {
        stop('newdata must contain impulse and response data frames')
    }
    matrices <- .cdrgam_distributional_block_lpmatrices(
        object,
        newdata$impulses,
        newdata$responses,
        chunk_size
    )
    if (identical(type, 'lpmatrix')) {
        return(lapply(matrices, `[[`, 'X'))
    }
    links <- matrix(
        NA_real_,
        nrow=nrow(newdata$responses),
        ncol=length(object$cdrgam$parameter_names),
        dimnames=list(NULL, object$cdrgam$parameter_names)
    )
    standard_errors <- links
    for (parameter in object$cdrgam$parameter_names) {
        range <- object$cdrgam$prediction$coefficient_ranges[[parameter]]
        assembled <- matrices[[parameter]]
        links[, parameter] <- drop(
            assembled$X %*% object$coefficients[range] + assembled$offset
        )
        if (isTRUE(se.fit)) {
            covariance <- object$Vp[range, range, drop=FALSE]
            standard_errors[, parameter] <- sqrt(pmax(
                0,
                rowSums((assembled$X %*% covariance) * assembled$X)
            ))
        }
    }
    fit <- links
    if (identical(type, 'response')) {
        fit[, 'scale'] <- object$family$linfo[[2L]]$linkinv(
            links[, 'scale']
        )
        if (isTRUE(se.fit)) {
            standard_errors[, 'scale'] <- standard_errors[, 'scale'] *
                abs(object$family$linfo[[2L]]$mu.eta(links[, 'scale']))
        }
    }
    if (!isTRUE(se.fit)) return(fit)
    list(fit=fit, se.fit=standard_errors)
}

#' @export
logLik.cdrgam_distributional_block <- function(object, ...) {
    value <- object$loglik
    attr(value, 'df') <- object$edf
    attr(value, 'nobs') <- length(object$y)
    class(value) <- 'logLik'
    value
}

#' @export
summary.cdrgam_distributional_block <- function(
        object,
        dispersion=NULL,
        freq=FALSE,
        re.test=TRUE,
        all.coefficients=FALSE,
        ...
) {
    smooth_edf <- vapply(object$smooth, function(smooth) {
        sum(object$coefficient.edf[smooth$first.para:smooth$last.para])
    }, numeric(1))
    covariance <- function(indices) {
        object$Vp[indices, indices, drop=FALSE]
    }
    parametric <- object$cdrgam$parametric_indices
    residual_df <- object$df.residual
    p_table <- if (length(parametric)) {
        .cdrgam_coefficient_table(
            object,
            parametric,
            covariance(parametric),
            residual_df,
            reference='z'
        )
    } else NULL
    s_table <- .cdrgam_smooth_table(
        object,
        smooth_edf,
        covariance,
        residual_df,
        chi_square=TRUE
    )
    formulas <- .cdrgam_summary_formulas(object)
    output <- list(
        call=object$call,
        family=object$family,
        formula=formulas$user,
        formulas=formulas,
        formula_strings=.cdrgam_formula_strings(formulas),
        p.coeff=if (is.null(p_table)) numeric() else p_table[, 'Estimate'],
        p.t=if (is.null(p_table)) numeric() else p_table[, 'z value'],
        p.pv=if (is.null(p_table)) numeric() else p_table[, 'Pr(>|z|)'],
        p.table=p_table,
        s.table=s_table,
        se=if (is.null(p_table)) numeric() else p_table[, 'Std. Error'],
        chi.sq=if (is.null(s_table)) numeric() else s_table[, 'Chi.sq'],
        s.pv=if (is.null(s_table)) numeric() else s_table[, 'p-value'],
        pTerms.pv=numeric(),
        pTerms.chi.sq=numeric(),
        pTerms.df=numeric(),
        m=length(object$smooth),
        edf=smooth_edf,
        residual.df=residual_df,
        scale=1,
        dispersion=1,
        r.sq=NA_real_,
        dev.expl=NA_real_,
        method='-REML',
        sp.criterion=object$reml,
        rank=length(object$coefficients),
        np=length(object$coefficients),
        n=length(object$y),
        sp=object$sp,
        backend=object$cdrgam$solver,
        rank_metadata=object$cdrgam$rank
    )
    if (isTRUE(all.coefficients)) {
        indices <- seq_along(object$coefficients)
        output$coefficients <- .cdrgam_coefficient_table(
            object,
            indices,
            covariance(indices),
            residual_df,
            reference='z'
        )
    }
    class(output) <- c(
        'summary.cdrgam_distributional_block',
        'summary.cdrgam_block',
        'summary.cdrgam',
        'summary.gam'
    )
    output
}
