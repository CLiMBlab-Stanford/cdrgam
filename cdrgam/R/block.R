.embed_penalties <- function(setup, smoothing_parameters) {
    p <- ncol(setup$X)
    total <- matrix(0, p, p)
    if (!length(setup$S)) {
        return(total)
    }
    for (i in seq_along(setup$S)) {
        indices <- setup$off[[i]] + seq_len(nrow(setup$S[[i]])) - 1L
        total[indices, indices] <- total[indices, indices] +
            smoothing_parameters[[i]] * setup$S[[i]]
    }
    total
}

.positive_log_determinant <- function(matrix) {
    values <- eigen(matrix, symmetric=TRUE, only.values=TRUE)$values
    # Use a machine-precision rank threshold. A sqrt(eps) threshold incorrectly
    # drops a valid penalty block when different smooths have very different
    # smoothing parameters, which also corrupts the REML log determinant.
    tolerance <- max(1, max(abs(values))) * .Machine$double.eps *
        max(100, nrow(matrix))
    positive <- values[values > tolerance]
    list(
        value=if (length(positive)) sum(log(positive)) else 0,
        rank=length(positive)
    )
}

.rank_tolerance <- function(value) {
    if (is.null(value)) return(sqrt(.Machine$double.eps))
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
            value <= 0) {
        stop('rank_tol must be a finite positive number')
    }
    value
}

.rank_penalty <- function(value) {
    if (is.null(value)) return(1e-6)
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
            value <= 0) {
        stop('rank_penalty must be a finite positive number')
    }
    value
}

.rank_error <- function(action=NULL) {
    suffix <- if (identical(action, 'drop')) {
        paste0(
            ' rank_action="drop" can only remove ordinary parametric aliases; ',
            'this deficiency spans smooth or IRF terms. Use ',
            'rank_action="minimum_norm" for a prediction-oriented fit or ',
            'rank_action="penalize" for regularized separation.'
        )
    } else {
        paste0(
            ' Revise the model, or refit with rank_action="minimum_norm" ',
            '(prediction-oriented), rank_action="drop" (ordinary aliases), ',
            'or rank_action="penalize" (regularized separation).'
        )
    }
    stop(
        'The penalized model matrix is rank deficient; one or more model ',
        'directions are not identifiable.', suffix,
        call.=FALSE
    )
}

.drop_parametric_aliases <- function(setup, tolerance) {
    count <- setup$nsdf
    info <- list(
        corrected=FALSE,
        original_dimension=ncol(setup$X),
        fitted_dimension=ncol(setup$X),
        dropped=character(),
        tolerance=tolerance,
        resolution='none'
    )
    if (count < 2L) return(list(setup=setup, info=info))
    parametric <- setup$X[, seq_len(count), drop=FALSE]
    decomposition <- qr(parametric, tol=tolerance, LAPACK=FALSE)
    if (decomposition$rank == count) return(list(setup=setup, info=info))
    if (!is.null(setup$paraPen) && length(setup$paraPen)) {
        stop(
            'Automatic correction of aliased penalized parametric terms is ',
            'not yet supported; revise paraPen or the model formula.',
            call.=FALSE
        )
    }
    keep_parametric <- sort(decomposition$pivot[seq_len(decomposition$rank)])
    dropped <- setdiff(seq_len(count), keep_parametric)
    smooth_columns <- if (count < ncol(setup$X)) {
        seq.int(count + 1L, ncol(setup$X))
    } else {
        integer()
    }
    keep <- c(keep_parametric, smooth_columns)
    setup$X <- setup$X[, keep, drop=FALSE]
    setup$nsdf <- length(keep_parametric)
    if (!is.null(setup$term.names)) {
        setup$term.names <- colnames(setup$X)
    }
    if (!is.null(setup$cmX) && length(setup$cmX) == length(keep)) {
        # This branch is retained for unusual setup objects whose centering
        # vector has already been reduced by their constructor.
        setup$cmX <- setup$cmX
    } else if (!is.null(setup$cmX) &&
            length(setup$cmX) == info$original_dimension) {
        setup$cmX <- setup$cmX[keep]
    }
    if (!is.null(setup$assign) && length(setup$assign) == count) {
        setup$assign <- setup$assign[keep_parametric]
    }
    shift <- length(dropped)
    if (length(setup$off)) setup$off <- setup$off - shift
    if (length(setup$smooth)) {
        setup$smooth <- lapply(setup$smooth, function(smooth) {
            smooth$first.para <- smooth$first.para - shift
            smooth$last.para <- smooth$last.para - shift
            smooth
        })
    }
    info$corrected <- TRUE
    info$fitted_dimension <- ncol(setup$X)
    info$dropped <- colnames(parametric)[dropped]
    info$resolution <- 'automatic-parametric-alias-drop'
    warning(
        'Removed aliased ordinary parametric coefficient(s): ',
        paste(info$dropped, collapse=', '),
        call.=FALSE
    )
    list(setup=setup, info=info)
}

.rank_regularization <- function(system, action, tolerance, penalty) {
    dense <- is.matrix(system)
    probe <- if (dense) {
        tryCatch(chol(system), error=function(e) NULL)
    } else {
        tryCatch(
            suppressWarnings(Matrix::Cholesky(
                Matrix::Matrix(system, sparse=TRUE),
                LDL=FALSE, perm=TRUE, super=FALSE
            )),
            error=function(e) NULL
        )
    }
    condition <- if (dense && !is.null(probe)) {
        tryCatch(rcond(system), error=function(e) 0)
    } else {
        NA_real_
    }
    if (!is.null(probe) && (!dense ||
            (is.finite(condition) && condition > tolerance))) {
        return(list(value=0, resolution='none'))
    }
    if (action %in% c('error', 'drop')) .rank_error(action)
    diagonal_scale <- max(abs(Matrix::diag(system)), 1)
    relative <- if (identical(action, 'minimum_norm')) tolerance else penalty
    value <- relative * diagonal_scale
    warning(
        if (identical(action, 'minimum_norm')) {
            paste0(
                'The model is rank deficient; applying a deterministic ',
                'minimum-norm ridge approximation (relative ridge ',
                format(relative, digits=4), '). Individual coefficients in ',
                'the unidentified subspace are convention-dependent.'
            )
        } else {
            paste0(
                'The model is rank deficient; applying a fixed ridge penalty ',
                '(relative penalty ', format(relative, digits=4), ').'
            )
        },
        call.=FALSE
    )
    list(
        value=value,
        resolution=if (identical(action, 'minimum_norm')) {
            'minimum-norm-ridge-approximation'
        } else {
            'fixed-ridge-penalty'
        }
    )
}

.setup_column_owners <- function(setup) {
    owners <- colnames(setup$X)
    if (setup$nsdf > 0L) {
        owners[seq_len(setup$nsdf)] <- paste0(
            'parametric:', owners[seq_len(setup$nsdf)]
        )
    }
    for (smooth in setup$smooth) {
        owners[smooth$first.para:smooth$last.para] <- smooth$label
    }
    owners
}

.audit_mgcv_setup <- function(setup, tolerance) {
    aliases <- .drop_parametric_aliases(setup, tolerance)
    setup <- aliases$setup
    weights <- setup$w
    if (is.null(weights)) weights <- rep.int(1, nrow(setup$X))
    weighted <- setup$X * sqrt(weights)
    system <- crossprod(weighted) +
        .embed_penalties(setup, rep.int(1, length(setup$S)))
    factor <- tryCatch(chol(system), error=function(e) NULL)
    condition <- if (is.null(factor)) 0 else
        tryCatch(rcond(system), error=function(e) 0)
    decomposition <- qr(system, tol=tolerance, LAPACK=FALSE)
    if (is.null(factor) || decomposition$rank < ncol(system)) {
        dependent <- if (decomposition$rank < ncol(system)) {
            decomposition$pivot[seq.int(decomposition$rank + 1L, ncol(system))]
        } else {
            seq_len(ncol(system))
        }
        implicated <- unique(.setup_column_owners(setup)[dependent])
        stop(
            'The globally assembled model remains rank deficient after ',
            'canonical constraints and benign alias removal. Implicated ',
            'terms: ', paste(implicated, collapse=', '), '. Revise the model ',
            'or use a custom backend with an explicit rank_action policy.',
            call.=FALSE
        )
    }
    list(
        setup=setup,
        info=list(
            dimension=ncol(system),
            rank=decomposition$rank,
            condition_indicator=condition,
            parametric=aliases$info,
            resolution='identified'
        )
    )
}

.fit_block_gaussian <- function(
        design,
        family,
        method,
        checkpoint=NULL,
        trace=FALSE,
        rank_action='error',
        rank_tol=NULL,
        rank_penalty=NULL,
        ...
) {
    family <- .as_family(family)
    reporter <- .new_solver_reporter(trace, 'block')
    reporter$phase('model setup')
    if (!identical(family$family, 'gaussian') ||
            !identical(family$link, 'identity')) {
        stop('The initial block backend supports only gaussian(identity)')
    }
    if (!is.null(method) && !(method %in% c('REML', 'fREML'))) {
        stop('The initial block backend supports only REML')
    }
    setup <- .fit_compressed_mgcv(
        y=design$responses[[design$response_name]],
        terms=design$terms,
        family=family,
        method='REML',
        engine='gam',
        base_formula=design$ordinary_formula,
        response_data=design$responses,
        user_formula=design$formula,
        preparation=list(
            configuration=design$configuration,
            plan=design$plan,
            stream=design$stream,
            specification=design$specification,
            simplifications=design$simplifications,
            scaling=design$scaling
        ),
        setup_only=TRUE,
        ...
    )
    tolerance <- .rank_tolerance(rank_tol)
    penalty_strength <- .rank_penalty(rank_penalty)
    alias_resolution <- .drop_parametric_aliases(setup, tolerance)
    setup <- alias_resolution$setup
    X <- setup$X
    y <- setup$y - setup$offset
    weights <- setup$w
    if (is.null(weights)) {
        weights <- rep.int(1, length(y))
    }
    if (any(!is.finite(weights)) || any(weights <= 0)) {
        stop('The initial block backend requires finite positive weights')
    }
    root_weights <- sqrt(weights)
    weighted_X <- X * root_weights
    weighted_y <- y * root_weights
    reporter$phase(
        'cross-product accumulation',
        observations=nrow(X),
        coefficients=ncol(X)
    )
    XtX <- crossprod(weighted_X)
    Xty <- crossprod(weighted_X, weighted_y)
    penalty_count <- length(setup$S)
    unit_penalty <- .embed_penalties(setup, rep.int(1, penalty_count))
    rank_resolution <- .rank_regularization(
        XtX + unit_penalty,
        rank_action,
        tolerance,
        penalty_strength
    )
    fixed_ridge <- rank_resolution$value
    signature_system <- XtX + unit_penalty + diag(fixed_ridge, ncol(X))
    checkpoint_signature <- .solver_checkpoint_signature(
        backend='block',
        formula=design$formula,
        observation_count=nrow(X),
        dimension=ncol(X),
        sp_names=names(setup$sp),
        Xty=Xty,
        system=signature_system
    )
    checkpoint_signature$optimizer <- list(
        gradient='finite',
        restarts=0L,
        lower=-25,
        upper=25
    )
    checkpoint_state <- .checkpoint_validate(
        .checkpoint_read(checkpoint),
        checkpoint_signature,
        checkpoint
    )
    resumed <- !is.null(checkpoint_state)
    if (is.null(checkpoint_state) || isTRUE(checkpoint_state$legacy)) {
        legacy_state <- checkpoint_state
        checkpoint_state <- .new_checkpoint_state(
            checkpoint_signature,
            'block',
            0L
        )
        if (!is.null(legacy_state$log_sp)) {
            checkpoint_state$current_log_sp <- legacy_state$log_sp
            checkpoint_state$best_log_sp <- legacy_state$log_sp
            checkpoint_state$current_criterion <- legacy_state$criterion
            checkpoint_state$best_criterion <- legacy_state$criterion
        }
    }
    checkpoint_best <- checkpoint_state$best_criterion
    if (is.null(checkpoint_best) || !is.finite(checkpoint_best)) {
        checkpoint_best <- Inf
    }
    evaluation_count <- if (is.null(checkpoint_state$evaluation_count)) {
        0L
    } else {
        as.integer(checkpoint_state$evaluation_count)
    }
    checkpoint_last_written <- evaluation_count
    checkpoint_last_written_time <- proc.time()[['elapsed']]
    last_level_one_report <- -Inf
    if (resumed) {
        reporter$emit(
            1L,
            'checkpoint resumed',
            path=checkpoint,
            stage=checkpoint_state$stage,
            evaluations=evaluation_count,
            best=checkpoint_best
        )
    }

    evaluate <- function(log_sp, retain=FALSE) {
        if (!retain) evaluation_count <<- evaluation_count + 1L
        evaluation_started <- proc.time()[['elapsed']]
        sp <- exp(log_sp)
        penalty <- .embed_penalties(setup, sp)
        penalty_det <- .positive_log_determinant(penalty)
        unpenalized_dimension <- ncol(X) - penalty_det$rank
        residual_df <- nrow(X) - unpenalized_dimension
        if (residual_df <= 0) {
            stop('Insufficient observations for Gaussian REML estimation')
        }
        system <- XtX + penalty + diag(fixed_ridge, ncol(X))
        factor <- tryCatch(chol(system), error=function(e) NULL)
        if (is.null(factor)) {
            return(if (retain) NULL else .Machine$double.xmax / 100)
        }
        coefficients <- backsolve(
            factor,
            forwardsolve(t(factor), Xty)
        )
        residual <- weighted_y - weighted_X %*% coefficients
        penalized_rss <- sum(residual^2) +
            drop(crossprod(coefficients, penalty %*% coefficients)) +
            fixed_ridge * sum(coefficients^2)
        if (!is.finite(penalized_rss) || penalized_rss <= 0) {
            return(if (retain) NULL else .Machine$double.xmax / 100)
        }
        log_det_system <- 2 * sum(log(diag(factor)))
        criterion <- residual_df * log(penalized_rss / residual_df) +
            log_det_system - penalty_det$value
        if (!is.finite(criterion)) {
            return(if (retain) NULL else .Machine$double.xmax / 100)
        }
        if (!retain) {
            new_best <- criterion < checkpoint_best
            now <- proc.time()[['elapsed']]
            level_one_due <- evaluation_count == 1L ||
                evaluation_count %% 10L == 0L ||
                (new_best && now - last_level_one_report >= 5)
            should_report <- reporter$level >= 2L || level_one_due
            if (should_report) {
                reporter$emit(
                    if (level_one_due) 1L else 2L,
                    if (new_best) 'new best' else 'evaluation',
                    evaluation=evaluation_count,
                    criterion=format(criterion, digits=10),
                    seconds=format(
                        proc.time()[['elapsed']] - evaluation_started,
                        digits=5
                    ),
                    sp=format(sp, digits=5)
                )
                if (level_one_due) last_level_one_report <<- now
            }
            checkpoint_state$current_log_sp <<- log_sp
            checkpoint_state$current_criterion <<- criterion
            checkpoint_state$evaluation_count <<- evaluation_count
            if (new_best) {
                checkpoint_best <<- criterion
                checkpoint_state$best_log_sp <<- log_sp
                checkpoint_state$best_criterion <<- criterion
            }
            if (!is.null(checkpoint) &&
                    (evaluation_count - checkpoint_last_written >= 10L ||
                     now - checkpoint_last_written_time >= 60)) {
                checkpoint_state$updated_at <<- as.character(Sys.time())
                .checkpoint_write(checkpoint_state, checkpoint)
                checkpoint_last_written <<- evaluation_count
                checkpoint_last_written_time <<- now
            }
        }
        if (!retain) {
            return(criterion)
        }
        list(
            criterion=criterion,
            sp=sp,
            penalty=penalty,
            factor=factor,
            coefficients=drop(coefficients),
            penalized_rss=penalized_rss,
            reml_df=residual_df
        )
    }

    # mgcv's variance-component intervals use the Hessian of the negative
    # restricted log likelihood with respect to log smoothing parameters and
    # (for Gaussian REML) log scale. `evaluate()` profiles scale out, so form
    # the equivalent unprofiled criterion at the solution after optimization.
    evaluate_unprofiled <- function(parameters) {
        log_sp <- parameters[seq_len(penalty_count)]
        log_scale <- parameters[[penalty_count + 1L]]
        sp <- exp(log_sp)
        scale <- exp(log_scale)
        penalty <- .embed_penalties(setup, sp)
        penalty_det <- .positive_log_determinant(penalty)
        unpenalized_dimension <- ncol(X) - penalty_det$rank
        residual_df <- nrow(X) - unpenalized_dimension
        system <- XtX + penalty + diag(fixed_ridge, ncol(X))
        factor <- tryCatch(chol(system), error=function(e) NULL)
        if (is.null(factor) || !is.finite(scale) || scale <= 0) {
            return(.Machine$double.xmax / 100)
        }
        coefficients <- backsolve(
            factor,
            forwardsolve(t(factor), Xty)
        )
        residual <- weighted_y - weighted_X %*% coefficients
        penalized_rss <- sum(residual^2) +
            drop(crossprod(coefficients, penalty %*% coefficients)) +
            fixed_ridge * sum(coefficients^2)
        log_det_system <- 2 * sum(log(diag(factor)))
        criterion <- penalized_rss / scale + residual_df * log_scale +
            log_det_system - penalty_det$value
        if (is.finite(criterion)) criterion else .Machine$double.xmax / 100
    }

    if (penalty_count) {
        canonical_initial <- rep.int(0, penalty_count)
        initial <- canonical_initial
        checkpoint_initial <- if (!is.null(checkpoint_state$current_log_sp)) {
            checkpoint_state$current_log_sp
        } else {
            checkpoint_state$best_log_sp
        }
        if (length(checkpoint_initial) == penalty_count &&
                all(is.finite(checkpoint_initial))) {
            initial <- checkpoint_initial
        }
        reporter$phase(
            'smoothing-parameter optimization',
            smoothing_parameters=penalty_count,
            gradient='finite'
        )
        if (checkpoint_state$stage %in% c('optimization_complete', 'complete') &&
                !is.null(checkpoint_state$optimization)) {
            optimization <- checkpoint_state$optimization
        } else {
            checkpoint_state$current_log_sp <- initial
            if (!is.null(checkpoint)) .checkpoint_write(checkpoint_state, checkpoint)
            optimization <- stats::optim(
                initial,
                evaluate,
                method='L-BFGS-B',
                lower=rep.int(-25, penalty_count),
                upper=rep.int(25, penalty_count),
                control=list(factr=1e7)
            )
            if (resumed && !isTRUE(all.equal(
                    initial,
                    canonical_initial,
                    tolerance=0
            ))) {
                reporter$emit(
                    1L,
                    'restart validation started',
                    warm_criterion=format(optimization$value, digits=10)
                )
                canonical_optimization <- stats::optim(
                    canonical_initial,
                    evaluate,
                    method='L-BFGS-B',
                    lower=rep.int(-25, penalty_count),
                    upper=rep.int(25, penalty_count),
                    control=list(factr=1e7)
                )
                if (canonical_optimization$value < optimization$value) {
                    optimization <- canonical_optimization
                }
            }
            checkpoint_state$stage <- 'optimization_complete'
            checkpoint_state$optimization <- optimization
            checkpoint_state$completed_runs <- list(optimization)
            checkpoint_state$best_log_sp <- optimization$par
            checkpoint_state$best_criterion <- optimization$value
            checkpoint_state$evaluation_count <- evaluation_count
            if (!is.null(checkpoint)) .checkpoint_write(checkpoint_state, checkpoint)
        }
        solution <- evaluate(optimization$par, retain=TRUE)
    } else {
        optimization <- list(
            par=numeric(),
            value=NA_real_,
            convergence=0L,
            message=NULL,
            counts=c(`function`=1L, gradient=NA_integer_)
        )
        solution <- evaluate(numeric(), retain=TRUE)
    }
    coefficients <- solution$coefficients
    names(coefficients) <- colnames(X)
    fitted_values <- drop(X %*% coefficients + setup$offset)
    raw_residuals <- setup$y - fitted_values
    inverse_system <- chol2inv(solution$factor)
    hat_trace <- sum(diag(XtX %*% inverse_system))
    scale <- solution$penalized_rss / solution$reml_df
    hessian_parameters <- c(log(solution$sp), log(scale))
    # evaluate_unprofiled is -2 REML, whereas mgcv stores the Hessian of the
    # negative log likelihood used by gam.vcomp().
    reporter$phase('outer Hessian', parameters=penalty_count + 1L)
    outer_hessian <- stats::optimHess(
        hessian_parameters,
        evaluate_unprofiled
    ) / 2
    covariance <- inverse_system * scale
    dimnames(covariance) <- list(names(coefficients), names(coefficients))
    rho_covariance <- tryCatch(
        solve(outer_hessian)[
            seq_len(penalty_count),
            seq_len(penalty_count),
            drop=FALSE
        ],
        error=function(e) NULL
    )
    unconditional_covariance <- NULL
    if (!is.null(rho_covariance) && all(is.finite(rho_covariance))) {
        coefficient_derivatives <- matrix(
            0,
            nrow=length(coefficients),
            ncol=penalty_count
        )
        for (i in seq_len(penalty_count)) {
            component_sp <- numeric(penalty_count)
            component_sp[[i]] <- solution$sp[[i]]
            derivative <- .embed_penalties(setup, component_sp)
            coefficient_derivatives[, i] <- -backsolve(
                solution$factor,
                forwardsolve(
                    t(solution$factor),
                    derivative %*% coefficients
                )
            )
        }
        unconditional_covariance <- covariance +
            coefficient_derivatives %*% rho_covariance %*%
            t(coefficient_derivatives)
        unconditional_covariance <-
            (unconditional_covariance + t(unconditional_covariance)) / 2
        dimnames(unconditional_covariance) <- dimnames(covariance)
    }
    term_metadata <- lapply(seq_along(design$terms), function(i) {
        term <- design$terms[[i]]
        smooth_index <- which(vapply(
            setup$smooth,
            function(smooth) paste(smooth$term, collapse=',') ==
                paste0('cdr_term_', i),
            logical(1)
        ))
        coefficient_index <- if (length(smooth_index) == 1L) {
            seq.int(
                setup$smooth[[smooth_index]]$first.para,
                setup$smooth[[smooth_index]]$last.para
            )
        } else {
            integer()
        }
        list(
            name=term$name,
            type=term$type,
            knots=term$knots,
            predictor_knots=term$predictor_knots,
            axis=term$axis,
            linear_predictors=term$linear_predictors,
            linear_predictor_summaries=term$linear_predictor_summaries,
            basis=term$basis,
            transform=term$transform,
            group=term$group,
            group_levels=term$group_levels,
            base_dimension=term$base_dimension,
            rank=term$rank,
            null.space.dim=term$null.space.dim,
            S.scale=term$S.scale,
            lag_scale=if (is.null(term$lag_scale)) 1 else term$lag_scale,
            predictor_scale=if (is.null(term$predictor_scale)) 1 else
                term$predictor_scale,
            amplitude_scale=if (is.null(term$amplitude_scale)) 1 else
                term$amplitude_scale,
            coefficient_index=coefficient_index
        )
    })
    identifiability <- design$identifiability
    identifiability$global <- list(
        dimension=ncol(X),
        rank=ncol(X),
        condition_indicator=tryCatch(
            rcond(XtX + unit_penalty + diag(fixed_ridge, ncol(X))),
            error=function(e) NA_real_
        ),
        parametric=alias_resolution$info,
        resolution=rank_resolution$resolution
    )

    reporter$phase('finalization')
    out <- list(
        coefficients=coefficients,
        fitted.values=fitted_values,
        residuals=raw_residuals,
        linear.predictors=fitted_values,
        family=family,
        sp=stats::setNames(solution$sp, names(setup$sp)),
        scale=scale,
        sig2=scale,
        reml.scale=scale,
        method='REML',
        smooth=setup$smooth,
        paraPen=setup$paraPen,
        full.sp=setup$full.sp,
        outer.info=list(hess=outer_hessian),
        Vp=covariance,
        Vc=unconditional_covariance,
        edf=hat_trace,
        df.residual=nrow(X) - hat_trace,
        y=setup$y,
        prior.weights=weights,
        X=X,
        offset=setup$offset,
        reml=solution$criterion,
        optimizer=optimization,
        cdrgam=list(
            schema_version=1L,
            engine='block',
            backend='block',
            formula=list(
                user=design$formula,
                normalized=design$normalized_formula,
                effective=design$effective_formula,
                mgcv=setup$formula
            ),
            preparation=list(
                configuration=design$configuration,
                plan=design$plan,
                stream=design$stream,
                specification=design$specification,
                simplifications=design$simplifications,
                scaling=design$scaling,
                identifiability=design$identifiability
            ),
            scaling=design$scaling,
            identifiability=identifiability,
            term_labels=names(design$terms),
            terms=term_metadata,
            prediction=list(setup=.cdr_prediction_setup(setup)),
            rank=list(
                action=rank_action,
                parametric=alias_resolution$info,
                resolution=rank_resolution$resolution,
                regularization=fixed_ridge,
                tolerance=tolerance
            ),
            solver='dense Gaussian REML reference solver'
        )
    )
    class(out) <- c('cdrgam_block', 'cdrgam')
    checkpoint_state$stage <- 'complete'
    checkpoint_state$optimization <- optimization
    checkpoint_state$best_log_sp <- optimization$par
    checkpoint_state$best_criterion <- solution$criterion
    checkpoint_state$evaluation_count <- evaluation_count
    checkpoint_state$converged <- identical(optimization$convergence, 0L)
    if (!is.null(checkpoint)) .checkpoint_write(checkpoint_state, checkpoint)
    reporter$emit(
        1L,
        'fit complete',
        criterion=format(solution$criterion, digits=10),
        evaluations=evaluation_count,
        converged=identical(optimization$convergence, 0L)
    )
    out
}

.fit_block_generalized <- function(
        design,
        family,
        method,
        checkpoint=NULL,
        trace=FALSE,
        rank_action='error',
        rank_tol=NULL,
        rank_penalty=NULL,
        ...
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
            'The generalized block backend currently supports ',
            'binomial(logit), poisson(log), and Gamma(log)'
        )
    }
    if (!is.null(method) && !(method %in% c('REML', 'fREML'))) {
        stop('The generalized block backend currently supports only REML')
    }
    if (!is.null(checkpoint)) {
        stop('Generalized block-backend checkpoints are not yet supported')
    }
    reporter <- .new_solver_reporter(trace, 'block')
    reporter$phase('model setup')
    setup <- .fit_compressed_mgcv(
        y=design$responses[[design$response_name]],
        terms=design$terms,
        family=family,
        method='REML',
        engine='gam',
        base_formula=design$ordinary_formula,
        response_data=design$responses,
        user_formula=design$formula,
        preparation=list(
            configuration=design$configuration,
            plan=design$plan,
            stream=design$stream,
            specification=design$specification,
            simplifications=design$simplifications,
            scaling=design$scaling
        ),
        setup_only=TRUE,
        ...
    )
    tolerance <- .rank_tolerance(rank_tol)
    penalty_strength <- .rank_penalty(rank_penalty)
    alias_resolution <- .drop_parametric_aliases(setup, tolerance)
    setup <- alias_resolution$setup
    unit_penalty <- .embed_penalties(setup, rep.int(1, length(setup$S)))
    rank_resolution <- .rank_regularization(
        crossprod(setup$X) + unit_penalty,
        rank_action,
        tolerance,
        penalty_strength
    )
    fixed_ridge <- rank_resolution$value
    reporter$phase(
        'smoothing-parameter optimization',
        smoothing_parameters=length(setup$S),
        gradient='finite',
        inner_solver='penalized IRLS'
    )
    result <- .cdrgam_dense_laml(
        setup,
        family,
        fixed_ridge=fixed_ridge
    )
    solution <- result$solution
    coefficients <- solution$coefficients
    names(coefficients) <- colnames(setup$X)
    inverse_system <- chol2inv(solution$factor)
    information <- crossprod(setup$X * sqrt(solution$working_weights))
    influence <- diag(inverse_system %*% information)
    edf <- sum(influence)
    reported_scale <- if (estimated_gamma) {
        .cdrgam_gamma_reported_scale(
            setup$y,
            solution$fitted_values,
            solution$prior_weights,
            edf
        )
    } else solution$scale
    covariance <- inverse_system * reported_scale
    dimnames(covariance) <- list(names(coefficients), names(coefficients))
    term_metadata <- lapply(seq_along(design$terms), function(i) {
        term <- design$terms[[i]]
        smooth_index <- which(vapply(
            setup$smooth,
            function(smooth) paste(smooth$term, collapse=',') ==
                paste0('cdr_term_', i),
            logical(1)
        ))
        coefficient_index <- if (length(smooth_index) == 1L) {
            seq.int(
                setup$smooth[[smooth_index]]$first.para,
                setup$smooth[[smooth_index]]$last.para
            )
        } else integer()
        list(
            name=term$name,
            type=term$type,
            knots=term$knots,
            predictor_knots=term$predictor_knots,
            axis=term$axis,
            linear_predictors=term$linear_predictors,
            linear_predictor_summaries=term$linear_predictor_summaries,
            basis=term$basis,
            transform=term$transform,
            group=term$group,
            group_levels=term$group_levels,
            base_dimension=term$base_dimension,
            rank=term$rank,
            null.space.dim=term$null.space.dim,
            S.scale=term$S.scale,
            lag_scale=if (is.null(term$lag_scale)) 1 else term$lag_scale,
            predictor_scale=if (is.null(term$predictor_scale)) 1 else
                term$predictor_scale,
            amplitude_scale=if (is.null(term$amplitude_scale)) 1 else
                term$amplitude_scale,
            coefficient_index=coefficient_index
        )
    })
    identifiability <- design$identifiability
    identifiability$global <- list(
        dimension=ncol(setup$X),
        rank=ncol(setup$X),
        condition_indicator=tryCatch(rcond(solution$system),
            error=function(error) NA_real_),
        parametric=alias_resolution$info,
        resolution=rank_resolution$resolution
    )
    output <- list(
        coefficients=coefficients,
        fitted.values=solution$fitted_values,
        residuals=solution$residuals,
        linear.predictors=solution$linear_predictors,
        family=family,
        sp=stats::setNames(solution$sp, names(setup$sp)),
        scale=reported_scale,
        sig2=reported_scale,
        reml.scale=solution$scale,
        method='REML',
        smooth=setup$smooth,
        paraPen=setup$paraPen,
        full.sp=setup$full.sp,
        outer.info=NULL,
        Vp=covariance,
        Vc=NULL,
        edf=edf,
        coefficient.edf=influence,
        df.residual=nrow(setup$X) - edf,
        y=setup$y,
        prior.weights=solution$prior_weights,
        working.weights=solution$working_weights,
        X=setup$X,
        offset=setup$offset,
        deviance=solution$deviance,
        reml=solution$criterion,
        optimizer=result$optimization,
        converged=isTRUE(solution$converged) &&
            identical(result$optimization$convergence, 0L),
        cdrgam=list(
            schema_version=1L,
            engine='block',
            backend='block',
            formula=list(
                user=design$formula,
                normalized=design$normalized_formula,
                effective=design$effective_formula,
                mgcv=setup$formula
            ),
            preparation=list(
                configuration=design$configuration,
                plan=design$plan,
                stream=design$stream,
                specification=design$specification,
                simplifications=design$simplifications,
                scaling=design$scaling,
                identifiability=design$identifiability
            ),
            scaling=design$scaling,
            identifiability=identifiability,
            term_labels=names(design$terms),
            terms=term_metadata,
            prediction=list(setup=.cdr_prediction_setup(setup)),
            rank=list(
                action=rank_action,
                parametric=alias_resolution$info,
                resolution=rank_resolution$resolution,
                regularization=fixed_ridge,
                tolerance=tolerance
            ),
            solver='dense generalized LAML reference solver'
        )
    )
    class(output) <- c('cdrgam_block', 'cdrgam')
    reporter$emit(
        1L,
        'fit complete',
        criterion=format(solution$criterion, digits=10),
        evaluations=result$evaluations,
        converged=output$converged
    )
    output
}

#' @export
coef.cdrgam_block <- function(object, ...) object$coefficients

#' @export
fitted.cdrgam_block <- function(object, ...) object$fitted.values

#' @export
residuals.cdrgam_block <- function(object, ...) object$residuals

#' @export
deviance.cdrgam_block <- function(object, ...) {
    if (!is.null(object$deviance)) return(object$deviance)
    sum(object$prior.weights * object$residuals^2)
}

#' @export
nobs.cdrgam_block <- function(object, ...) length(object$y)

#' @export
logLik.cdrgam_block <- function(object, ...) {
    n <- length(object$y)
    gaussian <- identical(object$family$family, 'gaussian') &&
        identical(object$family$link, 'identity')
    value <- if (gaussian) {
        -0.5 * (
            n * log(2 * pi * object$scale) +
                deviance.cdrgam_block(object) / object$scale
        )
    } else if (identical(object$family$family, 'Gamma')) {
        .cdrgam_gamma_loglik(
            object$y,
            object$fitted.values,
            object$prior.weights,
            object$scale
        )
    } else {
        -0.5 * object$family$aic(
            object$y,
            rep.int(1, n),
            object$fitted.values,
            object$prior.weights,
            deviance.cdrgam_block(object)
        )
    }
    estimated_dispersion <- gaussian ||
        identical(object$family$family, 'Gamma')
    attr(value, 'df') <- object$edf + as.integer(estimated_dispersion)
    attr(value, 'nobs') <- n
    class(value) <- 'logLik'
    value
}

#' @export
vcov.cdrgam_block <- function(object, unconditional=FALSE, ...) {
    if (isTRUE(unconditional)) {
        if (is.null(object$Vc)) {
            stop('Smoothing-parameter covariance is unavailable')
        }
        return(object$Vc)
    }
    object$Vp
}

#' @rdname predict.cdrgam
#' @export
predict.cdrgam_block <- function(object, newdata=NULL, ...) {
    if (is.null(newdata)) return(object$fitted.values)
    if (!is.list(newdata) ||
            !all(c('impulses', 'responses') %in% names(newdata))) {
        stop('newdata must contain impulse and response data frames')
    }
    .predict_cdrgam_streams(
        object,
        impulses=newdata$impulses,
        responses=newdata$responses,
        ...
    )
}

#' @export
print.cdrgam_block <- function(x, ...) {
    cat('Continuous-time deconvolutional GAM\n')
    cat('  backend:', x$cdrgam$solver, '\n')
    cat('  IRF terms:', paste(x$cdrgam$term_labels, collapse=', '), '\n')
    cat('  coefficients:', length(x$coefficients), '\n')
    cat('  REML criterion:', format(x$reml, digits=7), '\n')
    invisible(x)
}

#' @export
summary.cdrgam_block <- function(
        object,
        dispersion=NULL,
        freq=FALSE,
        re.test=TRUE,
        all.coefficients=FALSE,
        ...
) {
    out <- .cdrgam_custom_summary(
        object,
        covariance=function(indices) object$Vp[indices, indices, drop=FALSE],
        smooth_edf=.cdrgam_block_smooth_edf(object),
        dispersion=dispersion,
        all.coefficients=all.coefficients,
        all_variances=function() diag(object$Vp)
    )
    class(out) <- c('summary.cdrgam_block', 'summary.cdrgam', 'summary.gam')
    out
}

#' @export
print.summary.cdrgam_block <- function(x, ...) {
    .print_summary_cdrgam(x, ...)
    cat('Backend:', x$backend, '\n')
    if (isTRUE(x$rank_metadata$parametric$corrected)) {
        cat(
            'Aliased parametric coefficients removed:',
            paste(x$rank_metadata$parametric$dropped, collapse=', '), '\n'
        )
    }
    if (!is.null(x$rank_metadata$resolution) &&
            !identical(x$rank_metadata$resolution, 'none')) {
        cat('Rank resolution:', x$rank_metadata$resolution, '\n')
    }
    invisible(x)
}
