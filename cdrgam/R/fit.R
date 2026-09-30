#' Fit a prepared CDR-GAM design
#'
#' `cdrgam.fit()` fits the reusable design returned by [prepare_cdrgam()]. It
#' does not parse formulas or compile stream histories.
#'
#' @rdname cdrgam
#' @param design A `cdrgam_design` returned by [prepare_cdrgam()].
#' @inheritParams cdrgam
#' @return A fitted `cdrgam` object.
#' @export
cdrgam.fit <- function(
        design,
        family=stats::gaussian(),
        method=NULL,
        engine=c('bam', 'gam'),
        backend=c('mgcv', 'block', 'sparse'),
        checkpoint=NULL,
        solver_trace=FALSE,
        sparse_control=list(),
        rank_action=c('error', 'minimum_norm', 'drop', 'penalize'),
        rank_tol=NULL,
        rank_penalty=NULL,
        drop.unused.levels=NULL,
        ...
) {
    call <- .cdrgam_compact_call(
        match.call(),
        'cdrgam.fit',
        'design'
    )
    if (!inherits(design, 'cdrgam_design')) {
        stop('design must be a cdrgam_design returned by prepare_cdrgam()')
    }
    family <- .as_family(family)
    if (inherits(design, 'cdrgam_distributional_design')) {
        if (missing(engine)) engine <- 'gam'
        prepared_drop <- design$configuration$drop.unused.levels
        if (!is.null(drop.unused.levels) &&
                !identical(drop.unused.levels, prepared_drop)) {
            stop(
                'drop.unused.levels is fixed by prepare_cdrgam(); ',
                'reprepare the design with the requested value'
            )
        }
        backend <- match.arg(backend)
        engine <- match.arg(engine)
        if (!identical(family$family, 'gaulss')) {
            stop('The initial distributional design requires family="gaulss"')
        }
        if (!identical(match.arg(rank_action), 'error') ||
                !is.null(rank_tol) || !is.null(rank_penalty)) {
            stop(
                'Rank controls do not yet apply to distributional fits'
            )
        }
        if (!is.null(checkpoint) && !identical(backend, 'sparse')) {
            stop('Distributional checkpoints require backend="sparse"')
        }
        fit <- if (identical(backend, 'sparse')) {
            .fit_distributional_sparse(
                design,
                family=family,
                method=method,
                checkpoint=checkpoint,
                trace=solver_trace,
                sparse_control=sparse_control,
                ...
            )
        } else if (identical(backend, 'block')) {
            .fit_distributional_block(
                design,
                family=family,
                method=method,
                trace=solver_trace,
                ...
            )
        } else {
            if (length(sparse_control)) {
                stop('sparse_control applies only to backend="sparse"')
            }
            .fit_distributional_mgcv(
                design, family=family, method=method, engine=engine, ...
            )
        }
        fit$call <- call
        fit$cdrgam$call <- call
        return(fit)
    }
    prepared_drop <- design$configuration$drop.unused.levels
    if (!is.null(drop.unused.levels)) {
        if (length(drop.unused.levels) != 1L ||
                !is.logical(drop.unused.levels) || is.na(drop.unused.levels)) {
            stop('drop.unused.levels must be TRUE or FALSE')
        }
        if (!identical(drop.unused.levels, prepared_drop)) {
            stop(
                'drop.unused.levels is fixed by prepare_cdrgam(); ',
                'reprepare the design with the requested value'
            )
        }
    }
    backend <- match.arg(backend)
    rank_action <- match.arg(rank_action)
    trace_control <- .solver_trace_level(solver_trace)
    .rank_tolerance(rank_tol)
    .rank_penalty(rank_penalty)
    if (identical(backend, 'block')) {
        if (length(sparse_control)) {
            stop('sparse_control applies only to backend="sparse"')
        }
        design <- .materialize_cdr_design(design, sparse=FALSE)
        fitting_function <- if (identical(family$family, 'gaussian') &&
                identical(family$link, 'identity')) {
            .fit_block_gaussian
        } else {
            .fit_block_generalized
        }
        fit <- fitting_function(
            design=design, family=family, method=method,
            checkpoint=checkpoint, trace=solver_trace,
            rank_action=rank_action, rank_tol=rank_tol,
            rank_penalty=rank_penalty,
            drop.unused.levels=design$configuration$drop.unused.levels, ...
        )
        fit$call <- call
        fit$cdrgam$call <- call
        return(fit)
    }
    if (identical(backend, 'sparse')) {
        generalized <- !identical(family$family, 'gaussian') ||
            !identical(family$link, 'identity')
        if (generalized) {
            fit <- .fit_sparse_generalized(
                design=design, family=family, method=method,
                checkpoint=checkpoint, trace=solver_trace,
                sparse_control=sparse_control, rank_action=rank_action,
                rank_tol=rank_tol, rank_penalty=rank_penalty,
                drop.unused.levels=design$configuration$drop.unused.levels,
                ...
            )
            fit$call <- call
            fit$cdrgam$call <- call
            return(fit)
        }
        fit <- .fit_sparse_gaussian(
            design=design,
            family=family,
            method=method,
            checkpoint=checkpoint,
            trace=solver_trace,
            sparse_control=sparse_control,
            rank_action=rank_action,
            rank_tol=rank_tol,
            rank_penalty=rank_penalty,
            drop.unused.levels=design$configuration$drop.unused.levels,
            ...
        )
        boundary_action_value <- sparse_control[['boundary_action', exact=TRUE]]
        boundary_action <- if (is.null(boundary_action_value)) {
            'report'
        } else boundary_action_value
        boundary_reduced <- FALSE
        if (identical(boundary_action, 'reduce')) {
            reduction_design <- design
            reduction_tables <- list()
            reduction_pass <- 0L
            boundary_reporter <- .new_solver_reporter(
                solver_trace,
                'sparse'
            )
            repeat {
                planned <- .cdr_boundary_reduction_plan(
                    fit,
                    reduction_design,
                    log_sp_threshold=if (is.null(
                        sparse_control[['boundary_log_sp', exact=TRUE]]
                    )) 12 else sparse_control[[
                        'boundary_log_sp', exact=TRUE
                    ]],
                    score_tolerance=if (is.null(
                        sparse_control[[
                            'optimizer_gradient_tolerance', exact=TRUE
                        ]]
                    )) 1e-4 else sparse_control[[
                        'optimizer_gradient_tolerance', exact=TRUE
                    ]],
                    edf_tolerance=if (is.null(
                        sparse_control[[
                            'boundary_edf_tolerance', exact=TRUE
                        ]]
                    )) 0.01 else sparse_control[[
                        'boundary_edf_tolerance', exact=TRUE
                    ]]
                )
                if (!length(planned$entries)) break
                reduction_pass <- reduction_pass + 1L
                if (reduction_pass > 1000L) {
                    stop('Boundary reduction did not reach a fixed point')
                }
                boundary_reporter$phase(
                    'boundary reduction',
                    pass=reduction_pass,
                    reason='penalty curvature was empirically unsupported',
                    pilot_criterion=format(fit$reml, digits=10),
                    reductions=length(planned$entries),
                    pilot_coefficients=length(fit$coefficients),
                    pilot_smoothing_parameters=length(fit$sp)
                )
                for (entry in planned$entries) {
                    row <- planned$table[
                        planned$table$term_index == entry$term_index &
                            planned$table$status ==
                                'certified_boundary_candidate',
                        ,
                        drop=FALSE
                    ]
                    boundary_reporter$emit(
                        1L,
                        'boundary reduction selected',
                        pass=reduction_pass,
                        term=names(reduction_design$terms)[[
                            entry$term_index
                        ]],
                        components=if (nrow(row)) row$components[[1L]] else '',
                        dimension=paste0(
                            entry$original_dimension,
                            '->',
                            entry$effective_dimension
                        ),
                        penalties=paste(
                            entry$smoothing_parameters,
                            collapse=','
                        )
                    )
                }
                reduced_design <- .cdr_apply_boundary_reduction(
                    reduction_design,
                    planned
                )
                reduced_initial <- .cdr_boundary_warm_start(
                    fit,
                    reduced_design
                )
                reduced_control <- sparse_control
                reduced_control$boundary_action <- 'report'
                boundary_reporter$emit(
                    1L,
                    'reduced model refit started',
                    pass=reduction_pass,
                    reason='infinite-penalty limits require a reduced basis',
                    initialization='mapped current smoothing parameters',
                    checkpoint='disabled because the model structure changed'
                )
                fit <- .fit_sparse_gaussian(
                    design=reduced_design,
                    family=family,
                    method=method,
                    checkpoint=NULL,
                    trace=solver_trace,
                    sparse_control=reduced_control,
                    rank_action=rank_action,
                    rank_tol=rank_tol,
                    rank_penalty=rank_penalty,
                    drop.unused.levels=
                        reduced_design$configuration$drop.unused.levels,
                    initial_log_sp=reduced_initial,
                    initial_source='mapped boundary-fit smoothing parameters',
                    ...
                )
                boundary_reporter$emit(
                    1L,
                    'reduced model refit complete',
                    pass=reduction_pass,
                    criterion=format(fit$reml, digits=10),
                    converged=isTRUE(fit$converged),
                    coefficients=length(fit$coefficients),
                    smoothing_parameters=length(fit$sp)
                )
                applied <- planned$table
                selected <- applied$status == 'certified_boundary_candidate'
                applied$status[selected] <- 'applied_and_reoptimized'
                applied <- applied[selected, , drop=FALSE]
                applied$fit_stage <- if (reduction_pass == 1L) {
                    'pilot'
                } else {
                    'reoptimized'
                }
                reduction_tables[[length(reduction_tables) + 1L]] <- applied
                reduction_design <- reduced_design
                boundary_reduced <- TRUE
            }
            if (boundary_reduced) {
                remaining <- fit$cdrgam$boundary_reductions
                if (!is.null(remaining) && nrow(remaining)) {
                    remaining$fit_stage <- 'reoptimized'
                    reduction_tables[[length(reduction_tables) + 1L]] <-
                        remaining
                }
                reductions <- do.call(rbind, reduction_tables)
                rownames(reductions) <- NULL
                fit$cdrgam$boundary_reductions <- reductions
                fit$cdrgam$conditional_on_boundary_reduction <- TRUE
                fit$sparse$convergence$boundary_reductions <- reductions
                fit$sparse$convergence$conditional_on_boundary_reduction <- TRUE
            }
        }
        if (identical(boundary_action, 'reduce') && !boundary_reduced &&
                !isTRUE(fit$converged)) {
            warning(
                'Sparse REML optimizer did not converge (code ',
                fit$sparse$convergence$code, '): ',
                fit$sparse$convergence$message,
                call.=FALSE
            )
        }
        if (identical(boundary_action, 'reduce') && !boundary_reduced &&
                identical(
                    fit$sparse$convergence$hessian_positive_definite,
                    FALSE
                )) {
            warning(
                'Sparse REML outer Hessian is not positive definite; ',
                'variance-component intervals may be unreliable',
                call.=FALSE
            )
        }
        fit$call <- call
        fit$cdrgam$call <- call
        return(fit)
    }
    if (!is.null(checkpoint) || trace_control$level > 0L ||
            !is.null(trace_control$callback) || length(sparse_control) ||
            !identical(rank_action, 'error') || !is.null(rank_tol) ||
            !is.null(rank_penalty)) {
        stop(
            'checkpoint, solver_trace, sparse_control, and rank controls ',
            'apply only to custom backends'
        )
    }
    design <- .materialize_cdr_design(design, sparse=FALSE)
    y <- design$responses[[design$response_name]]
    fit <- .fit_compressed_mgcv(
        y=y,
        terms=design$terms,
        family=family,
        method=method,
        engine=engine,
        base_formula=design$ordinary_formula,
        response_data=design$responses,
        user_formula=design$formula,
        preparation=list(
            configuration=design$configuration,
            plan=design$plan,
            stream=design$stream,
            specification=design$specification,
            simplifications=design$simplifications,
            scaling=design$scaling,
            identifiability=design$identifiability,
            normalized_formula=design$normalized_formula,
            effective_formula=design$effective_formula
        ),
        drop.unused.levels=design$configuration$drop.unused.levels,
        ...
    )
    fit$call <- call
    fit$cdrgam$call <- call
    fit
}

