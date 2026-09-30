#' Report standardized fit convergence diagnostics
#'
#' @param object A fitted `cdrgam` model.
#' @return A named list containing at least `converged`, `code`, and `message`.
#'   Custom backends also report available objective counts, gradient norms,
#'   boundary information, and Hessian diagnostics.
#' @export
fit_diagnostics <- function(object) {
    if (!is_cdrgam(object)) stop('object must be a fitted cdrgam model')
    if (inherits(object, 'cdrgam_distributional_sparse')) {
        diagnostic_gradient <- if (length(
                object$optimizer$projected_gradient
            )) object$optimizer$projected_gradient else {
            object$optimizer$gradient
        }
        gradient_norm <- if (length(diagnostic_gradient)) {
            max(abs(diagnostic_gradient))
        } else NA_real_
        return(list(
            converged=isTRUE(object$converged),
            code=object$optimizer$convergence,
            message=if (isTRUE(object$converged)) {
                'full convergence'
            } else object$optimizer$message,
            total_objective_evaluations=object$distributional$evaluations,
            gradient_norm=gradient_norm,
            inner_gradient_norm=object$distributional$gradient_norm
        ))
    }
    if (inherits(object, 'cdrgam_sparse')) {
        return(object$sparse$convergence)
    }
    if (inherits(object, 'cdrgam_block')) {
        converged <- identical(object$optimizer$convergence, 0L)
        return(list(
            converged=converged,
            code=object$optimizer$convergence,
            message=if (converged) {
                'full convergence'
            } else if (is.null(object$optimizer$message)) {
                'optimizer convergence was not reached'
            } else object$optimizer$message
        ))
    }
    outer_message <- object$outer.info$conv
    converged <- if (!is.null(object$converged)) {
        isTRUE(object$converged)
    } else if (!is.null(outer_message)) {
        tolower(outer_message) %in% c('full convergence', 'converged')
    } else TRUE
    list(
        converged=converged,
        code=if (converged) 0L else NA_integer_,
        message=if (is.null(outer_message)) {
            'native mgcv convergence state'
        } else outer_message
    )
}

.cdrgam_simplification_sp_indices <- function(object) {
    cursor <- 0L
    lapply(object$smooth, function(smooth) {
        if (!is.null(smooth$first.sp) && !is.null(smooth$last.sp) &&
                is.finite(smooth$first.sp) && is.finite(smooth$last.sp)) {
            indices <- seq.int(smooth$first.sp, smooth$last.sp)
            cursor <<- max(cursor, indices)
            return(indices)
        }
        count <- length(smooth$sp)
        if (!count) return(integer())
        indices <- seq.int(cursor + 1L, cursor + count)
        cursor <<- cursor + count
        indices
    })
}

.cdrgam_simplification_edf <- function(object) {
    if (inherits(object, 'cdrgam_sparse')) {
        return(.cdrgam_sparse_smooth_edf(object))
    }
    if (inherits(object, 'cdrgam_block')) {
        return(.cdrgam_block_smooth_edf(object))
    }
    table <- summary(object)$s.table
    if (is.null(table)) numeric() else as.numeric(table[, 'edf'])
}

.cdrgam_simplification_gradient <- function(object) {
    gradient <- if (is.list(object$optimizer)) object$optimizer$gradient else NULL
    if (is.null(gradient)) gradient <- object$outer.info$grad
    if (!is.numeric(gradient) || length(gradient) != length(object$sp)) {
        return(rep.int(NA_real_, length(object$sp)))
    }
    as.numeric(gradient)
}

.cdrgam_simplification_tolerance <- function(object) {
    value <- object$sparse$control$optimizer_gradient_tolerance
    if (is.numeric(value) && length(value) == 1L && is.finite(value) &&
            value > 0) return(value)
    if (inherits(object, 'cdrgam_distributional_sparse')) 2e-3 else 1e-4
}

.cdrgam_simplification_term_info <- function(object, smooth) {
    indices <- seq.int(smooth$first.para, smooth$last.para)
    metadata <- object$cdrgam$terms
    if (!length(metadata)) return(NULL)
    matched <- which(vapply(metadata, function(term) {
        identical(as.integer(term$coefficient_index), as.integer(indices))
    }, logical(1)))
    if (length(matched) == 1L) metadata[[matched]] else NULL
}

.cdrgam_simplification_random_effect <- function(object, smooth) {
    if (inherits(smooth, 'random.effect') ||
            any(grepl('random.effect', class(smooth), fixed=TRUE))) return(TRUE)
    blocks <- object$cdrgam$prediction$random_effects
    if (!length(blocks)) return(FALSE)
    direct <- vapply(blocks, function(block) {
        is.list(block) && !is.null(block$coefficient_index)
    }, logical(1))
    if (!all(direct)) blocks <- unlist(blocks, recursive=FALSE)
    indices <- seq.int(smooth$first.para, smooth$last.para)
    any(vapply(blocks, function(block) {
        is.list(block) && identical(
            as.integer(block$coefficient_index), as.integer(indices)
        )
    }, logical(1)))
}

.cdrgam_simplification_recommendation <- function(
        object, label, smooth, metadata, collapsed, null_collapsed
) {
    random_effect <- .cdrgam_simplification_random_effect(object, smooth)
    grouped <- !is.null(metadata$group) || grepl('|', label, fixed=TRUE)
    type <- if (is.null(metadata$type)) '' else metadata$type
    if (isTRUE(collapsed)) {
        if (random_effect) return('remove the random effect')
        if (grouped) return('remove the grouped deviation')
        return('remove the term')
    }
    if (isTRUE(null_collapsed)) {
        return('replace the smooth with its unpenalized null-space form')
    }
    if (random_effect) return('remove the random effect')
    if (grouped) return('remove the grouped deviation before its population term')
    if (grepl('nonlinear', type, fixed=TRUE)) {
        return('make the predictor contribution linear or remove the interaction')
    }
    if (grepl('varying', type, fixed=TRUE)) {
        return('remove the time-varying component')
    }
    'reduce the basis dimension or remove the term'
}

.cdrgam_formula_text <- function(formula) {
    paste(deparse(formula, width.cutoff=500L), collapse=' ')
}

.cdrgam_simplification_spec_name <- function(spec) {
    name <- if (isTRUE(spec$constant)) 'irf(1)' else
        paste(spec$predictors, collapse=':')
    if (!is.null(spec$time)) name <- paste0(name, '~', spec$time)
    if (!is.null(spec$by)) name <- paste0(name, ':', spec$by)
    if (!is.null(spec$group)) name <- paste0(name, '|', spec$group)
    name
}

.cdrgam_drop_ordinary_smooth <- function(formula, smooth) {
    parsed <- .parse_cdr_formula(formula)
    terms_object <- stats::terms(parsed$ordinary_formula, keep.order=TRUE)
    labels <- attr(terms_object, 'term.labels')
    target_variables <- sort(as.character(smooth$term))
    matches <- which(vapply(labels, function(label) {
        expression <- tryCatch(str2lang(label), error=function(error) NULL)
        is.call(expression) && identical(as.character(expression[[1L]]), 's') &&
            identical(sort(all.vars(expression[[2L]])), target_variables)
    }, logical(1)))
    if (length(matches) != 1L) return(NULL)
    ordinary <- stats::formula(stats::drop.terms(
        terms_object, dropx=matches, keep.response=TRUE
    ))
    environment(ordinary) <- environment(formula)
    .compose_cdr_formula(ordinary, parsed$irfs)
}

.cdrgam_drop_irf <- function(formula, metadata) {
    parsed <- .parse_cdr_formula(formula)
    names <- vapply(
        parsed$irfs, .cdrgam_simplification_spec_name, character(1)
    )
    matches <- which(names == metadata$name)
    if (length(matches) != 1L) return(NULL)
    .compose_cdr_formula(parsed$ordinary_formula, parsed$irfs[-matches])
}

.cdrgam_simplification_patch <- function(
        object, parameter, smooth=NULL, metadata=NULL, action
) {
    formulas <- stats::formula(object, type='effective')
    distributional <- is.list(formulas)
    if (!distributional) formulas <- list(response=formulas)
    target <- if (distributional) parameter else names(formulas)[[1L]]
    if (is.na(target) || !(target %in% names(formulas))) return(NULL)
    before_formula <- formulas[[target]]
    after_formula <- if (identical(action, 'intercept_only_parameter')) {
        response <- if (length(before_formula) == 3L) {
            paste(deparse(before_formula[[2L]]), collapse='')
        } else NULL
        stats::as.formula(
            if (is.null(response)) '~ 1' else paste(response, '~ 1'),
            env=environment(before_formula)
        )
    } else if (!is.null(metadata)) {
        .cdrgam_drop_irf(before_formula, metadata)
    } else {
        .cdrgam_drop_ordinary_smooth(before_formula, smooth)
    }
    if (is.null(after_formula)) return(NULL)
    list(
        operation='replace_formula',
        path=if (distributional) paste0('formula.', target) else 'formula',
        parameter=if (distributional) target else NULL,
        before=.cdrgam_formula_text(before_formula),
        after=.cdrgam_formula_text(after_formula)
    )
}

.cdrgam_empty_simplification_candidates <- function() {
    data.frame(
        rank=integer(), priority=character(), parameter=character(),
        term=character(), action=character(), recommendation=character(),
        evidence=character(), automatable=logical(), patch_id=character(),
        edf=numeric(), coefficient_count=integer(),
        null_space_dimension=numeric(), smoothing_parameter=numeric(),
        outer_score=numeric(), evidence_score=numeric(), df_loss=numeric(),
        coefficient_loss=integer(), terms_touched=integer(),
        impact_fraction=numeric(), score=numeric(), stringsAsFactors=FALSE
    )
}

#' Suggest model simplifications from a fitted CDR-GAM
#'
#' The report ranks terms whose fitted effective degrees of freedom have
#' collapsed toward their penalized limit or whose unresolved outer derivative
#' favors stronger smoothing. For a nonconverged distributional fit, it also
#' identifies a parameter submodel when the simplification-directed score is
#' concentrated there. The function does not change or refit the model.
#'
#' @param object A fitted `cdrgam` model.
#' @param max_candidates Maximum number of candidates to return. Use `Inf` to
#'   retain every candidate.
#' @param conservatism Strength of the ranking penalty for model degrees of
#'   freedom removed by a candidate. Zero ranks only by diagnostic evidence.
#' @return A `cdrgam_simplification_report`. Use `as.data.frame()` to extract
#'   its ranked candidate table. The `patches` member contains exact,
#'   preconditioned formula replacements for candidates that are safe to
#'   automate mechanically; selecting and refitting them remains the caller's
#'   responsibility.
#' @export
suggest_simplifications <- function(object, max_candidates=10L, conservatism=1) {
    if (!is_cdrgam(object)) stop('object must be a fitted cdrgam model')
    if (length(max_candidates) != 1L || !is.numeric(max_candidates) ||
            is.na(max_candidates) || max_candidates < 1 ||
            !(is.infinite(max_candidates) || max_candidates == as.integer(max_candidates))) {
        stop('max_candidates must be a positive integer or Inf')
    }
    if (!is.numeric(conservatism) || length(conservatism) != 1L ||
            !is.finite(conservatism) || conservatism < 0) {
        stop('conservatism must be one finite nonnegative number')
    }
    labels <- .cdrgam_smooth_labels(object)
    if (!length(labels)) {
        output <- list(
            candidates=.cdrgam_empty_simplification_candidates(),
            patches=list(),
            converged=isTRUE(fit_diagnostics(object)$converged),
            gradient_tolerance=.cdrgam_simplification_tolerance(object),
            conservatism=conservatism
        )
        class(output) <- 'cdrgam_simplification_report'
        return(output)
    }
    edf <- .cdrgam_simplification_edf(object)
    if (length(edf) != length(labels)) {
        stop('Stored smooth metadata and effective degrees of freedom disagree')
    }
    sp_indices <- .cdrgam_simplification_sp_indices(object)
    gradient <- .cdrgam_simplification_gradient(object)
    tolerance <- .cdrgam_simplification_tolerance(object)
    converged <- isTRUE(fit_diagnostics(object)$converged)
    rows <- list()
    patches <- list()
    parameter_pressure <- numeric()
    total_edf <- max(1, sum(edf, na.rm=TRUE))
    for (i in seq_along(labels)) {
        smooth <- object$smooth[[i]]
        indices <- sp_indices[[i]]
        dimension <- smooth$last.para - smooth$first.para + 1L
        null_dimension <- if (is.null(smooth$null.space.dim)) {
            0
        } else smooth$null.space.dim
        penalized_dimension <- max(0, dimension - null_dimension)
        excess_edf <- max(0, edf[[i]] - null_dimension)
        collapsed <- null_dimension <= 0 &&
            edf[[i]] <= max(0.1, 0.02 * dimension)
        null_collapsed <- null_dimension > 0 && penalized_dimension > 0 &&
            excess_edf <= max(0.1, 0.02 * penalized_dimension)
        term_gradient <- gradient[indices]
        outer_score <- if (length(term_gradient) && any(is.finite(term_gradient))) {
            min(term_gradient, na.rm=TRUE)
        } else NA_real_
        pressure <- if (is.finite(outer_score)) max(0, -outer_score) else 0
        total_term_pressure <- if (length(term_gradient)) {
            sum(pmax(0, -term_gradient), na.rm=TRUE)
        } else 0
        parameter <- if (grepl(':', labels[[i]], fixed=TRUE) &&
                isTRUE(object$cdrgam$distributional)) {
            sub(':.*$', '', labels[[i]])
        } else NA_character_
        if (!is.na(parameter)) {
            existing_pressure <- if (parameter %in% names(parameter_pressure)) {
                parameter_pressure[[parameter]]
            } else 0
            parameter_pressure[[parameter]] <-
                existing_pressure + total_term_pressure
        }
        pressure_ratio <- pressure / tolerance
        if (!collapsed && !null_collapsed &&
                (converged || pressure_ratio < 2)) next
        metadata <- .cdrgam_simplification_term_info(object, smooth)
        evidence <- character()
        if (collapsed) evidence <- c(evidence, sprintf(
            'edf %.3g of %d coefficients', edf[[i]], dimension
        ))
        if (null_collapsed) evidence <- c(evidence, sprintf(
            'edf %.3g is near null-space dimension %.3g',
            edf[[i]], null_dimension
        ))
        if (pressure_ratio >= 2) evidence <- c(evidence, sprintf(
            'outer derivative %.3g favors stronger smoothing (%.1fx tolerance)',
            outer_score, pressure_ratio
        ))
        candidate_score <- max(
            if (collapsed) 100 * (1 - min(1, edf[[i]] /
                max(1, dimension))) else 0,
            if (null_collapsed) 80 * (1 - min(1, excess_edf /
                max(1, penalized_dimension))) else 0,
            if (pressure_ratio >= 2) min(90, 45 + 15 * log10(pressure_ratio)) else 0
        )
        smoothing_parameter <- if (length(indices)) {
            max(as.numeric(object$sp[indices]), na.rm=TRUE)
        } else NA_real_
        random_effect <- .cdrgam_simplification_random_effect(object, smooth)
        grouped <- !is.null(metadata$group) || grepl('|', labels[[i]], fixed=TRUE)
        action <- if (collapsed && random_effect) {
            'drop_random_effect'
        } else if (grouped) {
            'drop_grouped_deviation'
        } else if (collapsed) {
            'drop_term'
        } else if (null_collapsed) {
            'replace_with_null_space'
        } else 'review_term'
        automatable <- action %in% c(
            'drop_random_effect', 'drop_grouped_deviation', 'drop_term'
        ) && !identical(metadata$name, 'irf(1)')
        patch <- if (automatable) .cdrgam_simplification_patch(
            object, parameter, smooth, metadata, action
        ) else NULL
        automatable <- !is.null(patch)
        patch_id <- if (automatable) paste0('patch-', length(patches) + 1L) else NA_character_
        if (automatable) patches[[patch_id]] <- patch
        df_loss <- if (null_collapsed) excess_edf else edf[[i]]
        coefficient_loss <- if (null_collapsed) penalized_dimension else dimension
        impact_fraction <- min(1, max(0, df_loss / total_edf))
        adjusted_score <- candidate_score - conservatism *
            35 * sqrt(impact_fraction)
        rows[[length(rows) + 1L]] <- data.frame(
            rank=NA_integer_, priority=if (adjusted_score >= 70) 'high' else
                if (adjusted_score >= 50) 'moderate' else 'low',
            parameter=parameter,
            term=labels[[i]],
            action=action,
            recommendation=.cdrgam_simplification_recommendation(
                object, labels[[i]], smooth, metadata, collapsed, null_collapsed
            ),
            evidence=paste(evidence, collapse='; '),
            automatable=automatable, patch_id=patch_id,
            edf=edf[[i]], coefficient_count=as.integer(dimension),
            null_space_dimension=null_dimension,
            smoothing_parameter=smoothing_parameter,
            outer_score=outer_score, evidence_score=candidate_score,
            df_loss=df_loss, coefficient_loss=as.integer(coefficient_loss),
            terms_touched=1L, impact_fraction=impact_fraction,
            score=adjusted_score,
            stringsAsFactors=FALSE
        )
    }
    total_pressure <- sum(parameter_pressure)
    if (!converged && length(parameter_pressure) > 1L && total_pressure > 0) {
        dominant <- names(which.max(parameter_pressure))
        share <- unname(parameter_pressure[[dominant]] / total_pressure)
        if (share >= 0.67 && parameter_pressure[[dominant]] / tolerance >= 2) {
            evidence_score <- min(95, 70 + 25 * share)
            parameter_rows <- which(vapply(seq_along(labels), function(index) {
                grepl(paste0('^', dominant, ':'), labels[[index]])
            }, logical(1)))
            df_loss <- sum(edf[parameter_rows], na.rm=TRUE)
            coefficient_loss <- sum(vapply(parameter_rows, function(index) {
                object$smooth[[index]]$last.para -
                    object$smooth[[index]]$first.para + 1L
            }, integer(1)))
            impact_fraction <- min(1, max(0, df_loss / total_edf))
            score <- evidence_score - conservatism *
                (35 * sqrt(impact_fraction) + 5 * log1p(length(parameter_rows)))
            patch <- .cdrgam_simplification_patch(
                object, dominant, action='intercept_only_parameter'
            )
            automatable <- !is.null(patch)
            patch_id <- if (automatable) paste0('patch-', length(patches) + 1L) else NA_character_
            if (automatable) patches[[patch_id]] <- patch
            rows[[length(rows) + 1L]] <- data.frame(
                rank=NA_integer_, priority=if (score >= 70) 'high' else
                    if (score >= 50) 'moderate' else 'low', parameter=dominant,
                term=paste0('<', dominant, ' submodel>'),
                action='intercept_only_parameter',
                recommendation=paste0(
                    'start with an intercept-only ', dominant,
                    ' submodel, then restore terms incrementally'
                ),
                evidence=sprintf(
                    '%.0f%% of simplification-directed outer score is in this submodel',
                    100 * share
                ),
                automatable=automatable, patch_id=patch_id,
                edf=NA_real_, coefficient_count=NA_integer_,
                null_space_dimension=NA_real_, smoothing_parameter=NA_real_,
                outer_score=-parameter_pressure[[dominant]],
                evidence_score=evidence_score, df_loss=df_loss,
                coefficient_loss=as.integer(coefficient_loss),
                terms_touched=length(parameter_rows),
                impact_fraction=impact_fraction, score=score,
                stringsAsFactors=FALSE
            )
        }
    }
    candidates <- if (length(rows)) do.call(rbind, rows) else
        .cdrgam_empty_simplification_candidates()
    if (nrow(candidates)) {
        candidates <- candidates[order(
            -candidates$score, candidates$parameter, candidates$term,
            na.last=TRUE
        ), , drop=FALSE]
        if (is.finite(max_candidates)) {
            candidates <- utils::head(candidates, as.integer(max_candidates))
        }
        candidates$rank <- seq_len(nrow(candidates))
        rownames(candidates) <- NULL
    }
    output <- list(
        candidates=candidates,
        patches=patches,
        converged=converged,
        gradient_tolerance=tolerance,
        conservatism=conservatism
    )
    class(output) <- 'cdrgam_simplification_report'
    output
}

#' @export
as.data.frame.cdrgam_simplification_report <- function(x, ...) x$candidates

#' @export
print.cdrgam_simplification_report <- function(x, ...) {
    cat('CDR-GAM simplification guidance\n')
    cat('  fit converged:', if (isTRUE(x$converged)) 'yes' else 'no', '\n')
    if (!nrow(x$candidates)) {
        cat('  no strong simplification candidates identified\n')
        return(invisible(x))
    }
    cat('  candidates:', nrow(x$candidates), '\n\n')
    for (i in seq_len(nrow(x$candidates))) {
        candidate <- x$candidates[i, ]
        parameter <- if (is.na(candidate$parameter)) '' else
            paste0(' [', candidate$parameter, ']')
        cat(candidate$rank, '. ', candidate$term, parameter, ' (',
            candidate$priority, ')\n', sep='')
        cat('   Consider: ', candidate$recommendation, '.\n', sep='')
        cat('   Evidence: ', candidate$evidence, '.\n', sep='')
        cat('   Estimated df removed: ', format(candidate$df_loss, digits=3),
            '; adjusted score: ', format(candidate$score, digits=3), '.\n', sep='')
    }
    cat('\nRecommendations are diagnostics, not automatic model selection.\n')
    invisible(x)
}

#' Evaluate fitted impulse-response function terms
#'
#' @param object A fitted `cdrgam` model.
#' @param term Term name or one-based term number. Omit to evaluate every IRF.
#' @param lag Optional numeric lag grid. By default each term is evaluated over
#'   its knot range.
#' @param n Number of default grid points.
#' @param predictor Optional predictor-value grid for nonlinear IRFs.
#' @param n_predictor Number of default predictor grid points.
#' @param at Named list of evaluation values for smooth response-time and
#'   predictor axes. `predictor` remains an alias for the first non-lag axis;
#'   additional axes default to their fitted-grid medians.
#' @param group Optional grouping levels for grouped IRF deviations. By default
#'   all fitted levels are returned.
#' @param se Include pointwise standard errors when covariance is available.
#' @param unconditional Include smoothing-parameter uncertainty in standard
#'   errors when supported by the fitted backend.
#' @return A data frame with term, lag, estimate, and standard-error columns.
#' @export
estimate_irf <- function(
        object,
        term=NULL,
        lag=NULL,
        n=200,
        predictor=NULL,
        n_predictor=25,
        at=list(),
        group=NULL,
        se=TRUE,
        unconditional=FALSE
) {
    if (!is_cdrgam(object)) {
        stop('object must be a fitted cdrgam model')
    }
    metadata <- object$cdrgam$terms
    labels <- object$cdrgam$term_labels
    if (!length(metadata)) {
        stop('The fitted model does not contain stored IRF metadata')
    }
    indices <- if (is.null(term)) {
        seq_along(metadata)
    } else if (is.numeric(term) && length(term) == 1L) {
        as.integer(term)
    } else if (is.character(term) && length(term) == 1L) {
        which(labels == term)
    } else {
        integer()
    }
    if (!length(indices) || any(indices < 1L | indices > length(metadata))) {
        stop('term must uniquely identify a fitted IRF')
    }
    coefficients <- stats::coef(object)
    sparse_covariance <- isTRUE(se) && inherits(object, 'cdrgam_sparse')
    covariance <- if (isTRUE(se) && !sparse_covariance) {
        stats::vcov(object, unconditional=unconditional)
    } else {
        NULL
    }
    output <- list()
    for (j in seq_along(indices)) {
        i <- indices[[j]]
        info <- metadata[[i]]
        lag_scale <- if (is.null(info$lag_scale)) 1 else info$lag_scale
        predictor_scale <- if (is.null(info$predictor_scale)) {
            1
        } else info$predictor_scale
        amplitude_scale <- if (is.null(info$amplitude_scale)) {
            1
        } else info$amplitude_scale
        lag_grid <- lag
        if (is.null(lag_grid)) {
            lag_grid <- seq(
                min(info$knots) * lag_scale,
                max(info$knots) * lag_scale,
                length.out=n
            )
        }
        if (!is.numeric(lag_grid) || any(!is.finite(lag_grid))) {
            stop('lag must be a finite numeric vector')
        }
        is_surface <- (!is.null(info$axis) && length(info$axis) > 1L) ||
            startsWith(info$type, 'nonlinear') ||
            startsWith(info$type, 'varying')
        if (is_surface) {
            if (!is.null(info$axis)) {
                grid_values <- list(lag=lag_grid)
                for (axis_index in seq.int(2L, length(info$axis))) {
                    axis <- info$axis[[axis_index]]
                    requested <- at[[axis$variable]]
                    if (is.null(requested) && axis_index == 2L) {
                        requested <- predictor
                    }
                    if (is.null(requested)) {
                        requested <- if (axis_index == 2L) {
                            seq(
                                min(axis$grid), max(axis$grid),
                                length.out=n_predictor
                            ) * axis$scale
                        } else stats::median(axis$grid) * axis$scale
                    }
                    if (!is.numeric(requested) || any(!is.finite(requested))) {
                        stop('IRF axis values must be finite numeric vectors')
                    }
                    grid_values[[axis$variable]] <- requested
                }
                evaluation_grid <- do.call(expand.grid, c(
                    grid_values,
                    list(KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
                ))
                axis_data <- stats::setNames(lapply(seq_along(info$axis), function(index) {
                    axis <- info$axis[[index]]
                    value <- if (index == 1L) evaluation_grid$lag else
                        evaluation_grid[[axis$variable]]
                    value / axis$scale
                }), vapply(info$axis, `[[`, character(1), 'internal'))
                basis <- mgcv::PredictMat(
                    info$basis,
                    axis_data,
                    n=nrow(evaluation_grid)
                ) %*% info$transform
                evaluation_grid$predictor <- evaluation_grid[[
                    info$axis[[2L]]$variable
                ]]
            } else {
                predictor_grid <- predictor
                if (is.null(predictor_grid)) {
                    predictor_grid <- seq(
                        min(info$predictor_knots) * predictor_scale,
                        max(info$predictor_knots) * predictor_scale,
                        length.out=n_predictor
                    )
                }
                if (!is.numeric(predictor_grid) || any(!is.finite(predictor_grid))) {
                    stop('predictor must be a finite numeric vector')
                }
                evaluation_grid <- expand.grid(
                    lag=lag_grid,
                    predictor=predictor_grid
                )
                basis <- mgcv::PredictMat(
                    info$basis,
                    list(
                        cdr_delay=evaluation_grid$lag / lag_scale,
                        cdr_value=evaluation_grid$predictor / predictor_scale
                    ),
                    n=nrow(evaluation_grid)
                ) %*% info$transform
            }
        } else {
            evaluation_grid <- data.frame(
                lag=lag_grid,
                predictor=NA_real_
            )
            basis <- mgcv::PredictMat(
                info$basis,
                list(cdr_delta=lag_grid / lag_scale),
                n=length(lag_grid)
            )
            if (!is.null(info$transform)) {
                basis <- basis %*% info$transform
            }
        }
        basis <- basis / amplitude_scale
        coefficient_index <- info$coefficient_index
        grouped <- length(info$group_levels) > 0L
        levels_to_evaluate <- if (grouped) {
            if (is.null(group)) info$group_levels else as.character(group)
        } else {
            NA_character_
        }
        if (grouped && any(!(levels_to_evaluate %in% info$group_levels))) {
            stop('Unknown grouped-IRF deviation level requested')
        }
        for (group_level in levels_to_evaluate) {
            term_index <- coefficient_index
            if (grouped) {
                group_number <- match(group_level, info$group_levels)
                within_term <- (group_number - 1L) * info$base_dimension +
                    seq_len(info$base_dimension)
                term_index <- coefficient_index[within_term]
            }
            if (length(term_index) != ncol(basis)) {
                stop('Stored IRF coefficient metadata is inconsistent with its basis')
            }
            estimate <- drop(basis %*% coefficients[term_index])
            standard_error <- rep.int(NA_real_, nrow(evaluation_grid))
            if (!is.null(covariance) || sparse_covariance) {
                term_covariance <- if (sparse_covariance) {
                    .sparse_selected_vcov(
                        object,
                        term_index,
                        unconditional=unconditional
                    )
                } else covariance[
                    term_index,
                    term_index,
                    drop=FALSE
                ]
                standard_error <- sqrt(pmax(
                    0,
                    rowSums((basis %*% term_covariance) * basis)
                ))
            }
            result <- data.frame(
                term=labels[[i]],
                group=if (grouped) group_level else NA_character_,
                evaluation_grid,
                estimate=estimate,
                se=standard_error,
                stringsAsFactors=FALSE,
                check.names=FALSE
            )
            leading <- c('term', 'group', 'lag', 'predictor')
            output[[length(output) + 1L]] <- result[c(
                leading[leading %in% names(result)],
                setdiff(names(result), c(leading, 'estimate', 'se')),
                'estimate', 'se'
            )]
        }
    }
    columns <- unique(unlist(lapply(output, names)))
    output <- lapply(output, function(value) {
        missing <- setdiff(columns, names(value))
        for (name in missing) value[[name]] <- NA
        value[columns]
    })
    do.call(rbind, output)
}

#' Extract smoothing variance components
#'
#' Report the variance components implied by the smoothing parameters of the
#' expanded `mgcv` model. Native fits are passed directly to
#' [mgcv::gam.vcomp()]. For block-backend fits, the same calculation is applied
#' to the stored REML scale, penalty rescaling metadata, and outer Hessian.
#' This function does not estimate covariance parameters beyond those present
#' in the expanded model.
#'
#' @param object A fitted `cdrgam` model.
#' @param rescale Apply the penalty rescaling used in the original smooth
#'   specification, as in [mgcv::gam.vcomp()].
#' @param conf.lev Confidence level for intervals when an outer REML Hessian is
#'   available.
#' @return The value returned by [mgcv::gam.vcomp()]: normally a matrix of
#'   standard deviations and confidence intervals, or a named vector/list
#'   when intervals are unavailable.
#' @export
variance_components <- function(object, rescale=TRUE, conf.lev=0.95) {
    if (!is_cdrgam(object)) {
        stop('object must be a fitted cdrgam model')
    }
    if (inherits(object, 'gam')) {
        return(mgcv::gam.vcomp(
            object,
            rescale=rescale,
            conf.lev=conf.lev
        ))
    }
    if (!inherits(object, c('cdrgam_block', 'cdrgam_sparse'))) {
        stop('Unsupported cdrgam fitting backend')
    }
    if (inherits(object, 'cdrgam_sparse') &&
            is.null(object$outer.info$hess) &&
            is.environment(object$sparse$deferred_hessian)) {
        object$outer.info$hess <- .sparse_resolve_outer_hessian(object)
    }

    # gam.vcomp() only depends on this subset of a fitted gam object. Using it
    # as the single implementation also keeps linked penalties, parametric
    # penalties, rescaling, and interval conventions aligned with mgcv.
    surrogate <- list(
        sp=object$sp,
        full.sp=object$full.sp,
        smooth=object$smooth,
        paraPen=object$paraPen,
        reml.scale=object$reml.scale,
        sig2=object$sig2,
        method=object$method,
        outer.info=object$outer.info,
        family=object$family
    )
    class(surrogate) <- 'gam'
    mgcv::gam.vcomp(
        surrogate,
        rescale=rescale,
        conf.lev=conf.lev
    )
}

#' Simulate irregular impulse and response streams from fixed IRFs
#'
#' Predictor columns are constructed to be sample-orthogonal and standardized.
#' Responses are the additive convolution of every supplied IRF plus Gaussian
#' noise.
#'
#' @param irfs Named list of functions. One-argument functions define linear
#'   lag IRFs; two-argument functions define surfaces over lag and predictor.
#' @param n_impulses Number of impulses.
#' @param n_responses Number of responses.
#' @param duration Duration of the simulated series.
#' @param window Maximum causal lag included in the convolution.
#' @param intercept Response intercept.
#' @param noise_sd Gaussian response noise standard deviation.
#' @param seed Random seed.
#' @return A list containing impulse and response data, IRFs, and noiseless
#'   response components.
#' @export
simulate_cdr <- function(
        irfs,
        n_impulses=500,
        n_responses=600,
        duration=60,
        window=2,
        intercept=0,
        noise_sd=1,
        seed=NULL
) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    if (!is.list(irfs) || !length(irfs) || is.null(names(irfs)) ||
            any(!nzchar(names(irfs))) ||
            any(!vapply(irfs, is.function, logical(1)))) {
        stop('irfs must be a named, non-empty list of functions')
    }
    p <- length(irfs)
    if (n_impulses <= p || n_responses < 1L || duration <= 0 ||
            window <= 0 || noise_sd < 0) {
        stop('Invalid simulation dimensions, duration, window, or noise_sd')
    }
    impulse_time <- sort(stats::runif(n_impulses, 0, duration))
    response_time <- sort(stats::runif(n_responses, 0, duration))
    raw <- matrix(stats::rnorm(n_impulses * p), n_impulses, p)
    raw <- scale(raw, center=TRUE, scale=FALSE)
    predictors <- qr.Q(qr(raw)) * sqrt(n_impulses)
    colnames(predictors) <- names(irfs)
    impulses <- data.frame(time=impulse_time, predictors, check.names=FALSE)
    responses <- data.frame(time=response_time)

    links <- .build_history_links(
        impulses,
        responses,
        series=character(),
        impulse_time='time',
        response_time='time',
        window=c(0, window)
    )
    components <- matrix(0, nrow=n_responses, ncol=p)
    colnames(components) <- names(irfs)
    nonlinear <- vapply(irfs, function(fun) length(formals(fun)) >= 2L, logical(1))
    for (i in seq_along(irfs)) {
        predictor_value <- predictors[links$impulse_index, i]
        contribution <- if (nonlinear[[i]]) {
            irfs[[i]](links$delay, predictor_value)
        } else {
            predictor_value * irfs[[i]](links$delay)
        }
        if (length(contribution)) {
            summed <- rowsum(
                contribution,
                links$response_index,
                reorder=FALSE
            )
            components[as.integer(rownames(summed)), i] <- summed[, 1L]
        }
    }
    noiseless <- intercept + rowSums(components)
    responses$response <- noiseless + stats::rnorm(n_responses, sd=noise_sd)
    list(
        impulses=impulses,
        responses=responses,
        irfs=irfs,
        nonlinear=stats::setNames(nonlinear, names(irfs)),
        components=components,
        noiseless=noiseless,
        settings=list(
            window=window,
            intercept=intercept,
            noise_sd=noise_sd,
            seed=seed
        )
    )
}
