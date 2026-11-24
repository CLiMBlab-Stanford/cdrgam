.cdr_effect_axes <- function(info) {
    if (!is.null(info$axis)) return(info$axis)
    scale <- if (is.null(info$lag_scale)) 1 else info$lag_scale
    list(list(
        role='lag', variable=NA_character_, label='lag', internal='cdr_delta',
        grid=info$knots, scale=scale,
        summary=list(
            mean=mean(info$knots) * scale,
            sd=stats::sd(info$knots) * abs(scale),
            quantiles=stats::setNames(
                as.numeric(stats::quantile(
                    info$knots * scale,
                    c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), names=FALSE
                )),
                c('0', '0.1', '0.25', '0.5', '0.75', '0.9', '1')
            )
        )
    ))
}

.cdr_effect_axis_name <- function(axis) {
    if (identical(axis$role, 'lag')) 'lag' else axis$variable
}

# The term specifications aligned with the fitted terms, or NULL when they
# cannot be aligned. Boundary reduction can remove whole terms after the
# specifications were recorded; fits from before the specifications were
# subset with them hold one specification per original term. Those are
# aligned in order by predictors and grouping factor, since the fitted terms
# keep their original relative order.
.cdr_effect_specifications <- function(object) {
    specifications <- object$cdrgam$preparation$specification
    terms <- object$cdrgam$terms
    if (is.null(specifications)) return(NULL)
    if (length(specifications) == length(terms)) return(specifications)
    if (length(specifications) < length(terms)) return(NULL)
    key <- function(predictors, group) paste(
        paste(sort(unique(as.character(predictors))), collapse='\034'),
        if (is.null(group)) '' else as.character(group), sep='\035'
    )
    keys <- vapply(specifications, function(specification) {
        key(.cdr_effect_spec_predictors(specification), specification$group)
    }, character(1))
    aligned <- vector('list', length(terms))
    position <- 0L
    for (index in seq_along(terms)) {
        info <- terms[[index]]
        target <- key(unique(c(
            info$linear_predictors,
            vapply(Filter(function(axis) identical(axis$role, 'predictor'),
                .cdr_effect_axes(info)), `[[`, character(1), 'variable')
        )), info$group)
        found <- which(keys[seq_along(keys) > position] == target)
        if (!length(found)) return(NULL)
        position <- position + found[[1L]]
        aligned[[index]] <- specifications[[position]]
    }
    aligned
}

.cdr_effect_spec_predictors <- function(specification) {
    if (is.null(specification)) return(character())
    if (!is.null(specification$predictors)) {
        return(as.character(specification$predictors))
    }
    if (isTRUE(specification$constant) || is.null(specification$predictor) ||
            identical(specification$predictor, '1')) return(character())
    as.character(specification$predictor)
}

.cdr_effect_linear_predictors <- function(info, specification) {
    if (!is.null(info$linear_predictors)) return(info$linear_predictors)
    predictors <- .cdr_effect_spec_predictors(specification)
    if (!is.null(specification$k_p)) {
        smooth <- which(!vapply(specification$k_p, is.null, logical(1)))
        return(predictors[setdiff(seq_along(predictors), smooth)])
    }
    if (!isTRUE(specification$nonlinear)) predictors else character()
}

.cdr_effect_data_summary <- function(values, name) {
    if (!is.numeric(values) || any(!is.finite(values))) {
        stop('Effect axis ', name, ' must be finite and numeric')
    }
    probabilities <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
    list(
        mean=mean(values),
        sd=stats::sd(values),
        quantiles=stats::setNames(
            as.numeric(stats::quantile(values, probabilities, names=FALSE)),
            format(probabilities, trim=TRUE, scientific=FALSE)
        )
    )
}

.cdr_effect_restore_summaries <- function(object, impulses, responses) {
    if (is.null(impulses) && is.null(responses)) return(object)
    if (!is.null(impulses) && !is.data.frame(impulses)) {
        stop('impulses must be a data frame')
    }
    if (!is.null(responses) && !is.data.frame(responses)) {
        stop('responses must be a data frame')
    }
    specifications <- .cdr_effect_specifications(object)
    for (index in seq_along(object$cdrgam$terms)) {
        info <- object$cdrgam$terms[[index]]
        specification <- specifications[[index]]
        linear <- .cdr_effect_linear_predictors(info, specification)
        missing <- setdiff(linear, names(info$linear_predictor_summaries))
        for (variable in missing) {
            if (is.null(impulses) || !(variable %in% names(impulses))) next
            info$linear_predictor_summaries[[variable]] <-
                .cdr_effect_data_summary(impulses[[variable]], variable)
        }
        if (!is.null(info$axis)) {
            info$axis <- lapply(info$axis, function(axis) {
                if (!is.null(axis$summary) || identical(axis$role, 'lag')) {
                    return(axis)
                }
                stream <- if (identical(axis$role, 'time')) responses else impulses
                if (!is.null(stream) && axis$variable %in% names(stream)) {
                    axis$summary <- .cdr_effect_data_summary(
                        stream[[axis$variable]], axis$variable
                    )
                }
                axis
            })
        }
        object$cdrgam$terms[[index]] <- info
    }
    object
}

#' List the estimable impulse-response effects in a fitted model
#'
#' Returns the term and axis metadata used to validate and construct effect
#' evaluation grids. Axis ranges and summaries are reported in the source units
#' used by callers.
#'
#' @param object A fitted `cdrgam` model.
#' @param impulses,responses Optional training streams used to recover axis
#'   summaries from models fitted before summaries were stored.
#' @return A data frame with one row per impulse-response term. `predictors`,
#'   `axes`, `axis_summaries`, and `centers` are list columns. `centers` gives
#'   the source-unit value at which each centered linear predictor makes the
#'   term zero; it is empty for terms without centered linear predictors.
#' @export
effect_catalog <- function(object, impulses=NULL, responses=NULL) {
    if (!is_cdrgam(object)) stop('object must be a fitted cdrgam model')
    object <- .cdr_effect_restore_summaries(object, impulses, responses)
    terms <- object$cdrgam$terms
    labels <- object$cdrgam$term_labels
    specifications <- .cdr_effect_specifications(object)
    if (!length(terms) || length(terms) != length(labels)) {
        stop('The fitted model does not contain stored IRF metadata')
    }
    if (is.null(specifications)) {
        specifications <- rep.int(list(NULL), length(terms))
    }
    rows <- lapply(seq_along(terms), function(index) {
        info <- terms[[index]]
        specification <- specifications[[index]]
        axes <- .cdr_effect_axes(info)
        axis_names <- unique(c(
            vapply(axes, .cdr_effect_axis_name, character(1)),
            names(info$linear_predictor_summaries),
            .cdr_effect_linear_predictors(info, specification)
        ))
        summaries <- lapply(axes, `[[`, 'summary')
        names(summaries) <- vapply(axes, .cdr_effect_axis_name, character(1))
        summaries <- c(summaries, info$linear_predictor_summaries)
        row <- data.frame(
            index=index,
            term_id=paste0('irf-', index),
            label=labels[[index]],
            type=info$type,
            grouped=length(info$group_levels) > 0L,
            group=if (is.null(info$group)) NA_character_ else info$group,
            stringsAsFactors=FALSE
        )
        row$predictors <- I(list(.cdr_effect_spec_predictors(specification)))
        row$axes <- I(list(axis_names))
        row$axis_summaries <- I(list(summaries))
        row$centers <- I(list(if (is.null(info$linear_predictor_centers)) {
            stats::setNames(numeric(), character())
        } else info$linear_predictor_centers))
        row
    })
    do.call(rbind, rows)
}

.cdr_effect_select <- function(object, terms) {
    catalog <- effect_catalog(object)
    if (is.null(terms)) return(catalog$index)
    if (is.character(terms)) {
        selected <- match(terms, catalog$label)
        if (anyNA(selected)) stop(
            'Unknown effect terms: ', paste(terms[is.na(selected)], collapse=', ')
        )
        return(unique(selected))
    }
    if (!is.list(terms)) stop('terms must be labels or a selector list')
    allowed <- c('names', 'predictors', 'match', 'grouped')
    unknown <- setdiff(names(terms), allowed)
    if (length(unknown)) stop('Unknown term selector fields: ', paste(unknown, collapse=', '))
    keep <- rep.int(TRUE, nrow(catalog))
    if (!is.null(terms$names)) keep <- keep & catalog$label %in% terms$names
    if (!is.null(terms$grouped)) keep <- keep & catalog$grouped == terms$grouped
    if (!is.null(terms$predictors)) {
        requested <- as.character(terms$predictors)
        mode <- if (is.null(terms$match)) 'contains' else terms$match
        if (!(mode %in% c('contains', 'exact'))) stop(
            'terms.match must be contains or exact'
        )
        keep <- keep & vapply(catalog$predictors, function(value) {
            if (identical(requested, '*')) return(length(value) > 0L)
            if (identical(mode, 'exact')) setequal(value, requested) else
                all(requested %in% value)
        }, logical(1))
    }
    selected <- catalog$index[keep]
    if (!length(selected)) stop('The effect term selector matched nothing')
    selected
}

.cdr_effect_axis_request <- function(axes, axis) {
    name <- .cdr_effect_axis_name(axis)
    if (identical(axis$role, 'predictor')) {
        predictors <- axes$predictors
        if (!is.null(predictors[[name]])) return(predictors[[name]])
        if (!is.null(predictors[['*']])) return(predictors[['*']])
    }
    axes[[name]]
}

.cdr_effect_summary_quantile <- function(summary, probability) {
    values <- as.numeric(summary$quantiles)
    probabilities <- suppressWarnings(as.numeric(names(summary$quantiles)))
    if (!length(values) || any(!is.finite(values)) ||
            length(probabilities) != length(values) || any(!is.finite(probabilities))) {
        stop('Stored effect-axis quantiles are invalid')
    }
    stats::approx(
        probabilities, values, probability,
        rule=2, ties='ordered'
    )$y
}

.cdr_effect_axis_values <- function(axis, request, default_n) {
    scale <- if (is.null(axis$scale)) 1 else axis$scale
    fitted <- range(axis$grid * scale)
    summary <- axis$summary
    if (is.null(summary)) summary <- list(
        mean=mean(fitted), sd=stats::sd(fitted),
        quantiles=stats::setNames(fitted, c('0', '1'))
    )
    if (is.null(request)) {
        if (identical(axis$role, 'lag')) {
            return(seq(fitted[[1L]], fitted[[2L]], length.out=default_n))
        }
        return(.cdr_effect_summary_quantile(summary, 0.5))
    }
    if (is.numeric(request)) return(as.numeric(request))
    if (!is.list(request)) stop('Axis requests must be numeric or mappings')
    if (!is.null(request$at)) {
        at <- request$at
        if (is.numeric(at)) return(as.numeric(at))
        if (!is.list(at) || is.null(at$summary) ||
                !(at$summary %in% c('mean', 'median'))) {
            stop('Axis at must be numeric or name mean or median')
        }
        value <- if (identical(at$summary, 'mean')) summary$mean else
            .cdr_effect_summary_quantile(summary, 0.5)
        value + if (is.null(at$`offset-sd`)) 0 else at$`offset-sd` * summary$sd
    } else if (!is.null(request$values)) {
        as.numeric(unlist(request$values, use.names=FALSE))
    } else if (!is.null(request$quantiles)) {
        probabilities <- as.numeric(unlist(request$quantiles, use.names=FALSE))
        vapply(probabilities, function(probability) {
            .cdr_effect_summary_quantile(summary, probability)
        }, numeric(1))
    } else if (identical(request$grid, 'fitted')) {
        n <- if (is.null(request$n)) default_n else as.integer(request$n)
        seq(fitted[[1L]], fitted[[2L]], length.out=n)
    } else stop('Axis request must define at, values, quantiles, or grid=fitted')
}

# With `points`, its columns give paired values of the axes they name, one
# evaluation point per row, and every other axis must resolve to one value.
.cdr_effect_grid <- function(info, axes, n, points=NULL) {
    metadata <- .cdr_effect_axes(info)
    linear <- unique(c(
        names(info$linear_predictor_summaries), info$linear_predictors
    ))
    if (length(linear)) {
        metadata <- c(metadata, lapply(linear, function(variable) {
            summary <- info$linear_predictor_summaries[[variable]]
            if (is.null(summary)) summary <- list(
                mean=NA_real_, sd=NA_real_,
                quantiles=c('0'=NA_real_, '1'=NA_real_)
            )
            list(
                role='predictor', variable=variable, label=variable,
                internal=NA_character_, grid=range(summary$quantiles),
                scale=1, summary=summary
            )
        }))
    }
    axis_names <- vapply(metadata, .cdr_effect_axis_name, character(1))
    if (!is.null(points)) {
        unknown <- setdiff(names(points), axis_names)
        if (length(unknown)) stop(
            'points columns are not evaluation axes of the term: ',
            paste(unknown, collapse=', ')
        )
    }
    values <- lapply(seq_along(metadata), function(index) {
        name <- axis_names[[index]]
        if (name %in% names(points)) {
            value <- points[[name]]
            if (!is.numeric(value) || any(!is.finite(value))) {
                stop('points column ', name, ' must contain finite numbers')
            }
            return(as.numeric(value))
        }
        value <- .cdr_effect_axis_values(
            metadata[[index]], .cdr_effect_axis_request(axes, metadata[[index]]), n
        )
        if (!length(value) || any(!is.finite(value))) {
            stop('Effect axes must contain finite values')
        }
        if (!is.null(points) && length(value) != 1L) stop(
            'Axis ', name, ' is not a points column, so it must resolve to one value'
        )
        value
    })
    names(values) <- axis_names
    if (!is.null(points)) {
        return(data.frame(
            lapply(values, rep_len, length.out=nrow(points)), check.names=FALSE
        ))
    }
    do.call(expand.grid, c(
        values, list(KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
    ))
}

.cdr_effect_basis <- function(info, grid) {
    metadata <- .cdr_effect_axes(info)
    internal <- stats::setNames(lapply(metadata, function(axis) {
        grid[[.cdr_effect_axis_name(axis)]] / if (is.null(axis$scale)) 1 else axis$scale
    }), vapply(metadata, `[[`, character(1), 'internal'))
    basis <- mgcv::PredictMat(info$basis, internal, n=nrow(grid))
    if (!is.null(info$transform)) basis <- basis %*% info$transform
    amplitude <- if (is.null(info$amplitude_scale)) 1 else info$amplitude_scale
    basis <- basis / amplitude
    centers <- info$linear_predictor_centers
    for (variable in info$linear_predictors) {
        center <- if (variable %in% names(centers)) centers[[variable]] else 0
        basis <- basis * (grid[[variable]] - center)
    }
    basis
}


.cdr_effect_smooth_flags <- function(specification) {
    predictors <- .cdr_effect_spec_predictors(specification)
    flags <- if (is.null(specification$k_p)) {
        rep.int(isTRUE(specification$nonlinear), length(predictors))
    } else !vapply(specification$k_p, is.null, logical(1))
    stats::setNames(as.logical(flags), predictors)
}

# Population terms whose sum with the selected term gives the response to an
# impulse at the requested axis values: every ungrouped term on the same `by`
# covariate whose predictors are a subset of the selected term's predictors
# with matching smooth or linear treatment. A time-invariant term is included
# below a time-varying one; the reverse is not.
.cdr_effect_hierarchy <- function(object, index) {
    specifications <- .cdr_effect_specifications(object)
    if (is.null(specifications)) return(integer())
    selected <- specifications[[index]]
    predictors <- .cdr_effect_spec_predictors(selected)
    flags <- .cdr_effect_smooth_flags(selected)
    candidates <- vapply(seq_along(specifications), function(j) {
        if (j == index) return(FALSE)
        candidate <- specifications[[j]]
        if (!is.null(candidate$group) || !identical(candidate$by, selected$by)) {
            return(FALSE)
        }
        if (!is.null(candidate$k_t) &&
                (is.null(selected$k_t) || !identical(candidate$k_t, selected$k_t))) {
            return(FALSE)
        }
        shared <- .cdr_effect_spec_predictors(candidate)
        if (!all(shared %in% predictors)) return(FALSE)
        if (length(shared) == length(predictors) &&
                identical(is.null(candidate$k_t), is.null(selected$k_t))) {
            return(FALSE)
        }
        identical(.cdr_effect_smooth_flags(candidate)[shared], flags[shared])
    }, logical(1))
    which(candidates)
}

.cdr_effect_term_predictors <- function(object, index) {
    specifications <- .cdr_effect_specifications(object)
    if (is.null(specifications)) {
        info <- object$cdrgam$terms[[index]]
        return(unique(c(
            info$linear_predictors,
            vapply(Filter(function(axis) identical(axis$role, 'predictor'),
                .cdr_effect_axes(info)), `[[`, character(1), 'variable')
        )))
    }
    .cdr_effect_spec_predictors(specifications[[index]])
}

# `omit` names predictors left at the model default: terms that involve any
# of them contribute zero.
.cdr_effect_parts <- function(
        object, index, grid, composition, grouping, population, group_level,
        omit=character()
) {
    info <- object$cdrgam$terms[[index]]
    grouped <- length(info$group_levels) > 0L
    basis <- function(term) {
        .cdr_effect_basis(object$cdrgam$terms[[term]], grid)
    }
    parts <- list()
    add_part <- function(term, group_level=NULL) {
        if (length(omit) &&
                any(.cdr_effect_term_predictors(object, term) %in% omit)) return()
        parts[[length(parts) + 1L]] <<- .cdr_plot_term_part(
            object, term, basis(term), group_level
        )
    }
    if (!grouped || grouping %in% c('deviation', 'conditional')) {
        add_part(index, if (grouped) group_level else NULL)
    }
    base <- if (grouped) population else index
    if (grouped && grouping %in% c('population', 'conditional')) {
        add_part(population)
    }
    if (identical(composition, 'total') && (!grouped || grouping != 'deviation')) {
        group <- if (grouped) info$group else NULL
        for (parent in .cdr_effect_hierarchy(object, base)) {
            add_part(parent)
            if (!grouped || !identical(grouping, 'conditional')) next
            deviations <- which(vapply(seq_along(object$cdrgam$terms), function(j) {
                term <- object$cdrgam$terms[[j]]
                identical(term$group, group) && length(term$group_levels) > 0L &&
                    identical(.cdr_matching_term(object, j, 'population'), parent)
            }, logical(1)))
            for (deviation in deviations) {
                add_part(deviation, group_level)
            }
        }
    }
    parts
}

.cdr_effect_design <- function(parts, rows, indices=NULL) {
    if (is.null(indices)) {
        indices <- sort(unique(unlist(lapply(parts, `[[`, 'index'))))
    }
    design <- matrix(0, nrow=rows, ncol=length(indices))
    for (part in parts) {
        destination <- match(part$index, indices)
        design[, destination] <- design[, destination, drop=FALSE] + part$basis
    }
    list(index=indices, design=design)
}

.cdr_effect_boundary_object <- function(object, boundary, unconditional) {
    if (identical(boundary, 'estimate')) return(object)
    if (isTRUE(unconditional)) {
        stop(
            'boundary = "bound" is a conditional sensitivity analysis; ',
            'use unconditional = FALSE'
        )
    }
    if (!inherits(object, 'cdrgam_sparse')) {
        stop('boundary = "bound" is available only for sparse fits')
    }
    .sparse_boundary_bound_object(object)
}

.cdr_effect_covariance <- function(object, indices, unconditional) {
    if (inherits(object, 'cdrgam_sparse')) {
        return(.sparse_selected_vcov(object, indices, unconditional=unconditional))
    }
    stats::vcov(object, unconditional=unconditional)[indices, indices, drop=FALSE]
}

.cdr_effect_is_default <- function(request) {
    (is.atomic(request) && length(request) == 1L && is.na(request)) ||
        identical(request, 'default')
}

# Returns the reference grid, with attribute `omit` naming predictors left at
# the model default.
.cdr_effect_reference_grid <- function(info, grid, reference, n) {
    if (!length(reference)) return(NULL)
    metadata <- .cdr_effect_axes(info)
    available <- names(grid)
    unknown <- setdiff(names(reference), available)
    if (length(unknown)) stop(
        'Reference axes are not evaluation axes of the term: ',
        paste(unknown, collapse=', ')
    )
    output <- grid
    omit <- character()
    for (name in names(reference)) {
        request <- reference[[name]]
        if (.cdr_effect_is_default(request)) {
            axis <- Find(function(value) identical(.cdr_effect_axis_name(value), name),
                metadata)
            if (!is.null(axis) && !identical(axis$role, 'predictor')) stop(
                'Reference axis ', name, ' has no model default; give a value'
            )
            omit <- c(omit, name)
            output[[name]] <- NA_real_
            next
        }
        axis <- Find(function(value) identical(.cdr_effect_axis_name(value), name),
            metadata)
        if (is.null(axis)) {
            summary <- info$linear_predictor_summaries[[name]]
            axis <- list(
                role='predictor', variable=name, grid=range(summary$quantiles),
                scale=1, summary=summary
            )
        }
        value <- .cdr_effect_axis_values(axis, request, n)
        if (length(value) != 1L || !is.finite(value)) {
            stop('Reference axis ', name, ' must resolve to one finite value')
        }
        output[[name]] <- value
    }
    attr(output, 'omit') <- omit
    output
}

# Builds one linear functional per requested grouping level. Each element
# holds the evaluation grid, the reference grid, and the design matrix over
# the coefficient indices it uses.
.cdr_effect_functionals <- function(
        object, index, axes, composition, grouping, groups, reference, n,
        points=NULL
) {
    info <- object$cdrgam$terms[[index]]
    grid <- .cdr_effect_grid(info, axes, n, points)
    reference_grid <- .cdr_effect_reference_grid(info, grid, reference, n)
    grouped <- length(info$group_levels) > 0L
    if (!grouped && !identical(grouping, 'population')) stop(
        'grouping=', grouping, ' requires grouped IRF terms'
    )
    population <- if (grouped) .cdr_matching_term(object, index, 'population') else index
    if (grouped && grouping %in% c('population', 'conditional') &&
            length(population) != 1L) stop(
        'A unique population IRF could not be found for ',
        object$cdrgam$term_labels[[index]]
    )
    levels <- if (grouped && grouping %in% c('deviation', 'conditional')) {
        if (is.null(groups)) info$group_levels else as.character(groups)
    } else NA_character_
    if (grouped && any(!is.na(levels) & !(levels %in% info$group_levels))) {
        stop('Unknown grouped-IRF level requested')
    }
    lapply(levels, function(group_level) {
        parts <- .cdr_effect_parts(
            object, index, grid, composition, grouping, population, group_level
        )
        value <- .cdr_effect_design(parts, nrow(grid))
        if (!is.null(reference_grid)) {
            baseline <- .cdr_effect_design(.cdr_effect_parts(
                object, index, reference_grid, composition, grouping, population,
                group_level, omit=attr(reference_grid, 'omit')
            ), nrow(grid), value$index)
            value$design <- value$design - baseline$design
        }
        value$grid <- grid
        value$reference <- reference_grid
        value$group <- group_level
        value
    })
}

.cdr_effect_frame <- function(object, index, functional, composition, grouping) {
    output <- transform(
        functional$grid,
        term_id=paste0('irf-', index),
        term=object$cdrgam$term_labels[[index]],
        group=if (is.na(functional$group)) NA_character_ else functional$group,
        composition=composition,
        grouping=grouping
    )
    if (!is.null(functional$reference)) {
        for (name in attr(functional$reference, 'axes')) {
            output[[paste0(name, '_reference')]] <- functional$reference[[name]]
        }
    }
    output
}

.cdr_effect_validate <- function(object, axes, reference, level, n, points=NULL) {
    if (!is_cdrgam(object)) stop('object must be a fitted cdrgam model')
    if (!is.list(axes)) stop('axes must be a list')
    if (!is.null(points) && (!is.data.frame(points) || !nrow(points) ||
            !ncol(points) || is.null(names(points)) || any(!nzchar(names(points))) ||
            anyDuplicated(names(points)))) {
        stop('points must be a data frame with at least one row and uniquely named columns')
    }
    if (!is.list(reference) || (length(reference) && (is.null(names(reference)) ||
            any(!nzchar(names(reference))) || anyDuplicated(names(reference))))) {
        stop('reference must be a uniquely named list')
    }
    if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
            level <= 0 || level >= 1) stop('level must lie strictly between zero and one')
    n <- as.integer(n)
    if (is.na(n) || n < 2L) stop('n must be at least two')
    n
}

.cdr_effect_all_functionals <- function(
        object, terms, axes, composition, grouping, groups, reference, n,
        points=NULL
) {
    indices <- .cdr_effect_select(object, terms)
    output <- list()
    for (index in indices) {
        functionals <- .cdr_effect_functionals(
            object, index, axes, composition, grouping, groups, reference, n,
            points
        )
        for (functional in functionals) {
            if (!is.null(functional$reference)) {
                attr(functional$reference, 'axes') <- names(reference)
            }
            functional$term <- index
            output[[length(output) + 1L]] <- functional
        }
    }
    output
}

.cdr_effect_bind <- function(frames) {
    columns <- unique(unlist(lapply(frames, names), use.names=FALSE))
    frames <- lapply(frames, function(value) {
        for (name in setdiff(columns, names(value))) value[[name]] <- NA
        value[columns]
    })
    result <- do.call(rbind, frames)
    rownames(result) <- NULL
    result
}

#' Evaluate fitted effects on declarative axis grids
#'
#' Computes selected impulse-response effects and their pointwise uncertainty as
#' joint linear functionals of the fitted coefficients. Axes may vary over
#' numeric values, quantiles, or the fitted domain, or be fixed at a numeric
#' value or a training-data mean or median with an optional standard deviation
#' offset.
#'
#' With `composition = "total"`, the selected term is summed with every
#' ungrouped term in its hierarchy: terms on the same `by` covariate whose
#' predictors are a subset of the selected term's predictors, with the same
#' linear or smooth treatment of each shared predictor. For
#' `irf(x, z, k_p = list(NULL, 3))` this adds `irf(x)`, `irf(z, k_p = 3)`, and
#' `irf(1)` when they are present. A time-invariant term is included below a
#' time-varying term with the same predictors. Under
#' `grouping = "conditional"`, the grouped deviations of those terms at the
#' same level are also added.
#'
#' `reference` subtracts the same functional evaluated at reference values of
#' some axes. Axes not named in `reference` keep their evaluated values. The
#' standard errors account for the covariance of the two evaluations. A
#' predictor whose reference is `NA` or `"default"` is left at the model
#' default: every term that involves it contributes zero to the reference.
#' For a centered linear predictor this equals its center. For a smooth
#' predictor it is the training-average contribution implied by the term's
#' identifiability constraint. With every predictor at its default, the
#' reference under `composition = "total"` reduces to `irf(1)`, the response
#' to an average event.
#'
#' For example, with `composition = "total"`, `x` at its mean plus one
#' standard deviation, and `reference = list(x = list(at = list(summary =
#' "mean")))`, the result is the effect of a one-SD increase in `x` at each
#' requested value of `z`. A reference value of `z` gives differences between
#' the interaction curves.
#'
#' @param object A fitted `cdrgam` model.
#' @param terms Term labels, or a selector list with optional fields `names`
#'   (term labels), `predictors` (predictor names, or `"*"` for any term with
#'   predictors), `match` (`"contains"`, the default, or `"exact"`), and
#'   `grouped` (logical). `NULL` selects every impulse-response term.
#' @param impulses,responses Optional training streams used to recover axis
#'   summaries from models fitted before summaries were stored.
#' @param axes Named axis requests. Predictor requests are nested under
#'   `predictors`; `"*"` supplies a fallback for every predictor axis.
#' @param points Optional data frame of evaluation points, one per row, with
#'   columns named by axis, such as `lag` or a predictor. Its columns are
#'   paired rather than crossed, so the result has one row per point. Each
#'   axis that is not a column takes its value from `axes`, and that request
#'   must resolve to one value. Use it to evaluate along a path through the
#'   axes, such as each predictor value at its own lag.
#' @param composition Return the selected term alone or add its hierarchy of
#'   fitted marginal terms.
#' @param grouping Return population, grouped-deviation, or conditional effects.
#' @param groups Optional grouping levels.
#' @param reference Optional named list of reference values for evaluation
#'   axes, using the same request forms as `axes`. Each request must resolve
#'   to one value. A predictor may instead be `NA` or `"default"` to leave it
#'   at the model default. Names are axis names such as `lag` or a predictor.
#' @param se Include pointwise standard errors.
#' @param unconditional Include smoothing-parameter uncertainty when supported.
#'   Smoothing parameters classified as on the boundary of their range
#'   contribute no uncertainty, because the criterion has no local curvature
#'   along them. Sparse fits stop with an error when the outer Hessian is not
#'   positive definite after those parameters are excluded.
#' @param boundary `"estimate"` (the default) uses the fitted smoothing
#'   parameters. `"bound"` repeats the calculation with each boundary
#'   smoothing parameter set to its profile bound, the largest variance
#'   component the data still support, while the other smoothing parameters
#'   stay at their estimates. Compare the two to check whether conclusions
#'   depend on treating those variance components as zero. `"bound"` requires
#'   a sparse Gaussian REML fit with `sparse_control$boundary_sensitivity =
#'   TRUE` and `unconditional = FALSE`.
#' @param level Confidence level for `lower` and `upper`.
#' @param n Default number of points for fitted-domain grids.
#' @return A `cdrgam_effect_grid` data frame with evaluation axes, estimates,
#'   standard errors, and pointwise confidence limits on the linear-predictor
#'   scale. Each reference axis adds a
#'   `<axis>_reference` column, which is `NA` for a model default.
#' @seealso [test_effect()] for joint Wald tests of the same functionals.
#' @export
estimate_effect <- function(
        object, terms=NULL, axes=list(), impulses=NULL, responses=NULL,
        composition=c('term', 'total'),
        grouping=c('population', 'deviation', 'conditional'), groups=NULL,
        reference=list(), se=TRUE, unconditional=FALSE, level=0.95, n=200L,
        boundary=c('estimate', 'bound'), points=NULL
) {
    n <- .cdr_effect_validate(object, axes, reference, level, n, points)
    object <- .cdr_effect_restore_summaries(object, impulses, responses)
    object <- .cdr_effect_boundary_object(object, match.arg(boundary), unconditional)
    composition <- match.arg(composition)
    grouping <- match.arg(grouping)
    functionals <- .cdr_effect_all_functionals(
        object, terms, axes, composition, grouping, groups, reference, n, points
    )
    multiplier <- stats::qnorm((1 + level) / 2)
    frames <- lapply(functionals, function(functional) {
        estimate <- drop(functional$design %*% stats::coef(object)[functional$index])
        standard_error <- rep.int(NA_real_, length(estimate))
        if (isTRUE(se)) {
            covariance <- .cdr_effect_covariance(
                object, functional$index, unconditional
            )
            standard_error <- sqrt(pmax(0, rowSums(
                (functional$design %*% covariance) * functional$design
            )))
        }
        output <- .cdr_effect_frame(
            object, functional$term, functional, composition, grouping
        )
        output$estimate <- estimate
        output$se <- standard_error
        output$lower <- estimate - multiplier * standard_error
        output$upper <- estimate + multiplier * standard_error
        output
    })
    result <- .cdr_effect_bind(frames)
    class(result) <- c('cdrgam_effect_grid', class(result))
    attr(result, 'level') <- level
    attr(result, 'scale') <- 'linear-predictor'
    attr(result, 'reference') <- if (length(reference)) names(reference) else NULL
    result
}

# Draws of the coefficients at `indices` from N(beta-hat, V), one column per
# draw, through the eigendecomposition of V.
.cdr_effect_sample <- function(object, indices, unconditional, draws) {
    coefficients <- stats::coef(object)[indices]
    covariance <- .cdr_effect_covariance(object, indices, unconditional)
    decomposition <- eigen((covariance + t(covariance)) / 2, symmetric=TRUE)
    root <- decomposition$vectors %*%
        diag(sqrt(pmax(0, decomposition$values)), length(coefficients))
    normal <- matrix(stats::rnorm(length(coefficients) * draws),
        length(coefficients), draws)
    coefficients + root %*% normal
}

#' Sample coefficients shared by several effect simulations
#'
#' Draws once from the approximate coefficient posterior
#' \eqn{N(\hat\beta, V)}, restricted to the coefficients that the selected
#' terms use under `composition` and `grouping`. Passing the result to
#' [simulate_effect()] as `coefficients` evaluates every call on the same
#' draws, so summaries computed from different calls are joint: for
#' example, the lag at which one contrast peaks in each draw and another
#' effect at that lag in the same draw.
#'
#' @inheritParams simulate_effect
#' @return An object of class `cdrgam_coefficient_draws`: a list with
#'   `index`, the coefficient positions; `draws`, a matrix with one row per
#'   coefficient in `index` and one column per draw; `unconditional`; and
#'   `boundary`.
#' @seealso [simulate_effect()].
#' @export
sample_effect_coefficients <- function(
        object, terms=NULL, impulses=NULL, responses=NULL,
        composition=c('term', 'total'),
        grouping=c('population', 'deviation', 'conditional'), groups=NULL,
        unconditional=FALSE, draws=1000L, seed=NULL,
        boundary=c('estimate', 'bound')
) {
    .cdr_effect_validate(object, list(), list(), 0.95, 2L, NULL)
    draws <- suppressWarnings(as.integer(draws))
    if (length(draws) != 1L || is.na(draws) || draws < 1L) {
        stop('draws must be a positive integer')
    }
    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
            stop('seed must be one finite number')
        }
        exists <- exists('.Random.seed', envir=globalenv(), inherits=FALSE)
        if (exists) state <- get('.Random.seed', envir=globalenv())
        on.exit(if (exists) {
            assign('.Random.seed', state, envir=globalenv())
        } else rm('.Random.seed', envir=globalenv()), add=TRUE)
        set.seed(as.integer(seed))
    }
    boundary <- match.arg(boundary)
    object <- .cdr_effect_restore_summaries(object, impulses, responses)
    object <- .cdr_effect_boundary_object(object, boundary, unconditional)
    functionals <- .cdr_effect_all_functionals(
        object, terms, list(), match.arg(composition), match.arg(grouping), groups,
        list(), 2L
    )
    index <- sort(unique(unlist(lapply(functionals, `[[`, 'index'))))
    structure(list(
        index=index, draws=.cdr_effect_sample(object, index, unconditional, draws),
        unconditional=isTRUE(unconditional), boundary=boundary
    ), class='cdrgam_coefficient_draws')
}

#' Simulate fitted effects from the coefficient posterior
#'
#' Draws coefficient vectors from the approximate posterior
#' \eqn{N(\hat\beta, V)} and evaluates the functionals that
#' [estimate_effect()] would return for each draw. Use the draws for
#' uncertainty in nonlinear summaries of an effect, such as the lag at which
#' a response peaks, that pointwise standard errors do not determine. The
#' arguments select the same terms, axes, compositions, groupings, and
#' references.
#'
#' Draws for each term and grouping level are sampled separately, so draws
#' are joint across the rows of one functional but not across functionals.
#' `V` is the conditional Bayesian covariance, or the covariance corrected
#' for smoothing-parameter uncertainty when `unconditional = TRUE`. Draws are
#' on the linear-predictor scale. A rank-deficient `V` is handled through its
#' eigendecomposition, and directions with negative eigenvalues contribute
#' nothing.
#'
#' @inheritParams estimate_effect
#' @param draws Number of posterior draws.
#' @param seed Optional integer seed. When given, the global random-number
#'   state is restored on exit.
#' @param coefficients Optional coefficient draws from
#'   [sample_effect_coefficients()]. When given, every functional is
#'   evaluated on those draws, so draws are joint across functionals and
#'   across calls that share them; `draws`, `seed`, and `unconditional` are
#'   then taken from the sample.
#' @return A list with `grid`, a data frame like the one returned by
#'   [estimate_effect()] with `estimate` but no interval columns, and
#'   `draws`, a matrix with one row per row of `grid` and one column per
#'   draw.
#' @seealso [estimate_effect()] for pointwise intervals of the same
#'   functionals.
#' @export
simulate_effect <- function(
        object, terms=NULL, axes=list(), impulses=NULL, responses=NULL,
        composition=c('term', 'total'),
        grouping=c('population', 'deviation', 'conditional'), groups=NULL,
        reference=list(), unconditional=FALSE, draws=1000L, n=200L, seed=NULL,
        boundary=c('estimate', 'bound'), points=NULL, coefficients=NULL
) {
    n <- .cdr_effect_validate(object, axes, reference, 0.95, n, points)
    boundary <- match.arg(boundary)
    if (!is.null(coefficients)) {
        if (!inherits(coefficients, 'cdrgam_coefficient_draws')) {
            stop('coefficients must come from sample_effect_coefficients()')
        }
        if (!identical(coefficients$boundary, boundary)) stop(
            'coefficients were sampled with boundary = "', coefficients$boundary, '"'
        )
        unconditional <- coefficients$unconditional
        draws <- ncol(coefficients$draws)
        seed <- NULL
    }
    draws <- suppressWarnings(as.integer(draws))
    if (length(draws) != 1L || is.na(draws) || draws < 1L) {
        stop('draws must be a positive integer')
    }
    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
            stop('seed must be one finite number')
        }
        exists <- exists('.Random.seed', envir=globalenv(), inherits=FALSE)
        if (exists) state <- get('.Random.seed', envir=globalenv())
        on.exit(if (exists) {
            assign('.Random.seed', state, envir=globalenv())
        } else rm('.Random.seed', envir=globalenv()), add=TRUE)
        set.seed(as.integer(seed))
    }
    object <- .cdr_effect_restore_summaries(object, impulses, responses)
    object <- .cdr_effect_boundary_object(object, boundary, unconditional)
    composition <- match.arg(composition)
    grouping <- match.arg(grouping)
    functionals <- .cdr_effect_all_functionals(
        object, terms, axes, composition, grouping, groups, reference, n, points
    )
    shared <- coefficients
    simulated <- lapply(functionals, function(functional) {
        coefficients <- stats::coef(object)[functional$index]
        sample <- if (is.null(shared)) {
            .cdr_effect_sample(object, functional$index, unconditional, draws)
        } else {
            rows <- match(functional$index, shared$index)
            if (anyNA(rows)) stop(
                'coefficients do not cover term ', functional$term,
                '; sample them for every term this call evaluates'
            )
            shared$draws[rows, , drop=FALSE]
        }
        frame <- .cdr_effect_frame(
            object, functional$term, functional, composition, grouping
        )
        frame$estimate <- drop(functional$design %*% coefficients)
        list(frame=frame, draws=functional$design %*% sample)
    })
    grid <- .cdr_effect_bind(lapply(simulated, `[[`, 'frame'))
    attr(grid, 'scale') <- 'linear-predictor'
    attr(grid, 'reference') <- if (length(reference)) names(reference) else NULL
    values <- do.call(rbind, lapply(simulated, `[[`, 'draws'))
    dimnames(values) <- NULL
    list(grid=grid, draws=values)
}

#' Test fitted effects jointly over declarative axis grids
#'
#' Computes a Wald test of the null hypothesis that a functional returned by
#' [estimate_effect()] is zero at every grid point. The arguments select the
#' same terms, axes, compositions, groupings, and references. One test is
#' returned per selected term and grouping level.
#'
#' The statistic is \eqn{\hat\beta^\top L^\top (L V L^\top)^{-} L \hat\beta},
#' where `L` stacks the grid rows of the functional and `V` is the coefficient
#' covariance. The pseudo-inverse keeps eigenvalues above `tolerance` times the
#' largest, and the reference distribution is chi-squared with that rank as
#' its degrees of freedom. Use a grid dense enough to span the functional;
#' the rank stops growing once it does.
#'
#' Penalized coefficients make this test
#' approximate. Wood (2013) reports that such tests can reject too rarely when
#' the null holds on the boundary of the smoothing-parameter space.
#'
#' @inheritParams estimate_effect
#' @param tolerance Relative eigenvalue threshold for the rank of `L V L'`.
#' @return A `cdrgam_effect_test` data frame with `term_id`, `term`, `group`,
#'   `composition`, `grouping`, `statistic`, `df`, `p_value`, and `points`.
#'   `df` is zero and `statistic` and `p_value` are missing when the functional
#'   is identically zero.
#' @references Wood, S. N. (2013). On p-values for smooth components of an
#'   extended generalized additive model. *Biometrika*, 100(1), 221--228.
#' @seealso [estimate_effect()]
#' @export
test_effect <- function(
        object, terms=NULL, axes=list(), impulses=NULL, responses=NULL,
        composition=c('term', 'total'),
        grouping=c('population', 'deviation', 'conditional'), groups=NULL,
        reference=list(), unconditional=FALSE, n=50L, tolerance=1e-8,
        boundary=c('estimate', 'bound'), points=NULL
) {
    n <- .cdr_effect_validate(object, axes, reference, 0.95, n, points)
    if (!is.numeric(tolerance) || length(tolerance) != 1L ||
            !is.finite(tolerance) || tolerance <= 0 || tolerance >= 1) {
        stop('tolerance must lie strictly between zero and one')
    }
    object <- .cdr_effect_restore_summaries(object, impulses, responses)
    object <- .cdr_effect_boundary_object(object, match.arg(boundary), unconditional)
    composition <- match.arg(composition)
    grouping <- match.arg(grouping)
    functionals <- .cdr_effect_all_functionals(
        object, terms, axes, composition, grouping, groups, reference, n, points
    )
    rows <- lapply(functionals, function(functional) {
        estimate <- drop(functional$design %*% stats::coef(object)[functional$index])
        covariance <- .cdr_effect_covariance(object, functional$index, unconditional)
        projected <- functional$design %*% covariance %*% t(functional$design)
        decomposition <- eigen((projected + t(projected)) / 2, symmetric=TRUE)
        values <- decomposition$values
        keep <- values > tolerance * max(c(values, 0))
        rank <- sum(keep)
        statistic <- if (rank) {
            scores <- crossprod(decomposition$vectors[, keep, drop=FALSE], estimate)
            sum(scores^2 / values[keep])
        } else NA_real_
        data.frame(
            term_id=paste0('irf-', functional$term),
            term=object$cdrgam$term_labels[[functional$term]],
            group=if (is.na(functional$group)) NA_character_ else functional$group,
            composition=composition,
            grouping=grouping,
            statistic=statistic,
            df=rank,
            p_value=if (rank) stats::pchisq(statistic, rank, lower.tail=FALSE) else
                NA_real_,
            points=length(estimate),
            stringsAsFactors=FALSE
        )
    })
    result <- .cdr_effect_bind(rows)
    class(result) <- c('cdrgam_effect_test', class(result))
    attr(result, 'reference') <- if (length(reference)) names(reference) else NULL
    result
}
