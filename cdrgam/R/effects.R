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
    specifications <- object$cdrgam$preparation$specification
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
#' Returns the term metadata used to validate and construct effect grids. Axis
#' ranges and summaries are reported in the source units used by callers.
#'
#' @param object A fitted `cdrgam` model.
#' @param impulses,responses Optional training streams used to recover axis
#'   summaries from models fitted before summaries were stored.
#' @return A data frame with one row per impulse-response term. `predictors`,
#'   `axes`, and `axis_summaries` are list columns.
#' @export
effect_catalog <- function(object, impulses=NULL, responses=NULL) {
    if (!is_cdrgam(object)) stop('object must be a fitted cdrgam model')
    object <- .cdr_effect_restore_summaries(object, impulses, responses)
    terms <- object$cdrgam$terms
    labels <- object$cdrgam$term_labels
    specifications <- object$cdrgam$preparation$specification
    if (!length(terms) || length(terms) != length(labels)) {
        stop('The fitted model does not contain stored IRF metadata')
    }
    if (is.null(specifications) || length(specifications) != length(terms)) {
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

.cdr_effect_grid <- function(info, axes, n) {
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
    values <- lapply(metadata, function(axis) {
        value <- .cdr_effect_axis_values(
            axis, .cdr_effect_axis_request(axes, axis), n
        )
        if (!length(value) || any(!is.finite(value))) {
            stop('Effect axes must contain finite values')
        }
        value
    })
    names(values) <- vapply(metadata, .cdr_effect_axis_name, character(1))
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
    for (variable in info$linear_predictors) {
        basis <- basis * grid[[variable]]
    }
    basis
}

.cdr_effect_parent_chain <- function(object, index) {
    output <- integer()
    current <- index
    repeat {
        parent <- .cdr_matching_term(object, current, 'baseline')
        if (length(parent) != 1L || parent %in% c(current, output)) break
        output <- c(output, parent)
        current <- parent
    }
    output
}

.cdr_effect_evaluate <- function(
        object, index, axes, composition, grouping, groups, se, unconditional,
        level, n
) {
    info <- object$cdrgam$terms[[index]]
    grid <- .cdr_effect_grid(info, axes, n)
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
    output <- lapply(levels, function(group_level) {
        parts <- list()
        if (!grouped || grouping %in% c('deviation', 'conditional')) {
            parts[[length(parts) + 1L]] <- .cdr_plot_term_part(
                object, index, .cdr_effect_basis(info, grid),
                if (grouped) group_level else NULL
            )
        }
        base <- if (grouped) population else index
        if (grouped && grouping %in% c('population', 'conditional')) {
            parts[[length(parts) + 1L]] <- .cdr_plot_term_part(
                object, population,
                .cdr_effect_basis(object$cdrgam$terms[[population]], grid)
            )
        }
        if (identical(composition, 'total') &&
                (!grouped || grouping != 'deviation')) {
            for (parent in .cdr_effect_parent_chain(object, base)) {
                parts[[length(parts) + 1L]] <- .cdr_plot_term_part(
                    object, parent,
                    .cdr_effect_basis(object$cdrgam$terms[[parent]], grid)
                )
            }
        }
        values <- .cdr_plot_combine_parts(object, parts, se, unconditional)
        multiplier <- stats::qnorm((1 + level) / 2)
        transform(
            grid,
            term_id=paste0('irf-', index),
            term=object$cdrgam$term_labels[[index]],
            group=if (is.na(group_level)) NA_character_ else group_level,
            composition=composition,
            grouping=grouping,
            estimate=values$estimate,
            se=values$se,
            lower=values$estimate - multiplier * values$se,
            upper=values$estimate + multiplier * values$se
        )
    })
    do.call(rbind, output)
}

#' Evaluate fitted effects on declarative axis grids
#'
#' Computes selected impulse-response effects and their pointwise uncertainty as
#' joint linear functionals of the fitted coefficients. Tensor-product axes may
#' vary on grids or be fixed at numeric or training-summary values.
#'
#' @param object A fitted `cdrgam` model.
#' @param terms Term labels or a selector returned by the visualization schema.
#' @param impulses,responses Optional training streams used to recover axis
#'   summaries from models fitted before summaries were stored.
#' @param axes Named axis requests. Predictor requests are nested under
#'   `predictors`; `"*"` supplies a fallback for every predictor axis.
#' @param composition Return the selected term alone or add its hierarchy of
#'   fitted marginal terms.
#' @param grouping Return population, grouped-deviation, or conditional effects.
#' @param groups Optional grouping levels.
#' @param se Include pointwise standard errors.
#' @param unconditional Include smoothing-parameter uncertainty when supported.
#' @param level Confidence level for `lower` and `upper`.
#' @param n Default number of points for fitted-domain grids.
#' @return A data frame with evaluation axes, estimates, standard errors, and
#'   pointwise confidence limits.
#' @export
estimate_effect <- function(
        object, terms=NULL, axes=list(), impulses=NULL, responses=NULL,
        composition=c('term', 'deviation', 'total'),
        grouping=c('population', 'deviation', 'conditional'), groups=NULL,
        se=TRUE, unconditional=FALSE, level=0.95, n=200L
) {
    if (!is_cdrgam(object)) stop('object must be a fitted cdrgam model')
    object <- .cdr_effect_restore_summaries(object, impulses, responses)
    composition <- match.arg(composition)
    grouping <- match.arg(grouping)
    if (!is.list(axes)) stop('axes must be a list')
    if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
            level <= 0 || level >= 1) stop('level must lie strictly between zero and one')
    n <- as.integer(n)
    if (is.na(n) || n < 2L) stop('n must be at least two')
    indices <- .cdr_effect_select(object, terms)
    output <- lapply(indices, function(index) .cdr_effect_evaluate(
        object, index, axes, composition, grouping, groups, se,
        unconditional, level, n
    ))
    columns <- unique(unlist(lapply(output, names), use.names=FALSE))
    output <- lapply(output, function(value) {
        for (name in setdiff(columns, names(value))) value[[name]] <- NA
        value[columns]
    })
    result <- do.call(rbind, output)
    rownames(result) <- NULL
    class(result) <- c('cdrgam_effect_grid', class(result))
    attr(result, 'level') <- level
    attr(result, 'scale') <- 'linear-predictor'
    result
}
