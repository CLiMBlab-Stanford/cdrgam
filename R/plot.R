.cdr_plot_registry <- function(object) {
    metadata <- object$cdrgam$terms
    labels <- object$cdrgam$term_labels
    specifications <- object$cdrgam$preparation$specification
    if (!length(metadata) || length(metadata) != length(labels)) {
        stop('The fitted model does not contain stored IRF metadata')
    }
    if (is.null(specifications) || length(specifications) != length(metadata)) {
        specifications <- rep.int(list(NULL), length(metadata))
    }
    data.frame(
        index=seq_along(metadata),
        label=labels,
        type=vapply(metadata, `[[`, character(1), 'type'),
        grouped=vapply(metadata, function(info) {
            length(info$group_levels) > 0L
        }, logical(1)),
        surface=vapply(metadata, function(info) {
            (!is.null(info$axis) && length(info$axis) > 1L) ||
                startsWith(info$type, 'nonlinear') ||
                startsWith(info$type, 'varying')
        }, logical(1)),
        group=vapply(metadata, function(info) {
            if (is.null(info$group)) NA_character_ else info$group
        }, character(1)),
        stringsAsFactors=FALSE
    )
}

.cdr_source_scale_gam_plot <- function(object) {
    gam <- as_gam(object)
    scaling <- object$cdrgam$scaling
    if (is.null(scaling)) scaling <- object$cdrgam$preparation$scaling
    if (is.null(scaling) || !isTRUE(scaling$enabled)) return(gam)
    restored <- character()
    for (i in seq_along(gam$smooth)) {
        smooth <- gam$smooth[[i]]
        if (inherits(smooth, 'cdr.smooth')) next
        variables <- unique(c(smooth$term, smooth$by))
        variables <- variables[!is.na(variables) & variables != 'NA']
        divisors <- stats::setNames(vapply(variables, function(variable) {
            .cdr_scaling_divisor(scaling, 'responses', variable)
        }, numeric(1)), variables)
        divisors <- divisors[divisors != 1]
        if (!length(divisors)) next
        original <- smooth
        smooth$cdr_source_smooth <- original
        smooth$cdr_source_divisors <- divisors
        class(smooth) <- c('cdr.rescaled.smooth', class(smooth))
        gam$smooth[[i]] <- smooth
        restored <- union(restored, names(divisors))
    }
    if (!is.null(gam$model)) {
        for (variable in restored) {
            if (variable %in% names(gam$model) && is.numeric(gam$model[[variable]])) {
                gam$model[[variable]] <- gam$model[[variable]] *
                    .cdr_scaling_divisor(scaling, 'responses', variable)
            }
        }
    }
    gam
}

#' @export
Predict.matrix.cdr.rescaled.smooth <- function(object, data) {
    transformed <- data
    for (variable in names(object$cdr_source_divisors)) {
        if (variable %in% names(transformed)) {
            transformed[[variable]] <- transformed[[variable]] /
                object$cdr_source_divisors[[variable]]
        }
    }
    mgcv::Predict.matrix(object$cdr_source_smooth, transformed)
}

.cdr_plot_select <- function(registry, select) {
    if (is.null(select)) return(registry$index)
    if (is.numeric(select)) {
        indices <- as.integer(select)
        if (any(!is.finite(select)) || any(select != indices) ||
                any(indices < 1L | indices > nrow(registry))) {
            stop('Numeric select values must identify fitted IRF terms')
        }
        return(unique(indices))
    }
    if (!is.character(select) || !length(select)) {
        stop('select must contain IRF labels or one-based indices')
    }
    indices <- match(select, registry$label)
    if (anyNA(indices)) {
        stop('Unknown IRF terms: ', paste(select[is.na(indices)], collapse=', '))
    }
    unique(indices)
}

.cdr_plot_basis <- function(info, lag, predictor=NULL) {
    lag_scale <- if (is.null(info$lag_scale)) 1 else info$lag_scale
    predictor_scale <- if (is.null(info$predictor_scale)) {
        1
    } else info$predictor_scale
    amplitude_scale <- if (is.null(info$amplitude_scale)) {
        1
    } else info$amplitude_scale
    surface <- (!is.null(info$axis) && length(info$axis) > 1L) ||
        startsWith(info$type, 'nonlinear') ||
        startsWith(info$type, 'varying')
    if (surface) {
        if (is.null(predictor)) stop('A surface term requires predictor values')
        grid <- expand.grid(lag=lag, predictor=predictor)
        if (!is.null(info$axis)) {
            axis_data <- stats::setNames(lapply(seq_along(info$axis), function(index) {
                axis <- info$axis[[index]]
                if (index == 1L) return(grid$lag / axis$scale)
                if (index == 2L) return(grid$predictor / axis$scale)
                rep.int(stats::median(axis$grid), nrow(grid))
            }), vapply(info$axis, `[[`, character(1), 'internal'))
            basis <- mgcv::PredictMat(
                info$basis,
                axis_data,
                n=nrow(grid)
            ) %*% info$transform
        } else {
            basis <- mgcv::PredictMat(
                info$basis,
                list(
                    cdr_delay=grid$lag / lag_scale,
                    cdr_value=grid$predictor / predictor_scale
                ),
                n=nrow(grid)
            ) %*% info$transform
        }
    } else {
        grid <- data.frame(lag=lag, predictor=NA_real_)
        basis <- mgcv::PredictMat(
            info$basis,
            list(cdr_delta=lag / lag_scale),
            n=length(lag)
        )
        if (!is.null(info$transform)) basis <- basis %*% info$transform
    }
    basis <- basis / amplitude_scale
    list(grid=grid, basis=basis)
}

.cdr_spec_key <- function(specification, drop_group=FALSE) {
    if (is.null(specification)) return(NULL)
    specification <- specification[c(
        'predictors', 'constant', 'window', 'k_l', 'k_t', 'k_p',
        'bs_l', 'bs_t', 'bs_p', 'time', 'by', 'group'
    )]
    if (drop_group) specification$group <- NULL
    paste(utils::capture.output(dput(specification)), collapse='')
}

.cdr_matching_term <- function(object, index, kind=c('population', 'baseline')) {
    kind <- match.arg(kind)
    specifications <- object$cdrgam$preparation$specification
    if (is.null(specifications) || length(specifications) !=
            length(object$cdrgam$terms)) {
        return(integer())
    }
    selected <- specifications[[index]]
    if (identical(kind, 'baseline')) {
        if (!is.null(selected$k_t)) {
            baseline <- selected
            baseline$k_t <- NULL
            baseline$time <- NULL
            baseline$varying <- NULL
            baseline$group <- NULL
            candidates <- vapply(specifications, function(candidate) {
                is.null(candidate$group) && is.null(candidate$k_t) &&
                    identical(candidate$predictors, baseline$predictors) &&
                    identical(candidate$k_p, baseline$k_p) &&
                    identical(candidate$by, baseline$by)
            }, logical(1))
            return(which(candidates))
        }
        if (length(selected$predictors) == 1L &&
                any(!vapply(selected$k_p, is.null, logical(1)))) {
            candidates <- vapply(specifications, function(candidate) {
                isTRUE(candidate$constant) && is.null(candidate$group) &&
                    is.null(candidate$k_t) && is.null(candidate$by)
            }, logical(1))
            return(which(candidates))
        }
        return(integer())
    }
    key <- .cdr_spec_key(selected, drop_group=TRUE)
    candidates <- vapply(seq_along(specifications), function(i) {
        candidate <- specifications[[i]]
        is.null(candidate$group) &&
            identical(.cdr_spec_key(candidate, drop_group=TRUE), key)
    }, logical(1))
    which(candidates)
}

.cdr_plot_term_part <- function(object, index, basis, group_level=NULL) {
    info <- object$cdrgam$terms[[index]]
    indices <- info$coefficient_index
    if (length(info$group_levels)) {
        group_number <- match(group_level, info$group_levels)
        if (is.na(group_number)) {
            stop('Unknown grouped-IRF deviation level: ', group_level)
        }
        within <- (group_number - 1L) * info$base_dimension +
            seq_len(info$base_dimension)
        indices <- indices[within]
    }
    if (length(indices) != ncol(basis)) {
        stop('Stored IRF coefficient metadata is inconsistent with its basis')
    }
    list(index=indices, basis=basis)
}

.cdr_plot_combine_parts <- function(object, parts, se, unconditional) {
    indices <- sort(unique(unlist(lapply(parts, `[[`, 'index'))))
    design <- matrix(0, nrow=nrow(parts[[1L]]$basis), ncol=length(indices))
    for (part in parts) {
        destination <- match(part$index, indices)
        design[, destination] <- design[, destination, drop=FALSE] + part$basis
    }
    estimate <- drop(design %*% stats::coef(object)[indices])
    standard_error <- rep.int(NA_real_, nrow(design))
    if (isTRUE(se)) {
        covariance <- if (inherits(object, 'cdrgam_sparse')) {
            .sparse_selected_vcov(object, indices, unconditional=unconditional)
        } else {
            stats::vcov(object, unconditional=unconditional)[
                indices, indices, drop=FALSE
            ]
        }
        standard_error <- sqrt(pmax(
            0,
            rowSums((design %*% covariance) * design)
        ))
    }
    list(estimate=estimate, se=standard_error)
}

.cdr_plot_evaluate <- function(
        object,
        index,
        lag,
        predictor,
        component,
        group_levels,
        se,
        unconditional
) {
    info <- object$cdrgam$terms[[index]]
    evaluated <- .cdr_plot_basis(info, lag, predictor)
    grouped <- length(info$group_levels) > 0L
    if (grouped && identical(component, 'term')) component <- 'deviation'
    if (!grouped && component %in% c('population', 'conditional', 'deviation')) {
        if (identical(component, 'population')) component <- 'term' else {
            stop('component="', component, '" requires a grouped IRF term')
        }
    }
    levels <- if (grouped && component %in% c(
        'deviation', 'conditional', 'total'
    )) {
        if (is.null(group_levels)) info$group_levels else as.character(group_levels)
    } else {
        NA_character_
    }
    population_index <- if (grouped && component %in% c(
        'population', 'conditional', 'total'
    )) .cdr_matching_term(object, index, 'population') else integer()
    if (grouped && component %in% c('population', 'conditional', 'total') &&
            length(population_index) != 1L) {
        stop('A unique population IRF could not be found for ',
            object$cdrgam$term_labels[[index]])
    }
    baseline_source <- if (length(population_index)) population_index else index
    baseline_index <- if (identical(component, 'total')) {
        .cdr_matching_term(object, baseline_source, 'baseline')
    } else integer()
    if (identical(component, 'total') &&
            ((!is.null(info$axis) && length(info$axis) > 1L) ||
             startsWith(info$type, 'nonlinear') || startsWith(info$type, 'varying')) &&
            length(baseline_index) != 1L) {
        stop('A unique baseline IRF could not be found for ',
            object$cdrgam$term_labels[[index]])
    }
    output <- vector('list', length(levels))
    for (j in seq_along(levels)) {
        group_level <- levels[[j]]
        parts <- list()
        if (!grouped || component %in% c('term', 'deviation')) {
            parts[[length(parts) + 1L]] <- .cdr_plot_term_part(
                object, index, evaluated$basis,
                if (grouped) group_level else NULL
            )
        }
        if (grouped && component %in% c('population', 'conditional', 'total')) {
            population_basis <- .cdr_plot_basis(
                object$cdrgam$terms[[population_index]], lag, predictor
            )$basis
            parts[[length(parts) + 1L]] <- .cdr_plot_term_part(
                object, population_index, population_basis
            )
        }
        if (grouped && component %in% c('conditional', 'total')) {
            parts[[length(parts) + 1L]] <- .cdr_plot_term_part(
                object, index, evaluated$basis, group_level
            )
        }
        if (length(baseline_index)) {
            baseline_basis <- .cdr_plot_basis(
                object$cdrgam$terms[[baseline_index]], lag
            )$basis
            if (nrow(baseline_basis) != nrow(evaluated$basis)) {
                repeats <- nrow(evaluated$basis) / nrow(baseline_basis)
                baseline_basis <- baseline_basis[rep(
                    seq_len(nrow(baseline_basis)), times=repeats
                ), , drop=FALSE]
            }
            parts[[length(parts) + 1L]] <- .cdr_plot_term_part(
                object, baseline_index, baseline_basis
            )
        }
        values <- .cdr_plot_combine_parts(object, parts, se, unconditional)
        output[[j]] <- transform(
            evaluated$grid,
            term=object$cdrgam$term_labels[[index]],
            group=if (is.na(group_level)) NA_character_ else group_level,
            component=component,
            estimate=values$estimate,
            se=values$se
        )[, c('term', 'group', 'component', 'lag', 'predictor', 'estimate', 'se')]
    }
    do.call(rbind, output)
}

.cdr_plot_coefficient_groups <- function(object) {
    groups <- list()
    add_group <- function(label, indices, kind) {
        if (!length(indices)) return(invisible(NULL))
        groups[[length(groups) + 1L]] <<- list(
            label=label, indices=as.integer(indices), kind=kind
        )
    }
    for (i in seq_along(object$smooth)) {
        smooth <- object$smooth[[i]]
        if (inherits(smooth, 'cdr.smooth')) next
        add_group(
            smooth$label,
            seq.int(smooth$first.para, smooth$last.para),
            if (inherits(smooth, 'random.effect')) 'random' else 'smooth'
        )
    }
    random_effects <- object$cdrgam$prediction$random_effects
    if (length(random_effects)) {
        for (effect in random_effects) {
            add_group(effect$label, effect$coefficient_index, 'random')
        }
    }
    for (i in seq_along(object$cdrgam$terms)) {
        add_group(
            object$cdrgam$term_labels[[i]],
            object$cdrgam$terms[[i]]$coefficient_index,
            'irf'
        )
    }
    owned <- unique(unlist(lapply(groups, `[[`, 'indices')))
    parametric <- setdiff(seq_along(stats::coef(object)), owned)
    if (length(parametric)) {
        groups <- c(list(list(
            label='parametric', indices=parametric, kind='parametric'
        )), groups)
    }
    groups
}

.cdr_plot_coefficient_panels <- function(
        object, select, at, se, unconditional, ci_level, coefficient_limit
) {
    groups <- .cdr_plot_coefficient_groups(object)
    labels <- vapply(groups, `[[`, character(1), 'label')
    explicit_selection <- !is.null(select)
    selected <- if (!explicit_selection) seq_along(groups) else {
        if (is.numeric(select)) as.integer(select) else match(select, labels)
    }
    if (!length(selected) || anyNA(selected) ||
            any(selected < 1L | selected > length(groups))) {
        stop('select must identify coefficient groups: ', paste(labels, collapse=', '))
    }
    covariance <- if (isTRUE(se) && !inherits(object, 'cdrgam_sparse')) {
        stats::vcov(object, unconditional=unconditional)
    } else NULL
    coefficient_names <- names(stats::coef(object))
    if (is.null(coefficient_names)) {
        coefficient_names <- paste0('coefficient_', seq_along(stats::coef(object)))
    }
    coefficient_divisors <- .cdr_coefficient_divisors(object)
    panels <- list()
    for (i in selected) {
        group <- groups[[i]]
        indices <- group$indices
        requested <- at$coefficient
        if (!is.null(requested)) {
            keep <- coefficient_names[indices] %in% requested
            indices <- indices[keep]
        }
        if (!length(indices)) next
        if (length(indices) > coefficient_limit) {
            if (!explicit_selection && is.null(requested)) {
                warning(
                    'Omitting coefficient group ', group$label, ' with ',
                    length(indices), ' entries; select it explicitly and ',
                    'increase coefficient_limit to draw it',
                    call.=FALSE
                )
                next
            }
            stop(
                'Coefficient group ', group$label, ' has ', length(indices),
                ' entries; select at$coefficient or increase coefficient_limit'
            )
        }
        standard_error <- rep.int(NA_real_, length(indices))
        if (isTRUE(se)) {
            standard_error <- if (inherits(object, 'cdrgam_sparse')) {
                sqrt(pmax(0, diag(.sparse_selected_vcov(
                    object, indices, unconditional=unconditional
                ))))
            } else sqrt(pmax(0, diag(covariance[indices, indices, drop=FALSE])))
        }
        panels[[length(panels) + 1L]] <- list(
            view='coef',
            label=group$label,
            kind=group$kind,
            data=data.frame(
                coefficient=coefficient_names[indices],
                estimate=unname(
                    stats::coef(object)[indices] / coefficient_divisors[indices]
                ),
                se=standard_error / abs(coefficient_divisors[indices]),
                stringsAsFactors=FALSE
            ),
            ci_level=ci_level
        )
    }
    panels
}

.cdr_plot_panels <- function(
        object, view, select, at, component, se, unconditional,
        n, n_predictor, ci_level, coefficient_limit
) {
    if (identical(view, 'coef')) {
        return(.cdr_plot_coefficient_panels(
            object, select, at, se, unconditional, ci_level, coefficient_limit
        ))
    }
    registry <- .cdr_plot_registry(object)
    indices <- .cdr_plot_select(registry, select)
    panels <- list()
    for (index in indices) {
        row <- registry[registry$index == index, ]
        panel_view <- if (identical(view, 'auto')) {
            if (row$surface) 'surface' else 'irf'
        } else view
        if (identical(panel_view, 'surface') && !row$surface) {
            stop('view="surface" requires a nonlinear or varying IRF')
        }
        lag <- at$lag
        predictor <- at$predictor
        if (identical(panel_view, 'irf')) {
            if (is.null(lag)) {
                knots <- object$cdrgam$terms[[index]]$knots
                scale <- object$cdrgam$terms[[index]]$lag_scale
                if (is.null(scale)) scale <- 1
                lag <- seq(min(knots), max(knots), length.out=n) * scale
            }
            if (row$surface && is.null(predictor)) {
                knots <- object$cdrgam$terms[[index]]$predictor_knots
                scale <- object$cdrgam$terms[[index]]$predictor_scale
                if (is.null(scale)) scale <- 1
                predictor <- seq(min(knots), max(knots), length.out=3L) * scale
            }
            if (!row$surface && !is.null(predictor)) {
                base_predictor <- predictor
                predictor <- NULL
            } else base_predictor <- NULL
        } else if (identical(panel_view, 'predictor')) {
            if (is.null(lag)) {
                knots <- object$cdrgam$terms[[index]]$knots
                scale <- object$cdrgam$terms[[index]]$lag_scale
                if (is.null(scale)) scale <- 1
                lag <- seq(min(knots), max(knots), length.out=3L) * scale
            }
            if (is.null(predictor)) {
                knots <- object$cdrgam$terms[[index]]$predictor_knots
                if (is.null(knots)) predictor <- seq(-2, 2, length.out=n_predictor)
                else {
                    scale <- object$cdrgam$terms[[index]]$predictor_scale
                    if (is.null(scale)) scale <- 1
                    predictor <- seq(
                        min(knots), max(knots), length.out=n_predictor
                    ) * scale
                }
            }
            base_predictor <- if (row$surface) NULL else predictor
            if (!row$surface) predictor <- NULL
        } else if (identical(panel_view, 'surface')) {
            if (is.null(lag)) {
                knots <- object$cdrgam$terms[[index]]$knots
                scale <- object$cdrgam$terms[[index]]$lag_scale
                if (is.null(scale)) scale <- 1
                lag <- seq(min(knots), max(knots), length.out=n) * scale
            }
            if (is.null(predictor)) {
                knots <- object$cdrgam$terms[[index]]$predictor_knots
                scale <- object$cdrgam$terms[[index]]$predictor_scale
                if (is.null(scale)) scale <- 1
                predictor <- seq(
                    min(knots), max(knots), length.out=n_predictor
                ) * scale
            }
            base_predictor <- NULL
        } else stop('Unsupported CDR plot view: ', panel_view)
        if (!is.numeric(lag) || any(!is.finite(lag)) || !length(lag) ||
                (!is.null(predictor) && (!is.numeric(predictor) ||
                    any(!is.finite(predictor)) || !length(predictor)))) {
            stop('at$lag and at$predictor must be finite numeric vectors')
        }
        values <- .cdr_plot_evaluate(
            object, index, lag, predictor, component, at$group,
            se, unconditional
        )
        if (!is.null(base_predictor)) {
            expanded <- merge(
                values,
                data.frame(plot_predictor=base_predictor),
                by=NULL
            )
            expanded$estimate <- expanded$estimate * expanded$plot_predictor
            expanded$se <- expanded$se * abs(expanded$plot_predictor)
            expanded$predictor <- expanded$plot_predictor
            expanded$plot_predictor <- NULL
            values <- expanded
        }
        if (identical(panel_view, 'surface') && any(!is.na(values$group))) {
            for (group_level in unique(values$group)) {
                panels[[length(panels) + 1L]] <- list(
                    view=panel_view,
                    label=paste0(row$label, ': ', group_level),
                    data=values[values$group == group_level, ],
                    ci_level=ci_level
                )
            }
        } else {
            panels[[length(panels) + 1L]] <- list(
                view=panel_view,
                label=row$label,
                data=values,
                ci_level=ci_level
            )
        }
    }
    panels
}

.cdr_plot_curve <- function(panel, xlab, ylab, main, col, lwd, xlim, ylim) {
    data <- panel$data
    x_name <- if (identical(panel$view, 'predictor')) 'predictor' else 'lag'
    series_names <- if (identical(panel$view, 'predictor')) c('lag', 'group') else
        c('predictor', 'group')
    series_fields <- lapply(series_names, function(name) {
        value <- data[[name]]
        formatted <- format(value, digits=4, trim=TRUE)
        ifelse(is.na(value), '', paste0(name, '=', formatted))
    })
    series_labels <- apply(
        as.data.frame(series_fields, stringsAsFactors=FALSE),
        1L,
        function(value) paste(value[nzchar(value)], collapse=' | ')
    )
    series_labels[!nzchar(series_labels)] <- 'term'
    series_values <- factor(series_labels, levels=unique(series_labels))
    series <- split(seq_len(nrow(data)), series_values)
    colors <- if (is.null(col)) grDevices::hcl.colors(length(series), 'Dark 3') else
        rep(col, length.out=length(series))
    multiplier <- stats::qnorm((1 + panel$ci_level) / 2)
    limits <- data$estimate
    if (any(is.finite(data$se))) {
        limits <- c(limits, data$estimate - multiplier * data$se,
            data$estimate + multiplier * data$se)
    }
    if (is.null(xlim)) xlim <- range(data[[x_name]])
    if (is.null(ylim)) ylim <- range(limits[is.finite(limits)])
    graphics::plot(
        xlim, ylim, type='n',
        xlab=if (is.null(xlab)) {
            if (identical(x_name, 'lag')) 'Lag' else 'Predictor value'
        } else xlab,
        ylab=if (is.null(ylab)) 'Effect on linear predictor' else ylab,
        main=if (is.null(main)) panel$label else main
    )
    graphics::abline(h=0, col='grey80', lty=3)
    for (i in seq_along(series)) {
        rows <- series[[i]]
        rows <- rows[order(data[[x_name]][rows])]
        if (any(is.finite(data$se[rows]))) {
            graphics::polygon(
                c(data[[x_name]][rows], rev(data[[x_name]][rows])),
                c(
                    data$estimate[rows] - multiplier * data$se[rows],
                    rev(data$estimate[rows] + multiplier * data$se[rows])
                ),
                col=grDevices::adjustcolor(colors[[i]], alpha.f=0.18),
                border=NA
            )
        }
        graphics::lines(
            data[[x_name]][rows], data$estimate[rows],
            col=colors[[i]], lwd=lwd
        )
    }
    if (length(series) > 1L) {
        graphics::legend('topright', legend=names(series), col=colors, lwd=lwd,
            bty='n', cex=0.8)
    }
}

.cdr_plot_surface <- function(panel, scheme, xlab, ylab, main, col, ...) {
    data <- panel$data
    lag <- sort(unique(data$lag))
    predictor <- sort(unique(data$predictor))
    z <- matrix(NA_real_, nrow=length(lag), ncol=length(predictor))
    z[cbind(match(data$lag, lag), match(data$predictor, predictor))] <-
        data$estimate
    palette <- if (is.null(col)) grDevices::hcl.colors(64, 'Blue-Red 3') else col
    title <- if (is.null(main)) panel$label else main
    if (identical(scheme, 'persp')) {
        graphics::persp(
            lag, predictor, z,
            xlab=if (is.null(xlab)) 'Lag' else xlab,
            ylab=if (is.null(ylab)) 'Predictor value' else ylab,
            zlab='Effect', main=title, col=palette[[length(palette) %/% 2L]],
            ticktype='detailed', ...
        )
    } else if (identical(scheme, 'image')) {
        z_limit <- max(abs(z), na.rm=TRUE)
        if (!is.finite(z_limit) || z_limit == 0) z_limit <- 1
        graphics::image(
            lag, predictor, z,
            xlab=if (is.null(xlab)) 'Lag' else xlab,
            ylab=if (is.null(ylab)) 'Predictor value' else ylab,
            main=title, col=palette, zlim=c(-z_limit, z_limit), ...
        )
        graphics::contour(lag, predictor, z, add=TRUE, drawlabels=FALSE)
    } else {
        graphics::contour(
            lag, predictor, z,
            xlab=if (is.null(xlab)) 'Lag' else xlab,
            ylab=if (is.null(ylab)) 'Predictor value' else ylab,
            main=title, col=if (is.null(col)) 'black' else col[[1L]], ...
        )
    }
}

.cdr_plot_coefficients <- function(panel, xlab, main, col, ci_level) {
    data <- panel$data
    multiplier <- stats::qnorm((1 + ci_level) / 2)
    lower <- data$estimate - multiplier * data$se
    upper <- data$estimate + multiplier * data$se
    finite_interval <- is.finite(lower) & is.finite(upper)
    limits <- range(c(data$estimate, lower[finite_interval], upper[finite_interval]))
    graphics::plot(
        data$estimate, seq_len(nrow(data)),
        xlim=limits, yaxt='n', pch=19,
        xlab=if (is.null(xlab)) 'Coefficient estimate' else xlab,
        ylab='', main=if (is.null(main)) panel$label else main,
        col=if (is.null(col)) 'black' else col[[1L]]
    )
    graphics::axis(2, at=seq_len(nrow(data)), labels=data$coefficient,
        las=2, cex.axis=0.65)
    graphics::abline(v=0, col='grey80', lty=3)
    if (any(finite_interval)) graphics::segments(
        lower[finite_interval], which(finite_interval),
        upper[finite_interval], which(finite_interval),
        col=if (is.null(col)) 'black' else col[[1L]]
    )
}

#' Plot CDR effects and translated GAM components
#'
#' Evaluate and plot impulse-response functions by lag, predictor slices at
#' fixed delays, lag-by-predictor surfaces, or coefficient groups. All CDR
#' effects are shown on the additive linear-predictor scale. Use `view="gam"`
#' on a native `mgcv` fit to call [mgcv::plot.gam()] on the translated model.
#' When predictor rescaling was enabled, CDR and ordinary-smooth coordinates
#' are restored to native source units before display.
#'
#' @param x A fitted `cdrgam` model.
#' @param view Plot type: automatic CDR view, lag IRF, predictor slice,
#'   lag-by-predictor surface, translated GAM, or coefficients.
#' @param select IRF or coefficient-group labels or one-based indices. Native
#'   GAM views pass this argument to [mgcv::plot.gam()].
#' @param at Named conditioning values. Supported entries are `lag`,
#'   `predictor`, `group`, and `coefficient`.
#' @param component Effect composition. `"term"` shows the selected term.
#'   Grouped terms also support population, deviation, and conditional effects.
#'   `"total"` adds the corresponding baseline to tensor terms and adds the
#'   population effect to grouped deviations.
#' @param se Draw or return pointwise standard errors.
#' @param unconditional Include smoothing-parameter uncertainty where the
#'   fitted backend supports it.
#' @param ci_level Confidence level used for plotted intervals.
#' @param n Number of lag evaluation points.
#' @param n_predictor Number of predictor evaluation points.
#' @param surface Surface rendering style.
#' @param pages Number of pages. Zero places all panels on one page; a positive
#'   value distributes panels across that many pages.
#' @param ask Pause before advancing graphical pages on interactive devices.
#' @param draw Draw the panels. Set to false to return evaluated data only.
#' @param coefficient_limit Maximum coefficients allowed in one panel unless
#'   `at$coefficient` selects a smaller set. Oversized groups are omitted from
#'   the default coefficient view and error when explicitly selected.
#' @param xlab,ylab,main,col,lwd,xlim,ylim Graphical controls.
#' @param ... Additional arguments passed to surface plotting functions or to
#'   [mgcv::plot.gam()] for `view="gam"`.
#' @return Invisibly, a `cdrgam_plot_data` list containing evaluated panels.
#' @export
plot.cdrgam <- function(
        x,
        view=c('auto', 'irf', 'predictor', 'surface', 'gam', 'coef'),
        select=NULL,
        at=list(),
        component=c('term', 'total', 'population', 'deviation', 'conditional'),
        se=TRUE,
        unconditional=FALSE,
        ci_level=0.95,
        n=200,
        n_predictor=50,
        surface=c('persp', 'image', 'contour'),
        pages=0,
        ask=FALSE,
        draw=TRUE,
        coefficient_limit=200L,
        xlab=NULL,
        ylab=NULL,
        main=NULL,
        col=NULL,
        lwd=2,
        xlim=NULL,
        ylim=NULL,
        ...
) {
    if (!is_cdrgam(x)) stop('x must be a fitted cdrgam model')
    view <- match.arg(view)
    component <- match.arg(component)
    surface <- match.arg(surface)
    if (!is.list(at) || (length(at) && is.null(names(at)))) {
        stop('at must be a named list')
    }
    unknown_at <- setdiff(names(at), c('lag', 'predictor', 'group', 'coefficient'))
    if (length(unknown_at)) stop('Unknown at entries: ', paste(unknown_at, collapse=', '))
    if (length(ci_level) != 1L || !is.finite(ci_level) ||
            ci_level <= 0 || ci_level >= 1) stop('ci_level must lie between zero and one')
    numeric_controls <- c(n=n, n_predictor=n_predictor,
        coefficient_limit=coefficient_limit)
    if (any(!is.finite(numeric_controls)) || any(numeric_controls < 1) ||
            any(numeric_controls != as.integer(numeric_controls))) {
        stop('n, n_predictor, and coefficient_limit must be positive integers')
    }
    if (length(pages) != 1L || !is.numeric(pages) || !is.finite(pages) ||
            pages < 0 || pages != as.integer(pages)) {
        stop('pages must be a non-negative integer')
    }
    if (identical(view, 'gam')) {
        if (!inherits(x, 'gam')) {
            stop('view="gam" is available only for the native mgcv backend')
        }
        gam <- .cdr_source_scale_gam_plot(x)
        gam_select <- select
        if (is.character(gam_select)) {
            smooth_labels <- vapply(gam$smooth, `[[`, character(1), 'label')
            gam_select <- match(gam_select, smooth_labels)
            if (anyNA(gam_select)) {
                stop('Unknown translated GAM smooths: ',
                    paste(select[is.na(gam_select)], collapse=', '))
            }
        }
        result <- mgcv::plot.gam(
            gam, select=gam_select, pages=pages, ask=ask, se=se,
            unconditional=unconditional, xlab=xlab, ylab=ylab, main=main,
            xlim=xlim, ylim=ylim, ...
        )
        return(invisible(result))
    }
    panels <- .cdr_plot_panels(
        x, view, select, at, component, se, unconditional,
        as.integer(n), as.integer(n_predictor), ci_level,
        as.integer(coefficient_limit)
    )
    output <- structure(list(
        view=view,
        component=component,
        panels=panels
    ), class='cdrgam_plot_data')
    if (!isTRUE(draw) || !length(panels)) return(invisible(output))
    page_count <- if (pages <= 0L) 1L else min(as.integer(pages), length(panels))
    panels_per_page <- ceiling(length(panels) / page_count)
    rows <- floor(sqrt(panels_per_page))
    columns <- ceiling(panels_per_page / rows)
    old <- graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(old), add=TRUE)
    graphics::par(mfrow=c(rows, columns))
    old_ask <- grDevices::devAskNewPage(ask=ask)
    on.exit(grDevices::devAskNewPage(old_ask), add=TRUE)
    for (panel in panels) {
        if (panel$view %in% c('irf', 'predictor')) {
            .cdr_plot_curve(panel, xlab, ylab, main, col, lwd, xlim, ylim)
        } else if (identical(panel$view, 'surface')) {
            .cdr_plot_surface(panel, surface, xlab, ylab, main, col, ...)
        } else {
            .cdr_plot_coefficients(panel, xlab, main, col, ci_level)
        }
    }
    invisible(output)
}

#' Save CDR plots to PDF or raster images
#'
#' Open a graphics device, call [plot.cdrgam()], and close the device. PDF
#' output may contain multiple pages. Raster output inserts a numbered format
#' field before the extension when more than one page is requested.
#'
#' @param object A fitted `cdrgam` model.
#' @param file Output path.
#' @param device Output device. `"auto"` uses the file extension.
#' @param width,height Device dimensions in inches.
#' @param dpi Raster resolution.
#' @param pages Number of output pages passed to [plot.cdrgam()].
#' @param ... Additional arguments passed to [plot.cdrgam()].
#' @return Invisibly, the plot-data object returned by [plot.cdrgam()].
#' @export
save_cdrgam_plots <- function(
        object,
        file,
        device=c('auto', 'pdf', 'png', 'jpeg', 'tiff'),
        width=8,
        height=6,
        dpi=144,
        pages=1,
        ...
) {
    device <- match.arg(device)
    if (!is.character(file) || length(file) != 1L || !nzchar(file)) {
        stop('file must be one non-empty path')
    }
    if (identical(device, 'auto')) {
        extension <- tolower(tools::file_ext(file))
        device <- if (extension %in% c('jpg', 'jpeg')) 'jpeg' else extension
        if (!(device %in% c('pdf', 'png', 'jpeg', 'tiff'))) {
            stop('Cannot infer a supported device from the file extension')
        }
    }
    output_file <- file
    if (!identical(device, 'pdf') && pages > 1L && !grepl('%', file, fixed=TRUE)) {
        extension <- tools::file_ext(file)
        stem <- substr(file, 1L, nchar(file) - nchar(extension) - 1L)
        output_file <- paste0(stem, '-%03d.', extension)
    }
    if (identical(device, 'pdf')) {
        grDevices::pdf(output_file, width=width, height=height, onefile=TRUE)
    } else {
        opener <- getExportedValue('grDevices', device)
        opener(
            output_file,
            width=width * dpi,
            height=height * dpi,
            res=dpi
        )
    }
    on.exit(grDevices::dev.off(), add=TRUE)
    plot(object, pages=pages, ...)
}
