.cdr_empty_scaling <- function(enabled=FALSE) {
    list(
        enabled=isTRUE(enabled),
        time_divisor=1,
        variables=data.frame(
            stream=character(),
            variable=character(),
            roles=character(),
            divisor=numeric(),
            applied=logical(),
            reason=character(),
            stringsAsFactors=FALSE
        )
    )
}

.cdr_scaling_roles <- function(parsed, response_name, series) {
    impulse_roles <- list()
    response_roles <- list()
    add_role <- function(registry, variable, role) {
        if (is.null(variable) || !nzchar(variable)) return(registry)
        registry[[variable]] <- unique(c(registry[[variable]], role))
        registry
    }
    for (specification in parsed$irfs) {
        for (predictor in specification$predictors) {
            impulse_roles <- add_role(
                impulse_roles,
                predictor,
                'irf predictor'
            )
        }
        response_roles <- add_role(
            response_roles,
            specification$time,
            'IRF time axis'
        )
        response_roles <- add_role(
            response_roles,
            specification$by,
            'IRF by variable'
        )
        response_roles <- add_role(
            response_roles,
            specification$group,
            'IRF grouping variable'
        )
    }
    ordinary_variables <- setdiff(
        all.vars(parsed$ordinary_formula[[3L]]),
        response_name
    )
    for (variable in ordinary_variables) {
        response_roles <- add_role(
            response_roles,
            variable,
            'ordinary formula predictor'
        )
    }
    for (variable in series) {
        impulse_roles <- add_role(impulse_roles, variable, 'series identifier')
        response_roles <- add_role(response_roles, variable, 'series identifier')
    }
    list(impulses=impulse_roles, responses=response_roles)
}

.cdr_prepare_scaling <- function(
        parsed,
        impulses,
        responses,
        response_name,
        series,
        impulse_time,
        response_time,
        enabled
) {
    if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled)) {
        stop('rescale_predictors must be TRUE or FALSE')
    }
    if (!enabled) return(.cdr_empty_scaling(FALSE))
    impulse_times <- impulses[[impulse_time]]
    if (!is.numeric(impulse_times) || any(!is.finite(impulse_times))) {
        stop('impulse_time must name a finite numeric column in impulses')
    }
    time_divisor <- stats::sd(impulse_times)
    if (!is.finite(time_divisor) || time_divisor <= 0) {
        stop(
            'Cannot rescale time because the training impulse-time column ',
            'has zero or non-finite standard deviation',
            call.=FALSE
        )
    }
    roles <- .cdr_scaling_roles(parsed, response_name, series)
    grouping <- unique(c(
        series,
        vapply(parsed$irfs, function(specification) {
            if (is.null(specification$group)) NA_character_ else
                specification$group
        }, character(1))
    ))
    grouping <- grouping[!is.na(grouping)]
    make_rows <- function(data, registry, stream, time_name) {
        if (!length(registry)) return(NULL)
        lapply(names(registry), function(variable) {
            if (!(variable %in% names(data))) {
                stop(
                    'Scaling variable ', sQuote(variable), ' is missing from ',
                    stream,
                    call.=FALSE
                )
            }
            value <- data[[variable]]
            variable_roles <- registry[[variable]]
            applied <- FALSE
            divisor <- 1
            reason <- NULL
            if (variable %in% grouping) {
                reason <- 'grouping or series identifier'
            } else if (!is.numeric(value)) {
                reason <- 'non-numeric predictor'
            } else if (any(!is.finite(value))) {
                stop(
                    'Numeric predictor ', sQuote(variable), ' in ', stream,
                    ' must contain only finite values',
                    call.=FALSE
                )
            } else if (identical(variable, time_name)) {
                divisor <- time_divisor
                applied <- TRUE
                reason <- 'time coordinate'
            } else if (length(unique(value)) <= 2L) {
                reason <- if (length(unique(value)) == 1L) {
                    'constant numeric predictor'
                } else {
                    'binary numeric predictor'
                }
            } else {
                divisor <- stats::sd(value)
                if (!is.finite(divisor) || divisor <= 0) {
                    reason <- 'zero or non-finite standard deviation'
                    divisor <- 1
                } else {
                    applied <- TRUE
                    reason <- 'continuous numeric predictor'
                }
            }
            data.frame(
                stream=stream,
                variable=variable,
                roles=paste(variable_roles, collapse=', '),
                divisor=as.numeric(divisor),
                applied=applied,
                reason=reason,
                stringsAsFactors=FALSE
            )
        })
    }
    rows <- c(
        make_rows(impulses, roles$impulses, 'impulses', impulse_time),
        make_rows(responses, roles$responses, 'responses', response_time)
    )
    variables <- if (length(rows)) do.call(rbind, rows) else
        .cdr_empty_scaling(TRUE)$variables
    list(
        enabled=TRUE,
        time_divisor=as.numeric(time_divisor),
        variables=variables
    )
}

.cdr_apply_scaling <- function(data, scaling, stream) {
    if (is.null(scaling) || !isTRUE(scaling$enabled) ||
            !nrow(scaling$variables)) {
        return(data)
    }
    rows <- scaling$variables$stream == stream & scaling$variables$applied
    transformations <- scaling$variables[rows, , drop=FALSE]
    if (!nrow(transformations)) return(data)
    output <- data
    for (i in seq_len(nrow(transformations))) {
        variable <- transformations$variable[[i]]
        if (!(variable %in% names(output))) {
            stop(
                'New ', stream, ' data are missing scaled predictor: ',
                variable,
                call.=FALSE
            )
        }
        value <- output[[variable]]
        if (!is.numeric(value) || any(!is.finite(value))) {
            stop(
                'Scaled predictor ', sQuote(variable), ' in new ', stream,
                ' data must contain only finite numeric values',
                call.=FALSE
            )
        }
        output[[variable]] <- value / transformations$divisor[[i]]
    }
    output
}

.cdr_scaling_divisor <- function(scaling, stream, variable) {
    if (is.null(variable) || is.null(scaling) ||
            !isTRUE(scaling$enabled) || !nrow(scaling$variables)) {
        return(1)
    }
    rows <- scaling$variables$stream == stream &
        scaling$variables$variable == variable &
        scaling$variables$applied
    if (!any(rows)) 1 else scaling$variables$divisor[which(rows)[[1L]]]
}

.cdr_irf_scales <- function(specification, scaling) {
    if (is.null(scaling) || !isTRUE(scaling$enabled)) {
        return(list(
            lag=1,
            predictor=1,
            amplitude=1,
            axes=stats::setNames(rep.int(1, length(specification$predictors)),
                specification$predictors),
            time=1
        ))
    }
    impulse_divisors <- stats::setNames(vapply(
        specification$predictors,
        function(predictor) .cdr_scaling_divisor(
            scaling,
            'impulses',
            predictor
        ),
        numeric(1)
    ), specification$predictors)
    smooth <- !vapply(specification$k_p, is.null, logical(1))
    linear_divisors <- impulse_divisors[!smooth]
    predictor_divisor <- if (any(smooth)) {
        impulse_divisors[which(smooth)[[1L]]]
    } else if (!is.null(specification$k_t)) {
        .cdr_scaling_divisor(scaling, 'responses', specification$time)
    } else 1
    amplitude_divisor <- prod(linear_divisors)
    amplitude_divisor <- amplitude_divisor * .cdr_scaling_divisor(
        scaling,
        'responses',
        specification$by
    )
    list(
        lag=scaling$time_divisor,
        predictor=predictor_divisor,
        amplitude=amplitude_divisor,
        axes=impulse_divisors,
        time=.cdr_scaling_divisor(scaling, 'responses', specification$time)
    )
}

.cdr_expression_divisor <- function(expression, scaling) {
    if (is.symbol(expression)) {
        return(.cdr_scaling_divisor(
            scaling,
            'responses',
            as.character(expression)
        ))
    }
    if (!is.call(expression)) return(1)
    operator <- as.character(expression[[1L]])
    if (identical(operator, 'I') && length(expression) == 2L) {
        return(.cdr_expression_divisor(expression[[2L]], scaling))
    }
    if (operator %in% c(':', '*') && length(expression) == 3L) {
        return(
            .cdr_expression_divisor(expression[[2L]], scaling) *
                .cdr_expression_divisor(expression[[3L]], scaling)
        )
    }
    if (identical(operator, '^') && length(expression) == 3L &&
            is.numeric(expression[[3L]]) &&
            length(expression[[3L]]) == 1L) {
        return(
            .cdr_expression_divisor(expression[[2L]], scaling) ^
                as.numeric(expression[[3L]])
        )
    }
    # General transformations need not be homogeneous in their input units.
    # Their displayed coefficient therefore remains in fitted-coordinate
    # units rather than claiming an invalid multiplicative conversion.
    1
}

.cdr_coefficient_divisors <- function(object) {
    coefficients <- stats::coef(object)
    divisors <- rep.int(1, length(coefficients))
    names(divisors) <- names(coefficients)
    scaling <- object$cdrgam$scaling
    if (is.null(scaling)) scaling <- object$cdrgam$preparation$scaling
    if (is.null(scaling) || !isTRUE(scaling$enabled)) return(divisors)

    setup <- object$cdrgam$prediction$setup
    if (is.null(setup) && inherits(object, 'gam')) {
        setup <- .cdr_prediction_setup(object)
    }
    if (!is.null(setup$pterms) && !is.null(setup$nsdf) && setup$nsdf > 0L) {
        assignment <- setup$assign
        if (is.null(assignment) && !is.null(object$assign)) {
            assignment <- object$assign
        }
        labels <- attr(setup$pterms, 'term.labels')
        count <- min(as.integer(setup$nsdf), length(divisors))
        if (!is.null(assignment) && length(assignment) >= count) {
            for (i in seq_len(count)) {
                term_index <- assignment[[i]]
                if (is.na(term_index) || term_index < 1L ||
                        term_index > length(labels)) next
                expression <- tryCatch(
                    parse(text=labels[[term_index]], keep.source=FALSE)[[1L]],
                    error=function(e) NULL
                )
                if (!is.null(expression)) {
                    divisors[[i]] <- .cdr_expression_divisor(
                        expression,
                        scaling
                    )
                }
            }
        }
    }
    for (term in object$cdrgam$terms) {
        indices <- term$coefficient_index
        indices <- indices[indices >= 1L & indices <= length(divisors)]
        if (length(indices)) {
            amplitude <- if (is.null(term$amplitude_scale)) {
                1
            } else term$amplitude_scale
            divisors[indices] <- divisors[indices] * amplitude
        }
    }
    divisors
}
