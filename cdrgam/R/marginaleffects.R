.cdr_marginaleffects_enable <- function() {
    classes <- getOption('marginaleffects_model_classes')
    if (is.null(classes)) classes <- character()
    if (!('cdrgam_marginaleffects' %in% classes)) {
        options(marginaleffects_model_classes=c(
            classes, 'cdrgam_marginaleffects'
        ))
    }
    invisible(NULL)
}

.cdr_marginaleffects_specs <- function(object, terms) {
    specifications <- object$cdrgam$preparation$specification
    metadata <- object$cdrgam$terms
    if (is.null(specifications) || length(specifications) != length(metadata)) {
        stop('The fitted model does not contain complete IRF metadata')
    }
    indices <- .cdr_effect_select(object, terms)
    if (!length(indices)) {
        stop('No impulse-response terms were selected')
    }
    list(indices=indices, specifications=specifications, metadata=metadata)
}

.cdr_marginaleffects_pairs <- function(
        object, impulses, responses, term_indices, specifications
) {
    stream <- object$cdrgam$preparation$stream
    if (is.null(stream)) {
        stop('The fitted model does not contain stream metadata')
    }
    term_pairs <- vector('list', length(term_indices))
    names(term_pairs) <- as.character(term_indices)
    pair_keys <- character()
    for (position in seq_along(term_indices)) {
        index <- term_indices[[position]]
        links <- .build_history_links(
            impulses,
            responses,
            stream$series,
            stream$impulse_time,
            stream$response_time,
            specifications[[index]]$window
        )
        keys <- paste(links$response_index, links$impulse_index, sep='\034')
        term_pairs[[position]] <- keys
        pair_keys <- c(pair_keys, keys)
    }
    pair_keys <- unique(pair_keys)
    if (!length(pair_keys)) {
        stop('No impulse-response links occur in the selected term windows')
    }
    pieces <- strsplit(pair_keys, '\034', fixed=TRUE)
    response_index <- as.integer(vapply(pieces, `[[`, character(1), 1L))
    impulse_index <- as.integer(vapply(pieces, `[[`, character(1), 2L))
    membership <- lapply(term_pairs, match, table=pair_keys)
    membership <- lapply(membership, function(index) {
        index[!is.na(index)]
    })
    list(
        key=pair_keys,
        response_index=response_index,
        impulse_index=impulse_index,
        membership=membership
    )
}

.cdr_marginaleffects_event_design <- function(context, newdata, pair_map) {
    coefficient_count <- context$coefficient_count
    row_index <- integer()
    column_index <- integer()
    values <- numeric()
    response_index <- context$pair_response[pair_map]
    for (position in seq_along(context$term_indices)) {
        term_index <- context$term_indices[[position]]
        window <- context$specifications[[term_index]]$window
        selected <- which(
            newdata$.cdrgam_lag >= window[[1L]] &
                newdata$.cdrgam_lag <= window[[2L]]
        )
        if (!length(selected)) next
        info <- context$terms[[term_index]]
        specification <- context$specifications[[term_index]]
        grid <- newdata[selected, , drop=FALSE]
        grid$lag <- grid$.cdrgam_lag
        term_response <- response_index[selected]
        if (!is.null(info$axis)) {
            time_axes <- vapply(
                info$axis,
                function(axis) identical(axis$role, 'time'),
                logical(1)
            )
            for (axis in info$axis[time_axes]) {
                grid[[axis$variable]] <- context$responses[[axis$variable]][
                    term_response
                ]
            }
        }
        basis_info <- info
        basis_info$linear_predictors <- .cdr_effect_linear_predictors(
            info, specification
        )
        basis <- .cdr_effect_basis(basis_info, grid)
        if (!is.null(specification$by)) {
            basis <- basis * context$responses[[specification$by]][term_response]
        }
        if (length(info$group_levels)) {
            group <- as.character(
                context$responses[[specification$group]][term_response]
            )
            group_number <- match(group, info$group_levels)
            for (level in unique(group_number[!is.na(group_number)])) {
                rows <- which(group_number == level)
                within <- (level - 1L) * info$base_dimension +
                    seq_len(info$base_dimension)
                columns <- info$coefficient_index[within]
                block <- basis[rows, , drop=FALSE]
                row_index <- c(
                    row_index,
                    rep(selected[rows], times=ncol(block))
                )
                column_index <- c(
                    column_index,
                    rep(columns, each=nrow(block))
                )
                values <- c(values, as.vector(block))
            }
        } else {
            columns <- info$coefficient_index
            if (length(columns) != ncol(basis)) {
                stop('Stored IRF coefficient metadata is inconsistent')
            }
            row_index <- c(
                row_index,
                rep(selected, times=ncol(basis))
            )
            column_index <- c(
                column_index,
                rep(columns, each=nrow(basis))
            )
            values <- c(values, as.vector(basis))
        }
    }
    keep <- is.finite(values) & values != 0
    Matrix::sparseMatrix(
        i=row_index[keep],
        j=column_index[keep],
        x=values[keep],
        dims=c(nrow(newdata), coefficient_count)
    )
}

.cdr_marginaleffects_response_view <- function(
        object, impulses, responses, terms
) {
    if (!is.data.frame(impulses) || !is.data.frame(responses)) {
        stop('impulses and responses must be data frames')
    }
    if (isTRUE(object$cdrgam$distributional) ||
            inherits(object, c(
                'cdrgam_distributional_block',
                'cdrgam_distributional_sparse'
            ))) {
        stop('Distributional cdrgam models are not yet supported')
    }
    selected <- .cdr_marginaleffects_specs(object, terms)
    pairs <- .cdr_marginaleffects_pairs(
        object,
        impulses,
        responses,
        selected$indices,
        selected$specifications
    )
    predictors <- unique(unlist(lapply(selected$indices, function(index) {
        .cdr_effect_spec_predictors(selected$specifications[[index]])
    }), use.names=FALSE))
    reserved <- c(
        '.cdrgam_pair', '.cdrgam_response', '.cdrgam_impulse',
        '.cdrgam_lag', '.cdrgam_outcome', '.cdrgam_scenario', 'lag', 'rowid'
    )
    conflict <- intersect(predictors, reserved)
    if (length(conflict)) {
        stop('Impulse predictor names are reserved by the adapter: ',
            paste(conflict, collapse=', '))
    }
    missing_predictors <- setdiff(predictors, names(impulses))
    if (length(missing_predictors)) {
        stop('Impulse stream is missing predictors: ',
            paste(missing_predictors, collapse=', '))
    }
    stream <- object$cdrgam$preparation$stream
    response_time <- responses[[stream$response_time]][pairs$response_index]
    impulse_time <- impulses[[stream$impulse_time]][pairs$impulse_index]
    response_name <- all.vars(stats::formula(
        object, type='effective'
    )[[2L]])
    if (length(response_name) != 1L || !(response_name %in% names(responses))) {
        stop('Response stream is missing the fitted response variable')
    }
    data <- data.frame(
        .cdrgam_pair=seq_along(pairs$key),
        .cdrgam_response=pairs$response_index,
        .cdrgam_impulse=pairs$impulse_index,
        .cdrgam_lag=response_time - impulse_time,
        .cdrgam_outcome=responses[[response_name]][pairs$response_index],
        stringsAsFactors=FALSE
    )
    for (predictor in predictors) {
        value <- impulses[[predictor]][pairs$impulse_index]
        if (!is.numeric(value) || any(!is.finite(value))) {
            stop('Impulse predictor ', predictor, ' must be finite and numeric')
        }
        data[[predictor]] <- value
    }
    response_rows <- sort(unique(pairs$response_index))
    response_map <- match(pairs$response_index, response_rows)
    response_data <- responses[response_rows, , drop=FALSE]
    prediction_data <- list(impulses=impulses, responses=response_data)
    base_design <- stats::predict(
        object, newdata=prediction_data, type='lpmatrix'
    )
    if (is.list(base_design) || ncol(base_design) != length(stats::coef(object))) {
        stop('The fitted model did not return a compatible linear-predictor matrix')
    }
    base_link <- stats::predict(object, newdata=prediction_data, type='link')
    if (!is.numeric(base_link) || length(base_link) != nrow(response_data)) {
        stop('The fitted model did not return one linear predictor per response')
    }
    offset <- as.numeric(base_link - base_design %*% stats::coef(object))
    context <- list(
        coefficient_count=length(stats::coef(object)),
        terms=object$cdrgam$terms,
        term_indices=selected$indices,
        specifications=selected$specifications,
        membership=pairs$membership,
        pair_response=response_map,
        pair_impulse=pairs$impulse_index,
        responses=response_data,
        base_design=base_design,
        offset=offset,
        predictors=predictors
    )
    context$event_design <- .cdr_marginaleffects_event_design(
        context, data, seq_len(nrow(data))
    )
    context$coefficient_index <- which(
        Matrix::colSums(abs(Matrix::Matrix(base_design, sparse=TRUE))) > 0 |
            Matrix::colSums(abs(context$event_design)) > 0
    )
    if (!length(context$coefficient_index)) {
        stop('The selected response--impulse links have an empty design')
    }
    list(data=data, context=context)
}

#' Prepare a cdrgam model for marginaleffects
#'
#' Constructs an evaluation view with one row per response--impulse link.
#' Changing an impulse predictor in a row replaces that event's contribution
#' while holding the response's remaining impulse history and covariates fixed.
#'
#' The returned object implements the `marginaleffects` coefficient, covariance,
#' and prediction extension methods. Each row is one response--impulse
#' incidence. An ungrouped average over rows is incidence-weighted, so responses
#' with longer histories receive more weight. Use the internal
#' `.cdrgam_response` and `.cdrgam_impulse` columns for response- or event-level
#' aggregation.
#'
#' The view stores a response-level linear-predictor matrix for the supplied
#' streams. Covariance extraction is restricted to coefficients that affect
#' those responses, including selected sparse covariance solves for the sparse
#' backend. Distributional models are not yet supported.
#'
#' @param object A fitted single-parameter `cdrgam` model.
#' @param impulses,responses Impulse and response streams on the source scale.
#' @param terms Optional impulse-response term labels or selector list. The
#'   default includes every impulse-response term, including the rate term.
#' @param component Quantity represented by each evaluation row. The initial
#'   implementation supports complete response predictions.
#' @return A `cdrgam_marginaleffects` model view accepted by `marginaleffects`.
#'   Its attached model data contain response, impulse, and lag identifiers plus
#'   the scalar impulse predictors available for comparisons and slopes.
#' @examples
#' \dontrun{
#' view <- marginaleffects_view(fit, impulses, responses)
#' marginaleffects::slopes(view, variables="surprisal")
#' }
#' @export
marginaleffects_view <- function(
        object,
        impulses,
        responses,
        terms=NULL,
        component=c('response')
) {
    if (!is_cdrgam(object)) stop('object must be a fitted cdrgam model')
    component <- match.arg(component)
    prepared <- .cdr_marginaleffects_response_view(
        object, impulses, responses, terms
    )
    output <- list(
        fit=object,
        component=component,
        context=prepared$context
    )
    class(output) <- 'cdrgam_marginaleffects'
    attr(output, 'marginaleffects_modeldata') <- prepared$data
    .cdr_marginaleffects_enable()
    output
}

#' @export
formula.cdrgam_marginaleffects <- function(x, ...) {
    predictors <- x$context$predictors
    stats::reformulate(predictors, response='.cdrgam_outcome')
}

#' @export
coef.cdrgam_marginaleffects <- function(object, ...) {
    stats::coef(object$fit)[object$context$coefficient_index]
}

#' @export
vcov.cdrgam_marginaleffects <- function(object, unconditional=FALSE, ...) {
    indices <- object$context$coefficient_index
    if (inherits(object$fit, 'cdrgam_sparse')) {
        .sparse_selected_vcov(
            object$fit, indices, unconditional=unconditional
        )
    } else {
        stats::vcov(
            object$fit, unconditional=unconditional
        )[indices, indices, drop=FALSE]
    }
}

#' @export
predict.cdrgam_marginaleffects <- function(
        object,
        newdata=NULL,
        type=c('response', 'link', 'lpmatrix'),
        se.fit=FALSE,
        unconditional=FALSE,
        ...
) {
    type <- if (is.null(type)) 'response' else match.arg(type)
    if (is.null(newdata)) {
        newdata <- attr(object, 'marginaleffects_modeldata')
    }
    if (!is.data.frame(newdata)) stop('newdata must be a data frame')
    required <- c(
        '.cdrgam_pair', '.cdrgam_response', '.cdrgam_impulse', '.cdrgam_lag',
        object$context$predictors
    )
    missing <- setdiff(required, names(newdata))
    if (length(missing)) {
        stop('newdata is missing adapter columns: ', paste(missing, collapse=', '))
    }
    pair_map <- match(newdata$.cdrgam_pair,
        attr(object, 'marginaleffects_modeldata')$.cdrgam_pair)
    if (anyNA(pair_map)) {
        stop('newdata contains unknown response--impulse links')
    }
    current <- .cdr_marginaleffects_event_design(
        object$context, newdata, pair_map
    )
    delta <- current - object$context$event_design[pair_map, , drop=FALSE]
    response_index <- object$context$pair_response[pair_map]
    design <- object$context$base_design[response_index, , drop=FALSE] + delta
    if (identical(type, 'lpmatrix')) return(design)
    coefficients <- stats::coef(object$fit)
    link <- as.numeric(
        design %*% coefficients + object$context$offset[response_index]
    )
    estimate <- link
    if (identical(type, 'response')) {
        estimate <- object$fit$family$linkinv(estimate)
    }
    if (!isTRUE(se.fit)) return(estimate)
    indices <- object$context$coefficient_index
    covariance <- stats::vcov(object, unconditional=unconditional)
    selected_design <- design[, indices, drop=FALSE]
    standard_error <- sqrt(pmax(0, Matrix::rowSums(
        (selected_design %*% covariance) * selected_design
    )))
    if (identical(type, 'response')) {
        standard_error <- abs(object$fit$family$mu.eta(link)) * standard_error
    }
    list(fit=estimate, se.fit=standard_error)
}

#' Construct response-history scenarios
#'
#' Expands response--impulse incidences from a [marginaleffects_view()] over a
#' Cartesian grid of lag and impulse-predictor values. Each output row retains
#' the identity of its reference incidence, so prediction replaces one event's
#' contribution while holding the rest of that response's history fixed.
#'
#' @param object A `cdrgam_marginaleffects` view.
#' @param at Named list of values for `lag` or impulse predictors. `lag` is
#'   stored in the adapter column `.cdrgam_lag`.
#' @param rows Optional integer or logical subset of reference incidences.
#' @return A data frame suitable for `predict()`, [estimate_response_effect()],
#'   or direct use as `newdata` in `marginaleffects`.
#' @details Changing `lag` re-evaluates the selected IRF windows. A term stops
#'   contributing when the scenario lag lies outside its window.
#' @export
response_scenarios <- function(object, at=list(), rows=NULL) {
    if (!inherits(object, 'cdrgam_marginaleffects')) {
        stop('object must be a cdrgam_marginaleffects view')
    }
    if (!is.list(at) || (length(at) && (is.null(names(at)) ||
            any(!nzchar(names(at))) || anyDuplicated(names(at))))) {
        stop('at must be a uniquely named list')
    }
    allowed <- c('lag', object$context$predictors)
    unknown <- setdiff(names(at), allowed)
    if (length(unknown)) {
        stop('Unknown scenario axes: ', paste(unknown, collapse=', '))
    }
    reference <- attr(object, 'marginaleffects_modeldata')
    if (!is.null(rows)) {
        if (is.logical(rows)) {
            if (length(rows) != nrow(reference) || anyNA(rows)) {
                stop('Logical rows must have one non-missing value per incidence')
            }
        } else {
            rows <- as.integer(rows)
            if (!length(rows) || anyNA(rows) || any(rows < 1L) ||
                    any(rows > nrow(reference))) {
                stop('rows contains invalid incidence indices')
            }
        }
        reference <- reference[rows, , drop=FALSE]
    }
    axes <- lapply(at, function(value) {
        value <- as.numeric(value)
        if (!length(value) || any(!is.finite(value))) {
            stop('Scenario axes must contain finite numeric values')
        }
        value
    })
    grid <- if (length(axes)) {
        do.call(expand.grid, c(axes, KEEP.OUT.ATTRS=FALSE,
            stringsAsFactors=FALSE))
    } else data.frame(.placeholder=1L)
    scenario_count <- nrow(grid)
    incidence_count <- nrow(reference)
    output <- reference[rep(seq_len(incidence_count), times=scenario_count),
        , drop=FALSE]
    grid_rows <- rep(seq_len(scenario_count), each=incidence_count)
    for (name in names(axes)) {
        column <- if (identical(name, 'lag')) '.cdrgam_lag' else name
        output[[column]] <- grid[[name]][grid_rows]
    }
    output$.cdrgam_scenario <- grid_rows
    rownames(output) <- NULL
    output
}

#' Generate response-effect plot data
#'
#' Evaluates predictions, slopes, or comparisons on response-history scenarios
#' through the `marginaleffects` extension API. This preserves each response's
#' remaining impulse history and returns point estimates with pointwise
#' uncertainty in a tidy table. Grouping columns are passed to
#' `marginaleffects`; use `.cdrgam_scenario` to average each scenario over its
#' reference incidences, and add `.cdrgam_response` or `.cdrgam_impulse` when a
#' response- or event-weighted result is required.
#'
#' @param object A `cdrgam_marginaleffects` view.
#' @param newdata Scenario data, usually returned by [response_scenarios()].
#'   The default uses the unmodified incidences in the view.
#' @param estimand Quantity to estimate: fitted response values, numeric slopes,
#'   or finite comparisons.
#' @param variables Predictor specification passed to `marginaleffects`. It is
#'   required for slopes and comparisons.
#' @param by Optional columns defining averages. `NULL` retains one result per
#'   scenario row.
#' @param type Evaluate on the response or linear-predictor scale.
#' @param vcov Include fitted-coefficient uncertainty, or supply a covariance
#'   matrix or function accepted by [marginaleffects_view()].
#' @param level Confidence level for pointwise intervals.
#' @param ... Additional arguments passed to the selected `marginaleffects`
#'   function.
#' @return A `cdrgam_response_effect_grid` data frame. Adapter identity columns
#'   and scenario axes are retained when the requested aggregation permits it.
#' @export
estimate_response_effect <- function(
        object,
        newdata=NULL,
        estimand=c('prediction', 'slope', 'comparison'),
        variables=NULL,
        by=NULL,
        type=c('response', 'link'),
        vcov=TRUE,
        level=0.95,
        ...
) {
    if (!inherits(object, 'cdrgam_marginaleffects')) {
        stop('object must be a cdrgam_marginaleffects view')
    }
    if (!requireNamespace('marginaleffects', quietly=TRUE)) {
        stop('estimate_response_effect requires the marginaleffects package')
    }
    estimand <- match.arg(estimand)
    type <- match.arg(type)
    if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
            level <= 0 || level >= 1) {
        stop('level must lie strictly between zero and one')
    }
    if (is.null(newdata)) {
        newdata <- attr(object, 'marginaleffects_modeldata')
    }
    if (!is.data.frame(newdata)) stop('newdata must be a data frame')
    if (!is.null(by)) {
        by <- as.character(by)
        missing_by <- setdiff(by, names(newdata))
        if (length(missing_by)) {
            stop('by columns are missing from newdata: ',
                paste(missing_by, collapse=', '))
        }
    }
    if (!identical(estimand, 'prediction') && is.null(variables)) {
        stop('variables is required for slopes and comparisons')
    }
    arguments <- list(
        model=object, newdata=newdata, by=by, type=type, vcov=vcov,
        conf_level=level
    )
    if (!identical(estimand, 'prediction')) arguments$variables <- variables
    arguments <- c(arguments, list(...))
    estimate <- switch(
        estimand,
        prediction=do.call(marginaleffects::predictions, arguments),
        slope=do.call(marginaleffects::slopes, arguments),
        comparison=do.call(marginaleffects::comparisons, arguments)
    )
    output <- as.data.frame(estimate)
    class(output) <- c('cdrgam_response_effect_grid', class(output))
    attr(output, 'level') <- level
    attr(output, 'scale') <- type
    attr(output, 'estimand') <- estimand
    output
}

#' @exportS3Method marginaleffects::get_coef
get_coef.cdrgam_marginaleffects <- function(model, ...) {
    stats::coef(model)
}

#' @exportS3Method marginaleffects::set_coef
set_coef.cdrgam_marginaleffects <- function(model, coefs, ...) {
    indices <- model$context$coefficient_index
    expected <- names(stats::coef(model$fit))[indices]
    if (!is.null(names(coefs))) coefs <- coefs[expected]
    if (length(coefs) != length(expected) || anyNA(coefs)) {
        stop('Replacement coefficients do not match the fitted model')
    }
    model$fit$coefficients[indices] <- as.numeric(coefs)
    model
}

#' @exportS3Method marginaleffects::get_vcov
get_vcov.cdrgam_marginaleffects <- function(model, vcov=NULL, ...) {
    if (identical(vcov, FALSE)) return(NULL)
    indices <- model$context$coefficient_index
    select <- function(value) {
        if (!is.matrix(value)) stop('The covariance provider must return a matrix')
        if (nrow(value) == length(indices)) return(value)
        if (nrow(value) == length(stats::coef(model$fit))) {
            return(value[indices, indices, drop=FALSE])
        }
        stop('The covariance matrix does not match the adapter coefficients')
    }
    if (is.matrix(vcov)) return(select(vcov))
    if (is.function(vcov)) return(select(vcov(model$fit)))
    if (!is.null(vcov) && !isTRUE(vcov)) {
        stop('This cdrgam adapter supports the fitted covariance or a matrix/function')
    }
    stats::vcov(model, ...)
}

#' @exportS3Method marginaleffects::get_predict
get_predict.cdrgam_marginaleffects <- function(
        model, newdata, type=c('response', 'link'), ...
) {
    estimate <- stats::predict(model, newdata=newdata, type=type, ...)
    rowid <- if ('rowid' %in% names(newdata)) newdata$rowid else
        seq_len(nrow(newdata))
    data.frame(rowid=rowid, estimate=estimate)
}
