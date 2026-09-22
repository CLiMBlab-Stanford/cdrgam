.cdrgam_distributional_formula <- function(formula, response=NULL) {
    if (!inherits(formula, 'formula')) {
        stop('Every distributional predictor must be a formula')
    }
    if (length(formula) == 3L) {
        formula_response <- all.vars(formula[[2L]])
        if (length(formula_response) != 1L) {
            stop('A distributional response must name one response column')
        }
        if (!is.null(response) && !identical(formula_response, response)) {
            stop('Distributional formulas must use the same response')
        }
        return(formula)
    }
    if (length(formula) != 2L || is.null(response)) {
        stop('Only non-location distributional predictors may be one-sided')
    }
    output <- stats::as.formula(
        paste(response, '~', paste(deparse(formula[[2L]]), collapse='')),
        env=environment(formula)
    )
    output
}

.prepare_cdrgam_distributional <- function(formulas, arguments) {
    if (is.null(names(formulas)) || any(!nzchar(names(formulas))) ||
            anyDuplicated(names(formulas))) {
        stop('Distributional formulas must have unique parameter names')
    }
    if (!identical(names(formulas), c('location', 'scale'))) {
        stop(
            'The initial distributional API requires formulas named ',
            'location and scale, in that order'
        )
    }
    if (isTRUE(arguments$rescale_predictors)) {
        stop(
            'Distributional predictor rescaling is not yet supported; ',
            'prepare source variables explicitly or use ',
            'rescale_predictors=FALSE'
        )
    }
    location <- .cdrgam_distributional_formula(formulas$location)
    response <- all.vars(location[[2L]])
    if (length(response) != 1L) {
        stop('The location formula must name one response column')
    }
    compiled_formulas <- list(
        location=location,
        scale=.cdrgam_distributional_formula(formulas$scale, response)
    )
    designs <- lapply(compiled_formulas, function(formula) {
        do.call(prepare_cdrgam, c(list(formula=formula), arguments))
    })
    row_counts <- vapply(designs, function(design) {
        nrow(design$responses)
    }, integer(1))
    if (length(unique(row_counts)) != 1L) {
        stop('Distributional predictors retained different response rows')
    }
    output <- list(
        formula=formulas,
        normalized_formula=lapply(designs, `[[`, 'normalized_formula'),
        effective_formula=lapply(designs, `[[`, 'effective_formula'),
        response_name=response,
        responses=designs$location$responses,
        parameters=designs,
        simplifications=do.call(
            rbind,
            lapply(names(designs), function(parameter) {
                value <- designs[[parameter]]$simplifications
                value$parameter <- rep.int(parameter, nrow(value))
                value[c('parameter', setdiff(names(value), 'parameter'))]
            })
        ),
        configuration=list(
            drop.unused.levels=arguments$drop.unused.levels,
            rescale_predictors=FALSE
        )
    )
    class(output) <- c('cdrgam_distributional_design', 'cdrgam_design')
    output
}

.cdrgam_distributional_preparation <- function(design) {
    list(
        configuration=design$configuration,
        plan=design$plan,
        stream=design$stream,
        specification=design$specification,
        simplifications=design$simplifications,
        scaling=design$scaling,
        identifiability=design$identifiability,
        normalized_formula=design$normalized_formula,
        effective_formula=design$effective_formula
    )
}

.cdrgam_distributional_term_metadata <- function(
        term, fit, data_name, parameter
) {
    smooth_index <- which(vapply(fit$smooth, function(smooth) {
        identical(smooth$term, data_name) ||
            identical(paste(smooth$term, collapse=','), data_name)
    }, logical(1)))
    coefficient_index <- if (length(smooth_index) == 1L) {
        seq.int(
            fit$smooth[[smooth_index]]$first.para,
            fit$smooth[[smooth_index]]$last.para
        )
    } else integer()
    list(
        name=term$name,
        parameter=parameter,
        data_name=data_name,
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
}

.fit_distributional_mgcv <- function(
        design, family, method=NULL, engine=c('bam', 'gam'), ...
) {
    engine <- match.arg(engine)
    if (identical(engine, 'bam')) {
        stop('mgcv does not support general families such as gaulss in bam')
    }
    parameter_names <- names(design$parameters)
    formula_env <- new.env(parent=parent.frame())
    formula_env$s <- mgcv::s
    formula_env$te <- mgcv::te
    formula_env$ti <- mgcv::ti
    formula_env$t2 <- mgcv::t2
    data <- as.list(design$responses)
    translated <- vector('list', length(parameter_names))
    names(translated) <- parameter_names
    data_names <- vector('list', length(parameter_names))
    names(data_names) <- parameter_names
    for (parameter_index in seq_along(parameter_names)) {
        parameter <- parameter_names[[parameter_index]]
        parameter_design <- design$parameters[[parameter]]
        ordinary <- parameter_design$ordinary_formula
        response_text <- paste(deparse(ordinary[[2L]]), collapse='')
        base_rhs <- paste(deparse(ordinary[[3L]]), collapse='')
        smooth_terms <- character(length(parameter_design$terms))
        data_names[[parameter]] <- character(length(parameter_design$terms))
        for (term_index in seq_along(parameter_design$terms)) {
            term <- parameter_design$terms[[term_index]]
            data_name <- paste0(
                'cdr_', parameter, '_term_', term_index
            )
            spec_name <- paste0(data_name, '_spec')
            data_names[[parameter]][[term_index]] <- data_name
            data[[data_name]] <- seq_len(nrow(design$responses))
            formula_env[[spec_name]] <- list(
                S=term$S,
                S.scale=term$S.scale,
                rank=term$rank,
                null.space.dim=term$null.space.dim,
                p=ncol(term$X),
                X=term$X
            )
            smooth_terms[[term_index]] <- paste0(
                's(', data_name, ', bs="cdr", k=', ncol(term$X),
                ', xt=', spec_name, ')'
            )
        }
        rhs <- paste(c(base_rhs, smooth_terms), collapse=' + ')
        text <- if (parameter_index == 1L) {
            paste(response_text, '~', rhs)
        } else paste('~', rhs)
        translated[[parameter]] <- stats::as.formula(text, env=formula_env)
    }
    arguments <- c(
        list(formula=unname(translated), family=family, data=data),
        list(...)
    )
    if (!is.null(method)) arguments$method <- method
    fit <- if (identical(engine, 'bam')) {
        do.call(mgcv::bam, arguments)
    } else do.call(mgcv::gam, arguments)
    metadata <- list()
    labels <- character()
    parameter_terms <- vector('list', length(parameter_names))
    names(parameter_terms) <- parameter_names
    for (parameter in parameter_names) {
        parameter_design <- design$parameters[[parameter]]
        parameter_terms[[parameter]] <- integer(length(parameter_design$terms))
        for (term_index in seq_along(parameter_design$terms)) {
            term <- parameter_design$terms[[term_index]]
            data_name <- data_names[[parameter]][[term_index]]
            smooth_index <- which(vapply(fit$smooth, function(smooth) {
                identical(smooth$term, data_name) ||
                    identical(paste(smooth$term, collapse=','), data_name)
            }, logical(1)))
            if (length(smooth_index) == 1L) {
                fit$smooth[[smooth_index]]$S.scale <- term$S.scale
            }
            metadata[[length(metadata) + 1L]] <-
                .cdrgam_distributional_term_metadata(
                    term, fit, data_name, parameter
                )
            labels[[length(labels) + 1L]] <- paste0(
                parameter, ':', term$name
            )
            parameter_terms[[parameter]][[term_index]] <- length(metadata)
        }
    }
    preparations <- lapply(
        design$parameters,
        .cdrgam_distributional_preparation
    )
    fit$cdrgam <- list(
        schema_version=1L,
        engine=engine,
        backend='mgcv',
        distributional=TRUE,
        parameter_names=parameter_names,
        parameter_terms=parameter_terms,
        formula=list(
            user=design$formula,
            normalized=design$normalized_formula,
            effective=design$effective_formula,
            mgcv=translated
        ),
        preparation=list(parameters=preparations),
        scaling=NULL,
        identifiability=lapply(
            design$parameters,
            `[[`,
            'identifiability'
        ),
        term_labels=labels,
        terms=metadata,
        prediction=list(data_names=data_names)
    )
    class(fit) <- c('cdrgam', class(fit))
    fit
}

.predict_cdrgam_distributional_streams <- function(
        object,
        impulses,
        responses,
        type=c('response', 'link', 'lpmatrix'),
        se.fit=FALSE,
        unconditional=FALSE,
        chunk_size=10000,
        ...
) {
    type <- match.arg(type)
    data <- as.list(responses)
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
        names(matrices) <- object$cdrgam$prediction$data_names[[parameter]]
        data[names(matrices)] <- matrices
    }
    mgcv::predict.gam(
        object,
        newdata=data,
        type=type,
        se.fit=se.fit,
        unconditional=unconditional,
        ...
    )
}
