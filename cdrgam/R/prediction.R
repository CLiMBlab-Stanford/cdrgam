.cdr_prediction_setup <- function(setup) {
    training_frame <- if (!is.null(setup$mf)) setup$mf else setup$model
    variable_levels <- if (is.null(training_frame)) list() else lapply(
        training_frame[vapply(training_frame, is.factor, logical(1))],
        levels
    )
    list(
        pterms=setup$pterms,
        assign=setup$assign,
        xlevels=setup$xlevels,
        nsdf=setup$nsdf,
        smooth=setup$smooth,
        dimension=if (is.null(setup$X)) length(stats::coef(setup)) else
            ncol(setup$X),
        term_names=if (is.null(setup$X)) names(stats::coef(setup)) else
            colnames(setup$X),
        variable_levels=variable_levels
    )
}

.cdr_predict_irf_matrices <- function(
        object,
        impulses,
        responses,
        chunk_size,
        source_impulses=impulses,
        source_responses=responses
) {
    stream <- object$cdrgam$preparation$stream
    specifications <- object$cdrgam$preparation$specification
    metadata <- object$cdrgam$terms
    if (is.null(stream) || is.null(specifications) ||
            length(specifications) != length(metadata)) {
        stop('This fit does not contain the stream metadata required for prediction')
    }
    output <- vector('list', length(metadata))
    for (i in seq_along(metadata)) {
        info <- metadata[[i]]
        specification <- specifications[[i]]
        predictors <- if (is.null(specification$predictors)) {
            if (isTRUE(specification$constant)) character() else
                specification$predictor
        } else specification$predictors
        constant <- isTRUE(specification$constant)
        missing_predictors <- setdiff(predictors, names(impulses))
        if (length(missing_predictors)) {
            stop('Prediction impulse stream is missing: ',
                paste(missing_predictors, collapse=', '))
        }
        predictor_values <- lapply(predictors, function(predictor) {
            value <- impulses[[predictor]]
            if (!is.numeric(value) || any(!is.finite(value))) {
                stop('Prediction impulse predictors must be finite numeric columns')
            }
            value
        })
        names(predictor_values) <- predictors
        links <- .build_history_links(
            source_impulses,
            source_responses,
            stream$series,
            stream$impulse_time,
            stream$response_time,
            specification$window,
            history_length=if (is.null(stream$history_length)) Inf else
                stream$history_length
        )
        scaling <- object$cdrgam$scaling
        if (is.null(scaling)) scaling <- object$cdrgam$preparation$scaling
        if (!is.null(scaling) && isTRUE(scaling$enabled)) {
            links$delay <- links$delay / scaling$time_divisor
        }
        base_dimension <- if (is.null(info$base_dimension)) {
            length(info$coefficient_index)
        } else {
            info$base_dimension
        }
        base <- matrix(0, nrow=nrow(responses), ncol=base_dimension)
        if (length(links$delay)) {
            if (!is.null(info$axis)) {
                linear_predictors <- if (is.null(info$linear_predictors)) {
                    predictors[vapply(specification$k_p, is.null, logical(1))]
                } else info$linear_predictors
                linked_weights <- rep.int(1, length(links$delay))
                for (predictor in linear_predictors) {
                    linked_weights <- linked_weights * predictor_values[[predictor]][
                        links$impulse_index
                    ]
                }
                for (start in seq.int(1L, length(links$delay), by=chunk_size)) {
                    end <- min(length(links$delay), start + chunk_size - 1L)
                    rows <- start:end
                    axis_data <- setNames(lapply(info$axis, function(axis) {
                        switch(
                            axis$role,
                            lag=links$delay[rows],
                            time={
                                if (!(axis$variable %in% names(responses))) {
                                    stop('Prediction responses are missing time axis: ',
                                        axis$variable)
                                }
                                responses[[axis$variable]][links$response_index[rows]]
                            },
                            predictor=predictor_values[[axis$variable]][
                                links$impulse_index[rows]
                            ]
                        )
                    }), vapply(info$axis, `[[`, character(1), 'internal'))
                    basis <- mgcv::PredictMat(
                        info$basis,
                        axis_data,
                        n=length(rows)
                    ) %*% info$transform
                    basis <- basis * linked_weights[rows]
                    accumulated <- rowsum(
                        basis,
                        links$response_index[rows],
                        reorder=FALSE
                    )
                    response_rows <- as.integer(rownames(accumulated))
                    base[response_rows, ] <-
                        base[response_rows, , drop=FALSE] + accumulated
                }
            } else {
            values <- if (constant) rep.int(1, nrow(impulses)) else
                predictor_values[[1L]]
            linked_values <- values[links$impulse_index]
            tensor_values <- if (!is.null(specification$varying)) {
                if (!(specification$varying %in% names(responses))) {
                    stop(
                        'Prediction responses are missing varying covariate: ',
                        specification$varying
                    )
                }
                responses[[specification$varying]][links$response_index]
            } else {
                linked_values
            }
            for (start in seq.int(1L, length(links$delay), by=chunk_size)) {
                end <- min(length(links$delay), start + chunk_size - 1L)
                rows <- start:end
                tensor <- startsWith(info$type, 'nonlinear') ||
                    startsWith(info$type, 'varying')
                basis <- if (tensor) {
                    mgcv::PredictMat(
                        info$basis,
                        list(
                            cdr_delay=links$delay[rows],
                            cdr_value=tensor_values[rows]
                        ),
                        n=length(rows)
                    ) %*% info$transform
                } else {
                    lag_basis <- mgcv::PredictMat(
                        info$basis,
                        list(cdr_delta=links$delay[rows]),
                        n=length(rows)
                    )
                    if (!is.null(info$transform)) {
                        lag_basis <- lag_basis %*% info$transform
                    }
                    lag_basis * linked_values[rows]
                }
                if (!is.null(specification$varying)) {
                    basis <- basis * linked_values[rows]
                }
                accumulated <- rowsum(
                    basis,
                    links$response_index[rows],
                    reorder=FALSE
                )
                response_rows <- as.integer(rownames(accumulated))
                base[response_rows, ] <- base[response_rows, , drop=FALSE] +
                    accumulated
            }
            }
        }
        if (!is.null(specification$by)) {
            if (!(specification$by %in% names(responses))) {
                stop('Prediction responses are missing by covariate: ', specification$by)
            }
            by_values <- responses[[specification$by]]
            if (!is.numeric(by_values) || any(!is.finite(by_values))) {
                stop('Prediction by covariates must be finite numeric columns')
            }
            base <- base * by_values
        }
        if (!length(info$group_levels)) {
            output[[i]] <- base
            next
        }
        group <- info$group
        if (!(group %in% names(responses))) {
            stop('Prediction response stream is missing grouping factor: ', group)
        }
        group_index <- match(as.character(responses[[group]]), info$group_levels)
        unknown <- is.na(group_index)
        if (any(unknown)) {
            warning(
                'Factor levels ',
                paste(unique(as.character(responses[[group]])[unknown]), collapse=', '),
                ' for ', group, ' were not in the original fit; ',
                'their grouped deviations were set to zero',
                call.=FALSE
            )
        }
        known <- which(!unknown)
        if (inherits(object, 'cdrgam_sparse')) {
            expanded <- if (length(known)) {
                Matrix::sparseMatrix(
                    i=rep(known, each=base_dimension),
                    j=rep(
                        (group_index[known] - 1L) * base_dimension,
                        each=base_dimension
                    ) + rep.int(seq_len(base_dimension), times=length(known)),
                    x=as.vector(t(base[known, , drop=FALSE])),
                    dims=c(
                        nrow(responses),
                        length(info$group_levels) * base_dimension
                    )
                )
            } else {
                Matrix::Matrix(
                    0,
                    nrow=nrow(responses),
                    ncol=length(info$group_levels) * base_dimension,
                    sparse=TRUE
                )
            }
        } else {
            expanded <- matrix(
                0,
                nrow=nrow(responses),
                ncol=length(info$group_levels) * base_dimension
            )
        }
        if (length(known) && !inherits(object, 'cdrgam_sparse')) {
            for (level in unique(group_index[known])) {
                selected <- known[group_index[known] == level]
                columns <- (level - 1L) * base_dimension +
                    seq_len(base_dimension)
                expanded[selected, columns] <- base[selected, , drop=FALSE]
            }
        }
        output[[i]] <- expanded
    }
    output
}

.cdr_setup_lpmatrix <- function(
        prediction,
        responses,
        extra=list()
) {
    n <- nrow(responses)
    data <- as.list(responses)
    for (name in names(extra)) data[[name]] <- extra[[name]]
    output <- matrix(0, nrow=n, ncol=prediction$dimension)
    colnames(output) <- prediction$term_names
    offset <- numeric(n)
    if (prediction$nsdf > 0L) {
        parametric_terms <- stats::delete.response(prediction$pterms)
        intercept_only <- !length(attr(parametric_terms, 'term.labels')) &&
            identical(attr(parametric_terms, 'intercept'), 1L)
        if (intercept_only) {
            frame <- NULL
            parametric <- matrix(
                1,
                nrow=n,
                ncol=1L,
                dimnames=list(NULL, '(Intercept)')
            )
        } else {
            frame <- stats::model.frame(
                parametric_terms,
                data=data,
                xlev=prediction$xlevels,
                na.action=stats::na.pass
            )
            parametric <- stats::model.matrix(parametric_terms, frame)
        }
        required <- prediction$term_names[seq_len(prediction$nsdf)]
        positions <- match(required, colnames(parametric))
        if (anyNA(positions)) {
            stop('New response data do not reproduce the fitted parametric design')
        }
        output[, seq_len(prediction$nsdf)] <- parametric[, positions, drop=FALSE]
        frame_offset <- if (is.null(frame)) NULL else stats::model.offset(frame)
        if (!is.null(frame_offset)) offset <- frame_offset
    }
    for (smooth in prediction$smooth) {
        smooth_data <- data
        unknown <- rep.int(FALSE, n)
        if (inherits(smooth, 'random.effect')) {
            for (variable in smooth$term) {
                levels <- prediction$variable_levels[[variable]]
                if (is.null(levels)) next
                value <- as.character(smooth_data[[variable]])
                missing <- is.na(value) | !(value %in% levels)
                unknown <- unknown | missing
                if (any(missing)) {
                    warning(
                        'Factor levels ',
                        paste(unique(value[missing]), collapse=', '),
                        ' for ', variable, ' were not in the original fit; ',
                        'their random effects were set to zero',
                        call.=FALSE
                    )
                }
                value[missing] <- levels[[1L]]
                smooth_data[[variable]] <- factor(value, levels=levels)
            }
        }
        matrix <- if (inherits(smooth, 'cdr.smooth')) {
            Predict.matrix.cdr.smooth(smooth, smooth_data)
        } else {
            mgcv::PredictMat(smooth, smooth_data, n=n)
        }
        if (any(unknown)) matrix[unknown, ] <- 0
        if (!is.null(smooth$by) && !identical(smooth$by, 'NA')) {
            if (!(smooth$by %in% names(data))) {
                stop('New response data are missing by variable: ', smooth$by)
            }
            matrix <- matrix * data[[smooth$by]]
        }
        output[, smooth$first.para:smooth$last.para] <- matrix
    }
    list(X=output, offset=offset)
}

.cdr_predict_random_effect <- function(effect, responses) {
    data <- responses
    unknown <- rep.int(FALSE, nrow(data))
    for (variable in effect$variables) {
        value <- as.character(data[[variable]])
        missing <- is.na(value) | !(value %in% effect$levels[[variable]])
        unknown <- unknown | missing
        value[missing] <- effect$levels[[variable]][[1L]]
        data[[variable]] <- factor(value, levels=effect$levels[[variable]])
    }
    if (any(unknown)) {
        warning(
            'Factor levels absent from the original fit had their random ',
            'effects set to zero',
            call.=FALSE
        )
    }
    formula <- stats::as.formula(paste(
        '~', paste(effect$variables, collapse=':'), '- 1'
    ))
    raw <- Matrix::sparse.model.matrix(formula, data=data)
    output <- Matrix::Matrix(
        0,
        nrow=nrow(data),
        ncol=length(effect$column_names),
        sparse=TRUE
    )
    positions <- match(colnames(raw), effect$column_names)
    present <- which(!is.na(positions))
    if (length(present)) output[, positions[present]] <- raw[, present, drop=FALSE]
    if (any(unknown)) output[unknown, ] <- 0
    output
}

# Internal implementation shared by the fitted-class prediction methods.
.predict_cdrgam_streams <- function(
        object,
        impulses,
        responses,
        type=c('response', 'link', 'lpmatrix'),
        se.fit=FALSE,
        chunk_size=10000,
        unconditional=FALSE
) {
    type <- match.arg(type)
    if (!is_cdrgam(object) || !is.data.frame(impulses) ||
            !is.data.frame(responses)) {
        stop('object must be a cdrgam fit and streams must be data frames')
    }
    source_impulses <- impulses
    source_responses <- responses
    scaling <- object$cdrgam$scaling
    if (is.null(scaling)) scaling <- object$cdrgam$preparation$scaling
    impulses <- .cdr_apply_scaling(impulses, scaling, 'impulses')
    responses <- .cdr_apply_scaling(responses, scaling, 'responses')
    irf_matrices <- .cdr_predict_irf_matrices(
        object,
        impulses,
        responses,
        chunk_size,
        source_impulses=source_impulses,
        source_responses=source_responses
    )
    if (inherits(object, 'gam')) {
        extra <- stats::setNames(
            irf_matrices,
            paste0('cdr_term_', seq_along(irf_matrices))
        )
        assembled <- .cdr_setup_lpmatrix(
            .cdr_prediction_setup(object),
            responses,
            extra
        )
    } else {
        prediction <- object$cdrgam$prediction
        if (is.null(prediction$setup)) {
            stop('This fit predates stored prediction metadata')
        }
        if (inherits(object, 'cdrgam_block')) {
            extra <- stats::setNames(
                irf_matrices,
                paste0('cdr_term_', seq_along(irf_matrices))
            )
            assembled <- .cdr_setup_lpmatrix(
                prediction$setup,
                responses,
                extra
            )
        } else {
            ordinary <- .cdr_setup_lpmatrix(
                prediction$setup,
                responses
            )
            random <- lapply(
                prediction$random_effects,
                .cdr_predict_random_effect,
                responses=responses
            )
            prediction_parts <- c(list(ordinary$X), random, irf_matrices)
            prediction_parts <- lapply(
                prediction_parts,
                Matrix::Matrix,
                sparse=TRUE
            )
            assembled <- list(
                X=do.call(cbind, prediction_parts),
                offset=ordinary$offset
            )
        }
    }
    if (ncol(assembled$X) != length(stats::coef(object))) {
        stop('Prediction design does not match fitted coefficient dimension')
    }
    if (identical(type, 'lpmatrix')) return(assembled$X)
    linear_predictor <- as.numeric(
        assembled$X %*% stats::coef(object) + assembled$offset
    )
    fit <- if (identical(type, 'response')) {
        object$family$linkinv(linear_predictor)
    } else {
        linear_predictor
    }
    if (!isTRUE(se.fit)) return(fit)
    if (inherits(object, 'cdrgam_sparse') && !isTRUE(unconditional)) {
        variance <- numeric(nrow(assembled$X))
        for (start in seq.int(1L, nrow(assembled$X), by=256L)) {
            rows <- start:min(nrow(assembled$X), start + 255L)
            block <- assembled$X[rows, , drop=FALSE]
            solved <- .cdr_factor_solve(
                object$sparse$factor,
                Matrix::t(block)
            )
            variance[rows] <- Matrix::rowSums(
                block * t(as.matrix(solved))
            ) * object$scale
        }
        standard_error <- sqrt(pmax(0, variance))
    } else {
        covariance <- stats::vcov(object, unconditional=unconditional)
        standard_error <- sqrt(pmax(
            0,
            rowSums((assembled$X %*% covariance) * assembled$X)
        ))
    }
    if (identical(type, 'response')) {
        standard_error <- standard_error * abs(
            object$family$mu.eta(linear_predictor)
        )
    }
    list(fit=fit, se.fit=standard_error)
}
