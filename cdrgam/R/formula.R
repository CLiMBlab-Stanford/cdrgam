#' Specify an impulse-response term
#'
#' `irf()` is used inside a [cdrgam()] formula. Its predictor arguments name
#' numeric columns in the impulse stream, or the sole predictor is literal `1`
#' for the
#' deconvolutional intercept. Ordinary formula terms retain their usual
#' `mgcv` meaning and are evaluated against the response stream.
#'
#' @param predictor Unquoted numeric impulse-stream column, or literal `1`.
#' @param ... Additional unquoted numeric impulse-stream predictors. Multiple
#'   predictors define one interaction IRF.
#' @param window Inclusive lag window as `c(minimum, maximum)`, where lag is
#'   response time minus impulse time. Negative lags include impulses after the
#'   response time. `NULL` inherits the model-level `window`; without either,
#'   the default is `c(0, Inf)`. The minimum must be finite; the maximum may be
#'   infinite.
#' @param k_l Lag-basis dimension. This axis is always smooth and cannot be
#'   `NULL`.
#' @param k_t Optional response-time basis dimension. `NULL` makes the IRF
#'   stationary; a number models nonstationarity against `response_time`.
#' @param k_p Predictor-axis dimensions. Supply one entry per predictor;
#'   `NULL` means that predictor enters linearly, while a number gives a smooth
#'   marginal. Use a list, such as `list(NULL, 4)`, when mixing linear and
#'   smooth predictors because base R reduces `c(NULL, 4)` to `4`.
#' @param bs_l,bs_t Lag and response-time marginal basis names accepted by
#'   `mgcv`.
#' @param bs_p Predictor marginal basis names. Supply one value to recycle it
#'   or one entry per predictor.
#' @param k,bs,nonlinear,varying Compatibility controls for the original
#'   single-predictor API. They cannot be mixed with explicitly supplied
#'   axis-specific controls.
#' @param by Optional unquoted response-aligned numeric covariate multiplying
#'   the complete IRF contribution.
#' @param group Optional unquoted response-stream factor defining grouped IRF
#'   deviations.
#' @return An internal `cdrgam_irf_spec` when evaluated by the formula compiler.
#'   Before constraints, a tensor term's coefficient count is the product of
#'   its axis dimensions.
#' @export
irf <- function(
        predictor,
        ...,
        window=NULL,
        k_l=10,
        k_t=NULL,
        k_p=NULL,
        bs_l='cr',
        bs_t='cr',
        bs_p='cr',
        by=NULL,
        group=NULL,
        k=NULL,
        bs=NULL,
        nonlinear=NULL,
        varying=NULL
) {
    predictor_exprs <- as.list(substitute(list(predictor, ...)))[-1L]
    varying_expr <- substitute(varying)
    by_expr <- substitute(by)
    group_expr <- substitute(group)
    symbol_name <- function(expr, argument, allow_null=TRUE) {
        if (allow_null && identical(expr, quote(NULL))) {
            return(NULL)
        }
        if (!is.symbol(expr)) {
            stop(argument, ' must be an unquoted column name')
        }
        as.character(expr)
    }
    validate_window <- function(value) {
        if (is.null(value)) return(NULL)
        if (length(value) != 2L || !is.numeric(value) ||
                !is.finite(value[[1L]]) || is.na(value[[2L]]) ||
                value[[2L]] < value[[1L]]) {
            stop('window must satisfy finite min <= max; max may be Inf')
        }
        as.numeric(value)
    }
    validate_k <- function(value, argument, allow_null=FALSE) {
        if (allow_null && is.null(value)) return(NULL)
        if (!is.numeric(value) || length(value) != 1L ||
                !is.finite(value) || value < 3 || value != as.integer(value)) {
            stop(argument, ' must be NULL or one integer of at least 3')
        }
        as.integer(value)
    }
    validate_bs <- function(value, argument) {
        if (!is.character(value) || length(value) != 1L || is.na(value) ||
                !nzchar(value)) {
            stop(argument, ' must be one nonempty mgcv basis name')
        }
        value
    }
    window <- validate_window(window)
    by_name <- symbol_name(by_expr, 'by')
    group_name <- symbol_name(group_expr, 'group')
    if (!length(predictor_exprs)) stop('irf requires at least one predictor')
    constants <- vapply(predictor_exprs, function(expr) {
        is.numeric(expr) && length(expr) == 1L &&
            isTRUE(all.equal(as.numeric(expr), 1))
    }, logical(1))
    symbols <- vapply(predictor_exprs, is.symbol, logical(1))
    if (any(!(constants | symbols)) || (any(constants) &&
            (length(predictor_exprs) != 1L || !all(constants)))) {
        stop('IRF predictors must be unquoted column names or the sole literal 1')
    }
    constant <- all(constants)
    predictors <- if (constant) character() else
        vapply(predictor_exprs, as.character, character(1))
    if (anyDuplicated(predictors)) stop('IRF predictors must be unique')

    legacy <- !is.null(k) || !is.null(bs) || !is.null(nonlinear) ||
        !missing(varying)
    new_supplied <- !missing(k_l) || !missing(k_t) || !missing(k_p) ||
        !missing(bs_l) || !missing(bs_t) || !missing(bs_p)
    time_name <- NULL
    if (legacy) {
        if (new_supplied) {
            stop('Legacy k, bs, nonlinear, and varying cannot be mixed with axis-specific controls')
        }
        if (length(predictor_exprs) != 1L) {
            stop('Legacy IRF controls support one predictor only')
        }
        nonlinear <- if (is.null(nonlinear)) FALSE else nonlinear
        if (length(nonlinear) != 1L || !is.logical(nonlinear) ||
                is.na(nonlinear)) {
            stop('nonlinear must be TRUE or FALSE')
        }
        time_name <- symbol_name(varying_expr, 'varying')
        tensor <- nonlinear || !is.null(time_name)
        k <- if (is.null(k)) 10 else k
        bs <- if (is.null(bs)) 'cr' else bs
        dimensions <- if (tensor) 2L else 1L
        if (!is.numeric(k) || any(!is.finite(k)) || any(k < 3) ||
                !(length(k) %in% c(1L, dimensions))) {
            stop('k must supply dimensions of at least 3 for lag and predictor')
        }
        if (length(k) == 1L && tensor) k <- rep.int(k, 2L)
        if (!is.character(bs) || !(length(bs) %in% c(1L, dimensions))) {
            stop('bs must supply a basis name for lag and predictor')
        }
        if (length(bs) == 1L && tensor) bs <- rep.int(bs, 2L)
        k_l <- validate_k(k[[1L]], 'k')
        bs_l <- validate_bs(bs[[1L]], 'bs')
        k_t <- if (is.null(time_name)) NULL else
            validate_k(k[[2L]], 'k')
        bs_t <- if (is.null(time_name)) validate_bs(bs_l, 'bs') else
            validate_bs(bs[[2L]], 'bs')
        k_p <- if (isTRUE(nonlinear)) list(validate_k(k[[2L]], 'k')) else
            rep.int(list(NULL), length(predictors))
        bs_p <- if (isTRUE(nonlinear)) list(validate_bs(bs[[2L]], 'bs')) else
            rep.int(list(bs_l), length(predictors))
    } else {
        k_l <- validate_k(k_l, 'k_l')
        k_t <- validate_k(k_t, 'k_t', allow_null=TRUE)
        bs_l <- validate_bs(bs_l, 'bs_l')
        bs_t <- validate_bs(bs_t, 'bs_t')
        predictor_count <- length(predictors)
        if (!predictor_count) {
            if (!is.null(k_p)) stop('irf(1) cannot define predictor bases')
            k_p <- bs_p <- list()
        } else {
            if (is.null(k_p)) {
                k_p <- rep.int(list(NULL), predictor_count)
            } else if (is.list(k_p)) {
                if (length(k_p) == 1L && predictor_count > 1L) {
                    k_p <- rep(k_p, predictor_count)
                }
            } else {
                k_p <- as.list(k_p)
            }
            if (length(k_p) != predictor_count) {
                stop('k_p must have one entry per IRF predictor')
            }
            k_p <- lapply(k_p, validate_k, argument='k_p', allow_null=TRUE)
            if (is.list(bs_p)) {
                if (length(bs_p) == 1L && predictor_count > 1L) {
                    bs_p <- rep(bs_p, predictor_count)
                }
            } else if (length(bs_p) == 1L) {
                bs_p <- rep.int(list(bs_p), predictor_count)
            } else {
                bs_p <- as.list(bs_p)
            }
            if (length(bs_p) != predictor_count) {
                stop('bs_p must have one entry per IRF predictor')
            }
            bs_p <- lapply(bs_p, validate_bs, argument='bs_p')
        }
    }
    out <- list(
        predictors=predictors,
        predictor=if (constant) '1' else if (length(predictors) == 1L)
            predictors[[1L]] else paste(predictors, collapse=':'),
        constant=constant,
        window=window,
        k_l=k_l,
        k_t=k_t,
        k_p=k_p,
        bs_l=bs_l,
        bs_t=bs_t,
        bs_p=bs_p,
        time=time_name,
        k=c(k_l, if (!is.null(k_t)) k_t else
            unlist(k_p[!vapply(k_p, is.null, logical(1))])),
        bs=c(bs_l, if (!is.null(k_t)) bs_t else
            unlist(bs_p[!vapply(k_p, is.null, logical(1))])),
        nonlinear=any(!vapply(k_p, is.null, logical(1))),
        varying=time_name,
        by=by_name,
        group=group_name,
        api=if (legacy) 'legacy' else 'axis'
    )
    class(out) <- 'cdrgam_irf_spec'
    out
}

.format_irf_spec <- function(spec) {
    format_vector <- function(value) {
        paste0('c(', paste(format(value, trim=TRUE, scientific=FALSE), collapse=', '), ')')
    }
    format_nullable <- function(value) {
        if (is.null(value)) 'NULL' else format(value, trim=TRUE)
    }
    predictors <- if (isTRUE(spec$constant)) '1' else spec$predictors
    arguments <- c(
        predictors,
        paste0('window=', format_vector(spec$window)),
        paste0('k_l=', spec$k_l),
        paste0('k_t=', format_nullable(spec$k_t)),
        paste0('k_p=', if (!length(spec$k_p)) 'NULL' else paste0(
            'list(', paste(vapply(
                spec$k_p, format_nullable, character(1)
            ), collapse=', '), ')'
        )),
        paste0('bs_l=', encodeString(spec$bs_l, quote='"')),
        paste0('bs_t=', encodeString(spec$bs_t, quote='"')),
        paste0('bs_p=', if (!length(spec$bs_p)) {
            encodeString(spec$bs_l, quote='"')
        } else paste0('list(', paste(vapply(
            spec$bs_p, encodeString, character(1), quote='"'
        ), collapse=', '), ')'))
    )
    if (!is.null(spec$by)) arguments <- c(arguments, paste0('by=', spec$by))
    if (!is.null(spec$group)) arguments <- c(arguments, paste0('group=', spec$group))
    paste0('irf(', paste(arguments, collapse=', '), ')')
}

.compose_cdr_formula <- function(ordinary_formula, specs) {
    response <- paste(deparse(ordinary_formula[[2L]]), collapse='')
    rhs <- paste(deparse(ordinary_formula[[3L]]), collapse='')
    if (length(specs)) {
        rhs <- paste(c(rhs, vapply(specs, .format_irf_spec, character(1))), collapse=' + ')
    }
    stats::as.formula(
        paste(response, '~', rhs),
        env=environment(ordinary_formula)
    )
}

.parse_cdr_formula <- function(formula, window=NULL) {
    if (!inherits(formula, 'formula') || length(formula) != 3L) {
        stop('formula must be a two-sided formula')
    }
    terms_object <- stats::terms(formula, specials='irf', keep.order=TRUE)
    labels <- attr(terms_object, 'term.labels')
    factors <- attr(terms_object, 'factors')
    special_variables <- attr(terms_object, 'specials')$irf
    special_text <- rownames(factors)[special_variables]
    special_terms <- match(special_text, labels)
    included <- !is.na(special_terms)
    excluded_text <- special_text[!included]
    if (length(excluded_text) && any(excluded_text != 'irf(1)')) {
        stop('Only irf(1) may be removed as a formula term')
    }
    included_text <- special_text[included]
    included_terms <- special_terms[included]

    specs <- lapply(included_text, function(text) {
        call <- str2lang(text)
        call[[1L]] <- irf
        eval(call, envir=environment(formula), enclos=parent.frame())
    })
    if (!is.null(window) && (length(window) != 2L || !is.numeric(window) ||
            !is.finite(window[[1L]]) || is.na(window[[2L]]) ||
            window[[2L]] < window[[1L]])) {
        stop('window must satisfy finite min <= max; max may be Inf')
    }
    inherited_window <- if (is.null(window)) NULL else as.numeric(window)
    specs <- lapply(specs, function(spec) {
        if (is.null(spec$window)) {
            spec$window <- if (is.null(inherited_window)) c(0, Inf) else
                inherited_window
        }
        spec
    })
    is_fixed_rate <- function(spec) {
        isTRUE(spec$constant) && is.null(spec$k_t) &&
            is.null(spec$by) && is.null(spec$group)
    }
    explicit_rate <- which(vapply(
        specs,
        is_fixed_rate,
        logical(1)
    ))
    if (length(explicit_rate) > 1L) {
        stop('The formula may contain at most one fixed irf(1) term')
    }
    rate_suppressed <- 'irf(1)' %in% excluded_text
    if (!length(explicit_rate) && !rate_suppressed) {
        if (length(specs)) {
            lag_k <- vapply(specs, function(spec) spec$k_l, integer(1))
            rate_window <- if (!is.null(inherited_window)) {
                inherited_window
            } else c(
                min(vapply(specs, function(spec) spec$window[[1L]], numeric(1))),
                max(vapply(specs, function(spec) spec$window[[2L]], numeric(1)))
            )
            rate_spec <- irf(
                1,
                window=rate_window,
                k_l=max(lag_k),
                bs_l=specs[[1L]]$bs_l
            )
        } else {
            rate_spec <- irf(
                1,
                window=if (is.null(inherited_window)) c(0, Inf) else
                    inherited_window
            )
        }
        rate_spec$implicit <- TRUE
        specs <- c(list(rate_spec), specs)
    }
    specs <- lapply(specs, function(spec) {
        if (is.null(spec$implicit)) spec$implicit <- FALSE
        spec
    })
    ordinary_indices <- setdiff(seq_along(labels), included_terms)
    ordinary_formula <- if (!length(ordinary_indices)) {
        response_text <- paste(deparse(formula[[2L]]), collapse='')
        intercept <- attr(terms_object, 'intercept')
        stats::as.formula(
            paste(response_text, '~', if (intercept) '1' else '0')
        )
    } else {
        ordinary_terms <- stats::drop.terms(
            terms_object,
            dropx=included_terms,
            keep.response=TRUE
        )
        stats::formula(ordinary_terms)
    }
    environment(ordinary_formula) <- environment(formula)
    normalized_formula <- .compose_cdr_formula(ordinary_formula, specs)
    list(
        formula=formula,
        ordinary_formula=ordinary_formula,
        irfs=specs,
        rate_suppressed=rate_suppressed,
        normalized_formula=normalized_formula
    )
}

.stream_keys <- function(data, columns) {
    if (!length(columns)) {
        return(rep.int('__all__', nrow(data)))
    }
    missing <- setdiff(columns, names(data))
    if (length(missing)) {
        stop('Missing series columns: ', paste(missing, collapse=', '))
    }
    values <- data[columns]
    if (anyNA(values)) {
        stop('series columns must not contain missing values')
    }
    do.call(paste, c(lapply(values, as.character), sep='\034'))
}

.build_history_links <- function(
        impulses,
        responses,
        series,
        impulse_time,
        response_time,
        window,
        history_length=Inf
) {
    impulse_keys <- .stream_keys(impulses, series)
    response_keys <- .stream_keys(responses, series)
    impulse_times <- impulses[[impulse_time]]
    response_times <- responses[[response_time]]
    if (!is.numeric(impulse_times) || any(!is.finite(impulse_times))) {
        stop('impulse_time must name a finite numeric column in impulses')
    }
    if (!is.numeric(response_times) || any(!is.finite(response_times))) {
        stop('response_time must name a finite numeric column in responses')
    }

    impulse_groups <- split(seq_len(nrow(impulses)), impulse_keys)
    response_groups <- split(seq_len(nrow(responses)), response_keys)
    counts <- integer(nrow(responses))
    group_cache <- vector('list', length(response_groups))
    names(group_cache) <- names(response_groups)
    for (key in names(response_groups)) {
        response_rows <- response_groups[[key]]
        impulse_rows <- impulse_groups[[key]]
        if (is.null(impulse_rows) || !length(impulse_rows)) {
            next
        }
        impulse_rows <- impulse_rows[order(impulse_times[impulse_rows])]
        times <- impulse_times[impulse_rows]
        lower <- response_times[response_rows] - window[[2L]]
        upper <- response_times[response_rows] - window[[1L]]
        starts <- findInterval(lower, times, left.open=TRUE) + 1L
        ends <- findInterval(upper, times)
        if (is.finite(history_length)) {
            starts <- pmax.int(starts, ends - history_length + 1L)
        }
        counts[response_rows] <- pmax.int(0L, ends - starts + 1L)
        group_cache[[key]] <- list(
            response_rows=response_rows,
            impulse_rows=impulse_rows,
            starts=starts,
            ends=ends
        )
    }

    offsets <- c(0, cumsum(counts))
    link_count <- offsets[[length(offsets)]]
    impulse_index <- integer(link_count)
    for (group in group_cache) {
        if (is.null(group)) {
            next
        }
        for (j in seq_along(group$response_rows)) {
            response_row <- group$response_rows[[j]]
            if (!counts[[response_row]]) {
                next
            }
            destination <- seq.int(
                offsets[[response_row]] + 1L,
                offsets[[response_row + 1L]]
            )
            impulse_index[destination] <- group$impulse_rows[
                group$starts[[j]]:group$ends[[j]]
            ]
        }
    }
    response_index <- rep.int(seq_len(nrow(responses)), counts)
    delay <- response_times[response_index] - impulse_times[impulse_index]
    list(
        response_index=response_index,
        impulse_index=impulse_index,
        delay=delay,
        counts=counts,
        offsets=offsets
    )
}

.choose_history_layout <- function(counts, strategy) {
    strategy <- match.arg(strategy, c('auto', 'dense', 'ragged'))
    maximum <- if (length(counts)) max(counts) else 0L
    links <- sum(counts)
    slots <- length(counts) * maximum
    padding <- if (slots) 1 - links / slots else 0
    selected <- strategy
    if (identical(strategy, 'auto')) {
        selected <- if (maximum <= 64L && padding <= 0.5) 'dense' else 'ragged'
    }
    list(
        requested=strategy,
        selected=selected,
        responses=length(counts),
        links=links,
        maximum_history=maximum,
        median_history=stats::median(counts),
        dense_slots=slots,
        dense_padding=padding
    )
}

.compress_cdr_links <- function(
        links,
        weights,
        n,
        k,
        bs,
        chunk_size,
        name,
        response_multiplier=NULL
) {
    if (!length(weights) || length(weights) != length(links$delay)) {
        stop('An IRF has no impulse-response links or invalid weights')
    }
    by_sum <- numeric(n)
    summed <- rowsum(weights, links$response_index, reorder=FALSE)
    by_sum[as.integer(rownames(summed))] <- summed[, 1L]
    effective_sum <- if (is.null(response_multiplier)) {
        by_sum
    } else {
        by_sum * response_multiplier
    }
    centered <- .cdr_constant_sum(effective_sum)
    values <- unique(links$delay)
    if (length(values) < k) {
        stop('An IRF has fewer unique delays than its basis dimension k')
    }
    knots <- as.numeric(stats::quantile(
        values,
        seq(0, 1, length.out=k),
        names=FALSE
    ))
    if (any(diff(knots) <= 0)) {
        stop('Computed IRF knots are not strictly increasing')
    }
    cdr_delta <- NULL
    spec <- mgcv::s(cdr_delta, k=k, bs=bs)
    marginal <- mgcv::smoothCon(
        spec,
        data=list(cdr_delta=knots),
        knots=if (identical(bs, 'ps')) NULL else list(cdr_delta=knots),
        absorb.cons=FALSE,
        scale.penalty=FALSE,
        n=length(knots)
    )[[1L]]
    marginal_dimension <- ncol(marginal$X)
    transform <- if (centered) {
        constraint <- .cdr_basis_mean_constraint(
            marginal,
            'cdr_delta',
            links$delay,
            chunk_size
        )
        .cdr_constraint_transform(constraint, marginal_dimension)
    } else {
        diag(marginal_dimension)
    }
    basis_dimension <- ncol(transform)
    design <- matrix(0, nrow=n, ncol=basis_dimension)
    max_basis_row_norm <- 0
    for (start in seq.int(1L, length(weights), by=chunk_size)) {
        end <- min(length(weights), start + chunk_size - 1L)
        rows <- start:end
        raw_basis <- mgcv::PredictMat(
            marginal,
            list(cdr_delta=links$delay[rows]),
            n=length(rows)
        )
        max_basis_row_norm <- max(
            max_basis_row_norm,
            max(rowSums(abs(raw_basis)))
        )
        basis <- raw_basis %*% transform
        accumulated <- rowsum(
            basis * weights[rows],
            links$response_index[rows],
            reorder=FALSE
        )
        response_rows <- as.integer(rownames(accumulated))
        design[response_rows, ] <- design[response_rows, , drop=FALSE] +
            accumulated
    }
    basis_scale <- max_basis_row_norm^2
    penalties <- marginal$S
    penalty_scales <- numeric(length(penalties))
    for (i in seq_along(penalties)) {
        penalty_norm <- norm(penalties[[i]], type='I')
        penalty_scales[[i]] <- penalty_norm / basis_scale
        penalties[[i]] <- penalties[[i]] * basis_scale / penalty_norm
        penalties[[i]] <- crossprod(
            transform,
            penalties[[i]] %*% transform
        )
    }
    out <- list(
        X=design,
        S=penalties,
        rank=vapply(penalties, function(penalty) qr(penalty)$rank, integer(1)),
        null.space.dim=basis_dimension - qr(Reduce(`+`, penalties))$rank,
        knots=knots,
        basis=marginal,
        S.scale=penalty_scales,
        by_sum=effective_sum,
        name=name,
        type='linear',
        predictor_knots=NULL,
        transform=if (centered) transform else NULL,
        constraints=if (centered) 'linear-functional-centering' else character()
    )
    class(out) <- c('cdrgam_term', 'cdr_compressed_term')
    out
}

.compress_cdr_tensor_links <- function(
        links,
        predictor_values,
        n,
        k,
        bs,
        chunk_size,
        name,
        link_weights=NULL,
        type='nonlinear',
        response_multiplier=NULL
) {
    if (!length(links$delay) ||
            length(predictor_values) != length(links$delay)) {
        stop('A nonlinear IRF has no links or invalid predictor values')
    }
    if (!is.null(link_weights) &&
            (length(link_weights) != length(links$delay) ||
             any(!is.finite(link_weights)))) {
        stop('A tensor IRF has invalid convolution weights')
    }
    by_sum <- if (is.null(link_weights)) {
        tabulate(links$response_index, nbins=n)
    } else {
        output <- numeric(n)
        summed <- rowsum(link_weights, links$response_index, reorder=FALSE)
        output[as.integer(rownames(summed))] <- summed[, 1L]
        output
    }
    effective_sum <- if (is.null(response_multiplier)) {
        by_sum
    } else {
        by_sum * response_multiplier
    }
    centered <- .cdr_constant_sum(effective_sum)
    delay_values <- unique(links$delay)
    predictor_unique <- unique(predictor_values)
    if (length(delay_values) < k[[1L]] ||
            length(predictor_unique) < k[[2L]]) {
        stop('A nonlinear IRF has fewer unique values than a basis dimension')
    }
    delay_knots <- as.numeric(stats::quantile(
        delay_values,
        seq(0, 1, length.out=k[[1L]]),
        names=FALSE
    ))
    predictor_knots <- as.numeric(stats::quantile(
        predictor_unique,
        seq(0, 1, length.out=k[[2L]]),
        names=FALSE
    ))
    if (any(diff(delay_knots) <= 0) || any(diff(predictor_knots) <= 0)) {
        stop('Computed tensor knots are not strictly increasing')
    }
    cdr_delay <- NULL
    cdr_value <- NULL
    spec <- mgcv::te(cdr_delay, cdr_value, k=k, bs=bs)
    knot_grid <- expand.grid(
        cdr_delay=delay_knots,
        cdr_value=predictor_knots
    )
    tensor <- mgcv::smoothCon(
        spec,
        data=knot_grid,
        knots=list(
            cdr_delay=delay_knots,
            cdr_value=predictor_knots
        ),
        absorb.cons=FALSE,
        scale.penalty=FALSE,
        n=nrow(knot_grid)
    )[[1L]]
    full_basis_dimension <- prod(k)
    # A nonlinear predictor surface is defined as a deviation from the
    # predictor-independent rate IRF. Enforce zero mean over observed linked
    # predictor values separately across the lag marginal. This also prevents
    # multiple nonlinear IRFs from sharing identical lag-only columns.
    delay_constraint_basis <- mgcv::PredictMat(
        tensor$margin[[1L]],
        list(cdr_delay=delay_knots),
        n=length(delay_knots)
    )
    predictor_basis_mean <- numeric(k[[2L]])
    for (start in seq.int(1L, length(predictor_values), by=chunk_size)) {
        end <- min(length(predictor_values), start + chunk_size - 1L)
        rows <- start:end
        predictor_basis_mean <- predictor_basis_mean + colSums(
            mgcv::PredictMat(
                tensor$margin[[2L]],
                list(cdr_value=predictor_values[rows]),
                n=length(rows)
            )
        )
    }
    predictor_basis_mean <- predictor_basis_mean / length(predictor_values)
    constraint <- mgcv::tensor.prod.model.matrix(list(
        delay_constraint_basis,
        matrix(
            predictor_basis_mean,
            nrow=length(delay_knots),
            ncol=k[[2L]],
            byrow=TRUE
        )
    ))
    constraint_qr <- qr(t(constraint))
    constraint_rank <- constraint_qr$rank
    complete_q <- qr.Q(constraint_qr, complete=TRUE)
    transform <- complete_q[, seq.int(
        constraint_rank + 1L,
        full_basis_dimension
    ), drop=FALSE]
    constraints <- 'predictor-deviation-centering'
    if (centered) {
        centering_constraint <- tensor$C %*% transform
        centering_transform <- .cdr_constraint_transform(
            centering_constraint,
            ncol(transform)
        )
        if (ncol(centering_transform) < ncol(transform)) {
            transform <- transform %*% centering_transform
            constraints <- c(constraints, 'linear-functional-centering')
        }
    }
    basis_dimension <- ncol(transform)
    design <- matrix(0, nrow=n, ncol=basis_dimension)
    max_basis_row_norm <- 0
    for (start in seq.int(1L, length(predictor_values), by=chunk_size)) {
        end <- min(length(predictor_values), start + chunk_size - 1L)
        rows <- start:end
        basis <- mgcv::PredictMat(
            tensor,
            list(
                cdr_delay=links$delay[rows],
                cdr_value=predictor_values[rows]
            ),
            n=length(rows)
        ) %*% transform
        if (!is.null(link_weights)) basis <- basis * link_weights[rows]
        max_basis_row_norm <- max(
            max_basis_row_norm,
            max(rowSums(abs(basis)))
        )
        accumulated <- rowsum(
            basis,
            links$response_index[rows],
            reorder=FALSE
        )
        response_rows <- as.integer(rownames(accumulated))
        design[response_rows, ] <- design[response_rows, , drop=FALSE] +
            accumulated
    }
    basis_scale <- max_basis_row_norm^2
    penalties <- lapply(tensor$S, function(penalty) {
        crossprod(transform, penalty %*% transform)
    })
    penalty_scales <- numeric(length(penalties))
    for (i in seq_along(penalties)) {
        penalty_norm <- norm(penalties[[i]], type='I')
        penalty_scales[[i]] <- penalty_norm / basis_scale
        penalties[[i]] <- penalties[[i]] * basis_scale / penalty_norm
    }
    out <- list(
        X=design,
        S=penalties,
        rank=vapply(penalties, function(penalty) qr(penalty)$rank, integer(1)),
        null.space.dim=basis_dimension - qr(Reduce(`+`, penalties))$rank,
        knots=delay_knots,
        predictor_knots=predictor_knots,
        basis=tensor,
        transform=transform,
        S.scale=penalty_scales,
        by_sum=effective_sum,
        name=name,
        type=type,
        constraints=constraints
    )
    class(out) <- c('cdrgam_term', 'cdr_compressed_term')
    out
}

.compress_cdr_multi_tensor_links <- function(
        links,
        axes,
        linear_weights,
        n,
        chunk_size,
        name,
        response_multiplier=NULL
) {
    if (length(axes) < 2L || !length(links$delay) ||
            length(linear_weights) != length(links$delay)) {
        stop('A tensor IRF has no links or invalid axes')
    }
    if (any(!is.finite(linear_weights))) {
        stop('A tensor IRF has invalid linear predictor weights')
    }
    link_count <- length(links$delay)
    for (axis in axes) {
        if (length(axis$values) != link_count ||
                !is.numeric(axis$values) || any(!is.finite(axis$values))) {
            stop('Tensor IRF axes must contain one finite numeric value per link')
        }
    }
    by_sum <- numeric(n)
    summed <- rowsum(linear_weights, links$response_index, reorder=FALSE)
    by_sum[as.integer(rownames(summed))] <- summed[, 1L]
    effective_sum <- if (is.null(response_multiplier)) by_sum else
        by_sum * response_multiplier
    centered <- .cdr_constant_sum(effective_sum)

    internal_names <- paste0('cdr_axis_', seq_along(axes))
    representatives <- lapply(axes, function(axis) {
        unique_values <- unique(axis$values)
        if (length(unique_values) < axis$k) {
            stop(
                'IRF axis ', axis$label, ' has fewer unique values than ',
                'its basis dimension'
            )
        }
        values <- as.numeric(stats::quantile(
            unique_values,
            seq(0, 1, length.out=max(axis$k, min(20L, length(unique_values)))),
            names=FALSE
        ))
        unique(values)
    })
    names(representatives) <- internal_names
    basis_call <- as.call(c(
        list(quote(mgcv::te)),
        lapply(internal_names, as.name),
        list(
            k=vapply(axes, `[[`, integer(1), 'k'),
            bs=vapply(axes, `[[`, character(1), 'bs')
        )
    ))
    specification <- eval(basis_call)
    construction_data <- do.call(expand.grid, c(
        representatives,
        list(KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
    ))
    tensor <- mgcv::smoothCon(
        specification,
        data=construction_data,
        absorb.cons=FALSE,
        scale.penalty=FALSE,
        n=nrow(construction_data)
    )[[1L]]
    margin_dimensions <- vapply(tensor$margin, function(margin) {
        ncol(mgcv::PredictMat(
            margin,
            setNames(list(representatives[[margin$term]]), margin$term),
            n=length(representatives[[margin$term]])
        ))
    }, integer(1))
    full_dimension <- prod(margin_dimensions)

    constraints <- list()
    constraint_labels <- character()
    for (target in seq.int(2L, length(axes))) {
        mean_basis <- numeric(margin_dimensions[[target]])
        for (start in seq.int(1L, link_count, by=chunk_size)) {
            end <- min(link_count, start + chunk_size - 1L)
            rows <- start:end
            mean_basis <- mean_basis + colSums(mgcv::PredictMat(
                tensor$margin[[target]],
                setNames(list(axes[[target]]$values[rows]), internal_names[[target]]),
                n=length(rows)
            ))
        }
        mean_basis <- mean_basis / link_count
        other <- setdiff(seq_along(axes), target)
        index_grid <- do.call(expand.grid, c(
            lapply(other, function(index) seq_len(margin_dimensions[[index]])),
            list(KEEP.OUT.ATTRS=FALSE)
        ))
        matrices <- vector('list', length(axes))
        for (index in seq_along(axes)) {
            matrices[[index]] <- if (index == target) {
                matrix(mean_basis, nrow=nrow(index_grid),
                    ncol=length(mean_basis), byrow=TRUE)
            } else {
                position <- match(index, other)
                diag(margin_dimensions[[index]])[
                    index_grid[[position]], , drop=FALSE
                ]
            }
        }
        constraints[[length(constraints) + 1L]] <-
            mgcv::tensor.prod.model.matrix(matrices)
        constraint_labels <- c(
            constraint_labels,
            paste0('deviation-centering:', axes[[target]]$label)
        )
    }
    transform <- .cdr_constraint_transform(
        do.call(rbind, constraints),
        full_dimension
    )
    design <- matrix(0, nrow=n, ncol=ncol(transform))
    max_basis_row_norm <- 0
    for (start in seq.int(1L, link_count, by=chunk_size)) {
        end <- min(link_count, start + chunk_size - 1L)
        rows <- start:end
        newdata <- setNames(lapply(axes, function(axis) axis$values[rows]),
            internal_names)
        basis <- mgcv::PredictMat(tensor, newdata, n=length(rows)) %*%
            transform
        basis <- basis * linear_weights[rows]
        max_basis_row_norm <- max(
            max_basis_row_norm,
            max(rowSums(abs(basis)))
        )
        accumulated <- rowsum(
            basis,
            links$response_index[rows],
            reorder=FALSE
        )
        response_rows <- as.integer(rownames(accumulated))
        design[response_rows, ] <- design[response_rows, , drop=FALSE] +
            accumulated
    }
    if (centered && ncol(design)) {
        centering_transform <- .cdr_constraint_transform(
            matrix(colSums(design), nrow=1L),
            ncol(design)
        )
        if (ncol(centering_transform) < ncol(design)) {
            design <- design %*% centering_transform
            transform <- transform %*% centering_transform
            constraint_labels <- c(
                constraint_labels,
                'linear-functional-centering'
            )
        }
    }
    basis_scale <- max_basis_row_norm^2
    if (!is.finite(basis_scale) || basis_scale <= 0) {
        stop('A tensor IRF produced a zero or non-finite basis')
    }
    penalties <- lapply(tensor$S, function(penalty) {
        crossprod(transform, penalty %*% transform)
    })
    penalty_scales <- numeric(length(penalties))
    for (index in seq_along(penalties)) {
        penalty_norm <- norm(penalties[[index]], type='I')
        if (!is.finite(penalty_norm) || penalty_norm <= 0) {
            stop('A tensor IRF produced a zero or non-finite penalty')
        }
        penalty_scales[[index]] <- penalty_norm / basis_scale
        penalties[[index]] <- penalties[[index]] * basis_scale / penalty_norm
    }
    axis_metadata <- lapply(seq_along(axes), function(index) {
        axis <- axes[[index]]
        probabilities <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
        axis$summary <- list(
            mean=mean(axis$values),
            sd=stats::sd(axis$values),
            quantiles=stats::setNames(
                as.numeric(stats::quantile(
                    axis$values, probabilities, names=FALSE
                )),
                format(probabilities, trim=TRUE, scientific=FALSE)
            )
        )
        axis$values <- NULL
        axis$internal <- internal_names[[index]]
        axis$grid <- representatives[[index]]
        axis
    })
    out <- list(
        X=design,
        S=penalties,
        rank=vapply(penalties, function(penalty) qr(penalty)$rank, integer(1)),
        null.space.dim=ncol(design) - qr(Reduce(`+`, penalties))$rank,
        knots=representatives[[1L]],
        predictor_knots=representatives[[2L]],
        axis=axis_metadata,
        basis=tensor,
        transform=transform,
        S.scale=penalty_scales,
        by_sum=effective_sum,
        name=name,
        type='tensor',
        constraints=constraint_labels
    )
    class(out) <- c('cdrgam_term', 'cdr_compressed_term')
    out
}

.group_cdr_term <- function(term, group, group_name) {
    if (anyNA(group)) {
        stop('Grouped-IRF factors must not contain missing values')
    }
    group <- as.factor(group)
    levels <- levels(group)
    if (length(levels) < 2L) {
        stop('A grouped IRF deviation requires at least two observed levels')
    }
    base_dimension <- ncol(term$X)
    group_count <- length(levels)
    # Keep the response-by-basis matrix compact. Backends expand it either to
    # a dense matrix (native mgcv/reference block solver) or directly to a
    # sparse grouped matrix. Since each response belongs to exactly one group,
    # the infinity norm of the expanded design equals that of the base design.
    penalties <- term$S
    penalties[[length(penalties) + 1L]] <- diag(base_dimension)
    design_scale <- norm(term$X, type='I')^2
    penalty_scales <- numeric(length(penalties))
    for (i in seq_along(penalties)) {
        penalty_norm <- norm(penalties[[i]], type='I')
        penalty_scales[[i]] <- penalty_norm / design_scale
        penalties[[i]] <- penalties[[i]] * design_scale / penalty_norm
    }
    term$S <- penalties
    term$rank <- vapply(
        penalties,
        function(penalty) group_count * qr(penalty)$rank,
        integer(1)
    )
    term$null.space.dim <- 0L
    term$S.scale <- penalty_scales
    term$group <- group_name
    term$group_levels <- levels
    term$group_index <- as.integer(group)
    term$base_dimension <- base_dimension
    term$expanded_dimension <- base_dimension * group_count
    term$type <- paste0(term$type, '_group')
    term$name <- paste0(term$name, '|', group_name)
    term
}

.materialize_cdr_term <- function(term, sparse=FALSE) {
    if (is.null(term$group) || !length(term$group_levels)) {
        if (isTRUE(sparse)) {
            term$X <- Matrix::Matrix(term$X, sparse=TRUE)
            term$S <- lapply(term$S, Matrix::Matrix, sparse=TRUE)
        }
        return(term)
    }
    n <- nrow(term$X)
    k <- term$base_dimension
    groups <- length(term$group_levels)
    row_index <- rep(seq_len(n), each=k)
    column_index <- rep((term$group_index - 1L) * k, each=k) +
        rep.int(seq_len(k), times=n)
    values <- as.vector(t(term$X))
    if (isTRUE(sparse)) {
        term$X <- Matrix::sparseMatrix(
            i=row_index,
            j=column_index,
            x=values,
            dims=c(n, groups * k)
        )
        identity <- Matrix::Diagonal(groups)
        term$S <- lapply(term$S, function(penalty) {
            Matrix::kronecker(identity, Matrix::Matrix(penalty, sparse=TRUE))
        })
    } else {
        design <- matrix(0, nrow=n, ncol=groups * k)
        for (i in seq_len(groups)) {
            rows <- which(term$group_index == i)
            columns <- (i - 1L) * k + seq_len(k)
            design[rows, columns] <- term$X[rows, , drop=FALSE]
        }
        term$X <- design
        term$S <- lapply(term$S, function(penalty) {
            kronecker(diag(groups), penalty)
        })
    }
    term
}

.materialize_cdr_design <- function(design, sparse=FALSE) {
    design$terms <- lapply(
        design$terms,
        .materialize_cdr_term,
        sparse=sparse
    )
    names(design$terms) <- names(design$specification)
    # Specification names are not guaranteed to have been assigned; preserve
    # the compiler's term labels in that case.
    if (is.null(names(design$terms)) || any(!nzchar(names(design$terms)))) {
        names(design$terms) <- vapply(
            design$terms,
            `[[`,
            character(1),
            'name'
        )
    }
    design
}

.warn_unattested_random_effect_combinations <- function(formula, data) {
    labels <- attr(stats::terms(formula, keep.order=TRUE), 'term.labels')
    for (label in labels) {
        call <- tryCatch(str2lang(label), error=function(error) NULL)
        if (!is.call(call) || !identical(as.character(call[[1L]]), 's')) next
        arguments <- as.list(call)[-1L]
        argument_names <- names(arguments)
        bs_index <- which(argument_names == 'bs')
        if (length(bs_index) != 1L ||
                !identical(tryCatch(eval(arguments[[bs_index]]),
                    error=function(error) NULL), 're')) next
        unnamed <- which(is.na(argument_names) | argument_names == '')
        if (length(unnamed) < 2L ||
                !all(vapply(arguments[unnamed], is.symbol, logical(1)))) next
        variables <- vapply(arguments[unnamed], as.character, character(1))
        if (!all(variables %in% names(data)) ||
                !all(vapply(data[variables], is.factor, logical(1)))) next
        values <- data[variables]
        complete <- stats::complete.cases(values)
        marginal_levels <- vapply(
            values, function(value) nlevels(droplevels(value[complete])), integer(1)
        )
        cartesian <- prod(as.double(marginal_levels))
        attested <- if (any(complete)) {
            nrow(unique(values[complete, , drop=FALSE]))
        } else 0
        empty <- cartesian - attested
        if (empty > 0) warning(
            'Random-effect term ', label, ' retains ',
            format(empty, scientific=FALSE, big.mark=','),
            ' unattested factor combination(s) as zero design columns (',
            format(attested, scientific=FALSE, big.mark=','), ' attested of ',
            format(cartesian, scientific=FALSE, big.mark=','),
            '). Create a composite factor from the attested combinations and ',
            'use a one-factor bs="re" term to avoid those columns.',
            call.=FALSE
        )
    }
    invisible(NULL)
}

#' Compile untiled impulse and response streams for CDR-GAM fitting
#'
#' @param formula Extended model formula containing ordinary `mgcv` terms and
#'   optional [irf()] terms. An implicit `irf(1)` is added unless suppressed.
#' @param window Default inclusive lag window inherited by every [irf()] term
#'   that does not define its own `window`. `NULL` retains the term-level
#'   default `c(0, Inf)`.
#' @param impulses Data frame with one row per impulse.
#' @param responses Data frame with one row per response.
#' @param series Character vector of columns identifying independent series.
#' @param impulse_time Name of the impulse-time column.
#' @param response_time Name of the response-time column.
#' @param history One of `"auto"`, `"dense"`, or `"ragged"`.
#' @param history_length Maximum number of most-recent impulses retained per
#'   response and series inside each IRF's lag window. The default is
#'   unlimited.
#' @param chunk_size Maximum history links transformed at once in ragged mode,
#'   or responses transformed at once in dense mode.
#' @param rescale_predictors Divide continuous numeric predictors and internal
#'   lag coordinates by training-data standard deviations. Scaling is disabled
#'   by default. Binary, constant, categorical, grouping, and series variables
#'   are not changed.
#' @param drop.unused.levels Drop factor levels without training observations.
#'   The default matches [mgcv::gam()]. Retained levels remain part of the
#'   fitted factor vocabulary. For `cdrgam.fit()`, use `NULL` to inherit the
#'   choice already compiled into `design`; a supplied value must match it.
#' @param quiet Suppress the preparation plan message.
#' @return A reusable `cdrgam_design` object. Its `simplifications` data frame
#'   records automatic reductions of axis basis dimensions and
#'   lower-dimensional representations of discrete predictor axes.
#' @export
prepare_cdrgam <- function(
        formula,
        impulses,
        responses,
        window=NULL,
        series=character(),
        impulse_time='time',
        response_time='time',
        history=c('auto', 'dense', 'ragged'),
        history_length=Inf,
        chunk_size=10000,
        rescale_predictors=FALSE,
        quiet=FALSE,
        drop.unused.levels=TRUE
) {
    if (!is.data.frame(impulses) || !is.data.frame(responses)) {
        stop('impulses and responses must be data frames')
    }
    history <- match.arg(history)
    if (length(drop.unused.levels) != 1L ||
            !is.logical(drop.unused.levels) || is.na(drop.unused.levels)) {
        stop('drop.unused.levels must be TRUE or FALSE')
    }
    if (isTRUE(drop.unused.levels)) {
        impulses <- droplevels(impulses)
        responses <- droplevels(responses)
    }
    if (length(history_length) != 1L || !is.numeric(history_length) ||
            is.na(history_length) || history_length <= 0 ||
            (is.finite(history_length) && history_length != as.integer(history_length))) {
        stop('history_length must be a positive integer or Inf')
    }
    parsed <- .parse_cdr_formula(formula, window=window)
    parsed$irfs <- lapply(parsed$irfs, function(specification) {
        if (!is.null(specification$k_t) && is.null(specification$time)) {
            specification$time <- response_time
            specification$varying <- response_time
        }
        specification
    })
    parsed$normalized_formula <- .compose_cdr_formula(
        parsed$ordinary_formula,
        parsed$irfs
    )
    .warn_unattested_random_effect_combinations(
        parsed$ordinary_formula, responses
    )
    response_name <- all.vars(parsed$ordinary_formula[[2L]])
    if (length(response_name) != 1L || !(response_name %in% names(responses))) {
        stop('The formula response must name one column in responses')
    }
    scaling <- .cdr_prepare_scaling(
        parsed,
        impulses,
        responses,
        response_name,
        series,
        impulse_time,
        response_time,
        rescale_predictors
    )
    model_impulses <- .cdr_apply_scaling(impulses, scaling, 'impulses')
    model_responses <- .cdr_apply_scaling(responses, scaling, 'responses')
    terms <- vector('list', length(parsed$irfs))
    plans <- vector('list', length(parsed$irfs))
    simplifications <- list()
    record_simplification <- function(
            term,
            axis,
            action,
            requested,
            effective,
            reason
    ) {
        simplifications[[length(simplifications) + 1L]] <<- data.frame(
            term=term,
            axis=axis,
            action=action,
            requested=as.character(requested),
            effective=as.character(effective),
            reason=reason,
            stringsAsFactors=FALSE
        )
    }
    for (i in seq_along(parsed$irfs)) {
        spec <- parsed$irfs[[i]]
        predictors <- spec$predictors
        missing_predictors <- setdiff(predictors, names(impulses))
        if (length(missing_predictors)) {
            stop('IRF predictors not found in impulses: ',
                paste(missing_predictors, collapse=', '))
        }
        predictor_values <- lapply(predictors, function(predictor) {
            value <- model_impulses[[predictor]]
            if (!is.numeric(value) || any(!is.finite(value))) {
                stop('IRF predictors must be finite numeric columns')
            }
            value
        })
        names(predictor_values) <- predictors
        by_values <- NULL
        if (!is.null(spec$by)) {
            if (!(spec$by %in% names(responses))) {
                stop('By covariate not found in responses: ', spec$by)
            }
            by_values <- model_responses[[spec$by]]
            if (!is.numeric(by_values) || any(!is.finite(by_values))) {
                stop('By covariates must be finite numeric columns')
            }
        }
        links <- .build_history_links(
            impulses,
            responses,
            series,
            impulse_time,
            response_time,
            spec$window,
            history_length=history_length
        )
        model_links <- links
        model_links$delay <- links$delay / scaling$time_divisor
        term_name <- if (isTRUE(spec$constant)) 'irf(1)' else
            paste(predictors, collapse=':')
        delay_count <- length(unique(model_links$delay))
        if (delay_count >= 3L && spec$k_l > delay_count) {
            requested_k <- spec$k_l
            spec$k_l <- as.integer(delay_count)
            record_simplification(
                term_name,
                'lag',
                'basis_reduced',
                requested_k,
                spec$k_l,
                paste(delay_count, 'distinct linked delays')
            )
        }
        for (predictor_index in seq_along(predictors)) {
            if (is.null(spec$k_p[[predictor_index]])) next
            linked <- predictor_values[[predictor_index]][
                model_links$impulse_index
            ]
            value_count <- length(unique(linked))
            if (value_count <= 2L && value_count > 0L) {
                requested_k <- spec$k_p[[predictor_index]]
                spec$k_p[predictor_index] <- list(NULL)
                record_simplification(
                    term_name,
                    if (length(predictors) == 1L) 'predictor' else
                        paste0('predictor:', predictors[[predictor_index]]),
                    'nonlinear_to_linear',
                    requested_k,
                    'linear',
                    sprintf(
                        '%d distinct linked predictor value%s',
                        value_count,
                        if (value_count == 1L) '' else 's'
                    )
                )
            } else if (value_count >= 3L &&
                    spec$k_p[[predictor_index]] > value_count) {
                requested_k <- spec$k_p[[predictor_index]]
                spec$k_p[[predictor_index]] <- as.integer(value_count)
                record_simplification(
                    term_name,
                    if (length(predictors) == 1L) 'predictor' else
                        paste0('predictor:', predictors[[predictor_index]]),
                    'basis_reduced',
                    requested_k,
                    spec$k_p[[predictor_index]],
                    paste(value_count,
                        'distinct linked predictor values')
                )
            } else if (value_count < 1L) {
                stop(
                    'A nonlinear IRF for ', term_name,
                    ' has no linked predictor values',
                    call.=FALSE
                )
            }
        }
        time_values <- NULL
        if (!is.null(spec$k_t)) {
            if (!(spec$time %in% names(responses))) {
                stop('IRF time axis not found in responses: ', spec$time)
            }
            time_values <- model_responses[[spec$time]]
            if (!is.numeric(time_values) || any(!is.finite(time_values))) {
                stop('IRF time axes must be finite numeric columns')
            }
            time_count <- length(unique(time_values[model_links$response_index]))
            if (time_count < 3L) {
                stop('A nonstationary IRF requires at least three time values')
            }
            if (spec$k_t > time_count) {
                requested_k <- spec$k_t
                spec$k_t <- as.integer(time_count)
                record_simplification(
                    term_name, 'time', 'basis_reduced', requested_k,
                    spec$k_t, paste(time_count, 'distinct linked time values')
                )
            }
        }
        smooth_predictors <- which(!vapply(spec$k_p, is.null, logical(1)))
        linear_predictors <- setdiff(seq_along(predictors), smooth_predictors)
        link_weights <- rep.int(1, length(model_links$delay))
        for (predictor_index in linear_predictors) {
            link_weights <- link_weights * predictor_values[[predictor_index]][
                model_links$impulse_index
            ]
        }
        axes <- list(list(
            role='lag', variable=NA_character_, label='lag',
            values=model_links$delay, k=spec$k_l, bs=spec$bs_l
        ))
        if (!is.null(spec$k_t)) {
            axes[[length(axes) + 1L]] <- list(
                role='time', variable=spec$time, label=spec$time,
                values=time_values[model_links$response_index],
                k=spec$k_t, bs=spec$bs_t
            )
        }
        for (predictor_index in smooth_predictors) {
            axes[[length(axes) + 1L]] <- list(
                role='predictor', variable=predictors[[predictor_index]],
                label=predictors[[predictor_index]],
                values=predictor_values[[predictor_index]][
                    model_links$impulse_index
                ],
                k=spec$k_p[[predictor_index]],
                bs=spec$bs_p[[predictor_index]]
            )
        }
        spec$nonlinear <- length(smooth_predictors) > 0L
        spec$varying <- if (is.null(spec$k_t)) NULL else spec$time
        spec$k <- c(spec$k_l, if (!is.null(spec$k_t)) spec$k_t,
            unlist(spec$k_p[smooth_predictors]))
        spec$bs <- c(spec$bs_l, if (!is.null(spec$k_t)) spec$bs_t,
            unlist(spec$bs_p[smooth_predictors]))
        parsed$irfs[[i]] <- spec
        plan <- .choose_history_layout(links$counts, history)
        plans[[i]] <- c(list(term=term_name, window=spec$window), plan)
        if (length(axes) > 1L) {
            terms[[i]] <- .compress_cdr_multi_tensor_links(
                model_links,
                axes=axes,
                linear_weights=link_weights,
                n=nrow(responses),
                chunk_size=chunk_size,
                name=term_name,
                response_multiplier=by_values
            )
            if (identical(spec$api, 'legacy')) {
                terms[[i]]$type <- if (!is.null(spec$k_t)) 'varying' else
                    'nonlinear'
                if (length(axes) == 2L) {
                    terms[[i]]$constraints <- unique(c(
                        'predictor-deviation-centering',
                        terms[[i]]$constraints
                    ))
                }
            }
            terms[[i]]$linear_predictors <- predictors[linear_predictors]
            if (!is.null(spec$k_t)) {
                terms[[i]]$varying <- spec$time
                terms[[i]]$name <- paste0(term_name, '~', spec$time)
            }
        } else if (identical(plan$selected, 'dense')) {
            width <- plan$maximum_history
            if (!width) {
                stop('No impulses fall within the requested IRF window')
            }
            knots <- as.numeric(stats::quantile(
                unique(model_links$delay),
                seq(0, 1, length.out=spec$k_l),
                names=FALSE
            ))
            delay_matrix <- matrix(knots[[1L]], nrow=nrow(responses), ncol=width)
            weight_matrix <- matrix(0, nrow=nrow(responses), ncol=width)
            positions <- sequence(model_links$counts)
            index <- cbind(model_links$response_index, positions)
            delay_matrix[index] <- model_links$delay
            weight_matrix[index] <- link_weights
            terms[[i]] <- compress_cdr_smooth(
                delay_matrix,
                weight_matrix,
                k=spec$k_l,
                bs=spec$bs_l,
                knots=knots,
                chunk_size=chunk_size,
                name=term_name,
                response_multiplier=by_values,
                constraint_delays=model_links$delay
            )
        } else {
            terms[[i]] <- .compress_cdr_links(
                model_links,
                link_weights,
                n=nrow(responses),
                k=spec$k_l,
                bs=spec$bs_l,
                chunk_size=chunk_size,
                name=term_name,
                response_multiplier=by_values
            )
        }
        if (!is.null(spec$by)) {
            terms[[i]]$X <- terms[[i]]$X * by_values
            terms[[i]]$by <- spec$by
            terms[[i]]$name <- paste0(terms[[i]]$name, ':', spec$by)
        }
        if (!is.null(spec$group)) {
            if (!(spec$group %in% names(responses))) {
                stop('Grouped-IRF factor not found in responses: ', spec$group)
            }
            terms[[i]] <- .group_cdr_term(
                terms[[i]],
                responses[[spec$group]],
                spec$group
            )
        }
        term_scales <- .cdr_irf_scales(spec, scaling)
        terms[[i]]$lag_scale <- term_scales$lag
        terms[[i]]$predictor_scale <- term_scales$predictor
        terms[[i]]$amplitude_scale <- term_scales$amplitude
        if (length(linear_predictors)) {
            summaries <- lapply(linear_predictors, function(predictor_index) {
                values <- predictor_values[[predictor_index]][model_links$impulse_index]
                scale <- unname(term_scales$axes[[predictors[[predictor_index]]]])
                values <- values * scale
                probabilities <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
                list(
                    mean=mean(values), sd=stats::sd(values),
                    quantiles=stats::setNames(
                        as.numeric(stats::quantile(
                            values, probabilities, names=FALSE
                        )),
                        format(probabilities, trim=TRUE, scientific=FALSE)
                    )
                )
            })
            names(summaries) <- predictors[linear_predictors]
            terms[[i]]$linear_predictor_summaries <- summaries
        }
        if (!is.null(terms[[i]]$axis)) {
            terms[[i]]$axis <- lapply(terms[[i]]$axis, function(axis) {
                axis$scale <- switch(
                    axis$role,
                    lag=term_scales$lag,
                    time=term_scales$time,
                    predictor=unname(term_scales$axes[[axis$variable]])
                )
                axis$summary$mean <- axis$summary$mean * axis$scale
                axis$summary$sd <- axis$summary$sd * abs(axis$scale)
                axis$summary$quantiles <- axis$summary$quantiles * axis$scale
                axis
            })
            terms[[i]]$predictor_scale <- terms[[i]]$axis[[2L]]$scale
        }
    }
    identifiability <- list(
        constraints=lapply(terms, function(term) term$constraints),
        rate=list(requested=if (any(vapply(
            parsed$irfs,
            function(spec) {
                isTRUE(spec$constant) && !isTRUE(spec$implicit) &&
                    is.null(spec$k_t) && is.null(spec$by) &&
                    is.null(spec$group)
            },
            logical(1)
        ))) 'explicit' else if (parsed$rate_suppressed) 'suppressed' else 'implicit',
        status=if (parsed$rate_suppressed) 'suppressed' else 'retained',
        reason=NULL)
    )
    rate_index <- which(vapply(
        parsed$irfs,
        function(spec) {
            isTRUE(spec$constant) && is.null(spec$k_t) &&
                is.null(spec$by) && is.null(spec$group)
        },
        logical(1)
    ))
    if (length(rate_index)) {
        rate_term <- terms[[rate_index]]
        scale <- max(1, max(abs(rate_term$X)))
        degenerate <- !ncol(rate_term$X) ||
            max(abs(rate_term$X)) <= scale * .Machine$double.eps * 1000
        duplicate <- integer()
        if (isTRUE(parsed$irfs[[rate_index]]$implicit) && !degenerate) {
            candidates <- setdiff(seq_along(terms), rate_index)
            duplicate <- candidates[vapply(candidates, function(index) {
                candidate <- terms[[index]]
                if (!identical(dim(candidate$X), dim(rate_term$X))) {
                    return(FALSE)
                }
                candidate_scale <- max(abs(candidate$X))
                if (candidate_scale <= .Machine$double.eps * 1000) {
                    return(FALSE)
                }
                denominator <- sum(rate_term$X^2)
                if (!is.finite(denominator) || denominator <= 0) {
                    return(FALSE)
                }
                multiplier <- sum(candidate$X * rate_term$X) / denominator
                residual <- candidate$X - multiplier * rate_term$X
                max(abs(residual)) <= max(1, candidate_scale) *
                    sqrt(.Machine$double.eps)
            }, logical(1))]
        }
        if (degenerate) {
            if (!isTRUE(parsed$irfs[[rate_index]]$implicit)) {
                stop(
                    'Explicit irf(1) is unidentifiable: after the required ',
                    'centering constraint its design has no estimable ',
                    'variation. Remove irf(1) or revise the impulse history.',
                    call.=FALSE
                )
            }
        }
        if (degenerate || length(duplicate)) {
            reason <- if (degenerate) {
                'zero design after centering'
            } else {
                paste0('duplicates explicit term ', terms[[duplicate[[1L]]]]$name)
            }
            warning(
                'Removed the implicit irf(1): ', reason, '.',
                call.=FALSE
            )
            keep <- setdiff(seq_along(terms), rate_index)
            terms <- terms[keep]
            plans <- plans[keep]
            parsed$irfs <- parsed$irfs[keep]
            identifiability$constraints <- identifiability$constraints[keep]
            identifiability$rate$status <- 'removed'
            identifiability$rate$reason <- reason
            record_simplification(
                'irf(1)',
                'term',
                'term_removed',
                'retained',
                'removed',
                reason
            )
        }
    }
    if (!length(terms)) {
        stop('No identifiable IRF terms remain in the model')
    }
    names(terms) <- vapply(terms, `[[`, character(1), 'name')
    simplifications <- if (length(simplifications)) {
        do.call(rbind, simplifications)
    } else {
        data.frame(
            term=character(),
            axis=character(),
            action=character(),
            requested=character(),
            effective=character(),
            reason=character(),
            stringsAsFactors=FALSE
        )
    }
    out <- list(
        formula=formula,
        normalized_formula=parsed$normalized_formula,
        effective_formula=.compose_cdr_formula(parsed$ordinary_formula, parsed$irfs),
        ordinary_formula=parsed$ordinary_formula,
        response_name=response_name,
        responses=model_responses,
        terms=terms,
        specification=parsed$irfs,
        simplifications=simplifications,
        scaling=scaling,
        identifiability=identifiability,
        plan=plans,
        configuration=list(
            window=window,
            series=series,
            impulse_time=impulse_time,
            response_time=response_time,
            history=history,
            history_length=history_length,
            chunk_size=chunk_size,
            rescale_predictors=rescale_predictors,
            drop.unused.levels=drop.unused.levels
        ),
        stream=list(
            series=series,
            impulse_time=impulse_time,
            response_time=response_time,
            history_length=history_length
        )
    )
    class(out) <- 'cdrgam_design'
    if (!quiet) {
        if (isTRUE(scaling$enabled)) {
            message(sprintf(
                'Rescaling enabled: time and lag divisor %s; %d continuous predictor%s rescaled',
                format(scaling$time_divisor, digits=6),
                sum(scaling$variables$applied),
                if (sum(scaling$variables$applied) == 1L) '' else 's'
            ))
        }
        if (nrow(simplifications)) {
            for (i in seq_len(nrow(simplifications))) {
                change <- simplifications[i, ]
                message(sprintf(
                    'IRF %s: %s %s -> %s (%s)',
                    change$term,
                    change$axis,
                    change$requested,
                    change$effective,
                    change$reason
                ))
            }
        }
        for (plan in plans) {
            message(sprintf(
                'IRF %s: %s layout, %s responses, %s links, max history %s, %.1f%% dense padding',
                plan$term,
                plan$selected,
                format(plan$responses, big.mark=','),
                format(plan$links, big.mark=','),
                plan$maximum_history,
                100 * plan$dense_padding
            ))
        }
    }
    out
}

#' Fit a continuous-time deconvolutional GAM
#'
#' `cdrgam()` compiles untiled impulse and response streams and fits the
#' resulting model. Use [prepare_cdrgam()] followed by [cdrgam.fit()] when a
#' compiled design will be reused across fits.
#'
#' @inheritParams prepare_cdrgam
#' @param family An `mgcv` family.
#' @param method Smoothing-parameter estimation method.
#' @param engine Either `"bam"` or `"gam"` for the native backend.
#' @param backend Fitting backend. `"mgcv"` uses native fitting, `"block"`
#'   selects the dense Gaussian REML reference solver, and `"sparse"` uses
#'   sparse penalized normal equations. Select the trust-region optimizer with
#'   `sparse_control=list(gradient="exact", outer_optimizer="bfgs_trust")`.
#' @param checkpoint Optional checkpoint path for a custom backend. Checkpoints
#'   are written atomically after periodic REML evaluations and contain
#'   validated optimization/restart state. The trust optimizer checkpoints its
#'   dense curvature approximation, trust radius, counters, parameters, score,
#'   and criterion so interruption resumes the same optimization trajectory;
#'   older checkpoints without that state remain parameter-only warm starts.
#'   Reusing the path resumes an interrupted fit or skips an already completed
#'   outer optimization.
#' @param solver_trace Custom-backend progress reporting. `FALSE` or `0` is
#'   silent; `TRUE` or `1` reports phases, improving solutions, and
#'   trust-optimizer diagnostics; `2` reports every objective evaluation; and
#'   `3` adds low-level chunk/factorization events. Trust diagnostics include
#'   projected-gradient and tolerance ratios, step norms, trust radius,
#'   predicted and actual improvement, acceptance statistics, and expected
#'   post-fit Hessian work. At the numerical trust-radius floor, the sparse
#'   trust optimizer uses the resolved Hessian strategy to either certify
#'   practical convergence or perform a bounded positive-curvature reset.
#'   Automatic selection uses the same process memory budget as the post-fit
#'   Hessian calculation.
#'   A function receives the same progress events as named lists.
#' @param sparse_control Named control list for the sparse backend. The
#'   `gradient` entry may be `"auto"` (the default), `"finite"`, `"exact"`,
#'   `"stochastic"`, or `"hybrid"`; `gradient_probes` controls the fixed
#'   Rademacher trace probes used by the latter two methods; `"hybrid"` uses
#'   stochastic scores for a warm start and then refines with exact scores;
#'   `"auto"` compares predicted exact-score time and memory with the cost of
#'   finite differences after the first objective factorization;
#'   `gradient_cores` evaluates central finite-difference directions in
#'   parallel on non-Windows systems (use single-threaded BLAS when greater
#'   than one); `finite_difference_step` defaults to `1e-3`;
#'   `outer_optimizer` may be `"auto"` (the default), `"lbfgsb"`, or
#'   `"bfgs_trust"`; automatic selection uses safeguarded trust-region BFGS
#'   for exact gradients and L-BFGS-B otherwise; explicitly requesting
#'   `"bfgs_trust"` with `gradient="auto"` forces an exact gradient;
#'   `optimizer_maxit` defaults to `200`;
#'   `hessian` selects the post-fit outer-Hessian calculation: `"auto"`
#'   (the default) uses the analytic Hessian when its estimated peak memory fits
#'   a conservative fraction of the process's available memory and otherwise
#'   uses `"gradient"`; `"gradient"` differences the exact REML score;
#'   `"analytic"` differentiates the Gaussian REML criterion using the fitted
#'   factor and a selected inverse,
#'   `"profiled"` retains the slower reference calculation that differences the
#'   profiled scalar criterion, `"optimhess"` directly differences the criterion,
#'   `"defer"` retains the sparse workspace and computes the gradient-based
#'   Hessian on the first unconditional-inference request, and `"none"` skips
#'   smoothing-parameter uncertainty; the automatic memory check honors the
#'   effective Linux cgroup or SLURM allocation when available and the
#'   `cdrgam.memory_limit_bytes` option can impose a smaller process limit;
#'   `hessian_step` controls the
#'   finite-difference step and defaults to `1e-2`;
#'   `supernodal` optionally overrides the automatically selected CHOLMOD
#'   factorization form; `trace_method` may be `"auto"` (the default),
#'   `"solve"`, or `"schur_inverse"`; the latter contracts penalties with
#'   block-local inverse entries and requires `schur="always"`;
#'   `trace_chunk_size` bounds exact-score solves; and
#'   `boundary_action` controls empirically inactive penalty subspaces:
#'   `"report"` (the default) records them, `"reduce"` projects certified IRF
#'   curvature boundaries to their joint null space and reoptimizes, and
#'   `"error"` turns optimizer nonconvergence into an error.
#'   `boundary_log_sp` is the minimum log smoothing parameter considered for
#'   reduction and defaults to `12`. Reduced-model inference is conditional on
#'   the recorded data-selected boundary reduction. Full-rank boundaries that
#'   would remove an entire term are reported but require user confirmation.
#'   `schur = "always"` enables the experimental response-group Schur solver
#'   instead of the default `"never"`. It batches block elimination, exact
#'   traces, and solves through BLAS; a threaded BLAS implementation improves
#'   large-core runtime. `crossprod_chunk_size` bounds the
#'   number of response rows used while accumulating sparse normal equations;
#'   `restarts` requests additional deterministic smoothing-parameter starts.
#' @param rank_action Policy for scientifically meaningful nonidentifiability
#'   in a custom backend: `"error"` (the default), `"minimum_norm"`, `"drop"`,
#'   or `"penalize"`. Exact aliases among ordinary unpenalized parametric
#'   columns are always removed and reported. `"drop"` does not arbitrarily
#'   delete directions spanning smooth or IRF terms.
#' @param rank_tol Relative numerical rank tolerance. The default is
#'   `sqrt(.Machine$double.eps)`.
#' @param rank_penalty Relative fixed ridge used by
#'   `rank_action="penalize"`. The default is `1e-6`.
#' @param ... Additional arguments passed to `mgcv::bam()` or `mgcv::gam()`.
#' @return A fitted `cdrgam` object.
.cdrgam_compact_call <- function(call, function_name, data_arguments) {
    head <- call[[1L]]
    if (is.function(head) ||
            (is.call(head) && identical(head[[1L]], quote(`function`)))) {
        call[[1L]] <- as.name(function_name)
    }
    for (argument in intersect(data_arguments, names(call))) {
        value <- call[[argument]]
        if (!is.name(value) && !is.call(value)) {
            call[[argument]] <- as.name(argument)
        }
    }
    call
}

#' @export
cdrgam <- function(
        formula,
        impulses,
        responses,
        window=NULL,
        series=character(),
        impulse_time='time',
        response_time='time',
        history=c('auto', 'dense', 'ragged'),
        history_length=Inf,
        chunk_size=10000,
        rescale_predictors=FALSE,
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
        drop.unused.levels=TRUE,
        ...
) {
    call <- .cdrgam_compact_call(
        match.call(),
        'cdrgam',
        c('impulses', 'responses')
    )
    if (!inherits(formula, 'formula')) stop('formula must be a formula')
    design <- prepare_cdrgam(
        formula=formula,
        impulses=impulses,
        responses=responses,
        window=window,
        series=series,
        impulse_time=impulse_time,
        response_time=response_time,
        history=history,
        history_length=history_length,
        chunk_size=chunk_size,
        rescale_predictors=rescale_predictors,
        drop.unused.levels=drop.unused.levels
    )
    fit <- cdrgam.fit(
        design=design,
        family=family,
        method=method,
        engine=engine,
        backend=backend,
        checkpoint=checkpoint,
        solver_trace=solver_trace,
        sparse_control=sparse_control,
        rank_action=rank_action,
        rank_tol=rank_tol,
        rank_penalty=rank_penalty,
        ...
    )
    fit$call <- call
    fit$cdrgam$call <- call
    fit
}

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
        fit <- .fit_block_gaussian(
            design=design,
            family=family,
            method=method,
            checkpoint=checkpoint,
            trace=solver_trace,
            rank_action=rank_action,
            rank_tol=rank_tol,
            rank_penalty=rank_penalty,
            drop.unused.levels=design$configuration$drop.unused.levels,
            ...
        )
        fit$call <- call
        fit$cdrgam$call <- call
        return(fit)
    }
    if (identical(backend, 'sparse')) {
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
        boundary_plan <- fit$cdrgam$boundary_reductions
        if (identical(boundary_action, 'reduce') &&
                !is.null(boundary_plan) && nrow(boundary_plan)) {
            planned <- .cdr_boundary_reduction_plan(
                fit,
                design,
                log_sp_threshold=if (is.null(
                    sparse_control[['boundary_log_sp', exact=TRUE]]
                )) 12 else sparse_control[['boundary_log_sp', exact=TRUE]],
                score_tolerance=if (is.null(
                    sparse_control[['optimizer_gradient_tolerance', exact=TRUE]]
                )) 1e-4 else sparse_control[[
                    'optimizer_gradient_tolerance', exact=TRUE
                ]]
            )
            if (length(planned$entries)) {
                boundary_reporter <- .new_solver_reporter(
                    solver_trace,
                    'sparse'
                )
                boundary_reporter$phase(
                    'boundary reduction',
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
                        term=names(design$terms)[[entry$term_index]],
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
                reduced_design <- .cdr_apply_boundary_reduction(design, planned)
                reduced_initial <- .cdr_boundary_warm_start(
                    fit,
                    reduced_design
                )
                reduced_control <- sparse_control
                reduced_control$boundary_action <- 'report'
                boundary_reporter$emit(
                    1L,
                    'reduced model refit started',
                    reason='infinite-penalty limits require a reduced basis',
                    initialization='mapped converged pilot smoothing parameters',
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
                    initial_source=
                        'mapped converged boundary-pilot smoothing parameters',
                    ...
                )
                boundary_reporter$emit(
                    1L,
                    'reduced model refit complete',
                    criterion=format(fit$reml, digits=10),
                    converged=isTRUE(fit$converged),
                    coefficients=length(fit$coefficients),
                    smoothing_parameters=length(fit$sp)
                )
                remaining <- fit$cdrgam$boundary_reductions
                applied <- planned$table
                applied$status[applied$status ==
                    'certified_boundary_candidate'] <-
                    'applied_and_reoptimized'
                applied$fit_stage <- 'pilot'
                reductions <- if (is.null(remaining) || !nrow(remaining)) {
                    applied
                } else {
                    remaining$fit_stage <- 'reoptimized'
                    rbind(applied, remaining)
                }
                fit$cdrgam$boundary_reductions <- reductions
                fit$cdrgam$conditional_on_boundary_reduction <- TRUE
                fit$sparse$convergence$boundary_reductions <- reductions
                fit$sparse$convergence$conditional_on_boundary_reduction <- TRUE
                boundary_reduced <- TRUE
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

#' Extract a CDR-GAM formula
#'
#' @param x A fitted `cdrgam` object or prepared `cdrgam_design`.
#' @param type Formula representation: the formula supplied by the user, the
#'   normalized CDR formula including structural defaults, or the effective
#'   CDR formula after automatic simplification. Fitted models also support
#'   `"mgcv"`, the translated formula used by the fitting backend.
#' @param ... Unused.
#' @return A formula.
#' @rdname formula.cdrgam
#' @export
formula.cdrgam_design <- function(
        x,
        type=c('user', 'normalized', 'effective'),
        ...
) {
    type <- match.arg(type)
    switch(
        type,
        user=x$formula,
        normalized=x$normalized_formula,
        effective=x$effective_formula
    )
}

#' @rdname formula.cdrgam
#' @export
formula.cdrgam <- function(
        x,
        type=c('user', 'normalized', 'effective', 'mgcv'),
        ...
) {
    type <- match.arg(type)
    value <- x$cdrgam$formula[[type]]
    if (is.null(value)) x$cdrgam$formula$user else value
}
