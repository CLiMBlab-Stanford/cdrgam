#' Specify an impulse-response term
#'
#' `irf()` is used inside a [fit_cdrgam()] formula. Its first argument names a
#' numeric column in the impulse stream, or is the literal `1` for the
#' deconvolutional intercept. Ordinary formula terms retain their usual
#' `mgcv` meaning and are evaluated against the response stream.
#'
#' @param predictor Unquoted numeric impulse-stream column, or literal `1`.
#' @param window Inclusive lag window as `c(minimum, maximum)`.
#' @param k Lag-basis dimension. For nonlinear IRFs, one value or dimensions
#'   for lag and predictor value.
#' @param bs Basis name or, for nonlinear IRFs, lag and predictor basis names.
#'   Currently only `"cr"` is implemented.
#' @param nonlinear If true, model a tensor surface over predictor value and
#'   lag. In this case `k` and `bs` may each have two entries, ordered as lag
#'   then predictor value.
#' @param varying Optional unquoted response-aligned numeric covariate. This
#'   constructs a centered tensor surface over lag and the covariate.
#' @param by Optional unquoted response-aligned numeric covariate multiplying
#'   the complete IRF contribution.
#' @param group Optional unquoted response-stream factor defining random IRF
#'   deviations.
#' @return An internal `cdrgam_irf_spec` when evaluated by the formula compiler.
#' @export
irf <- function(
        predictor,
        window=c(0, Inf),
        k=10,
        bs='cr',
        nonlinear=FALSE,
        varying=NULL,
        by=NULL,
        group=NULL
) {
    predictor_expr <- substitute(predictor)
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
    if (length(window) != 2L || !is.numeric(window) ||
            !is.finite(window[[1L]]) || is.na(window[[2L]]) ||
            window[[1L]] < 0 || window[[2L]] < window[[1L]]) {
        stop('window must satisfy 0 <= min <= max; max may be Inf')
    }
    if (length(nonlinear) != 1L || !is.logical(nonlinear) || is.na(nonlinear)) {
        stop('nonlinear must be TRUE or FALSE')
    }
    varying_name <- symbol_name(varying_expr, 'varying')
    by_name <- symbol_name(by_expr, 'by')
    if (nonlinear && !is.null(varying_name)) {
        stop('nonlinear and varying cannot currently be combined in one IRF')
    }
    tensor_term <- nonlinear || !is.null(varying_name)
    required_dimensions <- if (tensor_term) 2L else 1L
    if (!is.numeric(k) || any(!is.finite(k)) || any(k < 3) ||
            !(length(k) %in% c(1L, required_dimensions))) {
        stop('k must supply dimensions of at least 3 for lag and predictor')
    }
    if (length(k) == 1L && tensor_term) {
        k <- rep.int(k, 2L)
    }
    if (!is.character(bs) || !(length(bs) %in% c(1L, required_dimensions))) {
        stop('bs must supply a basis name for lag and predictor')
    }
    if (length(bs) == 1L && tensor_term) {
        bs <- rep.int(bs, 2L)
    }
    if (any(bs != 'cr')) {
        stop('The initial formula compiler supports only bs="cr"')
    }
    constant <- is.numeric(predictor_expr) && length(predictor_expr) == 1L &&
        isTRUE(all.equal(as.numeric(predictor_expr), 1))
    if (!constant && !is.symbol(predictor_expr)) {
        stop('predictor must be an unquoted column name or the literal 1')
    }
    if (constant && nonlinear) {
        stop('irf(1) cannot use nonlinear=TRUE')
    }
    out <- list(
        predictor=if (constant) '1' else as.character(predictor_expr),
        constant=constant,
        window=as.numeric(window),
        k=as.integer(k),
        bs=bs,
        nonlinear=nonlinear,
        varying=varying_name,
        by=by_name,
        group=symbol_name(group_expr, 'group')
    )
    class(out) <- 'cdrgam_irf_spec'
    out
}

.format_irf_spec <- function(spec) {
    format_vector <- function(value) {
        paste0('c(', paste(format(value, trim=TRUE, scientific=FALSE), collapse=', '), ')')
    }
    arguments <- c(
        if (isTRUE(spec$constant)) '1' else spec$predictor,
        paste0('window=', format_vector(spec$window)),
        paste0('k=', format_vector(spec$k)),
        paste0('bs=', if (length(spec$bs) == 1L) {
            encodeString(spec$bs, quote='"')
        } else {
            paste0('c(', paste(encodeString(spec$bs, quote='"'), collapse=', '), ')')
        })
    )
    if (isTRUE(spec$nonlinear)) arguments <- c(arguments, 'nonlinear=TRUE')
    if (!is.null(spec$varying)) arguments <- c(arguments, paste0('varying=', spec$varying))
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

.parse_cdr_formula <- function(formula) {
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
    is_fixed_rate <- function(spec) {
        isTRUE(spec$constant) && is.null(spec$varying) &&
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
            lag_k <- vapply(specs, function(spec) spec$k[[1L]], integer(1))
            rate_spec <- irf(
                1,
                window=c(
                    min(vapply(specs, function(spec) spec$window[[1L]], numeric(1))),
                    max(vapply(specs, function(spec) spec$window[[2L]], numeric(1)))
                ),
                k=max(lag_k),
                bs=specs[[1L]]$bs[[1L]]
            )
        } else {
            rate_spec <- irf(1)
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
        knots=list(cdr_delta=knots),
        absorb.cons=FALSE,
        scale.penalty=FALSE,
        n=length(knots)
    )[[1L]]
    transform <- if (centered) {
        constraint <- .cdr_basis_mean_constraint(
            marginal,
            'cdr_delta',
            links$delay,
            chunk_size
        )
        .cdr_constraint_transform(constraint, k)
    } else {
        diag(k)
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

.group_cdr_term <- function(term, group, group_name) {
    if (anyNA(group)) {
        stop('Random-IRF grouping factors must not contain missing values')
    }
    group <- droplevels(as.factor(group))
    levels <- levels(group)
    if (length(levels) < 2L) {
        stop('A random IRF requires at least two observed grouping levels')
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

#' Compile untiled impulse and response streams for CDR-GAM fitting
#'
#' @param formula Extended model formula containing ordinary `mgcv` terms and
#'   optional [irf()] terms. An implicit `irf(1)` is added unless suppressed.
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
#' @param quiet Suppress the preparation plan message.
#' @return A reusable `cdrgam_design` object.
#' @export
prepare_cdrgam <- function(
        formula,
        impulses,
        responses,
        series=character(),
        impulse_time='time',
        response_time='time',
        history=c('auto', 'dense', 'ragged'),
        history_length=Inf,
        chunk_size=10000,
        quiet=FALSE
) {
    if (!is.data.frame(impulses) || !is.data.frame(responses)) {
        stop('impulses and responses must be data frames')
    }
    history <- match.arg(history)
    if (length(history_length) != 1L || !is.numeric(history_length) ||
            is.na(history_length) || history_length <= 0 ||
            (is.finite(history_length) && history_length != as.integer(history_length))) {
        stop('history_length must be a positive integer or Inf')
    }
    parsed <- .parse_cdr_formula(formula)
    response_name <- all.vars(parsed$ordinary_formula[[2L]])
    if (length(response_name) != 1L || !(response_name %in% names(responses))) {
        stop('The formula response must name one column in responses')
    }
    terms <- vector('list', length(parsed$irfs))
    plans <- vector('list', length(parsed$irfs))
    for (i in seq_along(parsed$irfs)) {
        spec <- parsed$irfs[[i]]
        if (!isTRUE(spec$constant) && !(spec$predictor %in% names(impulses))) {
            stop('IRF predictor not found in impulses: ', spec$predictor)
        }
        weights <- if (isTRUE(spec$constant)) {
            rep.int(1, nrow(impulses))
        } else {
            impulses[[spec$predictor]]
        }
        if (!is.numeric(weights) || any(!is.finite(weights))) {
            stop('IRF predictors must be finite numeric columns')
        }
        by_values <- NULL
        if (!is.null(spec$by)) {
            if (!(spec$by %in% names(responses))) {
                stop('By covariate not found in responses: ', spec$by)
            }
            by_values <- responses[[spec$by]]
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
        plan <- .choose_history_layout(links$counts, history)
        term_name <- if (isTRUE(spec$constant)) 'irf(1)' else spec$predictor
        plans[[i]] <- c(list(term=term_name, window=spec$window), plan)
        link_weights <- weights[links$impulse_index]
        if (!is.null(spec$varying)) {
            if (!(spec$varying %in% names(responses))) {
                stop('Varying covariate not found in responses: ', spec$varying)
            }
            varying_values <- responses[[spec$varying]]
            if (!is.numeric(varying_values) || any(!is.finite(varying_values))) {
                stop('Varying covariates must be finite numeric columns')
            }
            terms[[i]] <- .compress_cdr_tensor_links(
                links,
                varying_values[links$response_index],
                n=nrow(responses),
                k=spec$k,
                bs=spec$bs,
                chunk_size=chunk_size,
                name=term_name,
                link_weights=link_weights,
                type='varying',
                response_multiplier=by_values
            )
            terms[[i]]$varying <- spec$varying
            terms[[i]]$name <- paste0(term_name, '~', spec$varying)
        } else if (isTRUE(spec$nonlinear)) {
            terms[[i]] <- .compress_cdr_tensor_links(
                links,
                link_weights,
                n=nrow(responses),
                k=spec$k,
                bs=spec$bs,
                chunk_size=chunk_size,
                name=term_name,
                response_multiplier=by_values
            )
        } else if (identical(plan$selected, 'dense')) {
            width <- plan$maximum_history
            if (!width) {
                stop('No impulses fall within the requested IRF window')
            }
            knots <- as.numeric(stats::quantile(
                unique(links$delay),
                seq(0, 1, length.out=spec$k),
                names=FALSE
            ))
            delay_matrix <- matrix(knots[[1L]], nrow=nrow(responses), ncol=width)
            weight_matrix <- matrix(0, nrow=nrow(responses), ncol=width)
            positions <- sequence(links$counts)
            index <- cbind(links$response_index, positions)
            delay_matrix[index] <- links$delay
            weight_matrix[index] <- link_weights
            terms[[i]] <- compress_cdr_smooth(
                delay_matrix,
                weight_matrix,
                k=spec$k,
                bs=spec$bs,
                knots=knots,
                chunk_size=chunk_size,
                name=term_name,
                response_multiplier=by_values,
                constraint_delays=links$delay
            )
        } else {
            terms[[i]] <- .compress_cdr_links(
                links,
                link_weights,
                n=nrow(responses),
                k=spec$k,
                bs=spec$bs,
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
                stop('Random-IRF group not found in responses: ', spec$group)
            }
            terms[[i]] <- .group_cdr_term(
                terms[[i]],
                responses[[spec$group]],
                spec$group
            )
        }
    }
    identifiability <- list(
        constraints=lapply(terms, function(term) term$constraints),
        rate=list(requested=if (any(vapply(
            parsed$irfs,
            function(spec) {
                isTRUE(spec$constant) && !isTRUE(spec$implicit) &&
                    is.null(spec$varying) && is.null(spec$by) &&
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
            isTRUE(spec$constant) && is.null(spec$varying) &&
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
                identical(dim(candidate$X), dim(rate_term$X)) &&
                    max(abs(candidate$X - rate_term$X)) <=
                        max(1, max(abs(rate_term$X))) *
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
        }
    }
    if (!length(terms)) {
        stop('No identifiable IRF terms remain in the model')
    }
    names(terms) <- vapply(terms, `[[`, character(1), 'name')
    out <- list(
        formula=formula,
        normalized_formula=parsed$normalized_formula,
        effective_formula=.compose_cdr_formula(parsed$ordinary_formula, parsed$irfs),
        ordinary_formula=parsed$ordinary_formula,
        response_name=response_name,
        responses=responses,
        terms=terms,
        specification=parsed$irfs,
        identifiability=identifiability,
        plan=plans,
        stream=list(
            series=series,
            impulse_time=impulse_time,
            response_time=response_time,
            history_length=history_length
        )
    )
    class(out) <- 'cdrgam_design'
    if (!quiet) {
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
#' The primary interface accepts untiled impulse and response streams and an
#' extended formula containing [irf()] terms. A prepared [prepare_cdrgam()]
#' result can instead be supplied as the first argument for reuse across fits.
#'
#' @inheritParams prepare_cdrgam
#' @param family An `mgcv` family.
#' @param method Smoothing-parameter estimation method.
#' @param engine Either `"bam"` or `"gam"` for the native backend.
#' @param backend Fitting backend. `"mgcv"` uses native fitting, `"block"`
#'   selects the dense Gaussian REML reference solver, and `"sparse"` uses
#'   sparse penalized normal equations. `"sparse_trust"` selects the sparse
#'   engine with its exact REML score and safeguarded trust-region BFGS outer
#'   optimizer.
#' @param checkpoint Optional checkpoint path for a custom backend. Checkpoints
#'   are written atomically after periodic REML evaluations and
#'   contain validated optimization/restart state. Reusing the path resumes an
#'   interrupted fit or skips an already completed outer optimization.
#' @param solver_trace Custom-backend progress reporting. `FALSE` or `0` is
#'   silent; `TRUE` or `1` reports phases, improving solutions, and
#'   trust-optimizer diagnostics; `2` reports every objective evaluation; and
#'   `3` adds low-level chunk/factorization events. Trust diagnostics include
#'   projected-gradient and tolerance ratios, step norms, trust radius,
#'   predicted and actual improvement, acceptance statistics, and expected
#'   post-fit Hessian work. A function receives the same progress events as
#'   named lists.
#' @param sparse_control Named control list for the sparse backend. The
#'   `gradient` entry may be `"auto"` (the default), `"finite"`, `"exact"`,
#'   `"stochastic"`, or `"hybrid"`; `gradient_probes` controls the fixed
#'   Rademacher trace probes used by the latter two methods; `"hybrid"` always
#'   finishes with the exact finite-difference REML objective;
#'   `gradient_cores` evaluates central finite-difference directions in
#'   parallel on non-Windows systems (use single-threaded BLAS when greater
#'   than one); `finite_difference_step` defaults to `1e-3`;
#'   `hessian` selects the post-fit outer-Hessian calculation: `"profiled"`
#'   (the default) reconstructs the unprofiled Gaussian REML Hessian from the
#'   lower-dimensional profiled criterion, while `"optimhess"` directly
#'   differences the unprofiled criterion; `hessian_step` controls the
#'   finite-difference step for either method and defaults to `1e-2`;
#'   `supernodal` optionally overrides the automatically selected CHOLMOD
#'   factorization form; `trace_chunk_size` bounds exact-score solves; and
#'   `schur = "always"` enables the experimental response-group Schur solver
#'   instead of the default `"never"`. `crossprod_chunk_size` bounds the
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
#' @export
fit_cdrgam <- function(
        formula,
        impulses=NULL,
        responses=NULL,
        series=character(),
        impulse_time='time',
        response_time='time',
        history=c('auto', 'dense', 'ragged'),
        history_length=Inf,
        chunk_size=10000,
        family=stats::gaussian(),
        method=NULL,
        engine=c('bam', 'gam'),
        backend=c('mgcv', 'block', 'sparse', 'sparse_trust'),
        checkpoint=NULL,
        solver_trace=FALSE,
        sparse_control=list(),
        rank_action=c('error', 'minimum_norm', 'drop', 'penalize'),
        rank_tol=NULL,
        rank_penalty=NULL,
        ...
) {
    # Compatibility with the early low-level prototype: fit_cdrgam(y, terms).
    if (!inherits(formula, c('formula', 'cdrgam_design'))) {
        return(.fit_compressed_mgcv(
            y=formula,
            terms=impulses,
            family=family,
            method=method,
            engine=engine,
            ...
        ))
    }
    backend <- match.arg(backend)
    rank_action <- match.arg(rank_action)
    trace_control <- .solver_trace_level(solver_trace)
    .rank_tolerance(rank_tol)
    .rank_penalty(rank_penalty)
    design <- if (inherits(formula, 'cdrgam_design')) {
        formula
    } else {
        prepare_cdrgam(
            formula=formula,
            impulses=impulses,
            responses=responses,
            series=series,
            impulse_time=impulse_time,
            response_time=response_time,
            history=history,
            history_length=history_length,
            chunk_size=chunk_size
        )
    }
    if (identical(backend, 'block')) {
        if (length(sparse_control)) {
            stop('sparse_control applies only to backend="sparse"')
        }
        design <- .materialize_cdr_design(design, sparse=FALSE)
        return(.fit_block_gaussian(
            design=design,
            family=family,
            method=method,
            checkpoint=checkpoint,
            trace=solver_trace,
            rank_action=rank_action,
            rank_tol=rank_tol,
            rank_penalty=rank_penalty,
            ...
        ))
    }
    if (backend %in% c('sparse', 'sparse_trust')) {
        if (identical(backend, 'sparse_trust')) {
            if (!is.null(sparse_control$gradient) &&
                    !identical(sparse_control$gradient, 'exact')) {
                stop(
                    'backend="sparse_trust" requires ',
                    'sparse_control$gradient="exact"'
                )
            }
            if (!is.null(sparse_control$outer_optimizer) &&
                    !identical(
                        sparse_control$outer_optimizer,
                        'bfgs_trust'
                    )) {
                stop(
                    'backend="sparse_trust" requires ',
                    'sparse_control$outer_optimizer="bfgs_trust"'
                )
            }
            sparse_control$gradient <- 'exact'
            sparse_control$outer_optimizer <- 'bfgs_trust'
            if (is.null(sparse_control$optimizer_gradient_tolerance)) {
                sparse_control$optimizer_gradient_tolerance <- 2e-4
            }
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
            ...
        )
        if (identical(backend, 'sparse_trust')) {
            fit$cdrgam$backend <- 'sparse_trust'
            fit$cdrgam$solver <- paste(
                'sparse Gaussian REML solver',
                '(safeguarded trust-region BFGS)'
            )
        }
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
    .fit_compressed_mgcv(
        y=y,
        terms=design$terms,
        family=family,
        method=method,
        engine=engine,
        base_formula=design$ordinary_formula,
        response_data=design$responses,
        user_formula=design$formula,
        preparation=list(
            plan=design$plan,
            stream=design$stream,
            specification=design$specification,
            identifiability=design$identifiability,
            normalized_formula=design$normalized_formula,
            effective_formula=design$effective_formula
        ),
        ...
    )
}

#' Inspect the user-facing or expanded formula of a CDR-GAM
#'
#' @param object A fitted `cdrgam` object or prepared `cdrgam_design`.
#' @param type Formula representation: the formula supplied by the user, the
#'   normalized CDR formula including structural defaults, or the effective
#'   CDR formula after automatic simplification.
#' @return A formula.
#' @export
cdr_formula <- function(
        object,
        type=c('user', 'normalized', 'effective')
) {
    type <- match.arg(type)
    if (inherits(object, 'cdrgam_design')) {
        return(switch(
            type,
            user=object$formula,
            normalized=object$normalized_formula,
            effective=object$effective_formula
        ))
    }
    if (!is_cdrgam(object)) {
        stop('object must be a cdrgam fit or cdrgam_design')
    }
    formula <- object$cdrgam$formula[[type]]
    if (is.null(formula)) object$cdrgam$formula$user else formula
}

#' @rdname cdr_formula
#' @export
mgcv_formula <- function(object) {
    if (!is_cdrgam(object)) {
        stop('object must be a fitted cdrgam model')
    }
    object$cdrgam$formula$mgcv
}
