#' @importFrom mgcv Predict.matrix smooth.construct
NULL

.cdr_constraint_transform <- function(constraint, dimension, tolerance=1e-10) {
    if (is.null(constraint) || !length(constraint)) return(diag(dimension))
    constraint <- matrix(constraint, ncol=dimension)
    decomposition <- qr(t(constraint), tol=tolerance, LAPACK=FALSE)
    rank <- decomposition$rank
    if (!rank) return(diag(dimension))
    complete <- qr.Q(decomposition, complete=TRUE)
    if (rank >= dimension) {
        return(matrix(0, nrow=dimension, ncol=0L))
    }
    complete[, seq.int(rank + 1L, dimension), drop=FALSE]
}

.cdr_constant_sum <- function(values) {
    scale <- max(1, max(abs(values)))
    (max(values) - min(values)) <= scale * .Machine$double.eps * 1000
}

.cdr_basis_mean_constraint <- function(basis, variable, values, chunk_size) {
    # mgcv constructs its side constraint on the de-duplicated covariate
    # values retained by the smooth constructor, before the summation
    # convention expands them back to history positions.
    values <- unique(values)
    total <- numeric(ncol(basis$X))
    for (start in seq.int(1L, length(values), by=chunk_size)) {
        end <- min(length(values), start + chunk_size - 1L)
        block <- values[start:end]
        total <- total + colSums(mgcv::PredictMat(
            basis,
            stats::setNames(list(block), variable),
            n=length(block)
        ))
    }
    matrix(total / length(values), nrow=1L)
}

#' Build a compressed stationary CDR impulse-response term
#'
#' Construct the response-level design matrix for a stationary, one-dimensional
#' CDR impulse response without retaining the much larger basis matrix over all
#' response--impulse pairs.  For response `i`, the returned design row is
#'
#' `sum_j by[i, j] * b(t_delta[i, j])`,
#'
#' where `b()` is an `mgcv` smooth basis. The resulting term
#' can be fitted with [cdrgam()].
#'
#' The implementation supports univariate numeric marginal bases recognized
#' by `mgcv` and deliberately does not apply
#' a centering constraint unless the effective row sum is constant. This
#' matches `mgcv` linear-functional terms
#' when `rowSums(by)` is non-constant, which is the usual case for
#' predictor-weighted CDR terms.
#'
#' @param t_delta Numeric matrix of response-minus-impulse time differences.
#' @param by Numeric matrix of convolution weights with the same dimensions as
#'   `t_delta`. Invalid/padded history positions should have weight zero.
#' @param k Basis dimension.
#' @param bs Univariate numeric marginal basis recognized by `mgcv`.
#' @param knots Optional numeric knot vector of length `k`. If omitted, knots
#'   are placed at quantiles of the unique values in `t_delta`, following the
#'   `mgcv` cubic regression spline constructor.
#' @param chunk_size Number of responses whose history is transformed at once.
#' @param name Optional human-readable term name stored with the model.
#' @param response_multiplier Optional response-level multiplier applied after
#'   convolution. It is used when deciding whether mgcv centering is required.
#' @param constraint_delays Optional delay values defining the centering
#'   measure. The stream compiler uses actual, unpadded history links.
#' @return A `cdrgam_term` object containing the
#'   response-level design matrix, penalties, and underlying `mgcv` smooth.
#' @keywords internal
compress_cdr_smooth <- function(
        t_delta,
        by,
        k=10,
        bs='cr',
        knots=NULL,
        chunk_size=10000,
        name=NULL,
        response_multiplier=NULL,
        constraint_delays=NULL
) {
    if (!is.matrix(t_delta) || !is.numeric(t_delta)) {
        stop('t_delta must be a numeric matrix')
    }
    if (!is.matrix(by) || !is.numeric(by)) {
        stop('by must be a numeric matrix')
    }
    if (!identical(dim(t_delta), dim(by))) {
        stop('t_delta and by must have identical dimensions')
    }
    if (nrow(t_delta) < 1 || ncol(t_delta) < 1) {
        stop('t_delta and by must have at least one row and one column')
    }
    if (any(!is.finite(t_delta))) {
        stop('t_delta must contain only finite values')
    }
    if (any(!is.finite(by))) {
        stop('by must contain only finite values')
    }
    if (!is.character(bs) || length(bs) != 1L || is.na(bs) || !nzchar(bs)) {
        stop('bs must be one nonempty mgcv basis name')
    }
    if (length(k) != 1 || !is.numeric(k) || !is.finite(k) || k < 3) {
        stop('k must be a finite numeric scalar greater than or equal to 3')
    }
    k <- as.integer(k)
    if (length(chunk_size) != 1 || !is.numeric(chunk_size) ||
            !is.finite(chunk_size) || chunk_size < 1) {
        stop('chunk_size must be a positive finite numeric scalar')
    }
    chunk_size <- as.integer(chunk_size)
    if (!is.null(name) &&
            (length(name) != 1L || !is.character(name) || is.na(name) ||
             !nzchar(name))) {
        stop('name must be NULL or one non-empty character string')
    }

    if (!is.null(response_multiplier) &&
            (length(response_multiplier) != nrow(by) ||
             !is.numeric(response_multiplier) ||
             any(!is.finite(response_multiplier)))) {
        stop('response_multiplier must be finite with one value per response')
    }
    by_sum <- rowSums(by)
    effective_sum <- if (is.null(response_multiplier)) {
        by_sum
    } else {
        by_sum * response_multiplier
    }
    centered <- .cdr_constant_sum(effective_sum)

    if (is.null(knots)) {
        values <- unique(as.numeric(t_delta))
        if (length(values) < k) {
            stop('t_delta has fewer unique values than k')
        }
        knots <- as.numeric(stats::quantile(
            values,
            seq(0, 1, length.out=k),
            names=FALSE
        ))
    } else {
        knots <- as.numeric(knots)
        if (length(knots) != k || any(!is.finite(knots))) {
            stop('knots must contain exactly k finite numeric values')
        }
    }
    if (any(diff(knots) <= 0)) {
        stop('knots must be strictly increasing')
    }

    # Construct the mgcv marginal using only its knots. Predictions from this
    # object can then be evaluated in bounded response chunks. Penalty scaling
    # is reproduced below using the full covariate range, as smoothCon() does.
    cdr_delta <- NULL
    spec <- mgcv::s(cdr_delta, k=k, bs=bs)
    marginal <- mgcv::smoothCon(
        spec,
        data=list(cdr_delta=knots),
        knots=if (identical(bs, 'ps')) NULL else list(cdr_delta=knots),
        absorb.cons=FALSE,
        scale.penalty=FALSE,
        n=length(knots)
    )[[1]]
    marginal_dimension <- ncol(marginal$X)

    transform <- if (centered) {
        centering_values <- if (is.null(constraint_delays)) {
            as.numeric(t_delta)
        } else {
            as.numeric(constraint_delays)
        }
        if (!length(centering_values) || any(!is.finite(centering_values))) {
            stop('constraint_delays must contain finite delay values')
        }
        constraint <- .cdr_basis_mean_constraint(
            marginal,
            'cdr_delta',
            centering_values,
            chunk_size
        )
        .cdr_constraint_transform(constraint, marginal_dimension)
    } else {
        diag(marginal_dimension)
    }
    basis_dimension <- ncol(transform)

    n <- nrow(t_delta)
    width <- ncol(t_delta)
    design <- matrix(0, nrow=n, ncol=basis_dimension)
    max_basis_row_norm <- 0
    starts <- seq.int(1L, n, by=chunk_size)

    for (start in starts) {
        end <- min(n, start + chunk_size - 1L)
        rows <- start:end
        n_block <- length(rows)
        delta_block <- as.numeric(t_delta[rows, , drop=FALSE])
        by_block <- as.numeric(by[rows, , drop=FALSE])
        raw_basis_block <- mgcv::PredictMat(
            marginal,
            list(cdr_delta=delta_block),
            n=length(delta_block)
        )
        max_basis_row_norm <- max(
            max_basis_row_norm,
            max(rowSums(abs(raw_basis_block)))
        )
        basis_block <- raw_basis_block %*% transform
        weighted_basis <- basis_block * by_block
        for (column in seq_len(basis_dimension)) {
            design[rows, column] <- rowSums(matrix(
                weighted_basis[, column],
                nrow=n_block,
                ncol=width
            ))
        }
    }

    # mgcv scales each penalty by ||X||_inf^2 / ||S||_inf before fitting.
    # Here X refers to the unweighted history-level marginal basis, not the
    # response-level matrix after applying `by` and summing.
    basis_scale <- max_basis_row_norm^2
    penalties <- marginal$S
    penalty_scales <- numeric(length(penalties))
    if (length(penalties)) {
        for (i in seq_along(penalties)) {
            penalty_norm <- norm(penalties[[i]], type='I')
            if (!is.finite(penalty_norm) || penalty_norm <= 0) {
                stop('Encountered a degenerate mgcv penalty matrix')
            }
            penalty_scales[[i]] <- penalty_norm / basis_scale
            penalties[[i]] <- penalties[[i]] * basis_scale / penalty_norm
            penalties[[i]] <- crossprod(
                transform,
                penalties[[i]] %*% transform
            )
        }
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

.restore_cdr_penalty_scales <- function(object, terms) {
    for (i in seq_along(terms)) {
        smooth_index <- which(vapply(
            object$smooth,
            function(smooth) paste(smooth$term, collapse=',') ==
                paste0('cdr_term_', i),
            logical(1)
        ))
        if (length(smooth_index) == 1L) {
            object$smooth[[smooth_index]]$S.scale <- terms[[i]]$S.scale
        }
    }
    object
}

# Fit precomputed response-level CDR bases with the native mgcv backend.
.fit_compressed_mgcv <- function(
        y,
        terms,
        family=stats::gaussian(),
        method=NULL,
        engine=c('bam', 'gam'),
        intercept=TRUE,
        base_formula=NULL,
        response_data=NULL,
        user_formula=NULL,
        preparation=NULL,
        setup_only=FALSE,
        ...
) {
    engine <- match.arg(engine)
    if (inherits(terms, 'cdrgam_term')) {
        terms <- list(terms)
    }
    if (!is.list(terms) || !length(terms) ||
            any(!vapply(terms, inherits, logical(1), 'cdrgam_term'))) {
        stop('terms must be a cdrgam_term or a non-empty list of them')
    }
    n <- length(y)
    if (!n || any(!is.finite(y))) {
        stop('y must be a non-empty finite numeric vector')
    }
    if (any(vapply(terms, function(term) nrow(term$X) != n, logical(1)))) {
        stop('Every compressed term must have one row per response')
    }

    formula_env <- new.env(parent=parent.frame())
    # Ordinary mgcv terms in the user formula must work after `library(cdrgam)`
    # even when mgcv itself is not attached to the search path.
    formula_env$s <- mgcv::s
    formula_env$te <- mgcv::te
    formula_env$ti <- mgcv::ti
    formula_env$t2 <- mgcv::t2
    if (is.null(base_formula)) {
        data <- list(cdr_response=as.numeric(y))
    } else {
        if (!inherits(base_formula, 'formula')) {
            stop('base_formula must be a formula')
        }
        if (!is.data.frame(response_data)) {
            stop('response_data must be a data frame when base_formula is used')
        }
        data <- as.list(response_data)
    }
    formula_terms <- character(length(terms))
    for (i in seq_along(terms)) {
        data_name <- paste0('cdr_term_', i)
        spec_name <- paste0('cdr_spec_', i)
        # A unique row key prevents mgcv's covariate de-duplication from
        # collapsing the many structural zeros in grouped CDR designs. The
        # constructor receives the actual precomputed matrix through `xt`.
        data[[data_name]] <- seq_len(n)
        formula_env[[spec_name]] <- list(
            S=terms[[i]]$S,
            S.scale=terms[[i]]$S.scale,
            rank=terms[[i]]$rank,
            null.space.dim=terms[[i]]$null.space.dim,
            p=ncol(terms[[i]]$X),
            X=terms[[i]]$X
        )
        formula_terms[[i]] <- paste0(
            's(', data_name, ', bs="cdr", k=', ncol(terms[[i]]$X),
            ', xt=', spec_name, ')'
        )
    }
    cdr_rhs <- paste(formula_terms, collapse=' + ')
    if (is.null(base_formula)) {
        rhs <- cdr_rhs
        if (!isTRUE(intercept)) {
            rhs <- paste('0 +', rhs)
        }
        formula <- stats::as.formula(
            paste('cdr_response ~', rhs),
            env=formula_env
        )
    } else {
        response_text <- paste(deparse(base_formula[[2L]]), collapse='')
        base_rhs <- paste(deparse(base_formula[[3L]]), collapse='')
        rhs <- paste(base_rhs, cdr_rhs, sep=' + ')
        formula <- stats::as.formula(
            paste(response_text, '~', rhs),
            env=formula_env
        )
    }
    args <- list(formula=formula, family=family, data=data)
    if (!is.null(method)) {
        args$method <- method
    }
    args <- c(args, list(...))
    if (isTRUE(setup_only)) {
        args$fit <- FALSE
        setup <- do.call(mgcv::gam, args)
        return(.restore_cdr_penalty_scales(setup, terms))
    }
    setup_args <- args
    setup_args$fit <- FALSE
    if (identical(engine, 'bam')) {
        setup_args[c(
            'chunk.size', 'rho', 'AR.start', 'cluster', 'nthreads',
            'gc.level', 'use.chol', 'samfrac', 'coef'
        )] <- NULL
    }
    setup <- do.call(mgcv::gam, setup_args)
    setup <- .restore_cdr_penalty_scales(setup, terms)
    audit <- .audit_mgcv_setup(setup, .rank_tolerance(NULL))
    setup <- audit$setup
    if (!is.null(preparation)) {
        preparation$identifiability$global <- audit$info
    }
    fit <- if (identical(engine, 'bam')) {
        # bam's setup object retains only one fitting chunk, so it cannot be
        # used for a global audit. The full gam setup above performs the
        # backend-independent audit; bam then fits the same formula and lets
        # mgcv apply the already-verified benign parametric alias convention.
        do.call(mgcv::bam, args)
    } else {
        # Keep the audited setup separate from the object passed to mgcv.
        # Editing a gam setup's design columns is insufficient to rewrite all
        # formula/terms bookkeeping used later by summary.gam and prediction.
        # A direct fit lets mgcv apply its native alias handling coherently.
        do.call(mgcv::gam, args)
    }
    fit <- .restore_cdr_penalty_scales(fit, terms)
    labels <- names(terms)
    if (is.null(labels)) {
        labels <- rep('', length(terms))
    }
    for (i in seq_along(labels)) {
        if (!nzchar(labels[[i]]) && !is.null(terms[[i]]$name)) {
            labels[[i]] <- terms[[i]]$name
        }
        if (!nzchar(labels[[i]])) {
            labels[[i]] <- paste0('irf_', i)
        }
    }
    fit$cdrgam <- list(
        schema_version=1L,
        engine=engine,
        backend='mgcv',
        formula=list(
            user=if (is.null(user_formula)) formula else user_formula,
            normalized=if (is.null(preparation$normalized_formula)) {
                if (is.null(user_formula)) formula else user_formula
            } else preparation$normalized_formula,
            effective=if (is.null(preparation$effective_formula)) {
                if (is.null(user_formula)) formula else user_formula
            } else preparation$effective_formula,
            mgcv=formula
        ),
        preparation=preparation,
        scaling=if (is.null(preparation)) NULL else preparation$scaling,
        identifiability=if (is.null(preparation)) NULL else
            preparation$identifiability,
        term_labels=labels,
        terms=lapply(seq_along(terms), function(i) {
            term <- terms[[i]]
            smooth_index <- which(vapply(
                fit$smooth,
                function(smooth) paste(smooth$term, collapse=',') ==
                    paste0('cdr_term_', i),
                logical(1)
            ))
            coefficient_index <- if (length(smooth_index) == 1L) {
                seq.int(
                    fit$smooth[[smooth_index]]$first.para,
                    fit$smooth[[smooth_index]]$last.para
                )
            } else {
                integer()
            }
            list(
                name=term$name,
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
        })
    )
    class(fit) <- c('cdrgam', class(fit))
    fit
}

#' Test whether an object is a CDR-GAM fit
#'
#' @param x An R object.
#' @return A logical scalar.
#' @export
is_cdrgam <- function(x) inherits(x, 'cdrgam')

#' Drop the CDR-GAM subclass from a fitted model
#'
#' This is rarely necessary: standard `mgcv` methods already dispatch through
#' the inherited `gam` or `bam` class. It can be useful for code that checks
#' exact class names.
#'
#' @param x A `cdrgam` object.
#' @return The same fitted object with only its native `mgcv` classes.
#' @export
as_gam <- function(x) {
    if (!is_cdrgam(x)) {
        stop('x must be a cdrgam object')
    }
    if (!inherits(x, 'gam')) {
        stop('Only backend="mgcv" fits can be converted to a native gam object')
    }
    class(x) <- class(x)[class(x) != 'cdrgam']
    x
}

#' @export
print.cdrgam <- function(x, ...) {
    cat('Continuous-time deconvolutional GAM\n')
    cat('  engine:', x$cdrgam$engine, '\n')
    cat('  IRF terms:', paste(x$cdrgam$term_labels, collapse=', '), '\n\n')
    print(as_gam(x), ...)
    invisible(x)
}

#' @export
summary.cdrgam <- function(
        object,
        dispersion=NULL,
        freq=FALSE,
        re.test=TRUE,
        all.coefficients=FALSE,
        ...
) {
    output <- summary(
        as_gam(object), dispersion=dispersion, freq=freq,
        re.test=re.test, ...
    )
    divisors <- .cdr_coefficient_divisors(object)
    table <- output$p.table
    if (!is.null(table) && nrow(table)) {
        positions <- match(rownames(table), names(divisors))
        scale <- divisors[positions]
        scale[is.na(scale)] <- 1
        estimate <- match('Estimate', colnames(table))
        standard_error <- match('Std. Error', colnames(table))
        if (!is.na(estimate)) table[, estimate] <- table[, estimate] / scale
        if (!is.na(standard_error)) {
            table[, standard_error] <- table[, standard_error] / abs(scale)
        }
        output$p.table <- table
        output$p.coeff <- table[, 'Estimate']
    }
    formulas <- .cdrgam_summary_formulas(object)
    output$formula <- formulas$user
    output$formulas <- formulas
    output$formula_strings <- .cdrgam_formula_strings(formulas)
    output$cdrgam.scaling <- object$cdrgam$scaling
    if (!is.null(output$s.table)) {
        labels <- .cdrgam_smooth_labels(object)
        rownames(output$s.table) <- labels
        if (length(output$edf) == length(labels)) names(output$edf) <- labels
        if (length(output$s.pv) == length(labels)) names(output$s.pv) <- labels
        if (length(output$chi.sq) == length(labels)) {
            names(output$chi.sq) <- labels
        }
    }
    if (isTRUE(all.coefficients)) {
        indices <- seq_along(object$coefficients)
        output$coefficients <- .cdrgam_coefficient_table(
            object,
            indices,
            stats::vcov(object)[indices, indices, drop=FALSE],
            output$residual.df
        )
    }
    class(output) <- c('summary.cdrgam', class(output))
    output
}

#' Predict from a fitted CDR-GAM
#'
#' Supply new untiled streams as
#' `newdata=list(impulses=impulses, responses=responses)`. Stored training-data
#' scaling is applied before the response-level design is rebuilt. Native
#' `mgcv` fits also accept a response-side data frame. Factor labels are mapped
#' to the fitted training vocabulary. Unseen random-effect and grouped-IRF
#' levels produce a warning and contribute zero deviation.
#'
#' @param object A fitted `cdrgam` model.
#' @param newdata A named list containing `impulses` and `responses` data
#'   frames, or a response-side data frame for a native `mgcv` fit.
#' @param ... Prediction controls, including `type`, `se.fit`, `chunk_size`,
#'   and `unconditional`.
#' @return A prediction vector, linear-predictor matrix, or list containing
#'   `fit` and `se.fit`.
#' @export
predict.cdrgam <- function(object, newdata=NULL, ...) {
    if (is.list(newdata) && all(c('impulses', 'responses') %in% names(newdata))) {
        return(.predict_cdrgam_streams(
            object,
            impulses=newdata$impulses,
            responses=newdata$responses,
            ...
        ))
    }
    if (is.data.frame(newdata)) {
        scaling <- object$cdrgam$scaling
        if (is.null(scaling)) scaling <- object$cdrgam$preparation$scaling
        newdata <- .cdr_apply_scaling(newdata, scaling, 'responses')
    }
    NextMethod('predict', object=object, newdata=newdata)
}

# Internal mgcv extension: the covariate supplied to this smooth is an N by p
# response-level basis. mgcv flattens matrix covariates before dispatching the
# constructor, so reconstruct it using p stored in xt.
smooth.construct.cdr.smooth.spec <- function(object, data, knots) {
    p <- object$xt$p
    if (!is.null(object$xt$X)) {
        object$X <- object$xt$X
        object$cdr.X <- object$xt$X
        object$xt$X <- NULL
    } else {
        x <- data[[object$term]]
        if (length(x) %% p != 0) {
            stop(
                'Invalid precomputed CDR basis dimensions: ',
                length(x), ' values cannot be reshaped to ', p, ' columns'
            )
        }
        object$X <- matrix(x, ncol=p)
    }
    object$bs.dim <- p
    object$df <- p
    object$S <- object$xt$S
    object$S.scale <- object$xt$S.scale
    object$rank <- object$xt$rank
    object$null.space.dim <- object$xt$null.space.dim
    object$C <- matrix(0, 0, p)
    object$side.constrain <- FALSE
    object$plot.me <- FALSE
    object$no.rescale <- TRUE
    object$te.ok <- 0
    class(object) <- 'cdr.smooth'
    object
}

# Internal prediction method for the precomputed smooth.
Predict.matrix.cdr.smooth <- function(object, data) {
    x <- data[[object$term]]
    if (is.null(dim(x)) && !is.null(object$cdr.X)) {
        row_index <- as.integer(x)
        if (length(row_index) == length(x) &&
                all(is.finite(x)) && all(x == row_index) &&
                all(row_index >= 1L & row_index <= nrow(object$cdr.X))) {
            return(object$cdr.X[row_index, , drop=FALSE])
        }
    }
    if (length(x) %% object$bs.dim != 0) {
        stop('Invalid precomputed CDR basis dimensions in newdata')
    }
    matrix(x, ncol=object$bs.dim)
}
