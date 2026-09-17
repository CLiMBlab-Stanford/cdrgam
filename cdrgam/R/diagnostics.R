#' Evaluate fitted fixed-effect impulse-response functions
#'
#' @param object A fitted `cdrgam` model.
#' @param term Term name or one-based term number. Omit to evaluate every IRF.
#' @param lag Optional numeric lag grid. By default each term is evaluated over
#'   its knot range.
#' @param n Number of default grid points.
#' @param predictor Optional predictor-value grid for nonlinear IRFs.
#' @param n_predictor Number of default predictor grid points.
#' @param level Optional grouping levels for random IRFs. By default all fitted
#'   levels are returned.
#' @param se Include pointwise standard errors when covariance is available.
#' @param unconditional Include smoothing-parameter uncertainty in standard
#'   errors when supported by the fitted backend.
#' @return A data frame with term, lag, estimate, and standard-error columns.
#' @export
estimate_irf <- function(
        object,
        term=NULL,
        lag=NULL,
        n=200,
        predictor=NULL,
        n_predictor=25,
        level=NULL,
        se=TRUE,
        unconditional=FALSE
) {
    if (!is_cdrgam(object)) {
        stop('object must be a fitted cdrgam model')
    }
    metadata <- object$cdrgam$terms
    labels <- object$cdrgam$term_labels
    if (!length(metadata)) {
        stop('The fitted model does not contain stored IRF metadata')
    }
    indices <- if (is.null(term)) {
        seq_along(metadata)
    } else if (is.numeric(term) && length(term) == 1L) {
        as.integer(term)
    } else if (is.character(term) && length(term) == 1L) {
        which(labels == term)
    } else {
        integer()
    }
    if (!length(indices) || any(indices < 1L | indices > length(metadata))) {
        stop('term must uniquely identify a fitted IRF')
    }
    coefficients <- stats::coef(object)
    sparse_covariance <- isTRUE(se) && inherits(object, 'cdrgam_sparse')
    covariance <- if (isTRUE(se) && !sparse_covariance) {
        stats::vcov(object, unconditional=unconditional)
    } else {
        NULL
    }
    output <- list()
    for (j in seq_along(indices)) {
        i <- indices[[j]]
        info <- metadata[[i]]
        lag_grid <- lag
        if (is.null(lag_grid)) {
            lag_grid <- seq(min(info$knots), max(info$knots), length.out=n)
        }
        if (!is.numeric(lag_grid) || any(!is.finite(lag_grid))) {
            stop('lag must be a finite numeric vector')
        }
        is_surface <- startsWith(info$type, 'nonlinear') ||
            startsWith(info$type, 'varying')
        if (is_surface) {
            predictor_grid <- predictor
            if (is.null(predictor_grid)) {
                predictor_grid <- seq(
                    min(info$predictor_knots),
                    max(info$predictor_knots),
                    length.out=n_predictor
                )
            }
            if (!is.numeric(predictor_grid) || any(!is.finite(predictor_grid))) {
                stop('predictor must be a finite numeric vector')
            }
            evaluation_grid <- expand.grid(
                lag=lag_grid,
                predictor=predictor_grid
            )
            basis <- mgcv::PredictMat(
                info$basis,
                list(
                    cdr_delay=evaluation_grid$lag,
                    cdr_value=evaluation_grid$predictor
                ),
                n=nrow(evaluation_grid)
            )
            basis <- basis %*% info$transform
        } else {
            evaluation_grid <- data.frame(
                lag=lag_grid,
                predictor=NA_real_
            )
            basis <- mgcv::PredictMat(
                info$basis,
                list(cdr_delta=lag_grid),
                n=length(lag_grid)
            )
            if (!is.null(info$transform)) {
                basis <- basis %*% info$transform
            }
        }
        coefficient_index <- info$coefficient_index
        grouped <- length(info$group_levels) > 0L
        levels_to_evaluate <- if (grouped) {
            if (is.null(level)) info$group_levels else as.character(level)
        } else {
            NA_character_
        }
        if (grouped && any(!(levels_to_evaluate %in% info$group_levels))) {
            stop('Unknown random-IRF grouping level requested')
        }
        for (group_level in levels_to_evaluate) {
            term_index <- coefficient_index
            if (grouped) {
                group_number <- match(group_level, info$group_levels)
                within_term <- (group_number - 1L) * info$base_dimension +
                    seq_len(info$base_dimension)
                term_index <- coefficient_index[within_term]
            }
            if (length(term_index) != ncol(basis)) {
                stop('Stored IRF coefficient metadata is inconsistent with its basis')
            }
            estimate <- drop(basis %*% coefficients[term_index])
            standard_error <- rep.int(NA_real_, nrow(evaluation_grid))
            if (!is.null(covariance) || sparse_covariance) {
                term_covariance <- if (sparse_covariance) {
                    .sparse_selected_vcov(
                        object,
                        term_index,
                        unconditional=unconditional
                    )
                } else covariance[
                    term_index,
                    term_index,
                    drop=FALSE
                ]
                standard_error <- sqrt(pmax(
                    0,
                    rowSums((basis %*% term_covariance) * basis)
                ))
            }
            output[[length(output) + 1L]] <- data.frame(
                term=labels[[i]],
                group=if (grouped) group_level else NA_character_,
                lag=evaluation_grid$lag,
                predictor=evaluation_grid$predictor,
                estimate=estimate,
                se=standard_error,
                stringsAsFactors=FALSE
            )
        }
    }
    do.call(rbind, output)
}

#' Extract smoothing variance components
#'
#' Report the variance components implied by the smoothing parameters of the
#' expanded `mgcv` model. Native fits are passed directly to
#' [mgcv::gam.vcomp()]. For block-backend fits, the same calculation is applied
#' to the stored REML scale, penalty rescaling metadata, and outer Hessian.
#' This function does not estimate covariance parameters beyond those present
#' in the expanded model.
#'
#' @param object A fitted `cdrgam` model.
#' @param rescale Apply the penalty rescaling used in the original smooth
#'   specification, as in [mgcv::gam.vcomp()].
#' @param conf.lev Confidence level for intervals when an outer REML Hessian is
#'   available.
#' @return The value returned by [mgcv::gam.vcomp()]: normally a matrix of
#'   standard deviations and confidence intervals, or a named vector/list
#'   when intervals are unavailable.
#' @export
variance_components <- function(object, rescale=TRUE, conf.lev=0.95) {
    if (!is_cdrgam(object)) {
        stop('object must be a fitted cdrgam model')
    }
    if (inherits(object, 'gam')) {
        return(mgcv::gam.vcomp(
            object,
            rescale=rescale,
            conf.lev=conf.lev
        ))
    }
    if (!inherits(object, c('cdrgam_block', 'cdrgam_sparse'))) {
        stop('Unsupported cdrgam fitting backend')
    }

    # gam.vcomp() only depends on this subset of a fitted gam object. Using it
    # as the single implementation also keeps linked penalties, parametric
    # penalties, rescaling, and interval conventions aligned with mgcv.
    surrogate <- list(
        sp=object$sp,
        full.sp=object$full.sp,
        smooth=object$smooth,
        paraPen=object$paraPen,
        reml.scale=object$reml.scale,
        sig2=object$sig2,
        method=object$method,
        outer.info=object$outer.info,
        family=object$family
    )
    class(surrogate) <- 'gam'
    mgcv::gam.vcomp(
        surrogate,
        rescale=rescale,
        conf.lev=conf.lev
    )
}

#' Simulate irregular impulse and response streams from fixed IRFs
#'
#' Predictor columns are constructed to be sample-orthogonal and standardized.
#' Responses are the additive convolution of every supplied IRF plus Gaussian
#' noise.
#'
#' @param irfs Named list of functions. One-argument functions define linear
#'   lag IRFs; two-argument functions define surfaces over lag and predictor.
#' @param n_impulses Number of impulses.
#' @param n_responses Number of responses.
#' @param duration Duration of the simulated series.
#' @param window Maximum causal lag included in the convolution.
#' @param intercept Response intercept.
#' @param noise_sd Gaussian response noise standard deviation.
#' @param seed Random seed.
#' @return A list containing impulse and response data, IRFs, and noiseless
#'   response components.
#' @export
simulate_cdr <- function(
        irfs,
        n_impulses=500,
        n_responses=600,
        duration=60,
        window=2,
        intercept=0,
        noise_sd=1,
        seed=NULL
) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    if (!is.list(irfs) || !length(irfs) || is.null(names(irfs)) ||
            any(!nzchar(names(irfs))) ||
            any(!vapply(irfs, is.function, logical(1)))) {
        stop('irfs must be a named, non-empty list of functions')
    }
    p <- length(irfs)
    if (n_impulses <= p || n_responses < 1L || duration <= 0 ||
            window <= 0 || noise_sd < 0) {
        stop('Invalid simulation dimensions, duration, window, or noise_sd')
    }
    impulse_time <- sort(stats::runif(n_impulses, 0, duration))
    response_time <- sort(stats::runif(n_responses, 0, duration))
    raw <- matrix(stats::rnorm(n_impulses * p), n_impulses, p)
    raw <- scale(raw, center=TRUE, scale=FALSE)
    predictors <- qr.Q(qr(raw)) * sqrt(n_impulses)
    colnames(predictors) <- names(irfs)
    impulses <- data.frame(time=impulse_time, predictors, check.names=FALSE)
    responses <- data.frame(time=response_time)

    links <- .build_history_links(
        impulses,
        responses,
        series=character(),
        impulse_time='time',
        response_time='time',
        window=c(0, window)
    )
    components <- matrix(0, nrow=n_responses, ncol=p)
    colnames(components) <- names(irfs)
    nonlinear <- vapply(irfs, function(fun) length(formals(fun)) >= 2L, logical(1))
    for (i in seq_along(irfs)) {
        predictor_value <- predictors[links$impulse_index, i]
        contribution <- if (nonlinear[[i]]) {
            irfs[[i]](links$delay, predictor_value)
        } else {
            predictor_value * irfs[[i]](links$delay)
        }
        if (length(contribution)) {
            summed <- rowsum(
                contribution,
                links$response_index,
                reorder=FALSE
            )
            components[as.integer(rownames(summed)), i] <- summed[, 1L]
        }
    }
    noiseless <- intercept + rowSums(components)
    responses$response <- noiseless + stats::rnorm(n_responses, sd=noise_sd)
    list(
        impulses=impulses,
        responses=responses,
        irfs=irfs,
        nonlinear=stats::setNames(nonlinear, names(irfs)),
        components=components,
        noiseless=noiseless,
        settings=list(
            window=window,
            intercept=intercept,
            noise_sd=noise_sd,
            seed=seed
        )
    )
}
