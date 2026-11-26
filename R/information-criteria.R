# Information criteria for block and sparse fits, with mgcv's correction of
# the effective degrees of freedom for smoothing-parameter uncertainty
# (Wood, Pya, and Saefken 2016). mgcv's logLik.gam() uses edf2, computed in
# gam.fit3.post.proc(): the diagonal of Vc X'WX / scale, where Vc adds to Vp
# a first-order term for the dependence of the coefficients on the log
# smoothing parameters and a second-order term for the dependence of the
# Cholesky factor of the penalized system on them. The second term needs a
# dense p-by-p derivative per smoothing parameter, so it is computed only up
# to getOption('cdrgam.edf2_max_coefficients', 3000) coefficients.

.cdrgam_edf2_default_limit <- 3000L

# The Gaussian log-likelihood mgcv's logLik.gam() reports: one for the
# estimated scale minus half the family's AIC at the deviance, which
# evaluates the likelihood at the maximum-likelihood scale deviance / n
# rather than at the REML estimate. Other families keep their own
# log-likelihoods, which already agree with mgcv's (NULL here).
.cdrgam_mgcv_loglik <- function(object, deviance) {
    if (!identical(object$family$family, 'gaussian') || is.null(object$family$aic)) {
        return(NULL)
    }
    n <- length(object$y)
    1 - object$family$aic(
        object$y, rep.int(1, n), object$fitted.values, object$prior.weights,
        deviance
    ) / 2
}

# Derivative of the upper Cholesky factor R of A = R'R in the direction dA:
# Phi(R^-T dA R^-1) R, where Phi keeps the upper triangle and halves the
# diagonal. It solves R' dR + dR' R = dA with dR upper triangular, as
# mgcv's dchol() does.
.cdrgam_cholesky_derivative <- function(dA, R) {
    inner <- forwardsolve(t(R), t(forwardsolve(t(R), t(dA))))
    inner[lower.tri(inner)] <- 0
    diag(inner) <- diag(inner) / 2
    inner %*% R
}

# The pieces of the corrected covariance for a block or sparse fit: dense
# X'WX, the scaled penalty matrices with one smoothing parameter each, the
# coefficient covariance, per-coefficient edf1, and the outer Hessian.
.cdrgam_edf2_inputs <- function(object) {
    if (inherits(object, 'cdrgam_sparse')) {
        structure <- .cdrgam_sparse_test_structure(object)
        dimension <- object$sparse$dimension
        components <- object$sparse$penalty_components
        penalties <- lapply(seq_along(components), function(i) {
            as.matrix(object$sp[[i]] * methods::as(components[[i]], 'generalMatrix'))
        })
        penalty <- Reduce(`+`, penalties, matrix(0, dimension, dimension))
        gram <- as.matrix(.cdr_factor_system(object$sparse$factor)) - penalty
        covariance <- .sparse_selected_vcov(object, seq_len(dimension))
        hessian <- object$outer.info$hess
        if (is.null(hessian)) hessian <- .sparse_resolve_outer_hessian(object)
        edf1 <- structure$edf1_coefficients(seq_len(dimension))
    } else {
        structure <- .cdrgam_block_test_structure(object)
        gram <- as.matrix(structure$random$gram)
        penalties <- lapply(structure$random$penalties, function(penalty) {
            as.matrix(penalty$matrix)
        })
        covariance <- as.matrix(object$Vp)
        hessian <- object$outer.info$hess
        edf1 <- structure$edf1_coefficients(seq_len(ncol(gram)))
    }
    list(gram=(gram + t(gram)) / 2, penalties=penalties, covariance=covariance,
        hessian=hessian, edf1=edf1)
}

# Per-coefficient edf2 of a block or sparse fit, as mgcv computes it, or an
# error giving why it cannot be computed.
.cdrgam_edf2 <- function(object) {
    if (isTRUE(object$cdrgam$distributional)) stop(
        'The smoothing-parameter correction of the effective degrees of ',
        'freedom is not implemented for distributional fits', call.=FALSE
    )
    if (!(object$method %in% c('REML', 'ML'))) stop(
        'The smoothing-parameter correction of the effective degrees of ',
        'freedom needs a REML or ML fit; this fit used ', object$method,
        call.=FALSE
    )
    dimension <- length(object$coefficients)
    limit <- getOption('cdrgam.edf2_max_coefficients', .cdrgam_edf2_default_limit)
    if (dimension > limit) stop(
        'The smoothing-parameter correction of the effective degrees of ',
        'freedom needs dense ', dimension, '-by-', dimension, ' matrices, one ',
        'per smoothing parameter, which exceeds ',
        'cdrgam.edf2_max_coefficients = ', limit, '. Raise the option to ',
        'compute it, or compare models by held-out prediction.', call.=FALSE
    )
    inputs <- .cdrgam_edf2_inputs(object)
    count <- length(inputs$penalties)
    if (!count) return(inputs$edf1)
    if (count != length(object$sp)) stop(
        'The correction needs one smoothing parameter per penalty; this fit ',
        'links smoothing parameters', call.=FALSE
    )
    hessian <- inputs$hessian
    if (is.null(hessian) || nrow(hessian) < count || any(!is.finite(hessian))) stop(
        'The correction needs the outer Hessian of the smoothing-parameter ',
        'criterion, which this fit does not provide', call.=FALSE
    )
    scale <- object$scale
    coefficients <- object$coefficients
    penalized <- inputs$gram + Reduce(`+`, inputs$penalties)
    root <- tryCatch(chol(penalized), error=function(error) NULL)
    if (is.null(root)) stop(
        'The penalized system is not positive definite, so the correction ',
        'cannot be computed', call.=FALSE
    )
    # Derivatives of the coefficients with respect to each log smoothing
    # parameter: -H^{-1} lambda_i S_i beta.
    derivatives <- -backsolve(root, forwardsolve(t(root), vapply(
        inputs$penalties, function(penalty) drop(penalty %*% coefficients),
        numeric(dimension)
    )))
    derivatives <- matrix(derivatives, dimension, count)
    decomposition <- eigen((hessian + t(hessian)) / 2, symmetric=TRUE)
    values <- decomposition$values
    nonpositive <- values <= 0
    inverse_root <- ifelse(nonpositive, 0, 1 / sqrt(pmax(values, .Machine$double.xmin)))
    first <- (inverse_root * t(decomposition$vectors))[, seq_len(count), drop=FALSE]
    first_order <- crossprod(first %*% t(derivatives))
    # mgcv regularizes the smoothing-parameter covariance of the second-order
    # term by adding 1/10 to the Hessian eigenvalues, and drops a trailing
    # log-scale parameter.
    regularized <- ifelse(nonpositive, 0, values)
    smoothing <- crossprod((1 / sqrt(regularized + 1 / 10)) * t(decomposition$vectors))
    smoothing <- smoothing[seq_len(count), seq_len(count), drop=FALSE]
    inverse_derivatives <- lapply(inputs$penalties, function(penalty) {
        derivative <- .cdrgam_cholesky_derivative(penalty, root)
        -backsolve(root, t(backsolve(root, t(derivative), transpose=TRUE)), transpose=FALSE)
    })
    decomposed <- eigen((smoothing + t(smoothing)) / 2, symmetric=TRUE)
    second_order <- matrix(0, dimension, dimension)
    for (l in seq_along(decomposed$values)) {
        if (decomposed$values[[l]] <= 0) next
        combined <- Reduce(`+`, Map(`*`, inverse_derivatives, decomposed$vectors[, l]))
        second_order <- second_order + decomposed$values[[l]] * tcrossprod(combined)
    }
    corrected <- inputs$covariance + first_order + scale * second_order
    edf2 <- rowSums(corrected * inputs$gram) / scale
    if (sum(edf2) > sum(inputs$edf1)) edf2 <- inputs$edf1
    edf2
}

# Effective degrees of freedom for an information criterion: the corrected
# edf plus one for an estimated scale, capped at the number of coefficients
# plus that one, as in mgcv's logLik.gam().
.cdrgam_information_df <- function(object) {
    estimated_scale <- as.integer(
        object$family$family %in% c('gaussian', 'Gamma') ||
            isTRUE(object$scale.estimated)
    )
    min(sum(.cdrgam_edf2(object)) + estimated_scale,
        length(object$coefficients) + estimated_scale)
}

.cdrgam_information_criterion <- function(objects, penalty) {
    values <- vapply(objects, function(object) {
        log_likelihood <- stats::logLik(object)
        -2 * as.numeric(log_likelihood) +
            penalty(attr(log_likelihood, 'nobs')) * .cdrgam_information_df(object)
    }, numeric(1))
    if (length(objects) == 1L) return(values[[1L]])
    data.frame(
        df=vapply(objects, .cdrgam_information_df, numeric(1)),
        criterion=values
    )
}

.cdrgam_information_objects <- function(object, ..., call) {
    objects <- c(list(object), list(...))
    labels <- vapply(as.list(call)[-1L], function(value) {
        paste(deparse(value), collapse='')
    }, character(1))
    labels <- labels[!names(as.list(call)[-1L]) %in% 'k']
    names(objects) <- labels[seq_along(objects)]
    objects
}

#' Information criteria for block and sparse fits
#'
#' `AIC()` and `BIC()` of block and sparse fits correct the effective degrees
#' of freedom for smoothing-parameter uncertainty, as `mgcv` does for `gam`
#' fits through `logLik.gam()` (Wood, Pya, and Saefken 2016). Uncorrected,
#' the criterion favours overly complex models, especially ones with random
#' effects. The correction needs dense matrices of the size of the model,
#' one per smoothing parameter, so a fit with more coefficients than
#' `getOption("cdrgam.edf2_max_coefficients", 3000)` gives an error rather
#' than an uncorrected value. It needs a REML or ML fit.
#'
#' Criteria are comparable only between models fitted to the same responses,
#' with the same response variable on the same scale.
#'
#' @param object A block or sparse `cdrgam` fit.
#' @param ... Further fits to compare.
#' @param k The penalty per effective degree of freedom; 2 gives AIC.
#' @return A number for one fit; for several, a data frame with `df` and the
#'   criterion, one row per fit.
#' @references Wood, S. N., Pya, N., and Saefken, B. (2016). Smoothing
#'   parameter and model selection for general smooth models. *Journal of the
#'   American Statistical Association*, 111(516), 1548--1563.
#' @name information-criteria
#' @export
AIC.cdrgam_block <- function(object, ..., k=2) {
    objects <- .cdrgam_information_objects(object, ..., call=match.call())
    output <- .cdrgam_information_criterion(objects, function(n) k)
    if (is.data.frame(output)) {
        names(output)[[2L]] <- 'AIC'
        rownames(output) <- names(objects)
    }
    output
}

#' @rdname information-criteria
#' @export
AIC.cdrgam_sparse <- AIC.cdrgam_block

#' @rdname information-criteria
#' @export
BIC.cdrgam_block <- function(object, ...) {
    objects <- .cdrgam_information_objects(object, ..., call=match.call())
    output <- .cdrgam_information_criterion(objects, function(n) log(n))
    if (is.data.frame(output)) {
        names(output)[[2L]] <- 'BIC'
        rownames(output) <- names(objects)
    }
    output
}

#' @rdname information-criteria
#' @export
BIC.cdrgam_sparse <- BIC.cdrgam_block
