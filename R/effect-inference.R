# Sets the seed for the duration of the caller and returns a function that
# restores the previous random state.
.cdr_effect_seed <- function(seed) {
    if (is.null(seed)) return(function() invisible(NULL))
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
        stop('seed must be one finite number')
    }
    exists <- exists('.Random.seed', envir=globalenv(), inherits=FALSE)
    if (exists) state <- get('.Random.seed', envir=globalenv())
    set.seed(as.integer(seed))
    function() {
        if (exists) {
            assign('.Random.seed', state, envir=globalenv())
        } else if (exists('.Random.seed', envir=globalenv(), inherits=FALSE)) {
            rm('.Random.seed', envir=globalenv())
        }
        invisible(NULL)
    }
}

.cdr_inference_effect_fields <- c(
    'terms', 'axes', 'composition', 'grouping', 'groups', 'reference', 'n',
    'points'
)

# Resolves each named effect specification to its functionals and grid,
# with defaults as in estimate_effect().
.cdr_inference_effects <- function(object, effects) {
    if (!is.list(effects) || !length(effects) || is.null(names(effects)) ||
            any(!nzchar(names(effects))) || anyDuplicated(names(effects))) {
        stop('effects must be a non-empty, uniquely named list of effect specifications')
    }
    lapply(names(effects), function(name) {
        specification <- effects[[name]]
        if (!is.list(specification)) stop('Effect ', name, ' must be a list')
        unknown <- setdiff(names(specification), .cdr_inference_effect_fields)
        if (length(unknown)) stop(
            'Unknown fields in effect ', name, ': ', paste(unknown, collapse=', ')
        )
        axes <- if (is.null(specification$axes)) list() else specification$axes
        reference <- if (is.null(specification$reference)) list() else
            specification$reference
        n <- .cdr_effect_validate(
            object, axes, reference, 0.95,
            if (is.null(specification$n)) 200L else specification$n,
            specification$points
        )
        composition <- match.arg(specification$composition, c('term', 'total'))
        grouping <- match.arg(
            specification$grouping, c('population', 'deviation', 'conditional')
        )
        functionals <- .cdr_effect_all_functionals(
            object, specification$terms, axes, composition, grouping,
            specification$groups, reference, n, specification$points
        )
        frames <- lapply(functionals, function(functional) {
            .cdr_effect_frame(object, functional$term, functional, composition, grouping)
        })
        list(name=name, functionals=functionals, grid=.cdr_effect_bind(frames))
    })
}

# The stacked design of every effect row over the union of coefficients.
.cdr_inference_design <- function(resolved) {
    functionals <- unlist(lapply(resolved, `[[`, 'functionals'), recursive=FALSE)
    index <- sort(unique(unlist(lapply(functionals, `[[`, 'index'))))
    design <- do.call(rbind, lapply(functionals, function(functional) {
        rows <- matrix(0, nrow(functional$design), length(index))
        rows[, match(functional$index, index)] <- as.matrix(functional$design)
        rows
    }))
    list(index=index, design=design)
}

# A linear estimand as a matrix with one row per estimand and one column per
# stacked effect row. A list gives each estimand as weights by effect name;
# a weight of length one is recycled over that effect's rows.
.cdr_inference_linear <- function(estimand, sizes) {
    total <- sum(sizes)
    if (is.null(estimand)) {
        labels <- unlist(lapply(names(sizes), function(name) {
            paste0(name, '[', seq_len(sizes[[name]]), ']')
        }))
        return(structure(diag(total), dimnames=list(labels, NULL)))
    }
    if (is.matrix(estimand)) {
        if (!is.numeric(estimand) || ncol(estimand) != total) stop(
            'A matrix estimand needs one column per effect row (', total, ')'
        )
        if (is.null(rownames(estimand))) {
            rownames(estimand) <- paste0('estimand', seq_len(nrow(estimand)))
        }
        return(estimand)
    }
    if (is.null(names(estimand)) || any(!nzchar(names(estimand))) ||
            anyDuplicated(names(estimand))) {
        stop('A list estimand must be uniquely named')
    }
    offsets <- c(0L, cumsum(sizes))[seq_along(sizes)]
    names(offsets) <- names(sizes)
    rows <- lapply(names(estimand), function(label) {
        weights <- estimand[[label]]
        if (!is.list(weights) || is.null(names(weights)) ||
                !all(names(weights) %in% names(sizes))) stop(
            'Estimand ', label, ' must be a list of weights named by effect'
        )
        row <- numeric(total)
        for (name in names(weights)) {
            value <- weights[[name]]
            if (!is.numeric(value) || !(length(value) %in% c(1L, sizes[[name]]))) {
                stop('Weights for effect ', name, ' in estimand ', label,
                    ' must have length one or ', sizes[[name]])
            }
            columns <- offsets[[name]] + seq_len(sizes[[name]])
            row[columns] <- row[columns] + rep_len(value, sizes[[name]])
        }
        row
    })
    structure(do.call(rbind, rows), dimnames=list(names(estimand), NULL))
}

# Applies a function estimand to stacked effect values, split by effect.
.cdr_inference_apply <- function(estimand, values, sizes, grids) {
    pieces <- split(values, rep.int(names(sizes), sizes))[names(sizes)]
    result <- estimand(pieces, grids)
    if (!is.numeric(result) || !length(result)) {
        stop('estimand must return a non-empty numeric vector')
    }
    result
}

# A symmetric square root of a covariance matrix.
.cdr_inference_root <- function(covariance) {
    decomposition <- eigen((covariance + t(covariance)) / 2, symmetric=TRUE)
    decomposition$vectors %*%
        diag(sqrt(pmax(0, decomposition$values)), nrow(covariance))
}

# A Wald test that every estimand equals its null, using the pseudo-inverse
# of the estimand covariance.
.cdr_inference_joint <- function(difference, covariance, tolerance=1e-8) {
    keep <- is.finite(difference)
    if (!any(keep)) return(list(statistic=NA_real_, df=0L, p_value=NA_real_))
    covariance <- covariance[keep, keep, drop=FALSE]
    decomposition <- eigen((covariance + t(covariance)) / 2, symmetric=TRUE)
    values <- decomposition$values
    positive <- values > tolerance * max(c(values, 0))
    rank <- sum(positive)
    if (!rank) return(list(statistic=NA_real_, df=0L, p_value=NA_real_))
    scores <- crossprod(decomposition$vectors[, positive, drop=FALSE],
        difference[keep])
    statistic <- sum(scores^2 / values[positive])
    list(statistic=statistic, df=rank,
        p_value=stats::pchisq(statistic, rank, lower.tail=FALSE))
}

#' Inference for summary statistics of fitted effects
#'
#' Estimates, interval-estimates, and tests summary statistics of one or
#' more fitted effects, such as contrasts between groups, areas under an
#' impulse-response function, averages over observed covariates, or peak
#' lags. Every effect is evaluated as in [estimate_effect()], and all of
#' them are evaluated on one joint distribution of the coefficients they
#' use, so summaries that combine effects with shared coefficients have the
#' right variance.
#'
#' The effects are a linear map `D beta` of the coefficients. A linear
#' estimand `L D beta` has covariance `L D V D' L'`, and Wald inference is
#' exact given the Gaussian approximation to the coefficient posterior.
#' Give it as a matrix or a list of weights, or leave `estimand` `NULL` to
#' make every effect row an estimand.
#'
#' A function estimand is nonlinear in general. `method = "simulation"`
#' draws `beta ~ N(beta-hat, V)` and applies the function to each draw. The
#' estimate is the function at `beta-hat`, the standard error is the
#' standard deviation of the draws, the interval uses their quantiles, and
#' the p-value is the two-sided posterior tail area beyond `null`, whose
#' resolution is `1 / draws`. `method = "delta"` instead linearizes the
#' function at the estimate by central differences. It is invalid for
#' functions that are not differentiable there, such as maxima and
#' argmaxima.
#'
#' Simultaneous intervals and adjusted p-values cover all estimands of one
#' call. They use the maximum over estimands of the standardized absolute
#' deviation, from Gaussian draws with the estimand correlation for Wald and
#' delta inference, or from the coefficient draws for simulation. The
#' adjusted p-values are single-step max-modulus p-values, estimated from
#' those draws and never below the unadjusted p-value.
#'
#' @inheritParams estimate_effect
#' @param effects A named list of effect specifications. Each is a list with
#'   any of `terms`, `axes`, `composition`, `grouping`, `groups`,
#'   `reference`, `n`, and `points`, as in [estimate_effect()].
#' @param estimand `NULL`, a matrix with one row per estimand and one column
#'   per effect row (effects stacked in order), a named list whose elements
#'   are lists of weights named by effect (a weight of length one is
#'   recycled over that effect's rows), or a function. A function receives a
#'   named list of numeric effect values and a named list of effect grids,
#'   and returns a numeric vector, named to label the estimands.
#' @param null The value of each estimand under the null hypothesis, of
#'   length one or one per estimand.
#' @param method `"auto"` uses Wald inference for linear estimands and
#'   simulation for function estimands. `"delta"` applies only to function
#'   estimands.
#' @param draws Number of coefficient draws for simulation, and of Gaussian
#'   draws for simultaneous critical values.
#' @param seed Optional seed. The previous random state is restored.
#' @param keep_draws If `TRUE`, store the estimand draws of a simulation as
#'   attribute `draws`, with one row per estimand.
#' @return A `cdrgam_effect_inference` data frame with one row per estimand:
#'   `estimand`, `estimate`, `std_error`, `lower`, `upper`, `null`,
#'   `statistic` (the standardized difference from `null`), `p_value`,
#'   `simultaneous_lower`, `simultaneous_upper`, `adjusted_p_value`, and
#'   `method`. Attributes hold `covariance`, the estimand covariance;
#'   `joint_test`, a Wald test that every estimand equals its null; `grids`,
#'   the effect grids with their estimates; `level`; and `draws` when kept.
#' @seealso [estimate_effect()], [test_effect()], [simulate_effect()]
#' @export
effect_inference <- function(
        object, effects, estimand=NULL, null=0,
        method=c('auto', 'wald', 'simulation', 'delta'), level=0.95,
        draws=2000L, seed=NULL, unconditional=FALSE,
        boundary=c('estimate', 'bound'), impulses=NULL, responses=NULL,
        keep_draws=FALSE
) {
    .cdr_effect_validate(object, list(), list(), level, 2L, NULL)
    method <- match.arg(method)
    draws <- suppressWarnings(as.integer(draws))
    if (length(draws) != 1L || is.na(draws) || draws < 2L) {
        stop('draws must be an integer of at least two')
    }
    functional <- is.function(estimand)
    if (!functional && !is.null(estimand) && !is.matrix(estimand) &&
            !is.list(estimand)) {
        stop('estimand must be NULL, a matrix, a list of weights, or a function')
    }
    if (identical(method, 'auto')) method <- if (functional) 'simulation' else 'wald'
    if (!functional && method %in% c('simulation', 'delta')) stop(
        'method = "', method, '" needs a function estimand; ',
        'linear estimands use exact Wald inference'
    )
    if (functional && identical(method, 'wald')) stop(
        'method = "wald" needs a linear estimand; use "simulation" or "delta"'
    )
    restore <- .cdr_effect_seed(seed)
    on.exit(restore(), add=TRUE)
    object <- .cdr_effect_restore_summaries(object, impulses, responses)
    object <- .cdr_effect_boundary_object(object, match.arg(boundary), unconditional)
    resolved <- .cdr_inference_effects(object, effects)
    names(resolved) <- names(effects)
    stacked <- .cdr_inference_design(resolved)
    sizes <- vapply(resolved, function(effect) nrow(effect$grid), integer(1))
    coefficients <- stats::coef(object)[stacked$index]
    covariance <- as.matrix(.cdr_effect_covariance(object, stacked$index, unconditional))
    values <- drop(stacked$design %*% coefficients)
    grids <- lapply(seq_along(resolved), function(i) {
        grid <- resolved[[i]]$grid
        offset <- sum(sizes[seq_len(i - 1L)])
        grid$estimate <- values[offset + seq_len(sizes[[i]])]
        grid
    })
    names(grids) <- names(resolved)
    alpha <- 1 - level
    estimand_draws <- NULL
    if (identical(method, 'wald')) {
        weights <- .cdr_inference_linear(estimand, sizes)
        mapping <- weights %*% stacked$design
        estimate <- drop(mapping %*% coefficients)
        estimand_covariance <- mapping %*% covariance %*% t(mapping)
        labels <- rownames(weights)
    } else {
        estimate <- .cdr_inference_apply(estimand, values, sizes, grids)
        labels <- names(estimate)
        if (is.null(labels)) labels <- paste0('estimand', seq_along(estimate))
        if (identical(method, 'delta')) {
            jacobian <- matrix(0, length(estimate), length(values))
            for (k in seq_along(values)) {
                step <- 1e-5 * max(1, abs(values[[k]]))
                up <- down <- values
                up[[k]] <- up[[k]] + step
                down[[k]] <- down[[k]] - step
                jacobian[, k] <- (
                    .cdr_inference_apply(estimand, up, sizes, grids) -
                        .cdr_inference_apply(estimand, down, sizes, grids)
                ) / (2 * step)
            }
            mapping <- jacobian %*% stacked$design
            estimand_covariance <- mapping %*% covariance %*% t(mapping)
        } else {
            coefficient_draws <- coefficients + .cdr_inference_root(covariance) %*%
                matrix(stats::rnorm(length(coefficients) * draws),
                    length(coefficients), draws)
            value_draws <- stacked$design %*% coefficient_draws
            estimand_draws <- vapply(seq_len(draws), function(s) {
                result <- .cdr_inference_apply(estimand, value_draws[, s], sizes, grids)
                if (length(result) != length(estimate)) stop(
                    'estimand returned a different number of values for a draw'
                )
                as.numeric(result)
            }, numeric(length(estimate)))
            estimand_draws <- matrix(estimand_draws, nrow=length(estimate))
            estimand_covariance <- stats::cov(t(estimand_draws))
        }
    }
    estimand_covariance <- (estimand_covariance + t(estimand_covariance)) / 2
    dimnames(estimand_covariance) <- list(labels, labels)
    count <- length(estimate)
    if (!length(null) %in% c(1L, count) || !is.numeric(null) || anyNA(null)) {
        stop('null must be numeric with length one or one per estimand (', count, ')')
    }
    null <- rep_len(null, count)
    standard_error <- sqrt(pmax(0, diag(estimand_covariance)))
    informative <- standard_error > 0
    statistic <- ifelse(informative, (estimate - null) / standard_error, NA_real_)
    if (identical(method, 'simulation')) {
        lower <- apply(estimand_draws, 1L, stats::quantile, probs=alpha / 2,
            names=FALSE)
        upper <- apply(estimand_draws, 1L, stats::quantile, probs=1 - alpha / 2,
            names=FALSE)
        p_value <- vapply(seq_len(count), function(j) {
            if (!informative[[j]]) return(NA_real_)
            min(1, 2 * min(mean(estimand_draws[j, ] <= null[[j]]),
                mean(estimand_draws[j, ] >= null[[j]])))
        }, numeric(1))
        deviations <- abs(estimand_draws - estimate) / standard_error
    } else {
        multiplier <- stats::qnorm(1 - alpha / 2)
        lower <- estimate - multiplier * standard_error
        upper <- estimate + multiplier * standard_error
        p_value <- ifelse(informative, 2 * stats::pnorm(-abs(statistic)), NA_real_)
        correlation <- estimand_covariance[informative, informative, drop=FALSE] /
            tcrossprod(standard_error[informative])
        deviations <- matrix(NA_real_, count, draws)
        if (any(informative)) {
            deviations[informative, ] <- abs(.cdr_inference_root(correlation) %*%
                matrix(stats::rnorm(sum(informative) * draws), sum(informative), draws))
        }
    }
    maxima <- apply(deviations[informative, , drop=FALSE], 2L, max)
    critical <- if (any(informative)) {
        stats::quantile(maxima, level, names=FALSE)
    } else NA_real_
    adjusted <- vapply(seq_len(count), function(j) {
        if (!informative[[j]]) return(NA_real_)
        max(p_value[[j]], mean(maxima >= abs(statistic[[j]])))
    }, numeric(1))
    result <- data.frame(
        estimand=labels,
        estimate=as.numeric(estimate),
        std_error=standard_error,
        lower=lower,
        upper=upper,
        null=null,
        statistic=statistic,
        p_value=p_value,
        simultaneous_lower=estimate - critical * standard_error,
        simultaneous_upper=estimate + critical * standard_error,
        adjusted_p_value=adjusted,
        method=method,
        stringsAsFactors=FALSE
    )
    class(result) <- c('cdrgam_effect_inference', class(result))
    attr(result, 'covariance') <- estimand_covariance
    attr(result, 'joint_test') <- .cdr_inference_joint(
        estimate - null, estimand_covariance
    )
    attr(result, 'grids') <- grids
    attr(result, 'level') <- level
    if (isTRUE(keep_draws) && !is.null(estimand_draws)) {
        attr(result, 'draws') <- structure(estimand_draws, dimnames=list(labels, NULL))
    }
    result
}
