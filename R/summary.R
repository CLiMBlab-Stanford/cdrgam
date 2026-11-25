.cdrgam_summary_formulas <- function(object) {
    formulas <- object$cdrgam$formula
    fallback <- if (is.null(formulas$user)) stats::formula(object) else
        formulas$user
    output <- list(
        raw=fallback
    )
    if (!is.null(formulas$mgcv)) {
        output$mgcv <- .cdrgam_mgcv_formula_from_fit(object)
    }
    output
}

.cdrgam_summary_simplifications <- function(object) {
    empty <- data.frame(
        parameter=character(), term=character(), axis=character(),
        action=character(), requested=character(), effective=character(),
        reason=character(), stringsAsFactors=FALSE
    )
    preparation <- object$cdrgam$preparation
    if (is.null(preparation)) return(empty)
    if (!is.null(preparation$simplifications)) {
        value <- preparation$simplifications
        if (!nrow(value)) return(empty)
        value$parameter <- NA_character_
        return(value[c('parameter', setdiff(names(value), 'parameter'))])
    }
    parameters <- preparation$parameters
    if (is.null(parameters)) return(empty)
    values <- lapply(names(parameters), function(parameter) {
        value <- parameters[[parameter]]$simplifications
        if (is.null(value) || !nrow(value)) return(NULL)
        value$parameter <- parameter
        value[c('parameter', setdiff(names(value), 'parameter'))]
    })
    values <- Filter(Negate(is.null), values)
    if (!length(values)) empty else do.call(rbind, values)
}

.cdrgam_simplification_text <- function(change) {
    prefix <- if (is.na(change$parameter) || !nzchar(change$parameter)) '' else
        paste0('[', change$parameter, '] ')
    subject <- paste0(prefix, change$term)
    detail <- switch(
        change$action,
        basis_reduced=paste0(
            subject, ' ', change$axis, ' basis dimension: ',
            change$requested, ' -> ', change$effective
        ),
        nonlinear_to_linear=paste0(
            subject, ' ', change$axis, ': nonlinear basis (',
            change$requested, ') -> linear'
        ),
        term_removed=paste0(subject, ': removed'),
        boundary_nullspace_reduction=paste0(
            subject, ' penalty-space dimension: ', change$requested,
            ' -> ', change$effective
        ),
        paste0(
            subject, ' ', change$axis, ': ', change$action, ' (',
            change$requested, ' -> ', change$effective, ')'
        )
    )
    if (is.na(change$reason) || !nzchar(change$reason)) detail else
        paste0(detail, ' (', change$reason, ')')
}

.cdrgam_formula_strings <- function(formulas) {
    stringify <- function(formula) {
        if (inherits(formula, 'formula')) {
            return(paste(deparse(formula), collapse=' '))
        }
        if (is.list(formula)) return(lapply(formula, stringify))
        paste(deparse(formula), collapse=' ')
    }
    lapply(formulas, stringify)
}

.cdrgam_parametric_indices <- function(object) {
    if (!length(object$smooth)) return(seq_along(object$coefficients))
    count <- min(vapply(object$smooth, `[[`, numeric(1), 'first.para')) - 1L
    if (count > 0L) seq_len(count) else integer()
}

.cdrgam_coefficient_table <- function(
        object,
        indices,
        covariance,
        residual_df,
        include_p=TRUE,
        reference=c('t', 'z')
) {
    reference <- match.arg(reference)
    if (!length(indices)) return(NULL)
    coefficients <- object$coefficients[indices]
    variances <- if (is.matrix(covariance)) diag(covariance) else covariance
    standard_errors <- sqrt(pmax(0, variances))
    divisors <- .cdr_coefficient_divisors(object)[indices]
    statistic <- coefficients / standard_errors
    table <- cbind(
        Estimate=coefficients / divisors,
        `Std. Error`=standard_errors / abs(divisors),
        statistic
    )
    colnames(table)[[3L]] <- if (reference == 't') 't value' else 'z value'
    if (isTRUE(include_p)) {
        probability <- if (reference == 't') {
            2 * stats::pt(abs(statistic), df=residual_df, lower.tail=FALSE)
        } else {
            2 * stats::pnorm(abs(statistic), lower.tail=FALSE)
        }
        table <- cbind(
            table,
            probability
        )
        colnames(table)[[4L]] <- if (reference == 't') {
            'Pr(>|t|)'
        } else 'Pr(>|z|)'
    }
    rownames(table) <- names(coefficients)
    table
}

.cdrgam_smooth_labels <- function(object) {
    labels <- vapply(object$smooth, `[[`, character(1), 'label')
    metadata <- object$cdrgam$terms
    if (!length(metadata)) return(labels)
    for (i in seq_along(object$smooth)) {
        indices <- seq.int(
            object$smooth[[i]]$first.para,
            object$smooth[[i]]$last.para
        )
        match_index <- which(vapply(metadata, function(term) {
            identical(as.integer(term$coefficient_index), as.integer(indices))
        }, logical(1)))
        if (length(match_index) == 1L) {
            labels[[i]] <- object$cdrgam$term_labels[[match_index]]
        }
    }
    make.unique(labels)
}

.cdrgam_block_smooth_edf <- function(object) {
    source_weights <- if (is.null(object$working.weights)) {
        object$prior.weights
    } else object$working.weights
    weights <- sqrt(source_weights)
    weighted_design <- object$X * weights
    information <- crossprod(weighted_design)
    influence_covariance <- if (is.null(object$working.weights)) {
        object$Vp / object$scale
    } else object$Vp
    influence_diagonal <- diag(influence_covariance %*% information)
    vapply(object$smooth, function(smooth) {
        sum(influence_diagonal[smooth$first.para:smooth$last.para])
    }, numeric(1))
}

.cdrgam_sparse_smooth_edf <- function(object) {
    components <- object$sparse$penalty_components
    if (!length(components)) {
        return(vapply(object$smooth, function(smooth) {
            smooth$last.para - smooth$first.para + 1L
        }, numeric(1)))
    }
    supports <- .sparse_penalty_supports(components)
    schur_plan <- if (inherits(
            object$sparse$factor,
            'cdrgam_schur_factor'
    )) {
        .sparse_schur_trace_plan(
            list(
                core=object$sparse$factor$core,
                blocks=object$sparse$factor$blocks
            ),
            components
        )
    } else NULL
    traces <- object$sparse$penalty_trace
    if (!is.numeric(traces) || length(traces) != length(components) ||
            any(!is.finite(traces))) {
        traces <- .sparse_logdet_scores(
            object$sparse$factor,
            components,
            object$sp,
            supports=supports,
            chunk_size=if (is.null(object$sparse$control$trace_chunk_size)) {
                64L
            } else object$sparse$control$trace_chunk_size,
            schur_plan=schur_plan
        )
    }
    vapply(object$smooth, function(smooth) {
        indices <- smooth$first.para:smooth$last.para
        owned <- vapply(supports, function(support) {
            length(support) && all(support >= indices[[1L]]) &&
                all(support <= indices[[length(indices)]])
        }, logical(1))
        value <- length(indices) - sum(traces[owned])
        min(length(indices), max(0, value))
    }, numeric(1))
}

# Which smooths are random effects, whose coefficients mgcv's random-effect
# test integrates out when testing another term. Smooths carry a `random`
# flag. Sparse fits saved before the flag existed are recognized by their
# random-effect prediction entries, which record the coefficients they own.
.cdrgam_random_smooths <- function(object) {
    entries <- Filter(function(entry) identical(entry$kind, 'random_effect'),
        object$cdrgam$prediction$random_effects)
    vapply(object$smooth, function(smooth) {
        if (isTRUE(smooth$random)) return(TRUE)
        if (!is.null(smooth$random)) return(FALSE)
        any(vapply(entries, function(entry) {
            identical(as.numeric(range(entry$coefficient_index)),
                as.numeric(c(smooth$first.para, smooth$last.para)))
        }, logical(1)))
    }, logical(1))
}

# Solves M x = rhs for a symmetric positive definite M, sparse or dense.
# Returns NULL when M is not numerically positive definite.
.cdrgam_spd_solve <- function(system, rhs) {
    tryCatch({
        if (methods::is(system, 'sparseMatrix')) {
            factor <- Matrix::Cholesky(Matrix::forceSymmetric(system), LDL=FALSE)
            as.matrix(Matrix::solve(factor, rhs))
        } else {
            root <- chol(as.matrix(system))
            backsolve(root, forwardsolve(t(root), as.matrix(rhs)))
        }
    }, error=function(error) NULL)
}

# The structure Wood's test needs for a block fit: X'WX from the stored
# design and working weights, each smooth's edf1 from the influence matrix,
# and the penalties and covariance for random-effect tests.
.cdrgam_block_test_structure <- function(object) {
    source_weights <- if (is.null(object$working.weights)) {
        object$prior.weights
    } else object$working.weights
    information <- crossprod(object$X * sqrt(source_weights))
    influence_covariance <- if (is.null(object$working.weights)) {
        object$Vp / object$scale
    } else object$Vp
    edf1 <- .cdrgam_edf1(influence_covariance %*% information)
    sp <- if (is.null(object$full.sp)) object$sp else object$full.sp
    dimension <- ncol(information)
    penalties <- list()
    k <- 1L
    for (i in seq_along(object$smooth)) {
        smooth <- object$smooth[[i]]
        for (penalty in smooth$S) {
            columns <- smooth$first.para + seq_len(ncol(penalty)) - 1L
            placed <- Matrix::sparseMatrix(
                i=integer(), j=integer(), x=numeric(), dims=c(dimension, dimension)
            )
            placed[columns, columns] <- sp[[k]] * as.matrix(penalty)
            penalties[[length(penalties) + 1L]] <- list(smooth=i, matrix=placed)
            k <- k + 1L
        }
    }
    list(
        gram=function(indices) information[indices, indices, drop=FALSE],
        edf1=vapply(object$smooth, function(smooth) {
            sum(edf1[smooth$first.para:smooth$last.para])
        }, numeric(1)),
        random=list(
            gram=information,
            penalties=penalties,
            covariance_columns=function(indices) object$Vp[, indices, drop=FALSE],
            sig2=object$sig2
        )
    )
}

# The structure Wood's test needs for a sparse fit, whose factor holds
# H = X'WX + S(lambda). X'WX is recovered as H - S. Each smooth's edf1 uses
# F = I - H^{-1} S, whose rows and columns for the smooth's coefficients
# need one solve per coefficient. Everything stays sparse except blocks
# with one column per tested coefficient.
.cdrgam_sparse_test_structure <- function(object) {
    factor <- object$sparse$factor
    dimension <- object$sparse$dimension
    components <- object$sparse$penalty_components
    penalty <- Matrix::sparseMatrix(
        i=integer(), j=integer(), x=numeric(), dims=c(dimension, dimension)
    )
    for (i in seq_along(components)) {
        penalty <- penalty + object$sp[[i]] * methods::as(components[[i]], 'generalMatrix')
    }
    unit <- function(indices) as.matrix(Matrix::sparseMatrix(
        i=indices, j=seq_along(indices), x=1, dims=c(dimension, length(indices))
    ))
    edf1 <- vapply(object$smooth, function(smooth) {
        indices <- smooth$first.para:smooth$last.para
        identity <- unit(indices)
        # F[, indices] and t(F)[, indices] = (I - S H^{-1})[, indices].
        columns <- identity - as.matrix(.cdr_factor_solve(
            factor, as.matrix(penalty[, indices, drop=FALSE])
        ))
        rows <- identity - as.matrix(penalty %*% .cdr_factor_solve(factor, identity))
        sum(2 * diag(columns[indices, , drop=FALSE]) - colSums(columns * rows))
    }, numeric(1))
    supports <- .sparse_penalty_supports(components)
    owners <- vapply(supports, function(support) {
        owner <- which(vapply(object$smooth, function(smooth) {
            length(support) && all(support >= smooth$first.para) &&
                all(support <= smooth$last.para)
        }, logical(1)))
        if (length(owner) == 1L) owner else NA_integer_
    }, integer(1))
    random <- if (anyNA(owners)) {
        'not computed: a penalty spans several terms'
    } else list(
        gram=.cdr_factor_system(factor) - penalty,
        penalties=lapply(seq_along(owners), function(i) list(
            smooth=owners[[i]],
            matrix=object$sp[[i]] * methods::as(components[[i]], 'generalMatrix')
        )),
        covariance_columns=function(indices) {
            as.matrix(.cdr_factor_solve(factor, unit(indices))) * object$scale
        },
        solve_penalized=function(rhs) as.matrix(.cdr_factor_solve(factor, rhs)),
        sig2=object$sig2
    )
    list(
        gram=function(indices) {
            .cdr_factor_system_columns(factor, indices)[indices, , drop=FALSE] -
                as.matrix(penalty[indices, indices, drop=FALSE])
        },
        edf1=edf1,
        random=random
    )
}

# The structure Wood's test needs for a distributional fit: the stored
# likelihood Hessian plays the role of X'WX, and Vp is at unit scale.
.cdrgam_distributional_test_structure <- function(object) {
    information <- object$distributional$likelihood_hessian
    if (is.null(information)) {
        return(list(unavailable=paste(
            'not computed: the fit predates the stored likelihood Hessian;',
            'refit to test'
        )))
    }
    information <- as.matrix(information)
    edf1 <- .cdrgam_edf1(object$Vp %*% information)
    list(
        gram=function(indices) information[indices, indices, drop=FALSE],
        edf1=vapply(object$smooth, function(smooth) {
            sum(edf1[smooth$first.para:smooth$last.para])
        }, numeric(1)),
        random=paste(
            'not computed: random-effect tests are not implemented for',
            'distributional fits'
        )
    )
}

# Wood's (2013) test that a smooth term is zero, as computed by mgcv's
# summary.gam(). `gram` is the term's block of X'WX, `covariance` its
# coefficient covariance, and `rank` the reference rank, usually the term's
# alternative effective degrees of freedom edf1 capped at its dimension. A
# fractional rank uses Wood's modified pseudo-inverse, whose statistic is
# referred to a weighted sum of chi-squared variables. A positive
# `residual_df` gives the F-type reference for an estimated scale. The
# statistic depends on X only through X'WX, so any square root R with
# R'R = gram gives mgcv's value.
.cdrgam_wood_test <- function(coefficients, gram, covariance, rank,
        residual_df=-1) {
    gram <- (gram + t(gram)) / 2
    root <- eigen(gram, symmetric=TRUE)
    root <- sqrt(pmax(root$values, 0)) * t(root$vectors)
    projected <- root %*% covariance %*% t(root)
    decomposition <- eigen((projected + t(projected)) / 2, symmetric=TRUE)
    signs <- sign(decomposition$vectors[1L, ])
    signs[signs == 0] <- 1
    vectors <- sweep(decomposition$vectors, 2L, signs, '*')
    values <- decomposition$values
    whole <- max(0, floor(rank))
    fraction <- abs(rank - whole)
    upper <- if (fraction > 0) whole + 1L else whole
    estimated <- sum(values > max(values) * .Machine$double.eps^0.9)
    if (estimated < upper) {
        upper <- whole <- estimated
        fraction <- 0
        rank <- estimated
    }
    if (upper < ncol(vectors)) vectors <- vectors[, seq_len(upper), drop=FALSE]
    if (fraction > 0 && whole > 0) {
        if (whole > 1L) {
            leading <- seq_len(whole - 1L)
            vectors[, leading] <- t(t(vectors[, leading]) / sqrt(values[leading]))
        }
        tail <- whole:upper
        mixing <- sqrt(max(0, 0.5 * fraction * (1 - fraction)))
        scaling <- diag(values[tail]^-0.5, nrow=length(tail))
        mixed <- scaling %*% matrix(c(1, mixing, mixing, fraction), 2L, 2L) %*%
            scaling
        mixed <- eigen(mixed, symmetric=TRUE)
        mixed <- mixed$vectors %*% diag(sqrt(mixed$values)) %*% t(mixed$vectors)
        alternative <- vectors
        alternative[, tail] <- t(mixed %*% diag(c(-1, 1)) %*% t(vectors[, tail]))
        vectors[, tail] <- t(mixed %*% t(vectors[, tail]))
    } else {
        vectors <- if (whole == 0) {
            t(t(vectors) * sqrt(1 / values[[1L]]))
        } else t(t(vectors) / sqrt(values[seq_len(whole)]))
        alternative <- vectors
        if (whole == 1) rank <- 1
    }
    scores <- root %*% coefficients
    statistic <- sum(crossprod(vectors, scores)^2)
    alternative_statistic <- sum(crossprod(alternative, scores)^2)
    reference_rank <- rank
    p_value <- 2
    if (fraction > 0) {
        if (upper == 1L) {
            reference_rank <- weights <- 1
        } else {
            weights <- rep.int(1, upper)
            shifted <- fraction + 1
            weights[[whole]] <- (shifted + sqrt(shifted * (2 - shifted))) / 2
            weights[[upper]] <- shifted - weights[[whole]]
        }
        p_value <- if (residual_df <= 0) {
            (mgcv::psum.chisq(statistic, weights) +
                mgcv::psum.chisq(alternative_statistic, weights)) / 2
        } else {
            denominator <- max(1, round(residual_df))
            degrees <- c(rep.int(1, length(weights)), denominator)
            (mgcv::psum.chisq(0, c(weights, -statistic / denominator), df=degrees) +
                mgcv::psum.chisq(
                    0, c(weights, -alternative_statistic / denominator), df=degrees
                )) / 2
        }
    }
    if (p_value > 1) {
        p_value <- if (residual_df <= 0) {
            (stats::pchisq(statistic, reference_rank, lower.tail=FALSE) +
                stats::pchisq(alternative_statistic, reference_rank,
                    lower.tail=FALSE)) / 2
        } else {
            (stats::pf(statistic / reference_rank, reference_rank, residual_df,
                lower.tail=FALSE) +
                stats::pf(alternative_statistic / reference_rank, reference_rank,
                    residual_df, lower.tail=FALSE)) / 2
        }
    }
    list(statistic=statistic, p_value=min(1, p_value), rank=rank)
}

# Wood's (2013b) test that a fully penalized term is zero, giving the values
# of mgcv's reTest() and recov(). `random` supplies the whole-model X'WX as
# `gram`, the penalties as a list of `smooth` index and `matrix` (smoothing
# parameter times penalty, in full coefficient coordinates),
# `covariance_columns(indices)` returning those columns of Vp, `sig2`, and
# optionally `solve_penalized(rhs)` solving (X'WX + S) x = rhs. A positive
# `residual_df` gives the F-type reference for an estimated scale.
#
# Other random-effect terms I are integrated out of the null distribution.
# mgcv forms the triangular factor Rm of the tested columns T from a QR of
# the whole model. Rm only enters through Rm'Rm, which equals the inverse of
# the T block of M_K^{-1}, where K are the columns not integrated out and
#   M_K = G_KK + G_KI S_I^+ G_IK + S_KK
# for G = X'WX, S_I the penalty of the integrated terms, and S_KK the other
# penalties. With nothing integrated out, M is the penalized Hessian. The
# frequentist covariance needs only its T block,
#   Ve_TT = Vp_T. (G + G_.I S_I^+ G_I.) Vp_.T / sig2.
# Each test therefore costs one sparse Cholesky factorization and one solve
# per tested coefficient.
.cdrgam_random_effect_test <- function(object, random, m, residual_df,
        random_smooths=.cdrgam_random_smooths(object)) {
    gram <- random$gram
    dimension <- ncol(gram)
    others <- setdiff(which(random_smooths), m)
    integrated <- rep.int(FALSE, dimension)
    for (i in others) {
        integrated[object$smooth[[i]]$first.para:object$smooth[[i]]$last.para] <- TRUE
    }
    zero <- Matrix::sparseMatrix(
        i=integer(), j=integer(), x=numeric(), dims=c(dimension, dimension)
    )
    fixed_penalty <- zero
    random_penalty <- zero
    for (penalty in random$penalties) {
        if (penalty$smooth %in% others) {
            random_penalty <- random_penalty + penalty$matrix
        } else fixed_penalty <- fixed_penalty + penalty$matrix
    }
    target <- object$smooth[[m]]$first.para:object$smooth[[m]]$last.para
    covariance <- random$covariance_columns(target)
    unit <- function(indices, size) as.matrix(Matrix::sparseMatrix(
        i=indices, j=seq_along(indices), x=1, dims=c(size, length(indices))
    ))
    if (!any(integrated)) {
        solved <- if (is.null(random$solve_penalized)) {
            .cdrgam_spd_solve(gram + fixed_penalty, unit(target, dimension))
        } else random$solve_penalized(unit(target, dimension))
        block <- solved[target, , drop=FALSE]
        frequentist <- crossprod(covariance, as.matrix(gram %*% covariance))
    } else {
        kept <- which(!integrated)
        removed <- which(integrated)
        penalty2 <- as.matrix(random_penalty[removed, removed, drop=FALSE])
        if (max(abs(penalty2 - diag(diag(penalty2), nrow(penalty2)))) == 0) {
            values <- diag(penalty2)
            positive <- values > max(values) * .Machine$double.eps^0.8
            inverse <- Matrix::Diagonal(x=ifelse(positive, 1 / values, 0))
        } else {
            decomposition <- eigen((penalty2 + t(penalty2)) / 2, symmetric=TRUE)
            positive <- decomposition$values >
                max(decomposition$values) * .Machine$double.eps^0.8
            inverse <- decomposition$vectors %*%
                (ifelse(positive, 1 / decomposition$values, 0) *
                    t(decomposition$vectors))
        }
        coupling <- gram[, removed, drop=FALSE]
        marginal <- gram[kept, kept, drop=FALSE] +
            coupling[kept, , drop=FALSE] %*% inverse %*%
                Matrix::t(coupling[kept, , drop=FALSE]) +
            fixed_penalty[kept, kept, drop=FALSE]
        local <- match(target, kept)
        solved <- .cdrgam_spd_solve(marginal, unit(local, length(kept)))
        block <- if (is.null(solved)) NULL else solved[local, , drop=FALSE]
        projected <- as.matrix(Matrix::crossprod(coupling, covariance))
        frequentist <- crossprod(covariance, as.matrix(gram %*% covariance)) +
            crossprod(projected, as.matrix(inverse %*% projected))
    }
    if (is.null(block)) {
        stop('the penalized system for the random-effect test is not positive definite')
    }
    frequentist <- frequentist / random$sig2
    frequentist <- (frequentist + t(frequentist)) / 2
    precision <- solve((block + t(block)) / 2)
    coefficients <- object$coefficients[target]
    statistic <- drop(crossprod(coefficients, precision %*% coefficients)) /
        random$sig2
    covariance_root <- mgcv::mroot(frequentist)
    values <- eigen(crossprod(covariance_root, precision %*% covariance_root) /
        random$sig2, symmetric=TRUE, only.values=TRUE)$values
    values[values < 0] <- 0
    rank <- sum(values > max(values) * .Machine$double.eps^0.8)
    p_value <- if (residual_df > 0) {
        denominator <- max(1, round(residual_df))
        mgcv::psum.chisq(0, c(values, -statistic / denominator),
            df=c(rep.int(1, length(values)), denominator))
    } else mgcv::psum.chisq(statistic, values)
    list(statistic=statistic, p_value=min(1, max(0, p_value)), rank=rank)
}

# Per-coefficient alternative effective degrees of freedom 2 tr(F) - tr(FF),
# where F = V X'WX / scale is the influence matrix in coefficient space.
.cdrgam_edf1 <- function(influence) {
    2 * diag(influence) - rowSums(influence * t(influence))
}

# Tests every smooth with Wood's (2013) statistic, or fully penalized terms
# with Wood's (2013b) random-effect test. `structure` supplies
# `gram(indices)`, the X'WX block of those coefficients; `edf1`, the
# per-smooth alternative effective degrees of freedom; and `random`, the
# inputs to .cdrgam_random_effect_test() or a string giving why that test is
# unavailable. A structure with `unavailable` instead leaves every test
# uncomputed with that reason as its status.
.cdrgam_smooth_table <- function(
        object, edf, covariance, residual_df, chi_square=FALSE, structure=NULL
) {
    count <- length(object$smooth)
    if (!count) return(NULL)
    table <- matrix(
        NA_real_, nrow=count, ncol=4L,
        dimnames=list(
            .cdrgam_smooth_labels(object),
            c('edf', 'Ref.df', if (chi_square) 'Chi.sq' else 'F', 'p-value')
        )
    )
    table[, 'edf'] <- edf
    table[, 'Ref.df'] <- edf
    unavailable <- if (is.null(structure)) {
        'not computed: test structure unavailable'
    } else structure$unavailable
    test_status <- stats::setNames(
        rep.int(if (is.null(unavailable)) '' else unavailable, count),
        rownames(table)
    )
    if (!is.null(unavailable)) {
        attr(table, 'test.status') <- test_status
        return(table)
    }
    statistic_column <- if (chi_square) 'Chi.sq' else 'F'
    reference_df <- if (chi_square) -1 else residual_df
    random_smooths <- NULL
    for (i in seq_len(count)) {
        smooth <- object$smooth[[i]]
        if (is.null(smooth$null.space.dim)) {
            test_status[[i]] <- 'not computed: test structure unavailable'
            next
        }
        if (smooth$null.space.dim <= 0L) {
            if (is.character(structure$random) || is.null(structure$random)) {
                test_status[[i]] <- if (is.character(structure$random)) {
                    structure$random
                } else 'not computed: fully penalized term'
                next
            }
            if (is.null(random_smooths)) {
                random_smooths <- .cdrgam_random_smooths(object)
            }
            test <- .cdrgam_random_effect_test(
                object, structure$random, i, reference_df, random_smooths
            )
            status <- 'Wood (2013b) random effect'
        } else {
            indices <- smooth$first.para:smooth$last.para
            test <- .cdrgam_wood_test(
                object$coefficients[indices], structure$gram(indices),
                covariance(indices), min(length(indices), structure$edf1[[i]]),
                reference_df
            )
            status <- 'Wood (2013)'
        }
        table[i, 'Ref.df'] <- test$rank
        table[i, statistic_column] <- if (chi_square) test$statistic else
            test$statistic / test$rank
        table[i, 'p-value'] <- test$p_value
        test_status[[i]] <- status
    }
    attr(table, 'test.status') <- test_status
    table
}

.cdrgam_custom_summary <- function(
        object,
        covariance,
        smooth_edf,
        dispersion=NULL,
        all.coefficients=FALSE,
        all_variances=NULL,
        test_structure=NULL
) {
    scale <- if (is.null(dispersion)) object$scale else dispersion
    scale_ratio <- scale / object$scale
    total_edf <- if (inherits(object, 'cdrgam_sparse')) {
        .sparse_effective_df(object)
    } else object$edf
    residual_df <- length(object$y) - total_edf
    fixed_dispersion <- object$family$family %in% c('binomial', 'poisson')
    coefficient_reference <- if (fixed_dispersion) 'z' else 't'
    covariance_at_scale <- function(indices) covariance(indices) * scale_ratio
    parametric <- .cdrgam_parametric_indices(object)
    p_table <- if (length(parametric)) {
        .cdrgam_coefficient_table(
            object,
            parametric,
            covariance_at_scale(parametric),
            residual_df,
            reference=coefficient_reference
        )
    } else NULL
    s_table <- .cdrgam_smooth_table(
        object,
        smooth_edf,
        covariance_at_scale,
        residual_df,
        chi_square=fixed_dispersion,
        structure=test_structure
    )
    weights <- object$prior.weights
    mean_response <- sum(weights * object$y) / sum(weights)
    deviance <- if (is.null(object$deviance)) {
        sum(weights * object$residuals^2)
    } else object$deviance
    null_deviance <- if (identical(object$family$family, 'gaussian')) {
        sum(weights * (object$y - mean_response)^2)
    } else {
        sum(object$family$dev.resids(
            object$y,
            rep.int(mean_response, length(object$y)),
            weights
        ))
    }
    formulas <- .cdrgam_summary_formulas(object)
    gcv <- identical(object$method, 'GCV.Cp')
    prediction_error_method <- if (gcv &&
            !is.null(object$sparse$control$criterion)) {
        object$sparse$control$criterion
    } else if (gcv) 'GCV' else '-REML'
    output <- list(
        call=object$call,
        family=object$family,
        formula=formulas$raw,
        formulas=formulas,
        formula_strings=.cdrgam_formula_strings(formulas),
        simplifications=.cdrgam_summary_simplifications(object),
        p.coeff=if (length(parametric)) p_table[, 'Estimate'] else numeric(),
        p.t=if (length(parametric)) p_table[, 3L] else numeric(),
        p.pv=if (length(parametric)) p_table[, 4L] else numeric(),
        p.table=p_table,
        s.table=s_table,
        s.test=if (is.null(s_table)) character() else
            attr(s_table, 'test.status'),
        se=if (length(parametric)) p_table[, 'Std. Error'] else numeric(),
        chi.sq=if (is.null(s_table)) numeric() else
            if ('Chi.sq' %in% colnames(s_table)) s_table[, 'Chi.sq'] else
                s_table[, 'F'] * s_table[, 'Ref.df'],
        s.pv=if (is.null(s_table)) numeric() else s_table[, 'p-value'],
        pTerms.pv=numeric(),
        pTerms.chi.sq=numeric(),
        pTerms.df=numeric(),
        m=length(object$smooth),
        edf=smooth_edf,
        residual.df=residual_df,
        scale=scale,
        dispersion=scale,
        r.sq=if (identical(object$family$family, 'gaussian')) {
            1 - stats::var(sqrt(weights) * object$residuals) *
                (length(object$y) - 1) /
                (stats::var(sqrt(weights) *
                    (object$y - mean_response)) * residual_df)
        } else NULL,
        deviance=deviance,
        null.deviance=null_deviance,
        dev.expl=if (null_deviance > 0) 1 - deviance / null_deviance else NA_real_,
        method=prediction_error_method,
        sp.criterion=if (gcv) unname(object$gcv.ubre) else object$reml,
        rank=length(object$coefficients),
        np=length(object$coefficients),
        n=length(object$y),
        sp=object$sp,
        backend=object$cdrgam$solver,
        rank_metadata=object$cdrgam$rank
    )
    if (isTRUE(all.coefficients)) {
        indices <- seq_along(object$coefficients)
        coefficient_covariance <- if (is.null(all_variances)) {
            covariance_at_scale(indices)
        } else {
            all_variances() * scale_ratio
        }
        output$coefficients <- .cdrgam_coefficient_table(
            object,
            indices,
            coefficient_covariance,
            residual_df,
            reference=coefficient_reference
        )
    }
    output
}

.print_summary_cdrgam <- function(
        x,
        digits=max(3, getOption('digits') - 3),
        signif.stars=getOption('show.signif.stars'),
        ...
) {
    print(x$family)
    labels <- c(raw='Formula', mgcv='mgcv formula')
    for (name in names(labels)) {
        cat(labels[[name]], ':\n', sep='')
        value <- x$formulas[[name]]
        if (is.list(value) && !inherits(value, 'formula')) {
            for (parameter in names(value)) {
                cat('  ', parameter, ': ', sep='')
                cat(paste(deparse(value[[parameter]]), collapse='\n    '), '\n')
            }
        } else {
            cat(paste(deparse(value), collapse='\n'), '\n')
        }
    }
    if (nrow(x$simplifications)) {
        cat('\nAutomatic simplifications:\n')
        for (i in seq_len(nrow(x$simplifications))) {
            cat('- ', .cdrgam_simplification_text(x$simplifications[i, ]), '\n',
                sep='')
        }
    }
    if (length(x$p.coeff)) {
        cat('\nParametric coefficients:\n')
        stats::printCoefmat(
            x$p.table, digits=digits, signif.stars=signif.stars,
            na.print='NA', ...
        )
    }
    cat('\n')
    if (x$m > 0L) {
        cat('Approximate significance of smooth terms:\n')
        stats::printCoefmat(
            x$s.table, digits=digits, signif.stars=signif.stars,
            has.Pvalue=TRUE, na.print='not tested', cs.ind=1L, ...
        )
        if (anyNA(x$s.table[, ncol(x$s.table)])) {
            cat(
                'Tests marked "not tested" were not computed. For fully',
                'penalized terms, this does not indicate fit failure.\n',
                sep=' '
            )
        }
    }
    fit_statistics <- character()
    if (!is.null(x$r.sq) && length(x$r.sq) && is.finite(x$r.sq)) {
        fit_statistics <- c(
            fit_statistics,
            paste0('R-sq.(adj) = ', formatC(x$r.sq, digits=3, width=5))
        )
    }
    if (!is.null(x$dev.expl) && length(x$dev.expl) &&
            is.finite(x$dev.expl)) {
        fit_statistics <- c(fit_statistics, paste0(
            'Deviance explained = ',
            formatC(x$dev.expl * 100, digits=3, width=4), '%'
        ))
    }
    if (length(fit_statistics)) {
        cat('\n', paste(fit_statistics, collapse='  '), '\n', sep='')
    }
    if (!is.null(x$deviance) && length(x$deviance) &&
            is.finite(x$deviance) && !is.null(x$null.deviance) &&
            length(x$null.deviance) && is.finite(x$null.deviance)) {
        cat(
            'Deviance = ', formatC(x$deviance, digits=5),
            '  Null deviance = ', formatC(x$null.deviance, digits=5),
            '\n', sep=''
        )
    }
    cat(x$method, ' = ', formatC(x$sp.criterion, digits=5),
        '  Scale est. = ', formatC(x$scale, digits=5, width=8, flag='-'),
        '  n = ', x$n, '\n', sep='')
    if (!is.null(x$coefficients)) {
        cat('\nAll expanded coefficients (requested):\n')
        stats::printCoefmat(
            x$coefficients, digits=digits, signif.stars=signif.stars,
            na.print='NA', ...
        )
    }
    invisible(x)
}

#' Summarize a fitted CDR-GAM
#'
#' The default output follows [mgcv::summary.gam()]: it reports parametric
#' coefficients and one row per smooth term instead of expanded smooth and
#' random-effect coefficients. It prints the raw user formula, the native
#' `mgcv` export formula, and a direct account of any automatic
#' simplifications. Use [formula.cdrgam()] to inspect normalized and effective
#' formulas.
#'
#' Fully penalized random-effect and grouped IRF terms report effective degrees
#' of freedom and an explicit "not computed" test status. This status does not
#' indicate fit failure. Distributional Gaussian summaries report deviance,
#' null deviance, and deviance explained without an adjusted R-squared
#' statistic. Set `all.coefficients=TRUE` to request the expanded coefficient
#' table.
#'
#' @param object A fitted `cdrgam` model.
#' @param dispersion Optional known dispersion.
#' @param freq,re.test Controls passed to [mgcv::summary.gam()] by native fits.
#'   Custom backends retain these arguments for interface compatibility.
#' @param all.coefficients Include the expanded smooth and random-effect
#'   coefficient table. The default is `FALSE`.
#' @param x A CDR-GAM summary object.
#' @param ... Additional arguments passed to native [mgcv::summary.gam()] or
#'   to the print method.
#' @return A `summary.gam`-style object containing parametric and smooth-term
#'   tables, the raw and translated `mgcv` formulas, and recorded automatic
#'   simplifications. Native fits retain the complete `summary.gam` result.
#'   Custom backends provide its commonly used coefficient, smooth-term,
#'   fit-statistic, and formula fields. The `s.test` field states whether each
#'   smooth-term test was computed.
#' @name summary.cdrgam
#' @aliases summary.cdrgam_sparse
NULL

#' @rdname summary.cdrgam
#' @export
print.summary.cdrgam <- function(x, ...) .print_summary_cdrgam(x, ...)
