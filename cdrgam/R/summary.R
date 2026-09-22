.cdrgam_summary_formulas <- function(object) {
    formulas <- object$cdrgam$formula
    fallback <- if (is.null(formulas$user)) stats::formula(object) else
        formulas$user
    list(
        user=fallback,
        normalized=if (is.null(formulas$normalized)) fallback else
            formulas$normalized,
        effective=if (is.null(formulas$effective)) fallback else
            formulas$effective
    )
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

.cdrgam_smooth_table <- function(
        object, edf, covariance, residual_df, chi_square=FALSE
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
    for (i in seq_len(count)) {
        smooth <- object$smooth[[i]]
        if (is.null(smooth$null.space.dim) || smooth$null.space.dim <= 0L) next
        indices <- smooth$first.para:smooth$last.para
        selected <- covariance(indices)
        decomposition <- eigen(selected, symmetric=TRUE)
        values <- decomposition$values
        tolerance <- max(values) * sqrt(.Machine$double.eps)
        rank <- sum(values > tolerance)
        if (!rank) next
        retained <- values > tolerance
        precision <- decomposition$vectors[, retained, drop=FALSE] %*%
            ((1 / values[retained]) *
                t(decomposition$vectors[, retained, drop=FALSE]))
        coefficients <- object$coefficients[indices]
        statistic <- drop(crossprod(coefficients, precision %*% coefficients))
        reference_df <- max(edf[[i]], 1)
        table[i, 'Ref.df'] <- reference_df
        statistic_column <- if (chi_square) 'Chi.sq' else 'F'
        table[i, statistic_column] <- if (chi_square) {
            statistic
        } else statistic / reference_df
        table[i, 'p-value'] <- if (chi_square) {
            stats::pchisq(statistic, df=reference_df, lower.tail=FALSE)
        } else {
            stats::pf(
                table[i, 'F'],
                df1=reference_df,
                df2=residual_df,
                lower.tail=FALSE
            )
        }
    }
    table
}

.cdrgam_custom_summary <- function(
        object,
        covariance,
        smooth_edf,
        dispersion=NULL,
        all.coefficients=FALSE,
        all_variances=NULL
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
        chi_square=fixed_dispersion
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
    output <- list(
        call=object$call,
        family=object$family,
        formula=formulas$user,
        formulas=formulas,
        formula_strings=.cdrgam_formula_strings(formulas),
        p.coeff=if (length(parametric)) p_table[, 'Estimate'] else numeric(),
        p.t=if (length(parametric)) p_table[, 3L] else numeric(),
        p.pv=if (length(parametric)) p_table[, 4L] else numeric(),
        p.table=p_table,
        s.table=s_table,
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
        } else NA_real_,
        dev.expl=if (null_deviance > 0) 1 - deviance / null_deviance else NA_real_,
        method='-REML',
        sp.criterion=object$reml,
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
    labels <- c(user='Formula', normalized='Normalized formula',
        effective='Effective formula')
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
            has.Pvalue=TRUE, na.print='NA', cs.ind=1L, ...
        )
        if (anyNA(x$s.table[, 3L])) {
            cat(
                'Term-level tests are omitted for fully penalized',
                'random-effect and grouped-deviation terms.\n',
                sep=' '
            )
        }
    }
    cat('\nR-sq.(adj) = ', formatC(x$r.sq, digits=3, width=5),
        '  Deviance explained = ',
        formatC(x$dev.expl * 100, digits=3, width=4), '%\n', sep='')
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
#' random-effect coefficients. CDR formulas are shown in user, normalized, and
#' effective forms. Set `all.coefficients=TRUE` to request the expanded table.
#'
#' @param object A fitted `cdrgam` model.
#' @param dispersion Optional known dispersion.
#' @param freq,re.test Controls passed to [mgcv::summary.gam()] by native fits.
#'   Custom backends retain these arguments for interface compatibility.
#' @param all.coefficients Include every expanded coefficient.
#' @param ... Additional arguments passed to native [mgcv::summary.gam()].
#' @return A `summary.gam`-style object with CDR formula metadata. Native fits
#'   retain the complete `summary.gam` result. Custom backends provide its
#'   commonly used coefficient, smooth-term, fit-statistic, and formula fields.
#' @name summary.cdrgam
NULL

#' @export
print.summary.cdrgam <- function(x, ...) .print_summary_cdrgam(x, ...)
