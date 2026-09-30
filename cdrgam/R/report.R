#' Extract stable fit metadata and summary tables
#'
#' `fit_report()` returns the fitted-model information needed by artifact
#' writers and other downstream tools without requiring them to inspect the
#' internal representation of a `cdrgam` object.
#'
#' @param object A fitted `cdrgam` model.
#' @return `fit_metadata()` returns fitting metadata, formulas, distributional
#'   parameter names, and default plotting metadata. `fit_report()` adds
#'   printable summary text, parametric and smooth summary tables, and
#'   convergence diagnostics.
#' @export
fit_metadata <- function(object) {
    if (!is_cdrgam(object)) stop('object must be a fitted cdrgam model')
    rank <- object$cdrgam$rank
    fitting <- list(
        family=object$family$family,
        link=object$family$link,
        method=object$method,
        engine=object$cdrgam$engine,
        backend=object$cdrgam$backend
    )
    if (!is.null(rank)) fitting <- c(fitting, list(
        rank_action=rank$action,
        rank_tolerance=rank$tolerance,
        rank_penalty=rank$regularization
    ))
    if (inherits(object, 'cdrgam_sparse')) {
        fitting$sparse_control <- object$sparse$control
    }
    formula_text <- function(value) {
        if (inherits(value, 'formula')) {
            return(paste(deparse(value), collapse=' '))
        }
        if (is.list(value)) return(lapply(value, formula_text))
        as.character(value)
    }
    terms <- object$cdrgam$terms
    population <- which(vapply(terms, function(term) {
        !length(term$group_levels)
    }, logical(1)))
    list(
        fitting=fitting,
        formulas=lapply(
            object$cdrgam$formula[c('user', 'normalized', 'effective')],
            formula_text
        ),
        distributional=isTRUE(object$cdrgam$distributional),
        parameter_names=object$cdrgam$parameter_names,
        plotting=list(
            term_count=length(terms),
            population_terms=population
        )
    )
}

#' @rdname fit_report
#' @export
fit_report <- function(object) {
    metadata <- fit_metadata(object)
    model_summary <- summary(object)
    coefficients <- if (is.null(model_summary$p.table)) {
        data.frame()
    } else {
        data.frame(
            term=rownames(model_summary$p.table),
            model_summary$p.table,
            row.names=NULL,
            check.names=FALSE
        )
    }
    smooths <- if (is.null(model_summary$s.table)) {
        data.frame()
    } else {
        test_status <- model_summary$s.test
        if (is.null(test_status) ||
                length(test_status) != nrow(model_summary$s.table)) {
            test_status <- ifelse(
                is.finite(model_summary$s.table[, ncol(model_summary$s.table)]),
                'approximate',
                'not computed'
            )
        }
        data.frame(
            term=rownames(model_summary$s.table),
            model_summary$s.table,
            test_status=unname(test_status),
            row.names=NULL,
            check.names=FALSE
        )
    }
    c(metadata, list(
        summary_text=utils::capture.output(print(model_summary)),
        coefficients=coefficients,
        smooths=smooths,
        diagnostics=fit_diagnostics(object)
    ))
}
