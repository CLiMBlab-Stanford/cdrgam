#' Describe estimation procedures for methods reporting
#'
#' `methods_facts()` resolves the smoothing-parameter criterion that
#' [cdrgam.fit()] uses for a configuration and returns short descriptions of
#' the estimation procedure and the significance tests, each with the keys of
#' the references that support it. Report writers can use these facts instead
#' of restating the defaults and backend behavior of this package.
#'
#' The criterion follows the fitting code: with `method = NULL`, the native
#' `mgcv` backend uses the `mgcv` default (`"fREML"` for `bam()`, `"REML"` for
#' `gam()` with an extended or general family, and `"GCV.Cp"` otherwise), and
#' the custom backends use REML.
#'
#' @param backend,engine,method As in [cdrgam.fit()]. `method = NULL` selects
#'   the default criterion.
#' @param family A family object or family name, as in [cdrgam.fit()].
#' @return A list with elements `backend`, `engine` (`NULL` unless the backend
#'   is `"mgcv"`), `method` (the resolved criterion), `criterion` and
#'   `estimator` (sentence fragments), `citations` (a list of reference keys
#'   with elements `criterion` and `estimator`), and `term_tests`, a data frame with the columns `method` (the
#'   test labels used by [summary()] and [test_effect()]), `description`, and `citation`. The
#'   reference keys name entries in the bibliography that ships with
#'   `cdrgam.cli`.
#' @export
methods_facts <- function(
        backend=c('mgcv', 'block', 'sparse'),
        engine=c('bam', 'gam'),
        method=NULL,
        family=stats::gaussian()
) {
    backend <- match.arg(backend)
    engine <- match.arg(engine)
    family <- .as_family(family)
    special <- inherits(family, c('extended.family', 'general.family'))
    allowed <- switch(backend,
        mgcv=c('REML', 'fREML', 'ML', 'GCV.Cp', 'GACV.Cp', 'NCV', 'QNCV'),
        block=c('REML', 'fREML'),
        sparse=c('REML', 'fREML', 'GCV.Cp', 'QNCV')
    )
    if (!is.null(method) && (length(method) != 1L || !(method %in% allowed))) {
        stop('method must be one of ', paste(allowed, collapse=', '),
            ' for backend="', backend, '"')
    }
    resolved <- if (!is.null(method)) {
        if (!identical(backend, 'mgcv') && identical(method, 'fREML')) 'REML' else method
    } else if (identical(backend, 'mgcv')) {
        if (identical(engine, 'bam')) 'fREML' else if (special) 'REML' else 'GCV.Cp'
    } else 'REML'
    criterion <- switch(resolved,
        REML='restricted maximum likelihood',
        fREML='fast restricted maximum likelihood',
        ML='marginal likelihood',
        GCV.Cp=if (identical(family$family, 'gaussian')) {
            'generalized cross-validation'
        } else 'the unbiased risk estimator',
        GACV.Cp='generalized approximate cross-validation',
        NCV='neighbourhood cross-validation',
        QNCV='quasi-Newton neighbourhood cross-validation'
    )
    criterion_citations <- switch(resolved,
        REML=, ML=if (special) 'wood2016smoothing' else 'wood2011fast',
        fREML='wood2015bam',
        'wood2017generalized'
    )
    estimator <- switch(backend,
        mgcv=paste0('with the mgcv function ', engine, '()'),
        block=paste('with a dense penalized least-squares solver implemented in',
            'cdrgam for Gaussian responses'),
        sparse=paste('with a sparse penalized likelihood solver implemented in',
            'cdrgam that exploits the sparsity of the response-level design')
    )
    software <- if (identical(backend, 'mgcv')) {
        if (identical(engine, 'bam')) 'wood2015bam' else 'wood2017generalized'
    } else character()
    list(
        backend=backend,
        engine=if (identical(backend, 'mgcv')) engine else NULL,
        method=resolved,
        criterion=criterion,
        estimator=estimator,
        citations=list(criterion=criterion_citations, estimator=software),
        term_tests=data.frame(
            method=c('Wald', 'Wood (2013)', 'Wood (2013b) random effect'),
            description=c(
                'Wald tests of parametric terms',
                'tests of smooth terms based on a rank-truncated pseudo-inverse',
                'tests of fully penalized terms, such as random effects'
            ),
            citation=c('wood2017generalized', 'wood2013pvalues', 'wood2013random'),
            stringsAsFactors=FALSE
        )
    )
}
