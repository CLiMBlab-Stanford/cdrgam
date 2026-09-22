#' Construct a supported CDR-GAM response family
#'
#' `cdrgam_family()` is the validated family registry used by the command-line
#' harness. Family objects may still be passed directly to [cdrgam()] and
#' [cdrgam.fit()]. The native `mgcv` backend supports every combination
#' returned here. The dense block backend supports `gaussian(identity)` plus
#' `binomial(logit)`, `poisson(log)`, and estimated-dispersion `Gamma(log)`.
#' The sparse backend supports the same combinations with streamed PIRLS and
#' an exact Laplace score.
#'
#' @param family A supported family name or an existing family object.
#' @param link Optional link name. It may be supplied only when `family` is a
#'   name.
#' @return A standard R family object.
#' @export
cdrgam_family <- function(
        family=c('gaussian', 'binomial', 'poisson', 'Gamma'),
        link=NULL
) {
    if (missing(family)) family <- 'gaussian'
    if (inherits(family, 'family')) {
        if (!is.null(link)) {
            stop('link cannot be supplied with an existing family object')
        }
        return(family)
    }
    if (!is.character(family) || length(family) != 1L || is.na(family)) {
        stop('family must be a supported family name or family object')
    }
    constructors <- list(
        gaussian=stats::gaussian,
        binomial=stats::binomial,
        poisson=stats::poisson,
        Gamma=stats::Gamma
    )
    constructor <- constructors[[family, exact=TRUE]]
    if (is.null(constructor)) {
        stop(
            'Unsupported family: ', family, '. Supported families are: ',
            paste(names(constructors), collapse=', ')
        )
    }
    if (is.null(link)) return(constructor())
    if (!is.character(link) || length(link) != 1L || is.na(link) ||
            !nzchar(link)) {
        stop('link must be a nonempty link name')
    }
    constructor(link=link)
}

.as_family <- function(family) {
    if (is.character(family) || inherits(family, 'family')) {
        return(cdrgam_family(family))
    }
    if (is.function(family)) family <- family()
    if (!inherits(family, 'family')) {
        stop('family must resolve to a standard R family object')
    }
    family
}
