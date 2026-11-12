.cdrgam_mgcv_term_slug <- function(name) {
    if (identical(name, 'irf(1)')) return('irf_rate')
    name <- gsub('|', '_by_', name, fixed=TRUE)
    name <- gsub('[^[:alnum:]]+', '_', name)
    name <- gsub('^_+|_+$', '', name)
    if (!nzchar(name) || identical(name, '1')) name <- 'rate'
    paste0('irf_', name)
}

.cdrgam_mgcv_term_names <- function(terms, prefix=NULL) {
    labels <- names(terms)
    if (is.null(labels) || any(!nzchar(labels))) {
        labels <- vapply(terms, function(term) term$name, character(1))
    }
    labels <- vapply(labels, .cdrgam_mgcv_term_slug, character(1))
    if (!is.null(prefix)) labels <- paste(prefix, labels, sep='_')
    make.unique(labels, sep='_')
}

.cdrgam_mgcv_parameter_formula <- function(
        ordinary,
        term_names,
        response=TRUE,
        formula_environment=NULL
) {
    if (is.null(formula_environment)) {
        formula_environment <- environment(ordinary)
    }
    rhs <- paste(deparse(ordinary[[length(ordinary)]]), collapse='')
    rhs <- paste(c(rhs, term_names), collapse=' + ')
    text <- if (isTRUE(response)) {
        paste(
            paste(deparse(ordinary[[2L]]), collapse=''),
            '~', rhs
        )
    } else paste('~', rhs)
    stats::as.formula(text, env=formula_environment)
}

.cdrgam_mgcv_formula <- function(design) {
    if (inherits(design, 'cdrgam_distributional_design')) {
        parameters <- names(design$parameters)
        output <- lapply(seq_along(parameters), function(index) {
            parameter <- parameters[[index]]
            parameter_design <- design$parameters[[parameter]]
            .cdrgam_mgcv_parameter_formula(
                parameter_design$ordinary_formula,
                .cdrgam_mgcv_term_names(
                    parameter_design$terms,
                    prefix=parameter
                ),
                response=index == 1L
            )
        })
        names(output) <- parameters
        return(output)
    }
    .cdrgam_mgcv_parameter_formula(
        design$ordinary_formula,
        .cdrgam_mgcv_term_names(design$terms)
    )
}

.cdrgam_record_mgcv_formula <- function(fit, formula) {
    computational <- fit$cdrgam$formula$mgcv
    fit$cdrgam$formula$mgcv <- formula
    if (!is.null(computational)) {
        fit$cdrgam$formula$computational <- computational
    }
    fit
}

.cdrgam_private_mgcv_formula <- function(value) {
    if (is.list(value) && !inherits(value, 'formula')) {
        return(any(vapply(value, .cdrgam_private_mgcv_formula, logical(1))))
    }
    inherits(value, 'formula') && grepl(
        'cdr_term_|bs *= *["\']cdr["\']',
        paste(deparse(value), collapse='')
    )
}

.cdrgam_mgcv_formula_from_fit <- function(fit) {
    stored <- fit$cdrgam$formula$mgcv
    if (!is.null(stored) && !.cdrgam_private_mgcv_formula(stored)) {
        return(stored)
    }
    effective <- fit$cdrgam$formula$effective
    if (is.null(effective)) effective <- fit$cdrgam$formula$normalized
    if (is.null(effective)) effective <- fit$cdrgam$formula$user
    metadata <- fit$cdrgam$terms
    labels <- fit$cdrgam$term_labels
    for (index in seq_along(metadata)) {
        name <- metadata[[index]]$name
        if (is.null(name) || !length(name) || !nzchar(name[[1L]])) {
            metadata[[index]]$name <- if (index <= length(labels)) {
                labels[[index]]
            } else {
                paste0('term_', index)
            }
        }
    }
    if (is.list(effective) && !inherits(effective, 'formula')) {
        parameters <- names(effective)
        output <- lapply(seq_along(parameters), function(index) {
            parameter <- parameters[[index]]
            selected <- metadata[vapply(metadata, function(term) {
                identical(term$parameter, parameter)
            }, logical(1))]
            term_names <- vapply(selected, `[[`, character(1), 'name')
            ordinary <- .parse_cdr_formula(effective[[parameter]])$ordinary_formula
            .cdrgam_mgcv_parameter_formula(
                ordinary,
                .cdrgam_mgcv_term_names(
                    stats::setNames(selected, term_names),
                    prefix=parameter
                ),
                response=index == 1L
            )
        })
        names(output) <- parameters
        return(output)
    }
    term_names <- vapply(metadata, `[[`, character(1), 'name')
    ordinary <- .parse_cdr_formula(effective)$ordinary_formula
    .cdrgam_mgcv_parameter_formula(
        ordinary,
        .cdrgam_mgcv_term_names(stats::setNames(metadata, term_names))
    )
}

.cdrgam_mgcv_penalty <- function(term) {
    output <- lapply(term$S, function(penalty) as.matrix(penalty))
    output$rank <- as.integer(term$rank)
    output
}

.cdrgam_mgcv_export_parameter <- function(
        design,
        prefix=NULL,
        response=TRUE
) {
    terms <- lapply(design$terms, .materialize_cdr_term, sparse=FALSE)
    term_names <- .cdrgam_mgcv_term_names(terms, prefix=prefix)
    data <- as.list(design$responses)
    penalties <- list()
    for (index in seq_along(terms)) {
        term <- terms[[index]]
        name <- term_names[[index]]
        matrix <- as.matrix(term$X)
        colnames(matrix) <- paste0(name, '.', seq_len(ncol(matrix)))
        data[[name]] <- matrix
        penalties[[name]] <- .cdrgam_mgcv_penalty(term)
    }
    list(
        formula=.cdrgam_mgcv_parameter_formula(
            design$ordinary_formula,
            term_names,
            response=response
        ),
        data=data,
        paraPen=penalties,
        term_names=term_names
    )
}

.cdrgam_mgcv_prefit_penalties <- function(prefit, exported) {
    parameter_starts <- attr(prefit$nsdf, 'pstart')
    if (is.null(parameter_starts)) {
        parameter_starts <- cumsum(c(1L, head(prefit$nsdf, -1L)))
    }
    added_S <- list()
    added_off <- numeric()
    added_rank <- numeric()
    added_names <- character()
    for (parameter_index in seq_along(exported)) {
        parameter <- exported[[parameter_index]]
        labels <- attr(prefit$pterms[[parameter_index]], 'term.labels')
        assignment <- prefit$assign[[parameter_index]]
        for (term_name in parameter$term_names) {
            term_index <- match(term_name, labels)
            if (is.na(term_index)) {
                stop('Could not locate exported mgcv term ', term_name)
            }
            local <- which(assignment == term_index)
            penalty <- parameter$paraPen[[term_name]]
            ranks <- penalty$rank
            penalty$rank <- NULL
            if (!length(local) || any(vapply(
                    penalty,
                    function(value) nrow(value) != length(local) ||
                        ncol(value) != length(local),
                    logical(1)
                ))) {
                stop('Exported mgcv penalty dimensions do not match ', term_name)
            }
            count <- length(penalty)
            names <- if (count == 1L) term_name else
                paste0(term_name, seq_len(count))
            added_S <- c(added_S, penalty)
            added_off <- c(
                added_off,
                rep.int(parameter_starts[[parameter_index]] + min(local) - 1L,
                    count)
            )
            added_rank <- c(added_rank, ranks)
            added_names <- c(added_names, names)
        }
    }
    existing_count <- length(prefit$S)
    added_count <- length(added_S)
    prefit$S <- c(prefit$S, added_S)
    prefit$off <- c(prefit$off, added_off)
    prefit$rank <- c(prefit$rank, added_rank)
    prefit$sp <- stats::setNames(
        c(prefit$sp, rep.int(-1, added_count)),
        c(names(prefit$sp), added_names)
    )
    prefit$lsp0 <- stats::setNames(
        c(prefit$lsp0, rep.int(0, added_count)),
        c(names(prefit$lsp0), added_names)
    )
    if (!is.null(prefit$L)) {
        prefit$L <- rbind(
            cbind(prefit$L, matrix(0, nrow(prefit$L), added_count)),
            cbind(matrix(0, added_count, ncol(prefit$L)), diag(added_count))
        )
    }
    prefit$n.paraPen <- added_count
    prefit$pP <- list(
        S=added_S,
        off=added_off,
        sp=stats::setNames(rep.int(-1, added_count), added_names),
        L=diag(added_count),
        rank=added_rank,
        full.sp.names=added_names
    )
    if (existing_count + added_count != length(prefit$S)) {
        stop('Internal mgcv penalty export failed')
    }
    prefit
}

#' Export an exact native mgcv fitting specification
#'
#' Convert a prepared CDR-GAM design into arguments accepted by
#' [mgcv::gam()]. Impulse-response terms become meaningfully named matrix terms,
#' and native penalty interfaces carry their exact penalties. Scalar-response
#' specifications use `paraPen`. Distributional specifications use `mgcv`'s
#' prefit `G` interface because `mgcv::gam()` does not accept `paraPen` with a
#' formula list. The export does not use CDR-GAM's private compressed smooth
#' constructor.
#'
#' The returned data materialize grouped impulse-response designs. This is an
#' interoperability and validation interface, not the memory-efficient fitting
#' path used by the sparse backend.
#'
#' @param x A `cdrgam_design` returned by [prepare_cdrgam()]. Boundary-reduced
#'   fitted models are not supported.
#' @param family A family accepted by [mgcv::gam()]. Distributional designs
#'   require their corresponding general family.
#' @param method Optional smoothing-parameter estimation method.
#' @param ... Additional arguments to include in the `mgcv::gam()` call.
#' @return A named argument list suitable for `do.call(mgcv::gam, spec)`.
#' @export
as_mgcv_spec <- function(x, family=stats::gaussian(), method=NULL, ...) {
    if (!inherits(x, 'cdrgam_design')) {
        stop('x must be a cdrgam_design returned by prepare_cdrgam()')
    }
    family <- .as_family(family)
    dots <- list(...)
    if (inherits(x, 'cdrgam_distributional_design')) {
        parameters <- names(x$parameters)
        exported <- lapply(seq_along(parameters), function(index) {
            parameter <- parameters[[index]]
            .cdrgam_mgcv_export_parameter(
                x$parameters[[parameter]],
                prefix=parameter,
                response=index == 1L
            )
        })
        names(exported) <- parameters
        formulas <- lapply(exported, `[[`, 'formula')
        names(formulas) <- parameters
        data <- as.list(x$responses)
        penalties <- list()
        for (parameter in parameters) {
            data[names(exported[[parameter]]$data)] <-
                exported[[parameter]]$data
            penalties <- c(penalties, exported[[parameter]]$paraPen)
        }
        setup_arguments <- c(list(
            formula=unname(formulas),
            family=family,
            data=data,
            fit=FALSE
        ), dots)
        if (!is.null(method)) setup_arguments$method <- method
        prefit <- do.call(mgcv::gam, setup_arguments)
        prefit <- .cdrgam_mgcv_prefit_penalties(prefit, exported)
        final_dots <- dots[setdiff(
            names(dots),
            c('formula', 'family', 'data', 'fit', 'paraPen', 'G', 'method')
        )]
        arguments <- c(list(G=prefit), final_dots)
    } else {
        exported <- .cdrgam_mgcv_export_parameter(x)
        arguments <- c(list(
            formula=exported$formula,
            family=family,
            data=exported$data,
            paraPen=exported$paraPen
        ), dots)
    }
    if (!is.null(method)) arguments$method <- method
    arguments
}
