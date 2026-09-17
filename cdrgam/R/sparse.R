.sparse_global_penalty <- function(local, indices, dimension) {
    # Symmetric Matrix classes store only one triangle. Coerce to a general
    # sparse matrix before extracting triplets so the embedded global penalty
    # contains both off-diagonal halves.
    triplet <- methods::as(
        methods::as(Matrix::Matrix(local, sparse=TRUE), 'generalMatrix'),
        'TsparseMatrix'
    )
    embedded <- Matrix::sparseMatrix(
        i=indices[methods::slot(triplet, 'i') + 1L],
        j=indices[methods::slot(triplet, 'j') + 1L],
        x=methods::slot(triplet, 'x'),
        dims=c(dimension, dimension)
    )
    Matrix::forceSymmetric(embedded, uplo='U')
}

.sparse_grouped_component <- function(term) {
    structure(list(
        base=term$X,
        group_index=term$group_index,
        group_count=length(term$group_levels),
        base_dimension=term$base_dimension
    ), class='cdrgam_sparse_grouped_component')
}

.sparse_component_nrow <- function(component) {
    if (inherits(component, 'cdrgam_sparse_grouped_component')) {
        nrow(component$base)
    } else {
        nrow(component)
    }
}

.sparse_component_ncol <- function(component) {
    if (inherits(component, 'cdrgam_sparse_grouped_component')) {
        component$group_count * component$base_dimension
    } else {
        ncol(component)
    }
}

.sparse_component_nnzero <- function(component) {
    if (inherits(component, 'cdrgam_sparse_grouped_component')) {
        sum(component$base != 0)
    } else {
        Matrix::nnzero(component)
    }
}

.sparse_component_rows <- function(component, rows) {
    if (!inherits(component, 'cdrgam_sparse_grouped_component')) {
        return(component[rows, , drop=FALSE])
    }
    base <- component$base[rows, , drop=FALSE]
    k <- component$base_dimension
    output <- Matrix::sparseMatrix(
        i=rep(seq_along(rows), each=k),
        j=rep(
            (component$group_index[rows] - 1L) * k,
            each=k
        ) + rep.int(seq_len(k), times=length(rows)),
        x=as.vector(t(base)),
        dims=c(length(rows), component$group_count * k)
    )
    Matrix::drop0(output)
}

.sparse_component_fitted <- function(component, coefficients) {
    if (!inherits(component, 'cdrgam_sparse_grouped_component')) {
        return(drop(component %*% coefficients))
    }
    coefficient_matrix <- matrix(
        coefficients,
        nrow=component$base_dimension,
        ncol=component$group_count
    )
    rowSums(component$base * t(coefficient_matrix[
        , component$group_index, drop=FALSE
    ]))
}

.cdr_factor_solve <- function(factor, rhs) {
    if (!inherits(factor, 'cdrgam_schur_factor')) {
        return(Matrix::solve(factor, rhs))
    }
    rhs <- as.matrix(rhs)
    storage.mode(rhs) <- 'double'
    .Call(
        'cdrgam_schur_solve',
        as.integer(factor$core),
        lapply(factor$blocks, as.integer),
        factor$cross,
        factor$block_cholesky,
        factor$core_cholesky,
        rhs,
        PACKAGE='cdrgam'
    )
}

.cdr_factor_logdet <- function(factor) {
    if (!inherits(factor, 'cdrgam_schur_factor')) {
        return(as.numeric(Matrix::determinant(
            factor,
            logarithm=TRUE,
            sqrt=FALSE
        )$modulus))
    }
    factor$core_logdet + sum(vapply(
        factor$block_cholesky,
        function(cholesky) 2 * sum(log(diag(cholesky))),
        numeric(1)
    ))
}

.cdr_factor_nonzeros <- function(factor) {
    if (!inherits(factor, 'cdrgam_schur_factor')) {
        if (inherits(factor, 'CHMfactor')) {
            return(Matrix::nnzero(Matrix::expand(factor)$L))
        }
        return(sum(methods::slot(factor, 'x') != 0))
    }
    sum(factor$core_cholesky != 0) +
        sum(vapply(
            factor$block_cholesky,
            function(cholesky) sum(cholesky != 0),
            numeric(1)
        ))
}

.cdr_factor_condition_indicator <- function(factor) {
    diagonals <- if (inherits(factor, 'cdrgam_schur_factor')) {
        c(
            diag(factor$core_cholesky),
            unlist(lapply(factor$block_cholesky, diag), use.names=FALSE)
        )
    } else {
        Matrix::diag(Matrix::expand(factor)$L)
    }
    if (!length(diagonals) || any(!is.finite(diagonals))) return(NA_real_)
    (min(abs(diagonals)) / max(abs(diagonals)))^2
}

.factor_schur_system <- function(system, layout, supernodal) {
    core <- layout$core
    blocks <- layout$blocks
    factor <- tryCatch(
        .Call(
            'cdrgam_schur_factor_sparse',
            system,
            as.integer(core),
            lapply(blocks, as.integer),
            PACKAGE='cdrgam'
        ),
        error=function(e) NULL
    )
    if (is.null(factor)) return(NULL)
    core_cholesky <- factor[[1L]]
    block_cholesky <- factor[[2L]]
    cross <- factor[[3L]]
    out <- list(
        core=core,
        blocks=blocks,
        cross=cross,
        block_cholesky=block_cholesky,
        core_cholesky=core_cholesky,
        core_logdet=2 * sum(log(diag(core_cholesky))),
        dimension=nrow(system)
    )
    class(out) <- 'cdrgam_schur_factor'
    out
}

.sparse_penalty_logdet <- function(blocks, sp) {
    value <- 0
    rank <- 0L
    for (block in blocks) {
        if (isTRUE(block$identity)) {
            value <- value + block$dimension *
                log(sp[[block$sp_index[[1L]]]])
            rank <- rank + block$dimension
            next
        }
        penalty <- matrix(0, nrow=block$dimension, ncol=block$dimension)
        for (j in seq_along(block$sp_index)) {
            penalty <- penalty +
                sp[[block$sp_index[[j]]]] * block$S[[j]]
        }
        values <- eigen(penalty, symmetric=TRUE, only.values=TRUE)$values
        tolerance <- max(1, max(abs(values))) * .Machine$double.eps *
            max(100, block$dimension)
        positive <- values[values > tolerance]
        value <- value + block$repetitions *
            if (length(positive)) sum(log(positive)) else 0
        rank <- rank + block$repetitions * length(positive)
    }
    list(value=value, rank=rank)
}

.sparse_matrix_keys <- function(row, column, dimension) {
    # Matrix indices are 32-bit integers, but their column-major linear key
    # can exceed .Machine$integer.max for quite ordinary mixed models. IEEE
    # doubles represent every integer exactly up to 2^53, far beyond Matrix's
    # current dimensional limits.
    as.double(row) + as.double(dimension) * as.double(column)
}

.sparse_penalty_logdet_score <- function(blocks, sp, penalty_count) {
    score <- numeric(penalty_count)
    for (block in blocks) {
        if (isTRUE(block$identity)) {
            index <- block$sp_index[[1L]]
            score[[index]] <- score[[index]] + block$dimension
            next
        }
        penalty <- matrix(0, nrow=block$dimension, ncol=block$dimension)
        for (j in seq_along(block$sp_index)) {
            penalty <- penalty +
                sp[[block$sp_index[[j]]]] * block$S[[j]]
        }
        decomposition <- eigen(penalty, symmetric=TRUE)
        tolerance <- max(1, max(abs(decomposition$values))) *
            .Machine$double.eps * max(100, block$dimension)
        positive <- decomposition$values > tolerance
        inverse <- if (any(positive)) {
            vectors <- decomposition$vectors[, positive, drop=FALSE]
            vectors %*% (
                (1 / decomposition$values[positive]) * t(vectors)
            )
        } else {
            matrix(0, block$dimension, block$dimension)
        }
        for (j in seq_along(block$sp_index)) {
            index <- block$sp_index[[j]]
            derivative <- sp[[index]] * block$S[[j]]
            score[[index]] <- score[[index]] + block$repetitions *
                sum(inverse * derivative)
        }
    }
    score
}

.sparse_logdet_score <- function(factor, derivative, chunk_size=64L) {
    columns <- which(Matrix::colSums(abs(derivative)) != 0)
    if (!length(columns)) return(0)
    total <- 0
    for (start in seq.int(1L, length(columns), by=chunk_size)) {
        end <- min(length(columns), start + chunk_size - 1L)
        selected <- columns[start:end]
        solved <- .cdr_factor_solve(
            factor,
            derivative[, selected, drop=FALSE]
        )
        diagonal_block <- as.matrix(solved[selected, , drop=FALSE])
        total <- total + sum(diag(diagonal_block))
    }
    total
}

.sparse_penalty_supports <- function(components) {
    lapply(components, function(component) {
        which(Matrix::colSums(abs(component)) != 0)
    })
}

.sparse_trace_rhs_count <- function(supports) {
    if (!length(supports)) return(0L)
    unique_supports <- list()
    for (support in supports) {
        if (!any(vapply(
            unique_supports,
            identical,
            logical(1),
            y=support
        ))) {
            unique_supports[[length(unique_supports) + 1L]] <- support
        }
    }
    sum(lengths(unique_supports))
}

.sparse_logdet_scores <- function(
        factor,
        components,
        sp,
        supports=.sparse_penalty_supports(components),
        chunk_size=64L
) {
    count <- length(components)
    if (!count) return(numeric())
    scores <- numeric(count)
    groups <- list()
    for (i in seq_len(count)) {
        match_index <- which(vapply(
            groups,
            function(group) identical(supports[[group[[1L]]]], supports[[i]]),
            logical(1)
        ))
        if (length(match_index)) {
            groups[[match_index[[1L]]]] <- c(groups[[match_index[[1L]]]], i)
        } else {
            groups[[length(groups) + 1L]] <- i
        }
    }

    dimension <- nrow(components[[1L]])
    for (group in groups) {
        support <- supports[[group[[1L]]]]
        if (!length(support)) next
        if (length(group) == 1L) {
            i <- group[[1L]]
            scores[[i]] <- sp[[i]] * .sparse_logdet_score(
                factor,
                components[[i]],
                chunk_size=chunk_size
            )
            next
        }
        for (start in seq.int(1L, length(support), by=chunk_size)) {
            end <- min(length(support), start + chunk_size - 1L)
            selected <- support[start:end]
            selector <- Matrix::sparseMatrix(
                i=selected,
                j=seq_along(selected),
                x=1,
                dims=c(dimension, length(selected))
            )
            solved <- .cdr_factor_solve(factor, selector)
            inverse_columns <- as.matrix(solved[support, , drop=FALSE])
            for (i in group) {
                derivative_columns <- as.matrix(
                    components[[i]][support, selected, drop=FALSE]
                )
                scores[[i]] <- scores[[i]] + sp[[i]] *
                    sum(inverse_columns * derivative_columns)
            }
        }
    }
    scores
}

.split_sparse_random_effects <- function(formula, data) {
    terms_object <- stats::terms(formula, keep.order=TRUE)
    labels <- attr(terms_object, 'term.labels')
    parsed <- lapply(labels, function(label) {
        tryCatch(str2lang(label), error=function(e) NULL)
    })
    is_simple_re <- vapply(parsed, function(call) {
        if (!is.call(call) || !identical(as.character(call[[1L]]), 's')) {
            return(FALSE)
        }
        arguments <- as.list(call)[-1L]
        argument_names <- names(arguments)
        bs_index <- which(argument_names == 'bs')
        if (length(bs_index) != 1L ||
                !identical(eval(arguments[[bs_index]]), 're')) {
            return(FALSE)
        }
        unnamed <- which(is.na(argument_names) | argument_names == '')
        all(vapply(arguments[unnamed], is.symbol, logical(1))) &&
            all((argument_names[-unnamed] %in% 'bs'))
    }, logical(1))
    indices <- which(is_simple_re)
    random_effects <- lapply(indices, function(index) {
        call <- parsed[[index]]
        arguments <- as.list(call)[-1L]
        argument_names <- names(arguments)
        unnamed <- which(is.na(argument_names) | argument_names == '')
        variables <- vapply(arguments[unnamed], as.character, character(1))
        missing <- setdiff(variables, names(data))
        if (length(missing)) {
            stop('Random-effect columns not found: ', paste(missing, collapse=', '))
        }
        interaction_formula <- stats::as.formula(paste(
            '~', paste(variables, collapse=':'), '- 1'
        ))
        matrix <- Matrix::sparse.model.matrix(interaction_formula, data=data)
        list(
            label=paste0('s(', paste(variables, collapse=','), ')'),
            variables=variables,
            X=matrix
        )
    })
    remaining <- setdiff(seq_along(labels), indices)
    reduced <- if (!length(indices)) {
        formula
    } else if (!length(remaining)) {
        response_text <- paste(deparse(formula[[2L]]), collapse='')
        stats::as.formula(paste(
            response_text,
            '~',
            if (attr(terms_object, 'intercept')) '1' else '0'
        ))
    } else {
        stats::formula(stats::drop.terms(
            terms_object,
            dropx=indices,
            keep.response=TRUE
        ))
    }
    environment(reduced) <- environment(formula)
    list(formula=reduced, random_effects=random_effects)
}

.ordinary_sparse_setup <- function(design, family, ...) {
    split <- .split_sparse_random_effects(
        design$ordinary_formula,
        design$responses
    )
    formula <- split$formula
    formula_env <- new.env(parent=environment(formula))
    formula_env$s <- mgcv::s
    formula_env$te <- mgcv::te
    formula_env$ti <- mgcv::ti
    formula_env$t2 <- mgcv::t2
    environment(formula) <- formula_env
    args <- c(list(
        formula=formula,
        family=family,
        data=as.list(design$responses),
        method='REML',
        fit=FALSE
    ), list(...))
    list(
        setup=do.call(mgcv::gam, args),
        random_effects=split$random_effects
    )
}

.sparse_sp_names <- function(label, count) {
    if (count == 1L) label else paste0(label, seq_len(count))
}

.sparse_expanded_formula <- function(design) {
    response_text <- paste(deparse(design$ordinary_formula[[2L]]), collapse='')
    base_rhs <- paste(deparse(design$ordinary_formula[[3L]]), collapse='')
    cdr_terms <- vapply(seq_along(design$terms), function(i) {
        term <- design$terms[[i]]
        dimension <- if (is.null(term$expanded_dimension)) {
            ncol(term$X)
        } else {
            term$expanded_dimension
        }
        paste0(
            's(cdr_term_', i, ', bs="cdr", k=', dimension,
            ', xt=cdr_spec_', i, ')'
        )
    }, character(1))
    stats::as.formula(paste(
        response_text,
        '~',
        paste(c(base_rhs, cdr_terms), collapse=' + ')
    ))
}

.fit_sparse_gaussian <- function(
        design,
        family,
        method,
        checkpoint=NULL,
        trace=FALSE,
        sparse_control=list(),
        rank_action='error',
        rank_tol=NULL,
        rank_penalty=NULL,
        ...
) {
    family <- .as_family(family)
    if (!identical(family$family, 'gaussian') ||
            !identical(family$link, 'identity')) {
        stop('The sparse backend currently supports only gaussian(identity)')
    }
    if (!is.null(method) && !(method %in% c('REML', 'fREML'))) {
        stop('The sparse backend currently supports only REML')
    }
    if (!is.list(sparse_control)) {
        stop('sparse_control must be a named list')
    }
    reporter <- .new_solver_reporter(trace, 'sparse')
    reporter$phase('model setup')
    unknown_control <- setdiff(
        names(sparse_control),
        c(
            'gradient', 'supernodal', 'trace_chunk_size', 'schur',
            'crossprod_chunk_size', 'restarts', 'gradient_probes',
            'gradient_cores', 'finite_difference_step', 'hessian',
            'hessian_step', 'outer_optimizer', 'optimizer_maxit',
            'optimizer_gradient_tolerance', 'optimizer_trust_radius'
        )
    )
    if (length(unknown_control)) {
        stop('Unknown sparse_control entries: ', paste(unknown_control, collapse=', '))
    }
    gradient_method <- if (is.null(sparse_control$gradient)) {
        'auto'
    } else {
        match.arg(
            sparse_control$gradient,
            c('auto', 'finite', 'exact', 'stochastic', 'hybrid')
        )
    }
    gradient_probes <- if (is.null(sparse_control$gradient_probes)) {
        12L
    } else {
        value <- sparse_control$gradient_probes
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value < 1 || value != as.integer(value)) {
            stop('sparse_control$gradient_probes must be a positive integer')
        }
        as.integer(value)
    }
    gradient_cores <- if (is.null(sparse_control$gradient_cores)) {
        1L
    } else {
        value <- sparse_control$gradient_cores
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value < 1 || value != as.integer(value)) {
            stop('sparse_control$gradient_cores must be a positive integer')
        }
        as.integer(value)
    }
    if (.Platform$OS.type == 'windows' && gradient_cores > 1L) {
        stop('Parallel finite gradients are not supported on Windows')
    }
    finite_difference_step <- if (
        is.null(sparse_control$finite_difference_step)
    ) {
        1e-3
    } else {
        value <- sparse_control$finite_difference_step
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value <= 0) {
            stop('sparse_control$finite_difference_step must be positive')
        }
        value
    }
    hessian_method <- if (is.null(sparse_control$hessian)) {
        'profiled'
    } else {
        match.arg(sparse_control$hessian, c('optimhess', 'profiled'))
    }
    hessian_step <- if (is.null(sparse_control$hessian_step)) {
        1e-2
    } else {
        value <- sparse_control$hessian_step
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value <= 0) {
            stop('sparse_control$hessian_step must be positive')
        }
        value
    }
    outer_optimizer <- if (is.null(sparse_control$outer_optimizer)) {
        'lbfgsb'
    } else {
        match.arg(sparse_control$outer_optimizer, c('lbfgsb', 'bfgs_trust'))
    }
    optimizer_maxit <- if (is.null(sparse_control$optimizer_maxit)) {
        100L
    } else {
        value <- sparse_control$optimizer_maxit
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value < 1 || value != as.integer(value)) {
            stop('sparse_control$optimizer_maxit must be a positive integer')
        }
        as.integer(value)
    }
    optimizer_gradient_tolerance <- if (is.null(
            sparse_control$optimizer_gradient_tolerance
    )) {
        1e-4
    } else {
        value <- sparse_control$optimizer_gradient_tolerance
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value <= 0) {
            stop(
                'sparse_control$optimizer_gradient_tolerance must be positive'
            )
        }
        value
    }
    optimizer_trust_radius <- if (is.null(
            sparse_control$optimizer_trust_radius
    )) {
        2
    } else {
        value <- sparse_control$optimizer_trust_radius
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value <= 0) {
            stop('sparse_control$optimizer_trust_radius must be positive')
        }
        value
    }
    if (identical(outer_optimizer, 'bfgs_trust') &&
            !identical(gradient_method, 'exact')) {
        stop(
            'sparse_control$outer_optimizer="bfgs_trust" requires ',
            'gradient="exact"'
        )
    }
    supernodal <- if (is.null(sparse_control$supernodal)) {
        NULL
    } else {
        value <- sparse_control$supernodal
        if (length(value) != 1L || !is.logical(value)) {
            stop('sparse_control$supernodal must be TRUE, FALSE, or NA')
        }
        value
    }
    trace_chunk_size <- if (is.null(sparse_control$trace_chunk_size)) {
        256L
    } else {
        value <- sparse_control$trace_chunk_size
        if (length(value) != 1L || !is.numeric(value) ||
                !is.finite(value) || value < 1) {
            stop('sparse_control$trace_chunk_size must be a positive integer')
        }
        as.integer(value)
    }
    crossprod_chunk_size <- if (is.null(sparse_control$crossprod_chunk_size)) {
        10000L
    } else {
        value <- sparse_control$crossprod_chunk_size
        if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
                value < 1) {
            stop('sparse_control$crossprod_chunk_size must be positive')
        }
        if (is.infinite(value) || value > .Machine$integer.max) {
            .Machine$integer.max
        } else {
            as.integer(value)
        }
    }
    schur_method <- if (is.null(sparse_control$schur)) {
        'never'
    } else {
        match.arg(sparse_control$schur, c('never', 'always'))
    }
    restart_count <- if (is.null(sparse_control$restarts)) {
        0L
    } else {
        value <- sparse_control$restarts
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value < 0 || value != as.integer(value)) {
            stop('sparse_control$restarts must be a non-negative integer')
        }
        as.integer(value)
    }
    tolerance <- .rank_tolerance(rank_tol)
    penalty_strength <- .rank_penalty(rank_penalty)
    ordinary <- .ordinary_sparse_setup(design, family, ...)
    alias_resolution <- .drop_parametric_aliases(ordinary$setup, tolerance)
    setup <- alias_resolution$setup
    crossed_effects <- FALSE
    if (is.null(supernodal)) {
        # Simplicial numeric updates are substantially faster for the fill
        # pattern induced by crossed random effects; supernodal updates retain
        # a small advantage for a single grouping factor.
        ordinary_groups <- unlist(lapply(
            ordinary$random_effects,
            function(random_effect) {
                random_effect$variables[vapply(
                    design$responses[random_effect$variables],
                    is.factor,
                    logical(1)
                )]
            }
        ))
        irf_groups <- vapply(
            design$terms,
            function(term) if (is.null(term$group)) '' else term$group,
            character(1)
        )
        grouping_factors <- unique(c(
            ordinary_groups,
            irf_groups[nzchar(irf_groups)]
        ))
        crossed_effects <- length(grouping_factors) >= 2L
    }
    ordinary_dimension <- ncol(setup$X)
    matrices <- list(Matrix::Matrix(setup$X, sparse=TRUE))
    coefficient_names <- colnames(setup$X)
    smooths <- setup$smooth
    components <- list()
    blocks_by_key <- list()
    sp_names <- names(setup$sp)
    coefficient_groups <- list()
    random_prediction <- list()
    register_grouped_indices <- function(group, levels, indices_by_level) {
        if (is.null(coefficient_groups[[group]])) {
            coefficient_groups[[group]] <<- stats::setNames(
                replicate(length(levels), integer(), simplify=FALSE),
                levels
            )
        }
        for (i in seq_along(levels)) {
            coefficient_groups[[group]][[levels[[i]]]] <<- c(
                coefficient_groups[[group]][[levels[[i]]]],
                indices_by_level[[i]]
            )
        }
    }

    if (length(setup$S)) {
        for (i in seq_along(setup$S)) {
            indices <- setup$off[[i]] + seq_len(nrow(setup$S[[i]])) - 1L
            components[[length(components) + 1L]] <- list(
                indices=indices,
                S=setup$S[[i]]
            )
            key <- paste(setup$off[[i]], nrow(setup$S[[i]]), sep=':')
            if (is.null(blocks_by_key[[key]])) {
                blocks_by_key[[key]] <- list(
                    dimension=nrow(setup$S[[i]]),
                    repetitions=1L,
                    sp_index=integer(),
                    S=list()
                )
            }
            blocks_by_key[[key]]$sp_index <- c(
                blocks_by_key[[key]]$sp_index,
                length(components)
            )
            blocks_by_key[[key]]$S <- c(
                blocks_by_key[[key]]$S,
                list(setup$S[[i]])
            )
        }
    }

    coefficient_offset <- ordinary_dimension
    for (random_effect in ordinary$random_effects) {
        matrices[[length(matrices) + 1L]] <- random_effect$X
        random_dimension <- ncol(random_effect$X)
        indices <- coefficient_offset + seq_len(random_dimension)
        coefficient_names <- c(
            coefficient_names,
            colnames(random_effect$X)
        )
        components[[length(components) + 1L]] <- list(
            indices=indices,
            S=Matrix::Diagonal(random_dimension)
        )
        sp_index <- length(components)
        sp_names <- c(sp_names, random_effect$label)
        blocks_by_key[[paste0('re:', sp_index)]] <- list(
            dimension=random_dimension,
            repetitions=1L,
            sp_index=sp_index,
            S=list(),
            identity=TRUE
        )
        smooths[[length(smooths) + 1L]] <- list(
            label=random_effect$label,
            term=random_effect$variables,
            sp=stats::setNames(-1, random_effect$label),
            S=list(Matrix::Diagonal(random_dimension)),
            S.scale=1,
            rank=random_dimension,
            null.space.dim=0L,
            first.para=indices[[1L]],
            last.para=indices[[length(indices)]],
            id=NULL
        )
        random_prediction[[length(random_prediction) + 1L]] <- list(
            variables=random_effect$variables,
            column_names=colnames(random_effect$X),
            levels=stats::setNames(lapply(
                random_effect$variables,
                function(variable) {
                    value <- design$responses[[variable]]
                    if (is.factor(value)) levels(droplevels(value)) else
                        sort(unique(as.character(value)))
                }
            ), random_effect$variables),
            coefficient_index=indices
        )
        factor_variables <- random_effect$variables[vapply(
            design$responses[random_effect$variables],
            is.factor,
            logical(1)
        )]
        if (length(factor_variables) == 1L) {
            group <- factor_variables[[1L]]
            levels <- levels(droplevels(design$responses[[group]]))
            if (random_dimension == length(levels)) {
                register_grouped_indices(
                    group,
                    levels,
                    lapply(seq_along(levels), function(i) indices[[i]])
                )
            }
        }
        coefficient_offset <- coefficient_offset + random_dimension
    }

    term_metadata <- vector('list', length(design$terms))
    for (i in seq_along(design$terms)) {
        compact <- design$terms[[i]]
        base_dimension <- if (is.null(compact$base_dimension)) {
            ncol(compact$X)
        } else {
            compact$base_dimension
        }
        grouped <- !is.null(compact$group)
        if (grouped) {
            matrices[[length(matrices) + 1L]] <-
                .sparse_grouped_component(compact)
            identity <- Matrix::Diagonal(length(compact$group_levels))
            expanded_penalties <- lapply(compact$S, function(penalty) {
                Matrix::kronecker(
                    identity,
                    Matrix::Matrix(penalty, sparse=TRUE)
                )
            })
            term_dimension <- compact$expanded_dimension
        } else {
            matrices[[length(matrices) + 1L]] <-
                Matrix::Matrix(compact$X, sparse=TRUE)
            expanded_penalties <- lapply(
                compact$S,
                Matrix::Matrix,
                sparse=TRUE
            )
            term_dimension <- ncol(compact$X)
        }
        indices <- coefficient_offset + seq_len(term_dimension)
        coefficient_names <- c(
            coefficient_names,
            paste0('cdr_term_', i, '.', seq_len(term_dimension))
        )
        label <- paste0('s(cdr_term_', i, ')')
        term_sp_names <- .sparse_sp_names(label, length(compact$S))
        first_sp <- length(components) + 1L
        for (j in seq_along(expanded_penalties)) {
            components[[length(components) + 1L]] <- list(
                indices=indices,
                S=expanded_penalties[[j]]
            )
        }
        sp_indices <- seq.int(first_sp, length(components))
        sp_names <- c(sp_names, term_sp_names)
        repeats <- if (is.null(compact$group)) {
            1L
        } else {
            length(compact$group_levels)
        }
        blocks_by_key[[paste0('cdr:', i)]] <- list(
            dimension=base_dimension,
            repetitions=repeats,
            sp_index=sp_indices,
            S=compact$S
        )
        smooths[[length(smooths) + 1L]] <- list(
            label=label,
            term=paste0('cdr_term_', i),
            sp=stats::setNames(rep.int(-1, length(compact$S)), term_sp_names),
            S=compact$S,
            S.scale=compact$S.scale,
            rank=compact$rank,
            null.space.dim=compact$null.space.dim,
            first.para=indices[[1L]],
            last.para=indices[[length(indices)]],
            id=NULL
        )
        term_metadata[[i]] <- list(
            name=compact$name,
            type=compact$type,
            knots=compact$knots,
            predictor_knots=compact$predictor_knots,
            basis=compact$basis,
            transform=compact$transform,
            group=compact$group,
            group_levels=compact$group_levels,
            base_dimension=base_dimension,
            rank=compact$rank,
            null.space.dim=compact$null.space.dim,
            S.scale=compact$S.scale,
            coefficient_index=indices
        )
        if (!is.null(compact$group)) {
            register_grouped_indices(
                compact$group,
                compact$group_levels,
                lapply(seq_along(compact$group_levels), function(level) {
                    start <- (level - 1L) * base_dimension + 1L
                    indices[start:(start + base_dimension - 1L)]
                })
            )
        }
        coefficient_offset <- coefficient_offset + term_dimension
    }

    expanded_penalties <- NULL
    matrix_dimensions <- vapply(
        matrices,
        .sparse_component_ncol,
        integer(1)
    )
    matrix_offsets <- cumsum(c(0L, utils::head(matrix_dimensions, -1L)))
    dimension <- sum(matrix_dimensions)
    observation_count <- .sparse_component_nrow(matrices[[1L]])
    if (any(vapply(
        matrices,
        .sparse_component_nrow,
        integer(1)
    ) != observation_count)) {
        stop('Internal sparse design components have inconsistent row counts')
    }
    design_nonzeros <- sum(vapply(
        matrices,
        .sparse_component_nnzero,
        numeric(1)
    ))
    streamed_grouped_terms <- sum(vapply(
        matrices,
        inherits,
        logical(1),
        what='cdrgam_sparse_grouped_component'
    ))
    if (is.null(supernodal)) {
        # Simplicial updates win for small crossed systems, but supernodal
        # BLAS dominates after fill-in grows beyond this range.
        supernodal <- !crossed_effects || dimension > 750L
    }
    penalty_components <- lapply(components, function(component) {
        .sparse_global_penalty(component$S, component$indices, dimension)
    })
    names(penalty_components) <- sp_names
    penalty_supports <- .sparse_penalty_supports(penalty_components)
    exact_trace_rhs <- .sparse_trace_rhs_count(penalty_supports)
    gradient_penalty_count <- length(penalty_components)
    if (identical(gradient_method, 'auto')) {
        # Exact scores trade one factorization per finite-difference direction
        # for triangular solves over the union of penalty-supported columns.
        # Support per smoothing parameter is a substantially better predictor
        # than total coefficient dimension, especially for large models whose
        # smooths occupy small disjoint blocks.
        gradient_method <- if (
            exact_trace_rhs <= max(512L, 64L * gradient_penalty_count)
        ) 'exact' else 'finite'
    }
    blocks <- unname(blocks_by_key)
    y <- setup$y - setup$offset
    weights <- setup$w
    if (is.null(weights)) weights <- rep.int(1, length(y))
    if (any(!is.finite(weights)) || any(weights <= 0)) {
        stop('The sparse backend requires finite positive weights')
    }
    root_weights <- sqrt(weights)
    weighted_y <- y * root_weights
    crossprod_starts <- seq.int(
        1L,
        observation_count,
        by=min(crossprod_chunk_size, observation_count)
    )
    reporter$phase(
        'cross-product accumulation',
        observations=observation_count,
        coefficients=dimension,
        chunks=length(crossprod_starts)
    )
    XtX <- NULL
    Xty <- numeric(dimension)
    for (chunk_index in seq_along(crossprod_starts)) {
        start <- crossprod_starts[[chunk_index]]
        end <- min(observation_count, start + crossprod_chunk_size - 1L)
        rows <- start:end
        component_chunks <- lapply(
            matrices,
            .sparse_component_rows,
            rows=rows
        )
        design_chunk <- if (length(component_chunks) == 1L) {
            component_chunks[[1L]]
        } else {
            do.call(cbind, component_chunks)
        }
        weighted_chunk <- Matrix::Diagonal(x=root_weights[rows]) %*%
            design_chunk
        chunk_crossprod <- Matrix::crossprod(weighted_chunk)
        XtX <- if (is.null(XtX)) {
            chunk_crossprod
        } else {
            XtX + chunk_crossprod
        }
        Xty <- Xty + as.numeric(Matrix::crossprod(
            weighted_chunk,
            weighted_y[rows]
        ))
        if (reporter$level >= 3L) {
            reporter$emit(
                3L,
                'chunk complete',
                chunk=chunk_index,
                chunks=length(crossprod_starts),
                rows=end - start + 1L
            )
        }
    }
    XtX <- Matrix::forceSymmetric(XtX, uplo='U')
    penalty_count <- length(penalty_components)
    reporter$phase(
        'symbolic factorization',
        system_nonzeros=Matrix::nnzero(XtX),
        penalties=penalty_count
    )

    # Every REML evaluation has the same structural nonzero pattern. Build its
    # union once, map each penalty component into the compressed-column value
    # slot, and subsequently update only those numeric values. Absolute values
    # prevent cancellation from accidentally removing a structurally possible
    # entry from the symbolic factorization.
    system_pattern <- abs(Matrix::forceSymmetric(XtX, uplo='U'))
    for (component in penalty_components) {
        system_pattern <- system_pattern + abs(component)
    }
    system_template <- Matrix::forceSymmetric(system_pattern, uplo='U')
    template_columns <- rep.int(
        seq_len(dimension) - 1L,
        diff(methods::slot(system_template, 'p'))
    )
    template_keys <- .sparse_matrix_keys(
        methods::slot(system_template, 'i'),
        template_columns,
        dimension
    )
    map_system_values <- function(source) {
        source <- Matrix::forceSymmetric(source, uplo='U')
        triplet <- methods::as(source, 'TsparseMatrix')
        keys <- .sparse_matrix_keys(
            methods::slot(triplet, 'i'),
            methods::slot(triplet, 'j'),
            dimension
        )
        positions <- match(keys, template_keys)
        if (anyNA(positions)) {
            stop('Internal sparse-system pattern mapping failed')
        }
        list(positions=positions, values=methods::slot(triplet, 'x'))
    }
    base_mapping <- map_system_values(XtX)
    base_system_values <- numeric(length(template_keys))
    base_system_values[base_mapping$positions] <- base_mapping$values
    penalty_mappings <- lapply(penalty_components, map_system_values)
    unit_system_values <- base_system_values
    for (mapping in penalty_mappings) {
        unit_system_values[mapping$positions] <-
            unit_system_values[mapping$positions] + mapping$values
    }
    unit_system <- system_template
    methods::slot(unit_system, 'x') <- unit_system_values
    rank_resolution <- .rank_regularization(
        unit_system,
        rank_action,
        tolerance,
        penalty_strength
    )
    fixed_ridge <- rank_resolution$value
    if (fixed_ridge > 0) {
        diagonal_keys <- .sparse_matrix_keys(
            seq_len(dimension) - 1L,
            seq_len(dimension) - 1L,
            dimension
        )
        diagonal_positions <- match(diagonal_keys, template_keys)
        if (anyNA(diagonal_positions)) {
            stop('Internal sparse system is missing structural diagonal entries')
        }
        base_system_values[diagonal_positions] <-
            base_system_values[diagonal_positions] + fixed_ridge
        unit_system_values[diagonal_positions] <-
            unit_system_values[diagonal_positions] + fixed_ridge
        methods::slot(unit_system, 'x') <- unit_system_values
    }
    checkpoint_signature <- .solver_checkpoint_signature(
        backend='sparse',
        formula=design$formula,
        observation_count=observation_count,
        dimension=dimension,
        sp_names=sp_names,
        Xty=Xty,
        system=unit_system
    )
    checkpoint_signature$optimizer <- list(
        outer_optimizer=outer_optimizer,
        gradient=gradient_method,
        gradient_probes=gradient_probes,
        gradient_cores=gradient_cores,
        finite_difference_step=finite_difference_step,
        optimizer_maxit=optimizer_maxit,
        optimizer_gradient_tolerance=optimizer_gradient_tolerance,
        optimizer_trust_radius=optimizer_trust_radius,
        restarts=restart_count,
        lower=-25,
        upper=25
    )
    checkpoint_state <- .checkpoint_validate(
        .checkpoint_read(checkpoint),
        checkpoint_signature,
        checkpoint
    )
    resumed <- !is.null(checkpoint_state)
    if (is.null(checkpoint_state) || isTRUE(checkpoint_state$legacy)) {
        legacy_state <- checkpoint_state
        checkpoint_state <- .new_checkpoint_state(
            checkpoint_signature,
            'sparse',
            restart_count
        )
        if (!is.null(legacy_state$log_sp)) {
            checkpoint_state$current_log_sp <- legacy_state$log_sp
            checkpoint_state$best_log_sp <- legacy_state$log_sp
            checkpoint_state$current_criterion <- legacy_state$criterion
            checkpoint_state$best_criterion <- legacy_state$criterion
        }
    }
    checkpoint_best <- checkpoint_state$best_criterion
    if (is.null(checkpoint_best) || !is.finite(checkpoint_best)) {
        checkpoint_best <- Inf
    }
    evaluation_count <- if (is.null(checkpoint_state$evaluation_count)) {
        0L
    } else {
        as.integer(checkpoint_state$evaluation_count)
    }
    checkpoint_last_written <- evaluation_count
    checkpoint_last_written_time <- proc.time()[['elapsed']]
    last_level_one_report <- -Inf
    if (resumed) {
        reporter$emit(
            1L,
            'checkpoint resumed',
            path=checkpoint,
            stage=checkpoint_state$stage,
            evaluations=evaluation_count,
            best=checkpoint_best
        )
    }
    schur_layout <- NULL
    if (identical(schur_method, 'always') && length(coefficient_groups)) {
        candidate_order <- order(vapply(
            coefficient_groups,
            function(group) sum(lengths(group)),
            numeric(1)
        ), decreasing=TRUE)
        triplet <- methods::as(unit_system, 'TsparseMatrix')
        row_index <- methods::slot(triplet, 'i') + 1L
        column_index <- methods::slot(triplet, 'j') + 1L
        values <- methods::slot(triplet, 'x')
        for (candidate in candidate_order) {
            candidate_blocks <- unname(coefficient_groups[[candidate]])
            candidate_blocks <- candidate_blocks[lengths(candidate_blocks) > 0L]
            block_dimension <- sum(lengths(candidate_blocks))
            core <- setdiff(seq_len(dimension), unlist(candidate_blocks))
            if (length(candidate_blocks) < 2L || !length(core)) {
                next
            }
            labels <- integer(dimension)
            for (i in seq_along(candidate_blocks)) {
                labels[candidate_blocks[[i]]] <- i
            }
            left <- labels[row_index]
            right <- labels[column_index]
            invalid <- left > 0L & right > 0L & left != right & values != 0
            if (!any(invalid)) {
                schur_layout <- list(
                    group=names(coefficient_groups)[[candidate]],
                    core=core,
                    blocks=candidate_blocks
                )
                break
            }
        }
    }
    if (identical(schur_method, 'always') && is.null(schur_layout)) {
        stop('No valid response-aligned block partition was found for Schur fitting')
    }
    factor_template <- if (is.null(schur_layout)) {
        tryCatch(
            suppressWarnings(Matrix::Cholesky(
                unit_system,
                LDL=FALSE,
                perm=TRUE,
                super=supernodal
            )),
            error=function(e) NULL
        )
    } else {
        NULL
    }
    if (is.null(schur_layout) && is.null(factor_template)) .rank_error(rank_action)
    objective_cache <- new.env(hash=TRUE, parent=emptyenv())
    gradient_probe_matrix <- if (gradient_method %in% c('stochastic', 'hybrid')) {
        .deterministic_rademacher(dimension, gradient_probes)
    } else {
        NULL
    }
    numeric_update_count <- 0L
    parallel_gradient_factorizations <- 0L
    full_factorization_count <- 0L
    numeric_update_error <- NULL
    last_solution_key <- NULL
    last_solution <- NULL
    objective_seconds <- numeric()

    cache_key <- function(log_sp) {
        paste(formatC(log_sp, digits=17, format='fg'), collapse='\034')
    }

    factor_system <- function(sp) {
        system_values <- base_system_values
        for (i in seq_along(penalty_mappings)) {
            mapping <- penalty_mappings[[i]]
            system_values[mapping$positions] <-
                system_values[mapping$positions] + sp[[i]] * mapping$values
        }
        system <- system_template
        methods::slot(system, 'x') <- system_values
        if (!is.null(schur_layout)) {
            factor <- .factor_schur_system(system, schur_layout, supernodal)
            if (is.null(factor)) {
                full_factorization_count <<- full_factorization_count + 1L
                factor <- tryCatch(
                    Matrix::Cholesky(
                        system,
                        LDL=FALSE,
                        perm=TRUE,
                        super=supernodal
                    ),
                    error=function(e) NULL
                )
            } else {
                numeric_update_count <<- numeric_update_count + 1L
            }
            return(list(system=system, factor=factor))
        }
        factor <- tryCatch(
            Matrix::update(factor_template, system),
            error=function(e) {
                if (is.null(numeric_update_error)) {
                    numeric_update_error <<- conditionMessage(e)
                }
                NULL
            }
        )
        if (is.null(factor)) {
            full_factorization_count <<- full_factorization_count + 1L
            factor <- tryCatch(
                Matrix::Cholesky(
                    system,
                    LDL=FALSE,
                    perm=TRUE,
                    super=supernodal
                ),
                error=function(e) NULL
            )
        } else {
            numeric_update_count <<- numeric_update_count + 1L
        }
        list(system=system, factor=factor)
    }

    evaluate <- function(log_sp, retain=FALSE, record=TRUE) {
        if (!retain && record) evaluation_count <<- evaluation_count + 1L
        key <- cache_key(log_sp)
        if (!retain && record &&
                exists(key, envir=objective_cache, inherits=FALSE)) {
            if (reporter$level >= 3L) {
                reporter$emit(
                    3L,
                    'cached objective',
                    evaluation=evaluation_count
                )
            }
            return(get(key, envir=objective_cache, inherits=FALSE))
        }
        evaluation_started <- proc.time()[['elapsed']]
        sp <- exp(log_sp)
        factored <- factor_system(sp)
        if (is.null(factored$factor)) {
            if (!retain && record) {
                reporter$emit(
                    1L,
                    'factorization failed',
                    evaluation=evaluation_count
                )
            }
            return(if (retain) NULL else .Machine$double.xmax / 100)
        }
        penalty_det <- .sparse_penalty_logdet(blocks, sp)
        residual_df <- observation_count - (dimension - penalty_det$rank)
        if (residual_df <= 0) {
            stop('Insufficient observations for Gaussian REML estimation')
        }
        coefficients <- .cdr_factor_solve(factored$factor, Xty)
        # At the penalized least-squares solution,
        #   ||y-Xb||^2 + b'Sb = y'y - b'X'y.
        # This avoids an n-by-p matrix product on every REML evaluation.
        penalized_rss <- as.numeric(Matrix::crossprod(weighted_y) -
            Matrix::crossprod(coefficients, Xty))
        if (!is.finite(penalized_rss) || penalized_rss <= 0) {
            return(if (retain) NULL else .Machine$double.xmax / 100)
        }
        log_det_system <- .cdr_factor_logdet(factored$factor)
        criterion <- as.numeric(residual_df * log(penalized_rss / residual_df) +
            log_det_system - penalty_det$value
        )
        elapsed_objective <- proc.time()[['elapsed']] - evaluation_started
        objective_seconds <<- c(
            utils::tail(objective_seconds, 99L),
            elapsed_objective
        )
        if (!retain && record) assign(key, criterion, envir=objective_cache)
        if (!retain && record) {
            new_best <- criterion < checkpoint_best
            now <- proc.time()[['elapsed']]
            level_one_due <- evaluation_count == 1L ||
                evaluation_count %% 10L == 0L ||
                (new_best && now - last_level_one_report >= 5)
            should_report <- reporter$level >= 2L || level_one_due
            if (should_report) {
                reporter$emit(
                    if (level_one_due) 1L else 2L,
                    if (new_best) 'new best' else 'evaluation',
                    evaluation=evaluation_count,
                    criterion=format(criterion, digits=10),
                    seconds=format(elapsed_objective, digits=5),
                    sp=format(sp, digits=5)
                )
                if (level_one_due) last_level_one_report <<- now
            }
            checkpoint_state$current_log_sp <<- log_sp
            checkpoint_state$current_criterion <<- criterion
            checkpoint_state$evaluation_count <<- evaluation_count
            if (new_best) {
                checkpoint_best <<- criterion
                checkpoint_state$best_log_sp <<- log_sp
                checkpoint_state$best_criterion <<- criterion
            }
            if (!is.null(checkpoint) &&
                    (evaluation_count - checkpoint_last_written >= 10L ||
                     now - checkpoint_last_written_time >= 60)) {
                checkpoint_state$updated_at <<- as.character(Sys.time())
                .checkpoint_write(checkpoint_state, checkpoint)
                checkpoint_last_written <<- evaluation_count
                checkpoint_last_written_time <<- now
            }
        }
        retained <- c(factored, list(
            criterion=criterion,
            sp=sp,
            coefficients=as.numeric(coefficients),
            penalized_rss=penalized_rss,
            reml_df=residual_df,
            log_det_system=log_det_system,
            penalty_logdet=penalty_det$value
        ))
        last_solution_key <<- key
        last_solution <<- retained
        if (!retain) return(criterion)
        retained
    }

    retained_solution <- function(log_sp) {
        key <- cache_key(log_sp)
        if (!is.null(last_solution_key) &&
                identical(key, last_solution_key)) {
            last_solution
        } else {
            evaluate(log_sp, retain=TRUE)
        }
    }

    exact_gradient <- function(log_sp) {
        gradient_started <- proc.time()[['elapsed']]
        retained <- retained_solution(log_sp)
        if (is.null(retained)) return(rep.int(0, penalty_count))
        sp <- exp(log_sp)
        penalty_score <- .sparse_penalty_logdet_score(
            blocks,
            sp,
            penalty_count
        )
        determinant_scores <- .sparse_logdet_scores(
            retained$factor,
            penalty_components,
            sp,
            supports=penalty_supports,
            chunk_size=trace_chunk_size
        )
        score <- numeric(penalty_count)
        coefficients <- retained$coefficients
        for (i in seq_len(penalty_count)) {
            derivative <- sp[[i]] * penalty_components[[i]]
            rss_score <- as.numeric(Matrix::crossprod(
                coefficients,
                derivative %*% coefficients
            ))
            score[[i]] <- retained$reml_df * rss_score /
                retained$penalized_rss + determinant_scores[[i]] -
                penalty_score[[i]]
        }
        if (reporter$level >= 2L) {
            reporter$emit(
                2L,
                'exact gradient',
                seconds=format(
                    proc.time()[['elapsed']] - gradient_started,
                    digits=5
                ),
                maximum=format(max(abs(score)), digits=5)
            )
        }
        score
    }

    stochastic_gradient <- function(log_sp) {
        gradient_started <- proc.time()[['elapsed']]
        retained <- retained_solution(log_sp)
        if (is.null(retained)) return(rep.int(0, penalty_count))
        sp <- exp(log_sp)
        penalty_score <- .sparse_penalty_logdet_score(
            blocks,
            sp,
            penalty_count
        )
        inverse_probes <- .cdr_factor_solve(
            retained$factor,
            gradient_probe_matrix
        )
        inverse_probes <- as.matrix(inverse_probes)
        score <- numeric(penalty_count)
        coefficients <- retained$coefficients
        for (i in seq_len(penalty_count)) {
            derivative <- sp[[i]] * penalty_components[[i]]
            rss_score <- as.numeric(Matrix::crossprod(
                coefficients,
                derivative %*% coefficients
            ))
            derivative_probes <- as.matrix(
                derivative %*% gradient_probe_matrix
            )
            determinant_score <- mean(colSums(
                inverse_probes * derivative_probes
            ))
            score[[i]] <- retained$reml_df * rss_score /
                retained$penalized_rss + determinant_score -
                penalty_score[[i]]
        }
        if (reporter$level >= 2L) {
            reporter$emit(
                2L,
                'stochastic gradient',
                probes=gradient_probes,
                seconds=format(
                    proc.time()[['elapsed']] - gradient_started,
                    digits=5
                ),
                maximum=format(max(abs(score)), digits=5)
            )
        }
        score
    }

    finite_gradient <- function(log_sp) {
        gradient_started <- proc.time()[['elapsed']]
        points <- vector('list', 2L * penalty_count)
        widths <- numeric(penalty_count)
        for (i in seq_len(penalty_count)) {
            lower_point <- upper_point <- log_sp
            lower_point[[i]] <- max(
                lower_bound[[i]],
                log_sp[[i]] - finite_difference_step
            )
            upper_point[[i]] <- min(
                upper_bound[[i]],
                log_sp[[i]] + finite_difference_step
            )
            points[[2L * i - 1L]] <- lower_point
            points[[2L * i]] <- upper_point
            widths[[i]] <- upper_point[[i]] - lower_point[[i]]
        }
        worker <- function(point) evaluate(
            point,
            retain=FALSE,
            record=FALSE
        )
        values <- if (gradient_cores > 1L) {
            results <- parallel::mclapply(
                points,
                worker,
                mc.cores=min(gradient_cores, length(points)),
                mc.preschedule=TRUE,
                mc.set.seed=FALSE
            )
            valid <- vapply(
                results,
                function(value) is.numeric(value) && length(value) == 1L &&
                    is.finite(value),
                logical(1)
            )
            if (!all(valid)) {
                stop('Parallel finite-difference gradient evaluation failed')
            }
            vapply(results, as.numeric, numeric(1))
        } else {
            vapply(points, worker, numeric(1))
        }
        if (length(values) != length(points) || any(!is.finite(values))) {
            stop('Parallel finite-difference gradient evaluation failed')
        }
        evaluation_count <<- evaluation_count + length(points)
        parallel_gradient_factorizations <<-
            parallel_gradient_factorizations + length(points)
        checkpoint_state$evaluation_count <<- evaluation_count
        score <- (values[seq.int(2L, length(values), by=2L)] -
            values[seq.int(1L, length(values), by=2L)]) / widths
        reporter$emit(
            1L,
            'finite gradient',
            evaluations=length(points),
            cores=gradient_cores,
            seconds=format(
                proc.time()[['elapsed']] - gradient_started,
                digits=5
            ),
            maximum=format(max(abs(score)), digits=5)
        )
        score
    }

    canonical_initial <- rep.int(0, penalty_count)
    initial <- canonical_initial
    checkpoint_initial <- if (!is.null(checkpoint_state$current_log_sp)) {
        checkpoint_state$current_log_sp
    } else {
        checkpoint_state$best_log_sp
    }
    if (length(checkpoint_initial) == penalty_count &&
            all(is.finite(checkpoint_initial))) {
        initial <- checkpoint_initial
    }
    lower_bound <- rep.int(-25, penalty_count)
    upper_bound <- rep.int(25, penalty_count)
    optimization_arguments <- list(
        par=initial,
        fn=evaluate,
        method='L-BFGS-B',
        lower=lower_bound,
        upper=upper_bound,
        control=list(factr=1e7)
    )
    if (identical(gradient_method, 'exact')) {
        optimization_arguments$gr <- exact_gradient
    } else if (gradient_method %in% c('stochastic', 'hybrid')) {
        optimization_arguments$gr <- stochastic_gradient
    } else if (identical(gradient_method, 'finite') && gradient_cores > 1L) {
        optimization_arguments$gr <- finite_gradient
    }
    hessian_parameter_count <- if (identical(hessian_method, 'profiled')) {
        penalty_count
    } else {
        penalty_count + 1L
    }
    expected_hessian_factorizations <- 1L +
        2L * hessian_parameter_count^2L
    active_restart <- 0L
    optimizer_progress_callback <- function(record) {
        record$restart <- active_restart
        record$postfit_hessian_factorizations <-
            expected_hessian_factorizations
        record$median_objective_seconds <- if (length(objective_seconds)) {
            stats::median(objective_seconds)
        } else {
            NA_real_
        }
        record$estimated_postfit_hessian_seconds <-
            record$postfit_hessian_factorizations *
            record$median_objective_seconds
        checkpoint_state$current_log_sp <<- record$parameters
        checkpoint_state$current_criterion <<- record$criterion
        checkpoint_state$optimizer_progress <<- record
        optimizer_history <- checkpoint_state$optimizer_history
        if (is.null(optimizer_history)) optimizer_history <- list()
        optimizer_history[[length(optimizer_history) + 1L]] <- record
        if (length(optimizer_history) > 1000L) {
            optimizer_history <- utils::tail(optimizer_history, 1000L)
        }
        checkpoint_state$optimizer_history <<- optimizer_history
        checkpoint_state$updated_at <<- as.character(Sys.time())
        reporter$emit(
            1L,
            paste('outer', record$event),
            iteration=record$iteration,
            criterion=format(record$criterion, digits=10),
            projected_gradient=format(
                record$projected_gradient_max,
                digits=5
            ),
            gradient_ratio=format(record$gradient_ratio, digits=5),
            step_max=format(record$step_max, digits=5),
            trust_radius=format(record$trust_radius, digits=5),
            actual_improvement=format(
                record$actual_improvement,
                digits=5
            ),
            acceptance_ratio=format(record$acceptance_ratio, digits=5),
            accepted=record$accepted,
            rejected=record$rejected_steps,
            consecutive_rejections=record$consecutive_rejections,
            step_type=record$step_type,
            curvature_resets=record$curvature_resets,
            curvature_reset=record$curvature_reset,
            hessian_factorizations=record$postfit_hessian_factorizations,
            estimated_hessian_seconds=format(
                record$estimated_postfit_hessian_seconds,
                digits=5
            )
        )
        if (!is.null(checkpoint)) {
            .checkpoint_write(checkpoint_state, checkpoint)
            checkpoint_last_written <<- evaluation_count
            checkpoint_last_written_time <<- proc.time()[['elapsed']]
        }
        invisible(record)
    }
    run_optimizer <- function(arguments) {
        warm <- if (identical(outer_optimizer, 'bfgs_trust')) {
            .safeguarded_outer_bfgs(
                par=arguments$par,
                fn=arguments$fn,
                gr=exact_gradient,
                lower=arguments$lower,
                upper=arguments$upper,
                maxit=optimizer_maxit,
                gradient_tolerance=optimizer_gradient_tolerance,
                initial_radius=optimizer_trust_radius,
                progress=optimizer_progress_callback
            )
        } else {
            do.call(stats::optim, arguments)
        }
        if (!identical(gradient_method, 'hybrid')) return(warm)
        reporter$emit(
            1L,
            'finite-gradient refinement started',
            criterion=format(warm$value, digits=10)
        )
        refinement_arguments <- arguments
        refinement_arguments$par <- warm$par
        refinement_arguments$gr <- if (gradient_cores > 1L) {
            finite_gradient
        } else {
            NULL
        }
        refined <- do.call(stats::optim, refinement_arguments)
        refined$stochastic_warmup <- warm
        refined
    }
    start_points <- list(canonical_initial)
    if (restart_count > 0L) {
        restart_shifts <- c(-4, 4, -8, 8)
        for (i in seq_len(restart_count)) {
            if (i <= length(restart_shifts)) {
                candidate <- rep.int(restart_shifts[[i]], penalty_count)
            } else {
                candidate <- 6 * sin(
                    seq_len(penalty_count) * (i + sqrt(5))
                )
            }
            start_points[[length(start_points) + 1L]] <- candidate
        }
    }
    reporter$phase(
        'smoothing-parameter optimization',
        smoothing_parameters=penalty_count,
        runs=length(start_points),
        gradient=gradient_method,
        outer_optimizer=outer_optimizer
    )
    optimization_runs <- vector('list', length(start_points))
    saved_runs <- checkpoint_state$completed_runs
    if (length(saved_runs)) {
        for (i in seq_len(min(length(saved_runs), length(optimization_runs)))) {
            if (!is.null(saved_runs[[i]])) optimization_runs[[i]] <- saved_runs[[i]]
        }
    }
    if (checkpoint_state$stage %in% c('optimization_complete', 'complete') &&
            !is.null(checkpoint_state$optimization)) {
        optimization_runs <- checkpoint_state$completed_runs
    } else {
        for (run_index in seq_along(start_points)) {
            if (!is.null(optimization_runs[[run_index]])) {
                reporter$emit(
                    1L,
                    'restart reused',
                    restart=run_index - 1L,
                    criterion=optimization_runs[[run_index]]$value
                )
                next
            }
            canonical_start <- start_points[[run_index]]
            active_restart <- run_index - 1L
            start <- canonical_start
            resuming_interrupted_run <- FALSE
            if (identical(checkpoint_state$current_restart, run_index - 1L) &&
                    length(checkpoint_state$current_log_sp) == penalty_count &&
                    all(is.finite(checkpoint_state$current_log_sp))) {
                start <- checkpoint_state$current_log_sp
                resuming_interrupted_run <- resumed &&
                    !isTRUE(all.equal(start, canonical_start, tolerance=0))
            }
            checkpoint_state$current_restart <- run_index - 1L
            checkpoint_state$current_log_sp <- start
            checkpoint_state$stage <- 'optimization'
            if (!is.null(checkpoint)) .checkpoint_write(checkpoint_state, checkpoint)
            reporter$emit(1L, 'restart started', restart=run_index - 1L)
            arguments <- optimization_arguments
            arguments$par <- start
            candidate_run <- run_optimizer(arguments)
            # stats::optim() does not expose the L-BFGS-B curvature history.
            # A warm restart can consequently stop at a worse stationary
            # point. Re-running an interrupted L-BFGS-B start from its
            # canonical point preserves the deterministic uninterrupted
            # result. The safeguarded trust optimizer is explicitly designed
            # to resume from its saved parameters, so validating it this way
            # would discard the computational benefit of checkpointing.
            if (resuming_interrupted_run &&
                    identical(outer_optimizer, 'lbfgsb')) {
                reporter$emit(
                    1L,
                    'restart validation started',
                    restart=run_index - 1L,
                    warm_criterion=format(candidate_run$value, digits=10)
                )
                arguments$par <- canonical_start
                canonical_run <- run_optimizer(arguments)
                if (canonical_run$value < candidate_run$value) {
                    candidate_run <- canonical_run
                }
            }
            optimization_runs[[run_index]] <- candidate_run
            checkpoint_state$completed_runs <- optimization_runs
            checkpoint_state$current_log_sp <- optimization_runs[[run_index]]$par
            checkpoint_state$current_criterion <- optimization_runs[[run_index]]$value
            if (!is.null(checkpoint)) .checkpoint_write(checkpoint_state, checkpoint)
            reporter$emit(
                1L,
                'restart complete',
                restart=run_index - 1L,
                criterion=format(optimization_runs[[run_index]]$value, digits=10),
                convergence=optimization_runs[[run_index]]$convergence
            )
        }
    }
    run_values <- vapply(
        optimization_runs,
        function(run) run$value,
        numeric(1)
    )
    best_run <- which.min(run_values)
    optimization <- optimization_runs[[best_run]]
    checkpoint_state$stage <- 'optimization_complete'
    checkpoint_state$completed_runs <- optimization_runs
    checkpoint_state$optimization <- optimization
    checkpoint_state$best_log_sp <- optimization$par
    checkpoint_state$best_criterion <- optimization$value
    checkpoint_state$current_log_sp <- optimization$par
    checkpoint_state$current_criterion <- optimization$value
    checkpoint_state$evaluation_count <- evaluation_count
    if (!is.null(checkpoint)) .checkpoint_write(checkpoint_state, checkpoint)
    solution <- evaluate(optimization$par, retain=TRUE)
    if (is.null(solution) || !is.numeric(solution$sp)) {
        stop(
            'Sparse REML failed to retain a numeric solution; fields: ',
            paste(names(solution), collapse=', ')
        )
    }
    scale <- solution$penalized_rss / solution$reml_df

    evaluate_unprofiled <- function(parameters) {
        log_sp <- parameters[seq_len(penalty_count)]
        log_scale <- parameters[[penalty_count + 1L]]
        retained <- evaluate(log_sp, retain=TRUE)
        if (is.null(retained)) return(.Machine$double.xmax / 100)
        retained$penalized_rss / exp(log_scale) +
            retained$reml_df * log_scale +
            retained$log_det_system - retained$penalty_logdet
    }
    reporter$phase(
        'outer Hessian',
        parameters=penalty_count + 1L,
        method=hessian_method
    )
    checkpoint_hessian_step <- if (identical(hessian_method, 'profiled')) {
        hessian_step
    } else {
        NA_real_
    }
    reuse_hessian <- identical(checkpoint_state$stage, 'complete') &&
        is.matrix(checkpoint_state$outer_hessian) &&
        identical(checkpoint_state$hessian_method, hessian_method) &&
        identical(checkpoint_state$hessian_step, checkpoint_hessian_step)
    if (reuse_hessian) {
        outer_hessian <- checkpoint_state$outer_hessian
        hessian_evaluations <- checkpoint_state$hessian_evaluations
        reporter$emit(1L, 'completed Hessian reused')
    } else if (identical(hessian_method, 'profiled')) {
        # If U(rho, eta) is the unprofiled -2 REML criterion and eta is log
        # scale, profiling gives
        #   G''(rho) = U[rho,rho] - a a' / df,
        # where a_j = b' (lambda_j S_j) b / scale,
        # U[rho,eta] = -a, and U[eta,eta] = df. Reconstructing U'' from a
        # finite-difference Hessian of G avoids differencing the extra scale
        # dimension and therefore almost halves the required factorizations.
        profiled <- .central_difference_hessian(
            function(log_sp) evaluate(log_sp, retain=FALSE, record=FALSE),
            optimization$par,
            step=hessian_step
        )
        hessian_evaluations <- profiled$evaluations
        coefficients_at_solution <- solution$coefficients
        rss_scores <- vapply(seq_len(penalty_count), function(i) {
            derivative <- solution$sp[[i]] * penalty_components[[i]]
            as.numeric(Matrix::crossprod(
                coefficients_at_solution,
                derivative %*% coefficients_at_solution
            ))
        }, numeric(1))
        scale_scores <- rss_scores / scale
        unprofiled_rho <- profiled$hessian +
            tcrossprod(scale_scores) / solution$reml_df
        unprofiled_hessian <- rbind(
            cbind(unprofiled_rho, -scale_scores),
            c(-scale_scores, solution$reml_df)
        )
        outer_hessian <- unprofiled_hessian / 2
    } else {
        hessian_evaluations <- NA_integer_
        outer_hessian <- stats::optimHess(
            c(log(solution$sp), log(scale)),
            evaluate_unprofiled,
            control=list(ndeps=rep.int(
                hessian_step,
                penalty_count + 1L
            ))
        ) / 2
    }
    hessian_eigenvalues <- tryCatch(
        eigen(
            outer_hessian,
            symmetric=TRUE,
            only.values=TRUE
        )$values,
        error=function(e) rep.int(NA_real_, nrow(outer_hessian))
    )
    hessian_positive_definite <- all(is.finite(hessian_eigenvalues)) &&
        min(hessian_eigenvalues) > sqrt(.Machine$double.eps)
    boundary_tolerance <- 1e-6
    lower_boundary <- which(
        optimization$par - lower_bound <= boundary_tolerance
    )
    upper_boundary <- which(
        upper_bound - optimization$par <= boundary_tolerance
    )
    boundary <- c(
        stats::setNames(rep.int('lower', length(lower_boundary)),
            sp_names[lower_boundary]),
        stats::setNames(rep.int('upper', length(upper_boundary)),
            sp_names[upper_boundary])
    )
    gradient_norm <- if (identical(gradient_method, 'exact')) {
        max(abs(exact_gradient(optimization$par)))
    } else {
        NA_real_
    }
    convergence <- list(
        converged=identical(optimization$convergence, 0L),
        code=optimization$convergence,
        message=optimization$message,
        evaluations=optimization$counts,
        total_objective_evaluations=evaluation_count,
        gradient_norm=gradient_norm,
        boundary=boundary,
        hessian_positive_definite=hessian_positive_definite,
        hessian_min_eigenvalue=if (length(hessian_eigenvalues)) {
            min(hessian_eigenvalues)
        } else {
            NA_real_
        },
        hessian_method=hessian_method,
        hessian_evaluations=hessian_evaluations,
        restart_count=restart_count,
        restart_objectives=run_values,
        best_restart=best_run - 1L,
        resumed=resumed,
        checkpoint=checkpoint,
        global_optimum_certified=FALSE
    )
    if (!convergence$converged) {
        warning(
            'Sparse REML optimizer did not converge (code ',
            convergence$code, '): ', convergence$message,
            call.=FALSE
        )
    }
    if (!hessian_positive_definite) {
        warning(
            'Sparse REML outer Hessian is not positive definite; ',
            'variance-component intervals may be unreliable',
            call.=FALSE
        )
    }
    coefficients <- solution$coefficients
    names(coefficients) <- coefficient_names
    fitted_values <- setup$offset
    for (i in seq_along(matrices)) {
        indices <- matrix_offsets[[i]] + seq_len(matrix_dimensions[[i]])
        fitted_values <- fitted_values + drop(
            .sparse_component_fitted(matrices[[i]], coefficients[indices])
        )
    }
    fitted_values <- as.numeric(fitted_values)
    raw_residuals <- setup$y - fitted_values
    raw_residuals <- as.numeric(raw_residuals)
    sp <- stats::setNames(solution$sp, sp_names)
    identifiability <- design$identifiability
    identifiability$global <- list(
        dimension=dimension,
        rank=dimension,
        condition_indicator=.cdr_factor_condition_indicator(solution$factor),
        parametric=alias_resolution$info,
        resolution=rank_resolution$resolution
    )

    reporter$phase('finalization')
    out <- list(
        coefficients=coefficients,
        fitted.values=fitted_values,
        residuals=raw_residuals,
        linear.predictors=fitted_values,
        family=family,
        sp=sp,
        scale=scale,
        sig2=scale,
        reml.scale=scale,
        method='REML',
        smooth=smooths,
        paraPen=setup$paraPen,
        full.sp=setup$full.sp,
        outer.info=list(
            hess=outer_hessian,
            conv=convergence$converged,
            message=convergence$message
        ),
        converged=convergence$converged,
        df.residual=NA_real_,
        y=setup$y,
        prior.weights=weights,
        offset=setup$offset,
        reml=solution$criterion,
        optimizer=optimization,
        sparse=list(
            factor=solution$factor,
            dimension=dimension,
            nnzero=design_nonzeros,
            system_nnzero=Matrix::nnzero(solution$system),
            numeric_updates=numeric_update_count,
            parallel_gradient_factorizations=parallel_gradient_factorizations,
            full_factorizations=full_factorization_count,
            numeric_update_error=numeric_update_error,
            cached_objectives=length(ls(objective_cache, all.names=TRUE)),
            gradient=gradient_method,
            outer_optimizer=outer_optimizer,
            optimizer_progress=checkpoint_state$optimizer_progress,
            optimizer_history=checkpoint_state$optimizer_history,
            hessian=hessian_method,
            hessian_evaluations=hessian_evaluations,
            exact_trace_rhs=exact_trace_rhs,
            gradient_probes=if (gradient_method %in% c(
                'stochastic', 'hybrid'
            )) gradient_probes else 0L,
            factor_class=class(solution$factor)[[1L]],
            supernodal=supernodal,
            factor_nonzeros=.cdr_factor_nonzeros(solution$factor),
            condition_indicator=.cdr_factor_condition_indicator(solution$factor),
            crossprod_chunks=length(crossprod_starts),
            crossprod_chunk_size=crossprod_chunk_size,
            streamed_grouped_terms=streamed_grouped_terms,
            convergence=convergence,
            penalty_components=penalty_components,
            observation_count=observation_count,
            schur_group=if (!inherits(
                solution$factor,
                'cdrgam_schur_factor'
            )) NULL else schur_layout$group,
            schur_core_dimension=if (!inherits(
                solution$factor,
                'cdrgam_schur_factor'
            )) {
                NA_integer_
            } else {
                length(schur_layout$core)
            }
        ),
        cdrgam=list(
            schema_version=1L,
            engine='sparse',
            backend='sparse',
            formula=list(
                user=design$formula,
                normalized=design$normalized_formula,
                effective=design$effective_formula,
                mgcv=.sparse_expanded_formula(design)
            ),
            preparation=list(
                plan=design$plan,
                stream=design$stream,
                specification=design$specification,
                identifiability=design$identifiability
            ),
            identifiability=identifiability,
            rank=list(
                action=rank_action,
                parametric=alias_resolution$info,
                resolution=rank_resolution$resolution,
                regularization=fixed_ridge,
                tolerance=tolerance
            ),
            term_labels=names(design$terms),
            terms=term_metadata,
            prediction=list(
                setup=.cdr_prediction_setup(setup),
                random_effects=random_prediction
            ),
            solver='sparse Gaussian REML solver'
        )
    )
    class(out) <- c('cdrgam_sparse', 'cdrgam')
    checkpoint_state$stage <- 'complete'
    checkpoint_state$optimization <- optimization
    checkpoint_state$best_log_sp <- optimization$par
    checkpoint_state$best_criterion <- solution$criterion
    checkpoint_state$evaluation_count <- evaluation_count
    checkpoint_state$outer_hessian <- outer_hessian
    checkpoint_state$hessian_method <- hessian_method
    checkpoint_state$hessian_step <- checkpoint_hessian_step
    checkpoint_state$hessian_evaluations <- hessian_evaluations
    checkpoint_state$converged <- convergence$converged
    if (!is.null(checkpoint)) .checkpoint_write(checkpoint_state, checkpoint)
    reporter$emit(
        1L,
        'fit complete',
        criterion=format(solution$criterion, digits=10),
        evaluations=evaluation_count,
        converged=convergence$converged
    )
    out
}

#' @export
coef.cdrgam_sparse <- function(object, ...) object$coefficients

#' @export
fitted.cdrgam_sparse <- function(object, ...) object$fitted.values

#' @export
residuals.cdrgam_sparse <- function(object, ...) object$residuals

#' @export
deviance.cdrgam_sparse <- function(object, ...) {
    sum(object$prior.weights * object$residuals^2)
}

#' @export
nobs.cdrgam_sparse <- function(object, ...) length(object$y)

#' @export
logLik.cdrgam_sparse <- function(object, ...) {
    n <- length(object$y)
    value <- -0.5 * (
        n * log(2 * pi * object$scale) +
            deviance.cdrgam_sparse(object) / object$scale
    )
    attr(value, 'df') <- .sparse_effective_df(object) + 1
    attr(value, 'nobs') <- n
    class(value) <- 'logLik'
    value
}

.sparse_selected_vcov <- function(object, indices, unconditional=FALSE) {
    dimension <- object$sparse$dimension
    selector <- Matrix::sparseMatrix(
        i=indices,
        j=seq_along(indices),
        x=1,
        dims=c(dimension, length(indices))
    )
    columns <- .cdr_factor_solve(object$sparse$factor, selector)
    covariance <- as.matrix(columns[indices, , drop=FALSE]) * object$scale
    covariance <- (covariance + t(covariance)) / 2
    if (isTRUE(unconditional)) {
        penalty_count <- length(object$sp)
        hessian <- object$outer.info$hess
        rho_covariance <- tryCatch(
            solve(hessian)[seq_len(penalty_count), seq_len(penalty_count), drop=FALSE],
            error=function(e) NULL
        )
        if (is.null(rho_covariance) || any(!is.finite(rho_covariance))) {
            stop('Smoothing-parameter covariance is unavailable')
        }
        derivatives <- matrix(0, nrow=length(indices), ncol=penalty_count)
        for (i in seq_len(penalty_count)) {
            rhs <- object$sp[[i]] *
                object$sparse$penalty_components[[i]] %*%
                object$coefficients
            derivatives[, i] <- -.cdr_factor_solve(
                object$sparse$factor,
                rhs
            )[indices, 1L]
        }
        covariance <- covariance +
            derivatives %*% rho_covariance %*% t(derivatives)
        covariance <- (covariance + t(covariance)) / 2
    }
    coefficient_names <- names(object$coefficients)[indices]
    dimnames(covariance) <- list(coefficient_names, coefficient_names)
    covariance
}

.sparse_vcov_diagonal <- function(object, chunk_size=256L) {
    dimension <- object$sparse$dimension
    output <- numeric(dimension)
    for (start in seq.int(1L, dimension, by=chunk_size)) {
        indices <- start:min(dimension, start + chunk_size - 1L)
        output[indices] <- diag(.sparse_selected_vcov(object, indices))
    }
    stats::setNames(output, names(object$coefficients))
}

.sparse_effective_df <- function(object) {
    penalty_components <- object$sparse$penalty_components
    if (is.null(penalty_components)) return(NA_real_)
    penalty_trace <- 0
    for (i in seq_along(penalty_components)) {
        penalty_trace <- penalty_trace + .sparse_logdet_score(
            object$sparse$factor,
            object$sp[[i]] * penalty_components[[i]],
            chunk_size=256L
        )
    }
    object$sparse$dimension - penalty_trace
}

#' @export
vcov.cdrgam_sparse <- function(object, unconditional=FALSE, ...) {
    dimension <- object$sparse$dimension
    limit <- getOption('cdrgam.max_vcov_elements', 25e6)
    if (dimension^2 > limit) {
        stop(
            'The full covariance matrix would contain ',
            format(dimension^2, scientific=FALSE, big.mark=','),
            ' elements. Use estimate_irf() for selected covariance solves or ',
            'increase option cdrgam.max_vcov_elements explicitly.'
        )
    }
    .sparse_selected_vcov(
        object,
        seq_len(dimension),
        unconditional=unconditional
    )
}

#' @export
predict.cdrgam_sparse <- function(object, newdata=NULL, ...) {
    if (is.null(newdata)) return(object$fitted.values)
    if (!is.list(newdata) ||
            !all(c('impulses', 'responses') %in% names(newdata))) {
        stop('newdata must contain impulse and response data frames')
    }
    predict_cdrgam(
        object,
        impulses=newdata$impulses,
        responses=newdata$responses,
        ...
    )
}

#' @export
print.cdrgam_sparse <- function(x, ...) {
    cat('Continuous-time deconvolutional GAM\n')
    cat('  backend:', x$cdrgam$solver, '\n')
    cat('  IRF terms:', paste(x$cdrgam$term_labels, collapse=', '), '\n')
    cat('  coefficients:', length(x$coefficients), '\n')
    cat('  REML criterion:', format(x$reml, digits=7), '\n')
    cat('  converged:', if (isTRUE(x$converged)) 'yes' else 'no', '\n')
    if (length(x$sparse$convergence$boundary)) {
        cat(
            '  smoothing parameters at boundary:',
            paste(names(x$sparse$convergence$boundary), collapse=', '),
            '\n'
        )
    }
    invisible(x)
}

#' @export
summary.cdrgam_sparse <- function(object, ...) {
    covariance_diagonal <- .sparse_vcov_diagonal(object)
    standard_errors <- sqrt(pmax(0, covariance_diagonal))
    table <- cbind(
        Estimate=object$coefficients,
        `Std. Error`=standard_errors,
        `t value`=object$coefficients / standard_errors
    )
    edf <- .sparse_effective_df(object)
    df_residual <- object$sparse$observation_count - edf
    out <- list(
        call=NULL,
        family=object$family,
        coefficients=table,
        sp=object$sp,
        scale=object$scale,
        edf=edf,
        df.residual=df_residual,
        reml=object$reml,
        backend=object$cdrgam$solver,
        convergence=object$sparse$convergence,
        rank=object$cdrgam$rank
    )
    class(out) <- 'summary.cdrgam_sparse'
    out
}

#' @export
print.summary.cdrgam_sparse <- function(x, ...) {
    cat('CDR-GAM sparse fit summary\n')
    cat('Family:', x$family$family, '(', x$family$link, ')\n')
    cat('Backend:', x$backend, '\n')
    cat('Converged:', if (isTRUE(x$convergence$converged)) 'yes' else 'no', '\n\n')
    stats::printCoefmat(x$coefficients)
    cat('\nSmoothing parameters:\n')
    print(x$sp)
    cat('\nScale:', format(x$scale, digits=7), '\n')
    cat('Effective degrees of freedom:', format(x$edf, digits=7), '\n')
    cat('Residual degrees of freedom:', format(x$df.residual, digits=7), '\n')
    if (isTRUE(x$rank$parametric$corrected)) {
        cat(
            'Aliased parametric coefficients removed:',
            paste(x$rank$parametric$dropped, collapse=', '), '\n'
        )
    }
    if (!is.null(x$rank$resolution) &&
            !identical(x$rank$resolution, 'none')) {
        cat('Rank resolution:', x$rank$resolution, '\n')
    }
    if (length(x$convergence$boundary)) {
        cat(
            'Smoothing parameters at boundary:',
            paste(names(x$convergence$boundary), collapse=', '),
            '\n'
        )
    }
    invisible(x)
}
