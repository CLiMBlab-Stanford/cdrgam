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
    if (is.null(factor$transfer)) {
        return(.Call(
            'cdrgam_schur_solve',
            as.integer(factor$core),
            lapply(factor$blocks, as.integer),
            factor$cross,
            factor$block_cholesky,
            factor$core_cholesky,
            rhs,
            PACKAGE='cdrgam'
        ))
    }
    .Call(
        'cdrgam_schur_solve_batched',
        as.integer(factor$core),
        lapply(factor$blocks, as.integer),
        factor$transfer,
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
    out <- list(
        core=core,
        blocks=blocks,
        transfer=factor[[3L]],
        block_connections=factor[[4L]],
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

.sparse_penalty_logdet_derivatives <- function(blocks, sp, penalty_count) {
    score <- numeric(penalty_count)
    hessian <- matrix(0, penalty_count, penalty_count)
    for (block in blocks) {
        if (isTRUE(block$identity)) {
            index <- block$sp_index[[1L]]
            score[[index]] <- score[[index]] + block$dimension
            next
        }
        penalty <- matrix(0, nrow=block$dimension, ncol=block$dimension)
        derivatives <- vector('list', length(block$sp_index))
        for (j in seq_along(block$sp_index)) {
            index <- block$sp_index[[j]]
            derivatives[[j]] <- sp[[index]] * block$S[[j]]
            penalty <- penalty + derivatives[[j]]
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
        inverse_derivatives <- lapply(derivatives, function(derivative) {
            inverse %*% derivative
        })
        local_score <- vapply(
            inverse_derivatives,
            function(value) sum(diag(value)),
            numeric(1)
        )
        for (a in seq_along(block$sp_index)) {
            index_a <- block$sp_index[[a]]
            score[[index_a]] <- score[[index_a]] +
                block$repetitions * local_score[[a]]
            for (b in seq_len(a)) {
                index_b <- block$sp_index[[b]]
                cross <- sum(
                    inverse_derivatives[[b]] *
                        t(inverse_derivatives[[a]])
                )
                value <- block$repetitions * (
                    if (index_a == index_b) local_score[[a]] else 0
                ) - block$repetitions * cross
                hessian[index_a, index_b] <-
                    hessian[index_a, index_b] + value
                if (index_a != index_b) {
                    hessian[index_b, index_a] <-
                        hessian[index_b, index_a] + value
                }
            }
        }
    }
    list(score=score, hessian=hessian)
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

.cdrgam_memory_number <- function(value, default_multiplier=1) {
    if (!length(value) || is.na(value[[1L]])) return(NA_real_)
    value <- trimws(as.character(value[[1L]]))
    if (!nzchar(value) || value %in% c('max', 'unlimited')) return(Inf)
    matched <- regexec(
        '^([0-9]+(?:\\.[0-9]+)?)[[:space:]]*([KMGTPE]?)B?$',
        value,
        ignore.case=TRUE
    )
    pieces <- regmatches(value, matched)[[1L]]
    if (!length(pieces)) return(NA_real_)
    units <- c(K=1, M=2, G=3, T=4, P=5, E=6)
    suffix <- toupper(pieces[[3L]])
    multiplier <- if (nzchar(suffix)) 1024^units[[suffix]] else {
        default_multiplier
    }
    as.numeric(pieces[[2L]]) * multiplier
}

.cdrgam_read_memory_value <- function(path, default_multiplier=1) {
    value <- tryCatch(readLines(path, n=1L, warn=FALSE), error=function(e) '')
    .cdrgam_memory_number(value, default_multiplier=default_multiplier)
}

.cdrgam_proc_memory <- function(field) {
    lines <- tryCatch(
        readLines('/proc/self/status', warn=FALSE),
        error=function(e) character()
    )
    line <- grep(paste0('^', field, ':'), lines, value=TRUE)
    if (!length(line)) return(NA_real_)
    .cdrgam_memory_number(
        sub(paste0('^', field, ':[[:space:]]*'), '', line[[1L]]),
        default_multiplier=1024
    )
}

.cdrgam_r_memory <- function() {
    usage <- tryCatch(gc(), error=function(error) NULL)
    if (is.null(usage) || ncol(usage) < 2L) return(NA_real_)
    bytes <- sum(usage[, 2L]) * 1024^2
    if (is.finite(bytes) && bytes >= 0) bytes else NA_real_
}

.cdrgam_cgroup_memory_directories <- function() {
    directories <- '/sys/fs/cgroup'
    entries <- tryCatch(
        readLines('/proc/self/cgroup', warn=FALSE),
        error=function(e) character()
    )
    for (entry in entries) {
        fields <- strsplit(entry, ':', fixed=TRUE)[[1L]]
        if (length(fields) != 3L) next
        relative <- sub('^/+', '', fields[[3L]])
        if (!nzchar(fields[[2L]])) {
            directories <- c(
                directories,
                file.path('/sys/fs/cgroup', relative)
            )
        } else if ('memory' %in% strsplit(fields[[2L]], ',', fixed=TRUE)[[1L]]) {
            directories <- c(
                directories,
                file.path('/sys/fs/cgroup/memory', relative)
            )
        }
    }
    unique(normalizePath(directories, mustWork=FALSE))
}

.cdrgam_memory_availability <- function() {
    candidates <- list()
    add_candidate <- function(source, limit, used, available=NULL) {
        if (is.null(available)) available <- limit - used
        if (!is.finite(available) || available < 0) return(invisible(NULL))
        candidates[[length(candidates) + 1L]] <<- list(
            source=source,
            limit_bytes=limit,
            used_bytes=used,
            available_bytes=available
        )
        invisible(NULL)
    }

    process_rss <- .cdrgam_proc_memory('VmRSS')
    process_memory <- if (is.finite(process_rss)) {
        process_rss
    } else {
        .cdrgam_r_memory()
    }
    override <- getOption('cdrgam.memory_limit_bytes', NA_real_)
    if (is.numeric(override) && length(override) == 1L &&
            is.finite(override) && override > 0 && is.finite(process_memory)) {
        add_candidate('option', as.numeric(override), process_memory)
    }

    for (directory in .cdrgam_cgroup_memory_directories()) {
        current <- directory
        root <- normalizePath('/sys/fs/cgroup', mustWork=FALSE)
        repeat {
            v2_limit <- .cdrgam_read_memory_value(file.path(current, 'memory.max'))
            v2_used <- .cdrgam_read_memory_value(file.path(current, 'memory.current'))
            if (is.finite(v2_limit) && is.finite(v2_used)) {
                add_candidate('cgroup-v2', v2_limit, v2_used)
            }
            v1_limit <- .cdrgam_read_memory_value(
                file.path(current, 'memory.limit_in_bytes')
            )
            v1_used <- .cdrgam_read_memory_value(
                file.path(current, 'memory.usage_in_bytes')
            )
            if (is.finite(v1_limit) && is.finite(v1_used)) {
                add_candidate('cgroup-v1', v1_limit, v1_used)
            }
            parent <- dirname(current)
            if (identical(current, root) || identical(parent, current) ||
                    !startsWith(parent, root)) break
            current <- parent
        }
    }

    slurm_limit <- .cdrgam_memory_number(
        Sys.getenv('SLURM_MEM_PER_NODE'),
        default_multiplier=1024^2
    )
    if (!is.finite(slurm_limit)) {
        per_cpu <- .cdrgam_memory_number(
            Sys.getenv('SLURM_MEM_PER_CPU'),
            default_multiplier=1024^2
        )
        cpus <- suppressWarnings(as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE')))
        if (is.finite(per_cpu) && is.finite(cpus) && cpus > 0) {
            slurm_limit <- per_cpu * cpus
        }
    }
    if (is.finite(slurm_limit) && is.finite(process_rss)) {
        add_candidate('slurm', slurm_limit, process_rss)
    }

    meminfo <- tryCatch(
        readLines('/proc/meminfo', warn=FALSE),
        error=function(e) character()
    )
    available_line <- grep('^MemAvailable:', meminfo, value=TRUE)
    total_line <- grep('^MemTotal:', meminfo, value=TRUE)
    if (length(available_line)) {
        available <- .cdrgam_memory_number(
            sub('^MemAvailable:[[:space:]]*', '', available_line[[1L]]),
            default_multiplier=1024
        )
        total <- if (length(total_line)) {
            .cdrgam_memory_number(
                sub('^MemTotal:[[:space:]]*', '', total_line[[1L]]),
                default_multiplier=1024
            )
        } else {
            NA_real_
        }
        if (is.finite(available)) {
            add_candidate(
                'system',
                total,
                if (is.finite(total)) total - available else NA_real_,
                available=available
            )
        }
    }

    if (!length(candidates)) {
        return(list(
            source='unknown',
            limit_bytes=NA_real_,
            used_bytes=NA_real_,
            available_bytes=NA_real_
        ))
    }
    candidates[[which.min(vapply(
        candidates,
        `[[`,
        numeric(1),
        'available_bytes'
    ))]]
}

.sparse_analytic_hessian_memory <- function(supports, dimension, chunk_size) {
    union_count <- length(unique(unlist(supports, use.names=FALSE)))
    support_lengths <- lengths(supports)
    retained_elements <- as.double(union_count) * sum(support_lengths)
    largest_matrix_elements <- if (length(support_lengths)) {
        as.double(union_count) * max(support_lengths)
    } else {
        0
    }
    chunk_columns <- if (length(support_lengths)) {
        min(max(support_lengths), as.double(chunk_size))
    } else {
        0
    }
    chunk_elements <- as.double(dimension) * chunk_columns
    retained_bytes <- 8 * retained_elements
    transient_bytes <- 8 * (largest_matrix_elements + 3 * chunk_elements)
    list(
        retained_elements=retained_elements,
        retained_bytes=retained_bytes,
        transient_bytes=transient_bytes,
        estimated_peak_bytes=retained_bytes + transient_bytes
    )
}

.sparse_select_hessian_method <- function(
        requested,
        supports,
        dimension,
        chunk_size,
        memory=.cdrgam_memory_availability()
) {
    estimate <- .sparse_analytic_hessian_memory(
        supports,
        dimension,
        chunk_size
    )
    available <- memory$available_bytes
    reserve <- 1024^3
    budget <- if (is.finite(available)) {
        max(0, min(available / 2, available - reserve))
    } else {
        NA_real_
    }
    if (!identical(requested, 'auto')) {
        method <- requested
        reason <- 'explicit Hessian method'
    } else if (!is.finite(budget)) {
        method <- 'gradient'
        reason <- 'available memory could not be determined'
    } else if (estimate$estimated_peak_bytes <= budget) {
        method <- 'analytic'
        reason <- 'analytic Hessian fits the automatic memory budget'
    } else {
        method <- 'gradient'
        reason <- 'analytic Hessian exceeds the automatic memory budget'
    }
    c(
        list(
            policy_version=1L,
            requested=requested,
            method=method,
            reason=reason,
            memory_source=memory$source,
            memory_limit_bytes=memory$limit_bytes,
            memory_used_bytes=memory$used_bytes,
            memory_available_bytes=available,
            memory_budget_bytes=budget
        ),
        estimate
    )
}

.sparse_analytic_hessian_element_limit <- function(requested, selection) {
    if (identical(requested, 'auto') &&
            identical(selection$method, 'analytic')) {
        max(1, selection$retained_elements)
    } else {
        getOption('cdrgam.max_analytic_hessian_elements', 2e8)
    }
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

.sparse_schur_trace_plan <- function(layout, components) {
    dimension <- nrow(components[[1L]])
    block_count <- length(layout$blocks)
    block_label <- integer(dimension)
    block_position <- integer(dimension)
    core_position <- integer(dimension)
    core_position[layout$core] <- seq_along(layout$core)
    for (g in seq_len(block_count)) {
        indices <- layout$blocks[[g]]
        block_label[indices] <- g
        block_position[indices] <- seq_along(indices)
    }
    core <- vector('list', length(components))
    blocks <- vector('list', block_count)
    for (component_index in seq_along(components)) {
        triplet <- methods::as(
            methods::as(components[[component_index]], 'generalMatrix'),
            'TsparseMatrix'
        )
        rows <- methods::slot(triplet, 'i') + 1L
        columns <- methods::slot(triplet, 'j') + 1L
        values <- methods::slot(triplet, 'x')
        row_blocks <- block_label[rows]
        column_blocks <- block_label[columns]
        if (any(row_blocks != column_blocks & values != 0)) return(NULL)
        core_entries <- row_blocks == 0L
        core[[component_index]] <- list(
            i=core_position[rows[core_entries]],
            j=core_position[columns[core_entries]],
            x=values[core_entries]
        )
        block_entries <- which(row_blocks > 0L & values != 0)
        if (!length(block_entries)) next
        by_block <- split(block_entries, row_blocks[block_entries])
        for (block_name in names(by_block)) {
            selected <- by_block[[block_name]]
            g <- as.integer(block_name)
            blocks[[g]][[length(blocks[[g]]) + 1L]] <- list(
                component=component_index,
                i=block_position[rows[selected]],
                j=block_position[columns[selected]],
                x=values[selected]
            )
        }
    }
    list(core=core, blocks=blocks)
}

.sparse_schur_logdet_scores <- function(factor, plan, sp) {
    core_inverse <- chol2inv(factor$core_cholesky)
    scores <- numeric(length(plan$core))
    for (i in seq_along(plan$core)) {
        entries <- plan$core[[i]]
        if (length(entries$x)) {
            scores[[i]] <- sum(
                core_inverse[cbind(entries$i, entries$j)] * entries$x
            )
        }
    }
    batched_option <- getOption('cdrgam.batched_schur_trace')
    blas <- unname(extSoftVersion()[['BLAS']])
    accelerated_blas <- length(blas) == 1L && !is.na(blas) && grepl(
        'openblas|mkl|blis|accelerate|atlas',
        tolower(blas)
    )
    batched <- !is.null(factor$transfer) && if (is.null(batched_option)) {
        accelerated_blas
    } else {
        isTRUE(batched_option)
    }
    transformed <- if (batched) core_inverse %*% factor$transfer else NULL
    offset <- 0L
    for (g in seq_along(factor$blocks)) {
        entries <- plan$blocks[[g]]
        dimension <- length(factor$blocks[[g]])
        columns <- offset + seq_len(dimension)
        offset <- offset + dimension
        if (!length(entries)) next
        block_inverse <- chol2inv(factor$block_cholesky[[g]])
        if (batched) {
            local_transfer <- factor$transfer[, columns, drop=FALSE]
            block_inverse <- block_inverse + crossprod(
                local_transfer,
                transformed[, columns, drop=FALSE]
            )
        } else if (!is.null(factor$transfer)) {
            local_transfer <- t(factor$transfer[, columns, drop=FALSE])
            support <- which(colSums(abs(local_transfer)) != 0)
            if (length(support)) {
                local_transfer <- local_transfer[, support, drop=FALSE]
                block_inverse <- block_inverse + local_transfer %*%
                    core_inverse[support, support, drop=FALSE] %*%
                    t(local_transfer)
            }
        } else {
            local_transfer <- block_inverse %*% t(factor$cross[[g]])
            support <- which(colSums(abs(local_transfer)) != 0)
            if (length(support)) {
                local_transfer <- local_transfer[, support, drop=FALSE]
                block_inverse <- block_inverse + local_transfer %*%
                    core_inverse[support, support, drop=FALSE] %*%
                    t(local_transfer)
            }
        }
        for (record in entries) {
            scores[[record$component]] <- scores[[record$component]] + sum(
                block_inverse[cbind(record$i, record$j)] * record$x
            )
        }
    }
    sp * scores
}

.sparse_optimizer_selection <- function(
        gradient_requested,
        outer_optimizer_requested,
        trace_method,
        factor,
        objective_seconds,
        exact_trace_rhs,
        penalty_count,
        gradient_cores,
        dimension,
        saved=NULL
) {
    finite_evaluations <- 2 * penalty_count
    predicted_finite_seconds <- objective_seconds * finite_evaluations /
        max(1L, gradient_cores)
    schur_trace <- identical(trace_method, 'schur_inverse') &&
        inherits(factor, 'cdrgam_schur_factor')
    core_dimension <- if (schur_trace) length(factor$core) else NA_integer_
    block_count <- if (schur_trace) length(factor$blocks) else 0L
    block_dimensions <- if (schur_trace) lengths(factor$blocks) else integer()
    block_connections <- if (schur_trace) {
        if (!is.null(factor$block_connections)) {
            factor$block_connections
        } else {
            vapply(factor$cross, function(cross) {
                sum(rowSums(abs(cross)) != 0)
            }, integer(1))
        }
    } else {
        integer()
    }
    block_work <- if (schur_trace) sum(
        block_dimensions * block_connections^2 +
            block_dimensions^2 * block_connections + block_dimensions^3
    ) else 0
    core_factor_work <- if (schur_trace) {
        max(as.double(core_dimension)^3 / 3, 1)
    } else {
        NA_real_
    }
    exact_cost_ratio <- if (schur_trace) {
        max(1, 2 + block_work / core_factor_work)
    } else {
        max(1, exact_trace_rhs / 64)
    }
    predicted_exact_seconds <- objective_seconds * exact_cost_ratio
    exact_core_bytes <- if (schur_trace) {
        8 * as.double(core_dimension)^2
    } else {
        0
    }
    exact_memory_limit <- 4 * 1024^3
    exact_memory_ok <- exact_core_bytes <= exact_memory_limit
    compact_exact_trace <- exact_trace_rhs <= max(512L, 64L * penalty_count)

    saved_gradient <- if (is.list(saved) && identical(
        saved$policy_version, 1L
    )) saved$gradient else NULL
    saved_exact_seconds <- if (is.list(saved)) {
        saved$measured_exact_seconds
    } else {
        NULL
    }
    reuse_saved <- identical(gradient_requested, 'auto') &&
        length(saved_gradient) == 1L &&
        saved_gradient %in% c('exact', 'finite')
    if (reuse_saved) {
        gradient <- saved_gradient
        reason <- 'checkpoint resolution reused'
    } else if (!identical(gradient_requested, 'auto')) {
        gradient <- gradient_requested
        reason <- 'explicit gradient'
    } else if (identical(outer_optimizer_requested, 'bfgs_trust')) {
        gradient <- 'exact'
        reason <- 'explicit trust optimizer requires exact gradients'
    } else if (!exact_memory_ok) {
        gradient <- 'finite'
        reason <- 'predicted exact-gradient memory exceeds automatic limit'
    } else if (compact_exact_trace) {
        gradient <- 'exact'
        reason <- 'exact trace has compact penalty support'
    } else if (predicted_exact_seconds <= predicted_finite_seconds) {
        gradient <- 'exact'
        reason <- 'predicted exact gradient is no slower than finite differences'
    } else {
        gradient <- 'finite'
        reason <- 'predicted finite differences are faster'
    }
    outer_optimizer <- if (identical(outer_optimizer_requested, 'auto')) {
        if (identical(gradient, 'exact')) 'bfgs_trust' else 'lbfgsb'
    } else {
        outer_optimizer_requested
    }
    if (identical(outer_optimizer, 'bfgs_trust') &&
            !identical(gradient, 'exact')) {
        stop('Internal optimizer selection produced an invalid trust gradient')
    }
    list(
        policy_version=1L,
        gradient=gradient,
        outer_optimizer=outer_optimizer,
        reason=reason,
        objective_seconds=objective_seconds,
        predicted_exact_seconds=predicted_exact_seconds,
        predicted_finite_seconds=predicted_finite_seconds,
        exact_cost_ratio=exact_cost_ratio,
        exact_core_bytes=exact_core_bytes,
        exact_memory_limit=exact_memory_limit,
        exact_memory_ok=exact_memory_ok,
        exact_trace_rhs=exact_trace_rhs,
        compact_exact_trace=compact_exact_trace,
        system_dimension=dimension,
        finite_difference_evaluations=finite_evaluations,
        schur_trace=schur_trace,
        schur_core_dimension=core_dimension,
        schur_block_count=block_count,
        schur_block_dimension_max=if (block_count) {
            max(block_dimensions)
        } else {
            0L
        },
        schur_block_connection_mean=if (block_count) {
            mean(block_connections)
        } else {
            NA_real_
        },
        schur_block_connection_max=if (block_count) {
            max(block_connections)
        } else {
            NA_integer_
        },
        measured_exact_seconds=saved_exact_seconds,
        median_exact_seconds=if (length(saved_exact_seconds)) {
            stats::median(saved_exact_seconds)
        } else {
            NA_real_
        }
    )
}

.sparse_logdet_scores <- function(
        factor,
        components,
        sp,
        supports=.sparse_penalty_supports(components),
        chunk_size=64L,
        schur_plan=NULL
) {
    count <- length(components)
    if (!count) return(numeric())
    if (!is.null(schur_plan) && inherits(factor, 'cdrgam_schur_factor')) {
        return(.sparse_schur_logdet_scores(factor, schur_plan, sp))
    }
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

.sparse_analytic_profiled_hessian <- function(
        factor,
        blocks,
        sp,
        components,
        supports,
        coefficients,
        penalized_rss,
        reml_df,
        chunk_size=64L,
        element_limit=getOption('cdrgam.max_analytic_hessian_elements', 2e8)
) {
    penalty_count <- length(components)
    union_support <- sort(unique(unlist(supports, use.names=FALSE)))
    support_positions <- lapply(supports, match, table=union_support)
    stored_elements <- length(union_support) * sum(lengths(supports))
    if (!is.numeric(element_limit) || length(element_limit) != 1L ||
            !is.finite(element_limit) || element_limit < 1) {
        stop('The analytic Hessian element limit must be positive')
    }
    if (stored_elements > element_limit) {
        stop(
            'The analytic outer Hessian would retain ',
            format(stored_elements, scientific=FALSE, big.mark=','),
            ' selected-inverse elements; the current limit is ',
            format(element_limit, scientific=FALSE, big.mark=','),
            '. Use hessian="gradient" or raise option ',
            'cdrgam.max_analytic_hessian_elements explicitly.'
        )
    }
    inverse_derivatives <- vector('list', penalty_count)
    rss_scores <- numeric(penalty_count)
    coefficient_derivatives <- vector('list', penalty_count)
    system_scores <- numeric(penalty_count)
    for (i in seq_len(penalty_count)) {
        support <- supports[[i]]
        positions <- support_positions[[i]]
        derivative <- sp[[i]] * components[[i]]
        selected_inverse_derivative <- matrix(
            0,
            nrow=length(union_support),
            ncol=length(support)
        )
        if (length(support)) {
            for (start in seq.int(1L, length(support), by=chunk_size)) {
                selected <- start:min(
                    length(support),
                    start + chunk_size - 1L
                )
                columns <- support[selected]
                solved <- .cdr_factor_solve(
                    factor,
                    derivative[, columns, drop=FALSE]
                )
                selected_inverse_derivative[, selected] <- as.matrix(
                    solved[union_support, , drop=FALSE]
                )
            }
        }
        inverse_derivatives[[i]] <- selected_inverse_derivative
        local_rhs <- as.numeric(
            derivative[support, support, drop=FALSE] %*%
                coefficients[support]
        )
        rss_scores[[i]] <- sum(coefficients[support] * local_rhs)
        coefficient_derivatives[[i]] <- drop(
            selected_inverse_derivative %*% coefficients[support]
        )
        system_scores[[i]] <- if (length(support)) {
            sum(selected_inverse_derivative[cbind(
                positions,
                seq_along(support)
            )])
        } else {
            0
        }
    }
    rss_hessian <- matrix(0, penalty_count, penalty_count)
    system_hessian <- matrix(0, penalty_count, penalty_count)
    for (i in seq_len(penalty_count)) {
        support_i <- supports[[i]]
        positions_i <- support_positions[[i]]
        derivative_i <- sp[[i]] * components[[i]]
        rhs_i <- as.numeric(
            derivative_i[support_i, support_i, drop=FALSE] %*%
                coefficients[support_i]
        )
        for (j in seq_len(i)) {
            support_j <- supports[[j]]
            positions_j <- support_positions[[j]]
            rss_cross <- sum(
                rhs_i * coefficient_derivatives[[j]][positions_i]
            )
            second_rss <- if (i == j) rss_scores[[i]] else 0
            second_rss <- second_rss - 2 * rss_cross
            rss_value <- reml_df * (
                second_rss / penalized_rss -
                    rss_scores[[i]] * rss_scores[[j]] /
                        penalized_rss^2
            )
            trace_cross <- 0
            if (length(support_j)) {
                for (start in seq.int(
                        1L,
                        length(support_j),
                        by=chunk_size
                    )) {
                    selected <- start:min(
                        length(support_j),
                        start + chunk_size - 1L
                    )
                    left <- inverse_derivatives[[j]][
                        positions_i,
                        selected,
                        drop=FALSE
                    ]
                    right <- inverse_derivatives[[i]][
                        positions_j[selected],
                        ,
                        drop=FALSE
                    ]
                    trace_cross <- trace_cross + sum(left * t(right))
                }
            }
            system_value <- if (i == j) system_scores[[i]] else 0
            system_value <- system_value - trace_cross
            rss_hessian[i, j] <- rss_hessian[j, i] <- rss_value
            system_hessian[i, j] <- system_hessian[j, i] <- system_value
        }
    }
    penalty <- .sparse_penalty_logdet_derivatives(
        blocks,
        sp,
        penalty_count
    )
    hessian <- rss_hessian + system_hessian - penalty$hessian
    hessian <- (hessian + t(hessian)) / 2
    list(
        hessian=hessian,
        rhs=sum(lengths(supports)),
        retained_elements=stored_elements,
        rss_scores=rss_scores,
        system_scores=system_scores,
        penalty_scores=penalty$score
    )
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
    if (!count) character() else if (count == 1L) label else
        paste0(label, seq_len(count))
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

.cdr_penalty_null_basis <- function(penalties, tolerance=sqrt(.Machine$double.eps)) {
    if (!length(penalties)) return(NULL)
    normalized <- lapply(penalties, function(penalty) {
        penalty <- as.matrix(penalty)
        scale <- max(abs(penalty))
        if (!is.finite(scale) || scale <= 0) penalty else penalty / scale
    })
    combined <- Reduce(`+`, normalized)
    combined <- (combined + t(combined)) / 2
    decomposition <- eigen(combined, symmetric=TRUE)
    scale <- max(1, max(abs(decomposition$values)))
    keep <- abs(decomposition$values) <= scale * tolerance
    decomposition$vectors[, keep, drop=FALSE]
}

.cdr_penalty_component <- function(term, penalty_index) {
    grouped <- !is.null(term$group)
    if (grouped && penalty_index == length(term$S)) {
        return('overall group-deviation magnitude')
    }
    type <- sub('_group$', '', term$type)
    if (penalty_index == 1L) return('lag curvature')
    if (penalty_index == 2L && identical(type, 'nonlinear')) {
        return('predictor curvature')
    }
    if (penalty_index == 2L && identical(type, 'varying')) {
        return('varying-covariate curvature')
    }
    paste('penalty', penalty_index)
}

.cdr_boundary_reduction_plan <- function(
        fit,
        design,
        log_sp_threshold=12,
        score_tolerance=2e-4
) {
    empty <- list(entries=list(), table=data.frame(
        term=character(),
        term_index=integer(),
        penalties=character(),
        components=character(),
        original_dimension=integer(),
        effective_dimension=integer(),
        maximum_log_sp=numeric(),
        maximum_score=numeric(),
        status=character(),
        stringsAsFactors=FALSE
    ))
    if (!inherits(fit, 'cdrgam_sparse') || !length(fit$sp) ||
            is.null(fit$optimizer$gradient)) return(empty)
    exact_score <- identical(fit$sparse$gradient, 'exact')
    log_sp <- log(as.numeric(fit$sp))
    score <- as.numeric(fit$optimizer$gradient)
    names(log_sp) <- names(fit$sp)
    names(score) <- names(fit$sp)
    candidate <- which(
        is.finite(log_sp) & log_sp >= log_sp_threshold &
            is.finite(score) & abs(score) <= score_tolerance
    )
    if (!length(candidate)) return(empty)
    expression <- regexec(
        '^s\\(cdr_term_([0-9]+)\\)([0-9]*)$',
        names(log_sp)[candidate]
    )
    matches <- regmatches(names(log_sp)[candidate], expression)
    parsed <- lapply(seq_along(matches), function(i) {
        value <- matches[[i]]
        if (!length(value)) return(NULL)
        list(
            global=candidate[[i]],
            term=as.integer(value[[2L]]),
            penalty=if (nzchar(value[[3L]])) {
                as.integer(value[[3L]])
            } else 1L
        )
    })
    parsed <- parsed[!vapply(parsed, is.null, logical(1))]
    entries <- list()
    reports <- list()
    parsed_global <- if (length(parsed)) {
        vapply(parsed, `[[`, integer(1), 'global')
    } else integer()
    unmatched <- setdiff(candidate, parsed_global)
    for (global in unmatched) {
        reports[[length(reports) + 1L]] <- data.frame(
            term=names(log_sp)[[global]],
            term_index=NA_integer_,
            penalties=names(log_sp)[[global]],
            components='ordinary smooth penalty',
            original_dimension=NA_integer_,
            effective_dimension=NA_integer_,
            maximum_log_sp=log_sp[[global]],
            maximum_score=abs(score[[global]]),
            status='diagnostic_only_non_irf',
            stringsAsFactors=FALSE
        )
    }
    if (!length(parsed)) {
        if (!length(reports)) return(empty)
        return(list(entries=entries, table=do.call(rbind, reports)))
    }
    for (term_index in unique(vapply(parsed, `[[`, integer(1), 'term'))) {
        if (term_index < 1L || term_index > length(design$terms)) next
        rows <- parsed[vapply(parsed, `[[`, integer(1), 'term') == term_index]
        penalty_indices <- unique(vapply(rows, `[[`, integer(1), 'penalty'))
        term <- design$terms[[term_index]]
        penalty_indices <- penalty_indices[
            penalty_indices >= 1L & penalty_indices <= length(term$S)
        ]
        # Full-rank penalties remove the whole term. Keep that scientifically
        # consequential decision diagnostic-only unless a future interface
        # explicitly requests term deletion.
        rank_deficient <- penalty_indices[vapply(
            penalty_indices,
            function(index) qr(as.matrix(term$S[[index]]))$rank <
                ncol(term$S[[index]]),
            logical(1)
        )]
        full_rank <- setdiff(penalty_indices, rank_deficient)
        if (length(full_rank)) {
            full_rows <- rows[vapply(
                rows,
                function(row) row$penalty %in% full_rank,
                logical(1)
            )]
            full_global <- vapply(full_rows, `[[`, integer(1), 'global')
            reports[[length(reports) + 1L]] <- data.frame(
                term=names(design$terms)[[term_index]],
                term_index=term_index,
                penalties=paste(names(log_sp)[full_global], collapse=', '),
                components=paste(vapply(
                    full_rows,
                    function(row) .cdr_penalty_component(term, row$penalty),
                    character(1)
                ), collapse=', '),
                original_dimension=ncol(term$X),
                effective_dimension=0L,
                maximum_log_sp=max(log_sp[full_global]),
                maximum_score=max(abs(score[full_global])),
                status='full_term_boundary_requires_confirmation',
                stringsAsFactors=FALSE
            )
        }
        penalty_indices <- rank_deficient
        if (!length(penalty_indices)) next
        reduced_rows <- rows[vapply(
            rows,
            function(row) row$penalty %in% penalty_indices,
            logical(1)
        )]
        global <- vapply(reduced_rows, `[[`, integer(1), 'global')
        basis <- .cdr_penalty_null_basis(term$S[penalty_indices])
        if (is.null(basis) || !ncol(basis) ||
                ncol(basis) >= ncol(term$X)) next
        entry <- list(
            term_index=term_index,
            penalty_indices=penalty_indices,
            basis=basis,
            original_dimension=ncol(term$X),
            effective_dimension=ncol(basis),
            smoothing_parameters=names(log_sp)[global],
            log_sp=unname(log_sp[global]),
            score=unname(score[global])
        )
        if (exact_score) entries[[length(entries) + 1L]] <- entry
        reports[[length(reports) + 1L]] <- data.frame(
            term=names(design$terms)[[term_index]],
            term_index=term_index,
            penalties=paste(entry$smoothing_parameters, collapse=', '),
            components=paste(vapply(
                penalty_indices,
                function(index) .cdr_penalty_component(term, index),
                character(1)
            ), collapse=', '),
            original_dimension=entry$original_dimension,
            effective_dimension=entry$effective_dimension,
            maximum_log_sp=max(entry$log_sp),
            maximum_score=max(abs(entry$score)),
            status=if (exact_score) {
                'certified_boundary_candidate'
            } else {
                'requires_exact_score_for_reduction'
            },
            stringsAsFactors=FALSE
        )
    }
    if (!length(reports)) return(empty)
    table <- do.call(rbind, reports)
    list(entries=entries, table=table)
}

.cdr_apply_boundary_reduction <- function(design, plan) {
    if (!length(plan$entries)) return(design)
    for (entry in plan$entries) {
        index <- entry$term_index
        term <- design$terms[[index]]
        original_penalty_count <- length(term$S)
        basis <- entry$basis
        retained <- setdiff(seq_along(term$S), entry$penalty_indices)
        term$X <- term$X %*% basis
        term$transform <- if (is.null(term$transform)) {
            basis
        } else {
            term$transform %*% basis
        }
        term$S <- lapply(retained, function(j) {
            crossprod(basis, as.matrix(term$S[[j]]) %*% basis)
        })
        if (!is.null(term$S.scale)) term$S.scale <- term$S.scale[retained]
        if (length(term$S)) {
            nonzero <- vapply(term$S, function(penalty) {
                max(abs(penalty)) > .Machine$double.eps
            }, logical(1))
            term$S <- term$S[nonzero]
            if (!is.null(term$S.scale)) term$S.scale <- term$S.scale[nonzero]
            retained <- retained[nonzero]
        }
        term$boundary_penalty_map <- data.frame(
            reduced=seq_along(retained),
            original=retained,
            stringsAsFactors=FALSE
        )
        term$boundary_original_penalty_count <- original_penalty_count
        repetitions <- if (is.null(term$group_levels)) {
            1L
        } else length(term$group_levels)
        term$rank <- if (length(term$S)) {
            repetitions * vapply(
                term$S,
                function(penalty) qr(penalty)$rank,
                integer(1)
            )
        } else integer()
        term$null.space.dim <- repetitions * if (length(term$S)) {
            ncol(term$X) - qr(Reduce(`+`, term$S))$rank
        } else ncol(term$X)
        if (!is.null(term$base_dimension)) {
            term$base_dimension <- ncol(term$X)
            term$expanded_dimension <- ncol(term$X) *
                length(term$group_levels)
        }
        term$constraints <- unique(c(
            term$constraints,
            'empirical-boundary-nullspace'
        ))
        term$boundary_reduction <- plan$table[
            plan$table$term_index == index,
            ,
            drop=FALSE
        ]
        design$terms[[index]] <- term
        if (!is.null(design$identifiability$constraints) &&
                length(design$identifiability$constraints) >= index) {
            design$identifiability$constraints[[index]] <- term$constraints
        }
        design$simplifications <- rbind(
            design$simplifications,
            data.frame(
                term=names(design$terms)[[index]],
                axis='penalty',
                action='boundary_nullspace_reduction',
                requested=as.character(entry$original_dimension),
                effective=as.character(entry$effective_dimension),
                reason=paste(
                    'REML boundary penalties:',
                    paste(entry$smoothing_parameters, collapse=', ')
                ),
                stringsAsFactors=FALSE
            )
        )
    }
    design$boundary_reductions <- plan$table
    design
}

.cdr_boundary_warm_start <- function(fit, reduced_design) {
    values <- log(as.numeric(fit$sp))
    names(values) <- names(fit$sp)
    for (i in seq_along(reduced_design$terms)) {
        term <- reduced_design$terms[[i]]
        map <- term$boundary_penalty_map
        original_count <- term$boundary_original_penalty_count
        if (is.null(map) || is.null(original_count) || !nrow(map)) next
        label <- paste0('s(cdr_term_', i, ')')
        original_names <- .sparse_sp_names(label, original_count)
        reduced_names <- .sparse_sp_names(label, nrow(map))
        source <- original_names[map$original]
        if (anyNA(match(source, names(values)))) {
            stop('Internal boundary warm-start mapping failed')
        }
        values[reduced_names] <- values[source]
    }
    values
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
        initial_log_sp=NULL,
        initial_source=NULL,
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
            'gradient', 'supernodal', 'trace_method', 'trace_chunk_size',
            'schur',
            'crossprod_chunk_size', 'restarts', 'gradient_probes',
            'gradient_cores', 'finite_difference_step', 'hessian',
            'hessian_step', 'outer_optimizer', 'optimizer_maxit',
            'optimizer_gradient_tolerance', 'optimizer_trust_radius',
            'boundary_action', 'boundary_log_sp'
        )
    )
    if (length(unknown_control)) {
        stop('Unknown sparse_control entries: ', paste(unknown_control, collapse=', '))
    }
    control <- function(name) sparse_control[[name, exact=TRUE]]
    gradient_requested <- if (is.null(control('gradient'))) {
        'auto'
    } else {
        match.arg(
            control('gradient'),
            c('auto', 'finite', 'exact', 'stochastic', 'hybrid')
        )
    }
    gradient_probes <- if (is.null(control('gradient_probes'))) {
        12L
    } else {
        value <- control('gradient_probes')
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value < 1 || value != as.integer(value)) {
            stop('sparse_control$gradient_probes must be a positive integer')
        }
        as.integer(value)
    }
    gradient_cores <- if (is.null(control('gradient_cores'))) {
        1L
    } else {
        value <- control('gradient_cores')
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value < 1 || value != as.integer(value)) {
            stop('sparse_control$gradient_cores must be a positive integer')
        }
        as.integer(value)
    }
    if (.Platform$OS.type == 'windows' && gradient_cores > 1L) {
        warning(
            'Parallel finite gradients are not available on Windows; using one core',
            call.=FALSE
        )
        gradient_cores <- 1L
    }
    finite_difference_step <- if (
        is.null(control('finite_difference_step'))
    ) {
        1e-3
    } else {
        value <- control('finite_difference_step')
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value <= 0) {
            stop('sparse_control$finite_difference_step must be positive')
        }
        value
    }
    hessian_requested <- if (is.null(control('hessian'))) {
        'auto'
    } else {
        match.arg(
            control('hessian'),
            c(
                'auto', 'gradient', 'analytic', 'profiled', 'optimhess',
                'defer', 'none'
            )
        )
    }
    hessian_step <- if (is.null(control('hessian_step'))) {
        1e-2
    } else {
        value <- control('hessian_step')
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value <= 0) {
            stop('sparse_control$hessian_step must be positive')
        }
        value
    }
    outer_optimizer_requested <- if (is.null(control('outer_optimizer'))) {
        'auto'
    } else {
        match.arg(
            control('outer_optimizer'),
            c('auto', 'lbfgsb', 'bfgs_trust')
        )
    }
    optimizer_maxit <- if (is.null(control('optimizer_maxit'))) {
        200L
    } else {
        value <- control('optimizer_maxit')
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value < 1 || value != as.integer(value)) {
            stop('sparse_control$optimizer_maxit must be a positive integer')
        }
        as.integer(value)
    }
    optimizer_gradient_tolerance <- if (is.null(
            control('optimizer_gradient_tolerance')
    )) {
        1e-4
    } else {
        value <- control('optimizer_gradient_tolerance')
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value <= 0) {
            stop(
                'sparse_control$optimizer_gradient_tolerance must be positive'
            )
        }
        value
    }
    optimizer_trust_radius <- if (is.null(
            control('optimizer_trust_radius')
    )) {
        2
    } else {
        value <- control('optimizer_trust_radius')
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
                value <= 0) {
            stop('sparse_control$optimizer_trust_radius must be positive')
        }
        value
    }
    boundary_action <- if (is.null(control('boundary_action'))) {
        'report'
    } else {
        match.arg(control('boundary_action'), c('report', 'reduce', 'error'))
    }
    boundary_log_sp <- if (is.null(control('boundary_log_sp'))) {
        12
    } else {
        value <- control('boundary_log_sp')
        if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
            stop('sparse_control$boundary_log_sp must be finite')
        }
        value
    }
    if (identical(outer_optimizer_requested, 'bfgs_trust') &&
            !(gradient_requested %in% c('auto', 'exact'))) {
        stop(
            'sparse_control$outer_optimizer="bfgs_trust" requires ',
            'gradient="exact" (or gradient="auto", which resolves to exact)'
        )
    }
    supernodal <- if (is.null(control('supernodal'))) {
        NULL
    } else {
        value <- control('supernodal')
        if (length(value) != 1L || !is.logical(value)) {
            stop('sparse_control$supernodal must be TRUE, FALSE, or NA')
        }
        value
    }
    trace_method_requested <- if (is.null(control('trace_method'))) {
        'auto'
    } else {
        match.arg(
            control('trace_method'),
            c('auto', 'solve', 'schur_inverse')
        )
    }
    gradient_method <- gradient_requested
    outer_optimizer <- outer_optimizer_requested
    trace_method <- trace_method_requested
    trace_chunk_size <- if (is.null(control('trace_chunk_size'))) {
        256L
    } else {
        value <- control('trace_chunk_size')
        if (length(value) != 1L || !is.numeric(value) ||
                !is.finite(value) || value < 1) {
            stop('sparse_control$trace_chunk_size must be a positive integer')
        }
        as.integer(value)
    }
    crossprod_chunk_size <- if (is.null(control('crossprod_chunk_size'))) {
        10000L
    } else {
        value <- control('crossprod_chunk_size')
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
    schur_method <- if (is.null(control('schur'))) {
        'never'
    } else {
        match.arg(control('schur'), c('never', 'always'))
    }
    restart_count <- if (is.null(control('restarts'))) {
        0L
    } else {
        value <- control('restarts')
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
                    if (is.factor(value)) levels(value) else
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
            levels <- levels(design$responses[[group]])
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
        sp_indices <- if (length(expanded_penalties)) {
            seq.int(first_sp, length(components))
        } else integer()
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
            axis=compact$axis,
            linear_predictors=compact$linear_predictors,
            linear_predictor_summaries=compact$linear_predictor_summaries,
            basis=compact$basis,
            transform=compact$transform,
            group=compact$group,
            group_levels=compact$group_levels,
            base_dimension=base_dimension,
            rank=compact$rank,
            null.space.dim=compact$null.space.dim,
            S.scale=compact$S.scale,
            lag_scale=if (is.null(compact$lag_scale)) 1 else
                compact$lag_scale,
            predictor_scale=if (is.null(compact$predictor_scale)) 1 else
                compact$predictor_scale,
            amplitude_scale=if (is.null(compact$amplitude_scale)) 1 else
                compact$amplitude_scale,
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
        optimizer_selection_policy=1L,
        outer_optimizer=outer_optimizer_requested,
        gradient=gradient_requested,
        trace_method=trace_method_requested,
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
    schur_trace_plan <- if (is.null(schur_layout)) {
        NULL
    } else {
        .sparse_schur_trace_plan(schur_layout, penalty_components)
    }
    if (identical(trace_method, 'schur_inverse') &&
            is.null(schur_trace_plan)) {
        stop(
            'sparse_control$trace_method="schur_inverse" requires a Schur ',
            'factorization with block-diagonal penalties'
        )
    }
    if (identical(trace_method, 'auto')) {
        trace_method <- if (is.null(schur_trace_plan)) {
            'solve'
        } else {
            'schur_inverse'
        }
    }
    active_schur_trace_plan <- if (identical(trace_method, 'schur_inverse')) {
        schur_trace_plan
    } else {
        NULL
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
    gradient_probe_matrix <- NULL
    numeric_update_count <- 0L
    parallel_gradient_factorizations <- 0L
    full_factorization_count <- 0L
    numeric_update_error <- NULL
    last_solution_key <- NULL
    last_solution <- NULL
    objective_seconds <- numeric()
    exact_gradient_seconds <- numeric()
    optimizer_selection <- NULL

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
            chunk_size=trace_chunk_size,
            schur_plan=active_schur_trace_plan
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
        elapsed_gradient <- proc.time()[['elapsed']] - gradient_started
        exact_gradient_seconds <<- c(
            utils::tail(exact_gradient_seconds, 99L),
            elapsed_gradient
        )
        optimizer_selection$measured_exact_seconds <<-
            exact_gradient_seconds
        optimizer_selection$median_exact_seconds <<-
            stats::median(exact_gradient_seconds)
        checkpoint_state$optimizer_selection <<- optimizer_selection
        if (reporter$level >= 2L) {
            reporter$emit(
                2L,
                'exact gradient',
                seconds=format(elapsed_gradient, digits=5),
                maximum=format(max(abs(score)), digits=5)
            )
        }
        score
    }

    curvature_convergence_assessment <- function(
            log_sp,
            criterion,
            gradient,
            objective_noise=0
    ) {
        method <- if (hessian_method %in% c('analytic', 'gradient')) {
            hessian_method
        } else {
            'analytic'
        }
        reporter$emit(
            1L,
            'curvature convergence assessment started',
            method=method,
            criterion=format(criterion, digits=10),
            projected_gradient=format(max(abs(gradient)), digits=5)
        )
        retained <- retained_solution(log_sp)
        if (is.null(retained)) {
            stop('Could not retain the current sparse solution')
        }
        curvature <- if (identical(method, 'analytic')) {
            analytic <- .sparse_analytic_profiled_hessian(
                retained$factor,
                blocks,
                retained$sp,
                penalty_components,
                penalty_supports,
                retained$coefficients,
                retained$penalized_rss,
                retained$reml_df,
                chunk_size=trace_chunk_size,
                element_limit=.sparse_analytic_hessian_element_limit(
                    hessian_requested,
                    hessian_selection
                )
            )
            list(
                hessian=analytic$hessian,
                evaluations=0L,
                rhs=analytic$rhs,
                retained_elements=analytic$retained_elements
            )
        } else {
            differentiated <- .central_difference_jacobian(
                exact_gradient,
                log_sp,
                step=hessian_step
            )
            list(
                hessian=differentiated$hessian,
                evaluations=differentiated$evaluations,
                rhs=NA_integer_,
                retained_elements=0
            )
        }
        projected <- gradient
        at_lower <- log_sp <= lower_bound + 1e-10
        at_upper <- log_sp >= upper_bound - 1e-10
        projected[at_lower & projected > 0] <- 0
        projected[at_upper & projected < 0] <- 0
        assessment <- .analytic_outer_convergence_assessment(
            curvature$hessian,
            projected,
            criterion,
            optimizer_gradient_tolerance,
            optimizer_trust_radius,
            objective_noise
        )
        assessment$message <- sub('^analytic', method, assessment$message)
        assessment$diagnostics$method <- method
        assessment$diagnostics$evaluations <- curvature$evaluations
        assessment$diagnostics$rhs <- curvature$rhs
        assessment$diagnostics$retained_elements <-
            curvature$retained_elements
        reporter$emit(
            1L,
            if (assessment$converged) {
                'curvature convergence certified'
            } else {
                'curvature recovery requested'
            },
            method=method,
            predicted_improvement=format(
                assessment$diagnostics$predicted_improvement,
                digits=5
            ),
            newton_step_max=format(
                assessment$diagnostics$newton_step_max,
                digits=5
            ),
            unresolved_gradient=format(
                assessment$diagnostics$unresolved_gradient,
                digits=5
            ),
            minimum_eigenvalue=format(
                assessment$diagnostics$minimum_eigenvalue,
                digits=5
            )
        )
        assessment
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
    initialization <- 'canonical unit smoothing parameters'
    checkpoint_initial <- if (!is.null(checkpoint_state$current_log_sp)) {
        checkpoint_state$current_log_sp
    } else {
        checkpoint_state$best_log_sp
    }
    if (length(checkpoint_initial) == penalty_count &&
            all(is.finite(checkpoint_initial))) {
        initial <- checkpoint_initial
        initialization <- 'checkpoint'
    } else if (!is.null(initial_log_sp)) {
        if (!is.numeric(initial_log_sp) || is.null(names(initial_log_sp)) ||
                anyDuplicated(names(initial_log_sp)) ||
                any(!is.finite(initial_log_sp))) {
            stop('initial_log_sp must be a finite named numeric vector')
        }
        positions <- match(sp_names, names(initial_log_sp))
        if (anyNA(positions)) {
            stop(
                'Initial smoothing parameters do not cover the fitted model: ',
                paste(sp_names[is.na(positions)], collapse=', '),
                call.=FALSE
            )
        }
        initial <- as.numeric(initial_log_sp[positions])
        initialization <- if (is.null(initial_source)) {
            'supplied warm start'
        } else {
            initial_source
        }
    }
    checkpoint_complete <- checkpoint_state$stage %in% c(
        'optimization_complete', 'complete'
    ) && !is.null(checkpoint_state$optimization)
    probe_value <- evaluate(
        initial,
        retain=FALSE,
        record=!checkpoint_complete
    )
    if (!is.finite(probe_value) || is.null(last_solution)) {
        stop('Could not evaluate the initial sparse REML objective')
    }
    optimizer_selection <- .sparse_optimizer_selection(
        gradient_requested=gradient_requested,
        outer_optimizer_requested=outer_optimizer_requested,
        trace_method=trace_method,
        factor=last_solution$factor,
        objective_seconds=utils::tail(objective_seconds, 1L),
        exact_trace_rhs=exact_trace_rhs,
        penalty_count=penalty_count,
        gradient_cores=gradient_cores,
        dimension=dimension,
        saved=checkpoint_state$optimizer_selection
    )
    gradient_method <- optimizer_selection$gradient
    outer_optimizer <- optimizer_selection$outer_optimizer
    exact_gradient_seconds <- optimizer_selection$measured_exact_seconds
    if (is.null(exact_gradient_seconds)) exact_gradient_seconds <- numeric()
    gradient_probe_matrix <- if (gradient_method %in% c(
        'stochastic', 'hybrid'
    )) {
        .deterministic_rademacher(dimension, gradient_probes)
    } else {
        NULL
    }
    checkpoint_state$optimizer_selection <- optimizer_selection
    checkpoint_state$updated_at <- as.character(Sys.time())
    if (!is.null(checkpoint)) .checkpoint_write(checkpoint_state, checkpoint)
    reporter$emit(
        1L,
        'optimizer selected',
        requested_gradient=gradient_requested,
        gradient=gradient_method,
        requested_outer_optimizer=outer_optimizer_requested,
        outer_optimizer=outer_optimizer,
        trace_method=trace_method,
        reason=optimizer_selection$reason,
        objective_seconds=format(
            optimizer_selection$objective_seconds,
            digits=5
        ),
        predicted_exact_seconds=format(
            optimizer_selection$predicted_exact_seconds,
            digits=5
        ),
        predicted_finite_seconds=format(
            optimizer_selection$predicted_finite_seconds,
            digits=5
        ),
        exact_core_mib=format(
            optimizer_selection$exact_core_bytes / 1024^2,
            digits=5
        )
    )
    reporter$emit(
        1L,
        'optimizer initialization',
        source=initialization,
        smoothing_parameters=penalty_count,
        minimum_log_sp=format(min(initial), digits=5),
        maximum_log_sp=format(max(initial), digits=5)
    )
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
    hessian_selection <- .sparse_select_hessian_method(
        hessian_requested,
        penalty_supports,
        dimension,
        trace_chunk_size
    )
    hessian_method <- hessian_selection$method
    expected_hessian_factorizations <- switch(
        hessian_method,
        gradient=2L * penalty_count,
        analytic=0L,
        profiled=1L + 2L * penalty_count^2L,
        optimhess=1L + 2L * (penalty_count + 1L)^2L,
        defer=0L,
        none=0L
    )
    active_restart <- 0L
    optimizer_progress_callback <- function(record) {
        optimizer_state <- record$optimizer_state
        record$optimizer_state <- NULL
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
        checkpoint_state$optimizer_state <<- optimizer_state
        checkpoint_state$optimizer_progress <<- record
        optimizer_history <- checkpoint_state$optimizer_history
        if (is.null(optimizer_history)) optimizer_history <- list()
        optimizer_history[[length(optimizer_history) + 1L]] <- record
        if (length(optimizer_history) > 1000L) {
            optimizer_history <- utils::tail(optimizer_history, 1000L)
        }
        checkpoint_state$optimizer_history <<- optimizer_history
        checkpoint_state$updated_at <<- as.character(Sys.time())
        if (!is.null(checkpoint)) {
            .checkpoint_write(checkpoint_state, checkpoint)
            checkpoint_last_written <<- evaluation_count
            checkpoint_last_written_time <<- proc.time()[['elapsed']]
        }
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
                progress=optimizer_progress_callback,
                state=arguments$optimizer_state,
                convergence_assessment=curvature_convergence_assessment
            )
        } else {
            do.call(stats::optim, arguments)
        }
        if (!identical(gradient_method, 'hybrid')) return(warm)
        reporter$emit(
            1L,
            'exact-score refinement started',
            criterion=format(warm$value, digits=10)
        )
        refinement_arguments <- arguments
        refinement_arguments$par <- warm$par
        refinement_arguments$gr <- exact_gradient
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
        trace_method=trace_method,
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
            same_checkpoint_run <- identical(
                checkpoint_state$current_restart,
                run_index - 1L
            ) && identical(checkpoint_state$stage, 'optimization')
            if (same_checkpoint_run &&
                    length(checkpoint_state$current_log_sp) == penalty_count &&
                    all(is.finite(checkpoint_state$current_log_sp))) {
                start <- checkpoint_state$current_log_sp
                resuming_interrupted_run <- resumed &&
                    !isTRUE(all.equal(start, canonical_start, tolerance=0))
            }
            checkpoint_state$current_restart <- run_index - 1L
            checkpoint_state$current_log_sp <- start
            checkpoint_state$stage <- 'optimization'
            saved_optimizer_state <- checkpoint_state$optimizer_state
            use_optimizer_state <- identical(
                outer_optimizer,
                'bfgs_trust'
            ) && same_checkpoint_run && is.list(saved_optimizer_state) &&
                length(saved_optimizer_state$parameters) == penalty_count &&
                isTRUE(all.equal(
                    as.numeric(saved_optimizer_state$parameters),
                    as.numeric(start),
                    tolerance=0
                ))
            if (!use_optimizer_state) checkpoint_state$optimizer_state <- NULL
            if (!is.null(checkpoint)) .checkpoint_write(checkpoint_state, checkpoint)
            reporter$emit(1L, 'restart started', restart=run_index - 1L)
            arguments <- optimization_arguments
            arguments$par <- start
            if (identical(outer_optimizer, 'bfgs_trust')) {
                arguments$optimizer_state <- if (use_optimizer_state) {
                    saved_optimizer_state
                } else {
                    NULL
                }
            }
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
        requested=hessian_requested,
        method=hessian_method
    )
    reporter$emit(
        1L,
        'Hessian method selected',
        requested=hessian_requested,
        method=hessian_method,
        reason=hessian_selection$reason,
        memory_source=hessian_selection$memory_source,
        available_mib=format(
            hessian_selection$memory_available_bytes / 1024^2,
            digits=7
        ),
        budget_mib=format(
            hessian_selection$memory_budget_bytes / 1024^2,
            digits=7
        ),
        analytic_peak_mib=format(
            hessian_selection$estimated_peak_bytes / 1024^2,
            digits=7
        )
    )
    checkpoint_hessian_step <- if (hessian_method %in% c(
            'gradient', 'profiled'
        )) {
        hessian_step
    } else {
        NA_real_
    }
    analytic <- NULL
    reuse_hessian <- !(hessian_method %in% c('defer', 'none')) &&
        identical(checkpoint_state$stage, 'complete') &&
        is.matrix(checkpoint_state$outer_hessian) &&
        identical(checkpoint_state$hessian_method, hessian_method) &&
        identical(checkpoint_state$hessian_step, checkpoint_hessian_step)
    if (reuse_hessian) {
        outer_hessian <- checkpoint_state$outer_hessian
        hessian_evaluations <- checkpoint_state$hessian_evaluations
        hessian_rhs <- checkpoint_state$hessian_rhs
        hessian_retained_elements <-
            checkpoint_state$hessian_retained_elements
        reporter$emit(1L, 'completed Hessian reused')
    } else if (hessian_method %in% c(
            'gradient', 'analytic', 'profiled'
        )) {
        # If U(rho, eta) is the unprofiled -2 REML criterion and eta is log
        # scale, profiling gives
        #   G''(rho) = U[rho,rho] - a a' / df,
        # where a_j = b' (lambda_j S_j) b / scale,
        # U[rho,eta] = -a, and U[eta,eta] = df. Reconstructing U'' from a
        # finite-difference Hessian of G avoids differencing the extra scale
        # dimension and therefore almost halves the required factorizations.
        profiled <- if (identical(hessian_method, 'gradient')) {
            .central_difference_jacobian(
                exact_gradient,
                optimization$par,
                step=hessian_step
            )
        } else if (identical(hessian_method, 'analytic')) {
            analytic <- .sparse_analytic_profiled_hessian(
                solution$factor,
                blocks,
                solution$sp,
                penalty_components,
                penalty_supports,
                solution$coefficients,
                solution$penalized_rss,
                solution$reml_df,
                chunk_size=trace_chunk_size,
                element_limit=.sparse_analytic_hessian_element_limit(
                    hessian_requested,
                    hessian_selection
                )
            )
            list(hessian=analytic$hessian, evaluations=0L)
        } else {
            .central_difference_hessian(
                function(log_sp) evaluate(
                    log_sp,
                    retain=FALSE,
                    record=FALSE
                ),
                optimization$par,
                step=hessian_step
            )
        }
        hessian_evaluations <- profiled$evaluations
        hessian_rhs <- if (is.null(analytic)) NA_integer_ else analytic$rhs
        hessian_retained_elements <- if (is.null(analytic)) {
            NA_real_
        } else {
            analytic$retained_elements
        }
        coefficients_at_solution <- solution$coefficients
        rss_scores <- if (!is.null(analytic)) {
            analytic$rss_scores
        } else {
            vapply(seq_len(penalty_count), function(i) {
                derivative <- solution$sp[[i]] * penalty_components[[i]]
                as.numeric(Matrix::crossprod(
                    coefficients_at_solution,
                    derivative %*% coefficients_at_solution
                ))
            }, numeric(1))
        }
        scale_scores <- rss_scores / scale
        unprofiled_rho <- profiled$hessian +
            tcrossprod(scale_scores) / solution$reml_df
        unprofiled_hessian <- rbind(
            cbind(unprofiled_rho, -scale_scores),
            c(-scale_scores, solution$reml_df)
        )
        outer_hessian <- unprofiled_hessian / 2
    } else if (identical(hessian_method, 'optimhess')) {
        hessian_evaluations <- NA_integer_
        hessian_rhs <- NA_integer_
        hessian_retained_elements <- NA_real_
        outer_hessian <- stats::optimHess(
            c(log(solution$sp), log(scale)),
            evaluate_unprofiled,
            control=list(ndeps=rep.int(
                hessian_step,
                penalty_count + 1L
            ))
        ) / 2
    } else {
        hessian_evaluations <- 0L
        hessian_rhs <- 0L
        hessian_retained_elements <- 0
        outer_hessian <- NULL
    }
    effective_df_at_solution <- if (!is.null(analytic)) {
        dimension - sum(analytic$system_scores)
    } else {
        NA_real_
    }
    hessian_eigenvalues <- if (is.null(outer_hessian)) {
        numeric()
    } else {
        tryCatch(
            eigen(
                outer_hessian,
                symmetric=TRUE,
                only.values=TRUE
            )$values,
            error=function(e) rep.int(NA_real_, nrow(outer_hessian))
        )
    }
    hessian_positive_definite <- if (!length(hessian_eigenvalues)) {
        NA
    } else {
        all(is.finite(hessian_eigenvalues)) &&
            min(hessian_eigenvalues) > sqrt(.Machine$double.eps)
    }
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
    gradient_norm <- if (gradient_method %in% c('exact', 'hybrid')) {
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
        hessian_requested=hessian_requested,
        hessian_method=hessian_method,
        hessian_selection=hessian_selection,
        hessian_evaluations=hessian_evaluations,
        hessian_rhs=hessian_rhs,
        restart_count=restart_count,
        restart_objectives=run_values,
        best_restart=best_run - 1L,
        resumed=resumed,
        checkpoint=checkpoint,
        global_optimum_certified=FALSE
    )
    if (!convergence$converged && identical(boundary_action, 'error')) {
        stop(
            'Sparse REML optimizer did not converge (code ',
            convergence$code, '): ', convergence$message,
            call.=FALSE
        )
    }
    if (!convergence$converged && identical(boundary_action, 'report')) {
        warning(
            'Sparse REML optimizer did not converge (code ',
            convergence$code, '): ', convergence$message,
            call.=FALSE
        )
    }
    if (identical(hessian_positive_definite, FALSE) &&
            !identical(boundary_action, 'reduce')) {
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
    deferred_hessian <- NULL
    if (identical(hessian_method, 'defer')) {
        deferred_hessian <- new.env(parent=emptyenv())
        deferred_hessian$hessian <- NULL
        deferred_hessian$evaluations <- 0L
        deferred_hessian$state <- list(
            system_template=system_template,
            base_system_values=base_system_values,
            penalty_mappings=penalty_mappings,
            schur_layout=schur_layout,
            supernodal=supernodal,
            Xty=Xty,
            weighted_response_sum_squares=as.numeric(
                Matrix::crossprod(weighted_y)
            ),
            blocks=blocks,
            observation_count=observation_count,
            dimension=dimension,
            penalty_components=penalty_components,
            penalty_supports=penalty_supports,
            trace_chunk_size=trace_chunk_size,
            log_sp=optimization$par,
            scale=scale,
            reml_df=solution$reml_df,
            coefficients=solution$coefficients,
            step=hessian_step
        )
    }

    reporter$phase('finalization')
    solver_label <- if (identical(outer_optimizer, 'bfgs_trust')) {
        paste(
            'sparse Gaussian REML solver',
            '(safeguarded trust-region BFGS)'
        )
    } else {
        'sparse Gaussian REML solver'
    }
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
            control=list(
                gradient_requested=gradient_requested,
                gradient=gradient_method,
                gradient_probes=gradient_probes,
                gradient_cores=gradient_cores,
                finite_difference_step=finite_difference_step,
                hessian_requested=hessian_requested,
                hessian=hessian_method,
                hessian_step=hessian_step,
                outer_optimizer_requested=outer_optimizer_requested,
                outer_optimizer=outer_optimizer,
                optimizer_maxit=optimizer_maxit,
                optimizer_gradient_tolerance=optimizer_gradient_tolerance,
                optimizer_trust_radius=optimizer_trust_radius,
                boundary_action=boundary_action,
                boundary_log_sp=boundary_log_sp,
                supernodal=supernodal,
                trace_method_requested=trace_method_requested,
                trace_method=trace_method,
                trace_chunk_size=trace_chunk_size,
                schur=schur_method,
                crossprod_chunk_size=crossprod_chunk_size,
                restarts=restart_count
            ),
            factor=solution$factor,
            dimension=dimension,
            nnzero=design_nonzeros,
            system_nnzero=Matrix::nnzero(solution$system),
            numeric_updates=numeric_update_count,
            parallel_gradient_factorizations=parallel_gradient_factorizations,
            full_factorizations=full_factorization_count,
            numeric_update_error=numeric_update_error,
            cached_objectives=length(ls(objective_cache, all.names=TRUE)),
            gradient_requested=gradient_requested,
            gradient=gradient_method,
            trace_method_requested=trace_method_requested,
            trace_method=trace_method,
            outer_optimizer_requested=outer_optimizer_requested,
            outer_optimizer=outer_optimizer,
            optimizer_selection=optimizer_selection,
            optimizer_progress=checkpoint_state$optimizer_progress,
            optimizer_history=checkpoint_state$optimizer_history,
            hessian_requested=hessian_requested,
            hessian=hessian_method,
            hessian_selection=hessian_selection,
            hessian_evaluations=hessian_evaluations,
            hessian_rhs=hessian_rhs,
            hessian_retained_elements=hessian_retained_elements,
            effective_df=effective_df_at_solution,
            deferred_hessian=deferred_hessian,
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
            penalty_trace=if (is.null(analytic)) NULL else
                analytic$system_scores,
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
                configuration=design$configuration,
                plan=design$plan,
                stream=design$stream,
                specification=design$specification,
                simplifications=design$simplifications,
                scaling=design$scaling,
                identifiability=design$identifiability
            ),
            scaling=design$scaling,
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
            solver=solver_label
        )
    )
    class(out) <- c('cdrgam_sparse', 'cdrgam')
    boundary_plan <- .cdr_boundary_reduction_plan(
        out,
        design,
        log_sp_threshold=boundary_log_sp,
        score_tolerance=optimizer_gradient_tolerance
    )
    out$cdrgam$boundary_reductions <- boundary_plan$table
    out$sparse$convergence$boundary_reductions <- boundary_plan$table
    out$sparse$boundary_action <- boundary_action
    out$sparse$boundary_log_sp <- boundary_log_sp
    checkpoint_state$stage <- 'complete'
    checkpoint_state$optimization <- optimization
    checkpoint_state$best_log_sp <- optimization$par
    checkpoint_state$best_criterion <- solution$criterion
    checkpoint_state$evaluation_count <- evaluation_count
    checkpoint_state$outer_hessian <- outer_hessian
    checkpoint_state$hessian_requested <- hessian_requested
    checkpoint_state$hessian_method <- hessian_method
    checkpoint_state$hessian_selection <- hessian_selection
    checkpoint_state$hessian_step <- checkpoint_hessian_step
    checkpoint_state$hessian_evaluations <- hessian_evaluations
    checkpoint_state$hessian_rhs <- hessian_rhs
    checkpoint_state$hessian_retained_elements <- hessian_retained_elements
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
        if (is.null(hessian)) {
            hessian <- .sparse_resolve_outer_hessian(object)
        }
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

.sparse_schur_inverse_diagonal <- function(factor, chunk_size=256L) {
    core_cholesky <- factor$core_cholesky
    core_inverse_factor <- backsolve(
        core_cholesky,
        diag(nrow(core_cholesky))
    )
    core_diagonal <- rowSums(core_inverse_factor^2)
    block_diagonal <- unlist(lapply(
        factor$block_cholesky,
        function(cholesky) {
            if (nrow(cholesky) == 1L) return(1 / cholesky[[1L]]^2)
            inverse_factor <- backsolve(cholesky, diag(nrow(cholesky)))
            rowSums(inverse_factor^2)
        }
    ), use.names=FALSE)
    if (!is.null(factor$transfer) && ncol(factor$transfer)) {
        for (start in seq.int(1L, ncol(factor$transfer), by=chunk_size)) {
            columns <- start:min(
                ncol(factor$transfer),
                start + chunk_size - 1L
            )
            transformed <- backsolve(
                core_cholesky,
                factor$transfer[, columns, drop=FALSE],
                transpose=TRUE
            )
            block_diagonal[columns] <- block_diagonal[columns] +
                colSums(transformed^2)
        }
    } else if (is.null(factor$transfer)) {
        offset <- 0L
        for (g in seq_along(factor$blocks)) {
            block_size <- length(factor$blocks[[g]])
            columns <- offset + seq_len(block_size)
            offset <- offset + block_size
            block_inverse <- chol2inv(factor$block_cholesky[[g]])
            local_transfer <- block_inverse %*% t(factor$cross[[g]])
            transformed <- backsolve(
                core_cholesky,
                t(local_transfer),
                transpose=TRUE
            )
            block_diagonal[columns] <- block_diagonal[columns] +
                colSums(transformed^2)
        }
    }
    output <- numeric(factor$dimension)
    output[factor$core] <- core_diagonal
    output[unlist(factor$blocks, use.names=FALSE)] <- block_diagonal
    output
}

.sparse_deferred_gradient_hessian <- function(state) {
    system_template <- state$system_template
    unit_values <- state$base_system_values
    for (mapping in state$penalty_mappings) {
        unit_values[mapping$positions] <-
            unit_values[mapping$positions] + mapping$values
    }
    unit_system <- system_template
    methods::slot(unit_system, 'x') <- unit_values
    factor_template <- if (is.null(state$schur_layout)) {
        Matrix::Cholesky(
            unit_system,
            LDL=FALSE,
            perm=TRUE,
            super=state$supernodal
        )
    } else {
        NULL
    }
    factor_system <- function(sp) {
        values <- state$base_system_values
        for (i in seq_along(state$penalty_mappings)) {
            mapping <- state$penalty_mappings[[i]]
            values[mapping$positions] <- values[mapping$positions] +
                sp[[i]] * mapping$values
        }
        system <- system_template
        methods::slot(system, 'x') <- values
        factor <- if (is.null(state$schur_layout)) {
            tryCatch(
                Matrix::update(factor_template, system),
                error=function(e) NULL
            )
        } else {
            .factor_schur_system(
                system,
                state$schur_layout,
                state$supernodal
            )
        }
        if (is.null(factor)) {
            factor <- Matrix::Cholesky(
                system,
                LDL=FALSE,
                perm=TRUE,
                super=state$supernodal
            )
        }
        list(system=system, factor=factor)
    }
    score <- function(log_sp) {
        sp <- exp(log_sp)
        factored <- factor_system(sp)
        penalty_det <- .sparse_penalty_logdet(state$blocks, sp)
        reml_df <- state$observation_count -
            (state$dimension - penalty_det$rank)
        coefficients <- as.numeric(.cdr_factor_solve(
            factored$factor,
            state$Xty
        ))
        penalized_rss <- state$weighted_response_sum_squares -
            as.numeric(Matrix::crossprod(coefficients, state$Xty))
        penalty_score <- .sparse_penalty_logdet_score(
            state$blocks,
            sp,
            length(sp)
        )
        determinant_scores <- .sparse_logdet_scores(
            factored$factor,
            state$penalty_components,
            sp,
            supports=state$penalty_supports,
            chunk_size=state$trace_chunk_size
        )
        vapply(seq_along(sp), function(i) {
            derivative <- sp[[i]] * state$penalty_components[[i]]
            rss_score <- as.numeric(Matrix::crossprod(
                coefficients,
                derivative %*% coefficients
            ))
            reml_df * rss_score / penalized_rss +
                determinant_scores[[i]] - penalty_score[[i]]
        }, numeric(1))
    }
    profiled <- .central_difference_jacobian(
        score,
        state$log_sp,
        step=state$step
    )
    rss_scores <- vapply(seq_along(state$penalty_components), function(i) {
        derivative <- exp(state$log_sp[[i]]) *
            state$penalty_components[[i]]
        as.numeric(Matrix::crossprod(
            state$coefficients,
            derivative %*% state$coefficients
        ))
    }, numeric(1))
    scale_scores <- rss_scores / state$scale
    unprofiled_rho <- profiled$hessian +
        tcrossprod(scale_scores) / state$reml_df
    hessian <- rbind(
        cbind(unprofiled_rho, -scale_scores),
        c(-scale_scores, state$reml_df)
    ) / 2
    list(hessian=hessian, evaluations=profiled$evaluations)
}

.sparse_resolve_outer_hessian <- function(object) {
    if (!is.null(object$outer.info$hess)) return(object$outer.info$hess)
    deferred <- object$sparse$deferred_hessian
    if (!is.environment(deferred)) {
        stop(
            'Smoothing-parameter covariance was not computed because ',
            'the model was fitted with sparse_control$hessian="none". ',
            'Refit with hessian="gradient" or hessian="defer" for ',
            'unconditional inference.'
        )
    }
    if (is.matrix(deferred$hessian)) return(deferred$hessian)
    if (is.null(deferred$state)) {
        stop('Deferred smoothing-parameter Hessian state is unavailable')
    }
    computed <- .sparse_deferred_gradient_hessian(deferred$state)
    deferred$hessian <- computed$hessian
    deferred$evaluations <- computed$evaluations
    deferred$state <- NULL
    deferred$hessian
}

.sparse_vcov_diagonal <- function(object, chunk_size=256L) {
    dimension <- object$sparse$dimension
    if (inherits(object$sparse$factor, 'cdrgam_schur_factor')) {
        output <- .sparse_schur_inverse_diagonal(
            object$sparse$factor,
            chunk_size=chunk_size
        ) * object$scale
    } else {
        output <- numeric(dimension)
        for (start in seq.int(1L, dimension, by=chunk_size)) {
            indices <- start:min(dimension, start + chunk_size - 1L)
            output[indices] <- diag(.sparse_selected_vcov(object, indices))
        }
    }
    stats::setNames(output, names(object$coefficients))
}

.sparse_effective_df <- function(object) {
    if (is.numeric(object$sparse$effective_df) &&
            length(object$sparse$effective_df) == 1L &&
            is.finite(object$sparse$effective_df)) {
        return(object$sparse$effective_df)
    }
    penalty_components <- object$sparse$penalty_components
    if (is.null(penalty_components)) return(NA_real_)
    if (inherits(object$sparse$factor, 'cdrgam_schur_factor')) {
        layout <- list(
            core=object$sparse$factor$core,
            blocks=object$sparse$factor$blocks
        )
        plan <- .sparse_schur_trace_plan(layout, penalty_components)
        if (!is.null(plan)) {
            scores <- .sparse_schur_logdet_scores(
                object$sparse$factor,
                plan,
                object$sp
            )
            return(object$sparse$dimension - sum(scores))
        }
    }
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

#' @rdname predict.cdrgam
#' @export
predict.cdrgam_sparse <- function(object, newdata=NULL, ...) {
    if (is.null(newdata)) return(object$fitted.values)
    if (!is.list(newdata) ||
            !all(c('impulses', 'responses') %in% names(newdata))) {
        stop('newdata must contain impulse and response data frames')
    }
    .predict_cdrgam_streams(
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
    reductions <- x$cdrgam$boundary_reductions
    if (!is.null(reductions) && nrow(reductions)) {
        applied <- sum(reductions$status == 'applied_and_reoptimized')
        cat('  boundary reductions:', nrow(reductions),
            'identified,', applied, 'applied\n')
    }
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
summary.cdrgam_sparse <- function(
        object,
        dispersion=NULL,
        freq=FALSE,
        re.test=TRUE,
        all.coefficients=FALSE,
        ...
) {
    out <- .cdrgam_custom_summary(
        object,
        covariance=function(indices) {
            .sparse_selected_vcov(object, indices)
        },
        smooth_edf=.cdrgam_sparse_smooth_edf(object),
        dispersion=dispersion,
        all.coefficients=all.coefficients,
        all_variances=function() .sparse_vcov_diagonal(object)
    )
    out$convergence <- object$sparse$convergence
    out$boundary_reductions <- object$cdrgam$boundary_reductions
    out$conditional_on_boundary_reduction <- isTRUE(
        object$cdrgam$conditional_on_boundary_reduction
    )
    class(out) <- c('summary.cdrgam_sparse', 'summary.cdrgam', 'summary.gam')
    out
}

#' @export
print.summary.cdrgam_sparse <- function(x, ...) {
    .print_summary_cdrgam(x, ...)
    cat('Backend:', x$backend, '\n')
    cat('Converged:', if (isTRUE(x$convergence$converged)) 'yes' else 'no', '\n')
    if (isTRUE(x$rank_metadata$parametric$corrected)) {
        cat(
            'Aliased parametric coefficients removed:',
            paste(x$rank_metadata$parametric$dropped, collapse=', '), '\n'
        )
    }
    if (!is.null(x$rank_metadata$resolution) &&
            !identical(x$rank_metadata$resolution, 'none')) {
        cat('Rank resolution:', x$rank_metadata$resolution, '\n')
    }
    if (!is.null(x$boundary_reductions) && nrow(x$boundary_reductions)) {
        cat('\nBoundary penalty reductions:\n')
        print(x$boundary_reductions, row.names=FALSE)
        if (isTRUE(x$conditional_on_boundary_reduction)) {
            cat('Inference is conditional on the applied boundary reductions.\n')
        }
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
