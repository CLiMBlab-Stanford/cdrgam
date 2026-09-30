.cdrgam_positive_integer <- function(value, field) {
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
            value < 1 || value != as.integer(value)) {
        stop(field, ' must be a positive integer')
    }
    as.integer(value)
}

.cdrgam_available_cores <- function(value=NULL) {
    if (!is.null(value)) {
        return(.cdrgam_positive_integer(value, 'sparse_control$cores'))
    }
    detected <- suppressWarnings(parallelly::availableCores(omit=0L))
    if (!length(detected) || !is.finite(detected[[1L]]) ||
            detected[[1L]] < 1) return(1L)
    as.integer(detected[[1L]])
}

.cdrgam_parallel_plan <- function(cores, tasks, workers=NULL) {
    cores <- .cdrgam_positive_integer(cores, 'core budget')
    tasks <- .cdrgam_positive_integer(tasks, 'parallel task count')
    forked <- .Platform$OS.type != 'windows'
    requested_workers <- if (is.null(workers)) {
        if (forked) min(cores, tasks) else 1L
    } else {
        .cdrgam_positive_integer(workers, 'parallel worker count')
    }
    resolved_workers <- min(requested_workers, cores, tasks)
    if (!forked) resolved_workers <- 1L
    list(
        cores=cores,
        workers=resolved_workers,
        blas_threads=max(1L, cores %/% resolved_workers),
        forked=forked,
        worker_source=if (is.null(workers)) 'automatic' else 'explicit'
    )
}

.cdrgam_score_parallel_plan <- function(cores, tasks, workers=NULL) {
    automatic <- is.null(workers)
    if (automatic) {
        workers <- max(1L, floor(sqrt(as.double(cores))))
    }
    plan <- .cdrgam_parallel_plan(cores, tasks, workers)
    plan$worker_source <- if (automatic) 'automatic' else 'explicit'
    plan
}

.cdrgam_memory_parallel_plan <- function(
        cores, tasks, workers=NULL, per_worker_bytes,
        memory_fraction=0.5, reserve_bytes=512 * 1024^2,
        memory=.cdrgam_memory_availability()
) {
    plan <- .cdrgam_parallel_plan(cores, tasks, workers)
    requested_workers <- plan$workers
    if (length(per_worker_bytes) != 1L || !is.numeric(per_worker_bytes) ||
            !is.finite(per_worker_bytes) || per_worker_bytes <= 0) {
        stop('per_worker_bytes must be one positive finite number')
    }
    if (length(memory_fraction) != 1L || !is.numeric(memory_fraction) ||
            !is.finite(memory_fraction) || memory_fraction <= 0 ||
            memory_fraction > 1) {
        stop('memory_fraction must lie in (0, 1]')
    }
    if (length(reserve_bytes) != 1L || !is.numeric(reserve_bytes) ||
            !is.finite(reserve_bytes) || reserve_bytes < 0) {
        stop('reserve_bytes must be one nonnegative finite number')
    }
    available <- memory$available_bytes
    budget <- if (is.finite(available)) {
        max(0, min(memory_fraction * available, available - reserve_bytes))
    } else NA_real_
    memory_workers <- if (is.finite(budget)) {
        max(1L, floor(budget / per_worker_bytes))
    } else plan$workers
    resolved_workers <- min(plan$workers, memory_workers)
    plan$workers <- as.integer(resolved_workers)
    plan$blas_threads <- max(1L, plan$cores %/% plan$workers)
    plan$requested_workers <- as.integer(requested_workers)
    plan$per_worker_bytes <- as.numeric(per_worker_bytes)
    plan$memory_workers <- as.integer(memory_workers)
    plan$memory_limited <- plan$workers < requested_workers
    plan$memory_fraction <- memory_fraction
    plan$memory_reserve_bytes <- reserve_bytes
    plan$memory_budget_bytes <- budget
    plan$memory_source <- memory$source
    plan$memory_available_bytes <- available
    plan
}

.cdrgam_sparse_worker_bytes <- function(
        system, factor, transient_bytes=0
) {
    system_nonzeros <- as.double(Matrix::nnzero(system))
    factor_nonzeros <- as.double(.cdr_factor_nonzeros(factor))
    if (length(transient_bytes) != 1L || !is.numeric(transient_bytes) ||
            !is.finite(transient_bytes) || transient_bytes < 0) {
        stop('transient_bytes must be one nonnegative finite number')
    }
    max(
        1,
        48 * system_nonzeros + 64 * factor_nonzeros + transient_bytes
    )
}

.cdrgam_sparse_score_batch_plan <- function(
        task_count, workers, derivative_bytes, batch_size=NULL
) {
    task_count <- .cdrgam_positive_integer(task_count, 'score task count')
    workers <- min(
        .cdrgam_positive_integer(workers, 'score worker count'),
        task_count
    )
    maximum_batch <- task_count
    derivative_bytes <- max(1, as.double(derivative_bytes))
    memory <- .cdrgam_memory_availability()
    memory_batch <- if (is.finite(memory$available_bytes)) {
        max(1L, floor(
            0.55 * memory$available_bytes / derivative_bytes
        ))
    } else 1L
    resolved_batch_size <- if (is.null(batch_size)) {
        min(maximum_batch, memory_batch)
    } else {
        min(maximum_batch, .cdrgam_positive_integer(
            batch_size, 'sparse_control$score_batch_size'
        ))
    }
    groups <- split(
        seq_len(task_count),
        ceiling(seq_len(task_count) / resolved_batch_size)
    )
    list(
        workers=min(workers, task_count),
        groups=groups,
        batch_size=resolved_batch_size,
        derivative_bytes=derivative_bytes,
        memory=memory
    )
}

.cdrgam_sparse_trace_chunk_size <- function(
        dimension, workers, memory=.cdrgam_memory_availability()
) {
    dimension <- .cdrgam_positive_integer(dimension, 'trace dimension')
    workers <- .cdrgam_positive_integer(workers, 'trace worker count')
    if (!is.finite(memory$available_bytes)) return(256L)
    candidate <- floor(
        0.1 * memory$available_bytes / (16 * dimension * workers)
    )
    candidate <- min(1024L, max(64L, as.integer(candidate)))
    as.integer(2^floor(log(candidate, base=2)))
}

.cdrgam_blas_threads <- function() {
    value <- suppressWarnings(RhpcBLASctl::blas_get_num_procs())
    if (!length(value) || !is.finite(value[[1L]]) || value[[1L]] < 1) {
        return(1L)
    }
    as.integer(value[[1L]])
}

.cdrgam_with_blas_threads <- function(threads, code) {
    threads <- .cdrgam_positive_integer(threads, 'BLAS thread count')
    previous <- .cdrgam_blas_threads()
    on.exit(RhpcBLASctl::blas_set_num_threads(previous), add=TRUE)
    RhpcBLASctl::blas_set_num_threads(threads)
    force(code)
}
