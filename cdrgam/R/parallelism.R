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
        workers=min(workers, length(groups)),
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
