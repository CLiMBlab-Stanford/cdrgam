library(cdrgam)

positive_integer <- cdrgam:::.cdrgam_positive_integer
available_cores <- cdrgam:::.cdrgam_available_cores
parallel_plan <- cdrgam:::.cdrgam_parallel_plan
score_parallel_plan <- cdrgam:::.cdrgam_score_parallel_plan
batch_plan <- cdrgam:::.cdrgam_sparse_score_batch_plan
trace_chunk_size <- cdrgam:::.cdrgam_sparse_trace_chunk_size
with_blas_threads <- cdrgam:::.cdrgam_with_blas_threads
blas_threads <- cdrgam:::.cdrgam_blas_threads

stopifnot(
    identical(positive_integer(3, 'test value'), 3L),
    identical(available_cores(5), 5L),
    available_cores() >= 1L
)

batched <- batch_plan(10L, 4L, derivative_bytes=1024, batch_size=3L)
single_batch <- batch_plan(10L, 4L, derivative_bytes=1024, batch_size=10L)
bounded_trace_chunk <- trace_chunk_size(
    1000L, 2L,
    memory=list(available_bytes=64 * 1024^2)
)
automatic_score <- score_parallel_plan(8L, 40L)
stopifnot(
    batched$workers == 4L,
    single_batch$workers == 4L,
    length(single_batch$groups) == 1L,
    identical(unname(lengths(batched$groups)), c(3L, 3L, 3L, 1L)),
    batched$batch_size == 3L,
    batched$derivative_bytes == 1024,
    bounded_trace_chunk >= 64L,
    bounded_trace_chunk <= 1024L,
    log(bounded_trace_chunk, base=2) %% 1 == 0,
    automatic_score$workers == if (.Platform$OS.type == 'windows') 1L else 2L,
    automatic_score$blas_threads == if (
        .Platform$OS.type == 'windows'
    ) 8L else 4L,
    identical(automatic_score$worker_source, 'automatic')
)

many_tasks <- parallel_plan(8L, 40L)
few_tasks <- parallel_plan(8L, 2L)
explicit <- parallel_plan(8L, 40L, 2L)
if (.Platform$OS.type == 'windows') {
    stopifnot(
        many_tasks$workers == 1L,
        many_tasks$blas_threads == 8L,
        few_tasks$workers == 1L,
        few_tasks$blas_threads == 8L
    )
} else {
    stopifnot(
        many_tasks$workers == 8L,
        many_tasks$blas_threads == 1L,
        few_tasks$workers == 2L,
        few_tasks$blas_threads == 4L,
        explicit$workers == 2L,
        explicit$blas_threads == 4L
    )
}
stopifnot(
    many_tasks$workers * many_tasks$blas_threads <= many_tasks$cores,
    few_tasks$workers * few_tasks$blas_threads <= few_tasks$cores,
    identical(explicit$worker_source, 'explicit')
)

initial_threads <- blas_threads()
observed_threads <- with_blas_threads(1L, blas_threads())
stopifnot(observed_threads == 1L, blas_threads() == initial_threads)

if (.Platform$OS.type != 'windows') {
    fork_plan <- parallel_plan(4L, 2L)
    worker_threads <- with_blas_threads(
        fork_plan$blas_threads,
        unlist(parallel::mclapply(
            seq_len(fork_plan$workers),
            function(index) cdrgam:::.cdrgam_blas_threads(),
            mc.cores=fork_plan$workers,
            mc.set.seed=FALSE
        ))
    )
    stopifnot(
        all(worker_threads == fork_plan$blas_threads),
        blas_threads() == initial_threads
    )
}

cat('Parallel phase-planning checks passed.\n')
