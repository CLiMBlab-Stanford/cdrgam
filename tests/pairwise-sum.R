library(cdrgam)

# The pairwise accumulator returns the same sum as adding terms one at a time,
# for counts that exercise every carry pattern of the binary counter.
set.seed(7311)
term <- function() Matrix::rsparsematrix(40, 40, density=0.05)
for (count in c(0L, 1L, 2L, 3L, 7L, 8L, 13L)) {
    terms <- replicate(count, term(), simplify=FALSE)
    accumulator <- cdrgam:::.sparse_pairwise_sum()
    for (value in terms) accumulator$add(value)
    total <- accumulator$result()
    if (!count) {
        stopifnot(is.null(total))
        next
    }
    stopifnot(
        inherits(total, 'Matrix'),
        max(abs(as.matrix(total - Reduce(`+`, terms)))) < 1e-12
    )
}

# Chunked cross products summed pairwise equal the full cross product.
X <- Matrix::rsparsematrix(500, 60, density=0.1)
accumulator <- cdrgam:::.sparse_pairwise_sum()
for (start in seq.int(1L, 500L, by=37L)) {
    rows <- start:min(500L, start + 36L)
    accumulator$add(Matrix::crossprod(X[rows, , drop=FALSE]))
}
stopifnot(
    max(abs(as.matrix(accumulator$result() - Matrix::crossprod(X)))) < 1e-10
)
