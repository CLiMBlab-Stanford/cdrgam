library(cdrgam)

chunk_crossprod <- cdrgam:::.sparse_chunk_crossprod
agree <- function(X) {
    fast <- chunk_crossprod(X)
    reference <- as.matrix(Matrix::crossprod(X))
    isTRUE(Matrix::isSymmetric(fast)) &&
        max(abs(as.matrix(fast) - reference)) <=
            1e-10 * max(1, max(abs(reference)))
}

set.seed(2027)
rows <- 400L
dense <- Matrix::Matrix(matrix(stats::rnorm(rows * 60L), rows), sparse=TRUE)
sparse <- Matrix::rsparsematrix(rows, 300L, density=0.01)
empty <- Matrix::Matrix(0, rows, 5L, sparse=TRUE)

stopifnot(
    # Dense population-style columns mixed with sparse ones, in either order
    # and with empty columns between them.
    agree(cbind(dense, sparse)),
    agree(cbind(sparse, empty, dense)),
    # Too few dense columns falls back to the sparse product.
    agree(cbind(dense[, 1:10], sparse)),
    # Every column dense, and a chunk with no dense columns at all.
    agree(dense),
    agree(sparse)
)
