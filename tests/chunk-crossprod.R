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

# Columns that share no row give exact zeros in the dense blocks; these are
# not stored, so the pattern matches a plain sparse cross product. Here two
# groups of dense columns occupy disjoint rows, as for two subjects in one
# chunk.
first <- Matrix::Matrix(rbind(
    matrix(stats::rnorm(200L * 40L), 200L),
    matrix(0, 200L, 40L)
), sparse=TRUE)
second <- Matrix::Matrix(rbind(
    matrix(0, 200L, 40L),
    matrix(stats::rnorm(200L * 40L), 200L)
), sparse=TRUE)
disjoint <- cbind(first, second)
disjoint_product <- chunk_crossprod(disjoint)
stopifnot(
    agree(disjoint),
    all(methods::slot(disjoint_product, 'x') != 0),
    length(methods::slot(disjoint_product, 'x')) ==
        length(methods::slot(
            Matrix::forceSymmetric(Matrix::crossprod(disjoint), uplo='U'),
            'x'
        ))
)

# Row extraction through a transposed source matches direct subscripting for
# ungrouped sparse components, dense matrices, and grouped components.
row_source <- cdrgam:::.sparse_component_row_source
component_rows <- cdrgam:::.sparse_component_rows
ungrouped <- Matrix::rsparsematrix(500L, 40L, density=0.05)
dense_matrix <- matrix(stats::rnorm(500L * 3L), 500L)
grouped <- function(base) structure(list(
    base=base,
    group_index=rep(1:5, each=100L),
    group_count=5L,
    base_dimension=ncol(base)
), class='cdrgam_sparse_grouped_component')
same_rows <- function(component, rows) {
    expected <- as.matrix(component_rows(component, rows))
    observed <- as.matrix(component_rows(row_source(component), rows))
    identical(dim(expected), dim(observed)) && identical(expected, observed)
}
for (rows in list(1:120, 121:500, 250L)) {
    stopifnot(
        same_rows(ungrouped, rows),
        same_rows(dense_matrix, rows),
        same_rows(grouped(dense_matrix), rows)
    )
}

# A chunk assembled over its nonempty columns, with its cross product placed
# back at the global indices, matches the cross product of the full chunk.
set.seed(11)
components <- list(
    Matrix::rsparsematrix(500L, 30L, density=0.05),
    grouped(matrix(stats::rnorm(500L * 3L), 500L)),
    Matrix::Matrix(cbind(1, stats::rnorm(500L)), sparse=TRUE),
    grouped(matrix(stats::rnorm(500L * 2L), 500L))
)
offsets <- cumsum(c(0, vapply(
    components, cdrgam:::.sparse_component_ncol, numeric(1)
)))[seq_along(components)]
total_columns <- sum(vapply(components, cdrgam:::.sparse_component_ncol, numeric(1)))
sources <- lapply(components, row_source)
weights <- stats::runif(500L)
for (range in list(c(1L, 120L), c(180L, 260L), c(401L, 500L))) {
    rows <- range[[1L]]:range[[2L]]
    full <- do.call(cbind, lapply(components, component_rows, rows=rows))
    full <- Matrix::Diagonal(x=weights[rows]) %*% full
    reference <- as.matrix(Matrix::crossprod(full))
    chunk <- cdrgam:::.sparse_compact_chunk(
        sources, offsets, range[[1L]], range[[2L]], weights[rows]
    )
    expanded <- cdrgam:::.sparse_expand_symmetric(
        chunk_crossprod(chunk$X), chunk$columns, total_columns
    )
    stopifnot(
        ncol(chunk$X) < total_columns,
        max(abs(as.matrix(expanded) - reference)) < 1e-10 * max(1, abs(reference))
    )
}
