#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>
#include <string.h>

static int matrix_nrow(SEXP x) {
    return INTEGER(getAttrib(x, R_DimSymbol))[0];
}

static int matrix_ncol(SEXP x) {
    return INTEGER(getAttrib(x, R_DimSymbol))[1];
}

static void check_info(int info, const char *operation) {
    if (info != 0) error("LAPACK %s failed with code %d", operation, info);
}

SEXP cdrgam_schur_factor(SEXP core, SEXP blocks, SEXP crosses) {
    int q = matrix_nrow(core);
    int group_count = length(blocks);
    int info;
    const double one = 1.0;
    const double minus_one = -1.0;

    if (matrix_ncol(core) != q || length(crosses) != group_count)
        error("Invalid Schur factor dimensions");

    SEXP schur = PROTECT(duplicate(core));
    SEXP cholesky = PROTECT(allocVector(VECSXP, group_count));

    for (int g = 0; g < group_count; ++g) {
        SEXP block = VECTOR_ELT(blocks, g);
        SEXP cross = VECTOR_ELT(crosses, g);
        int p = matrix_nrow(block);
        if (matrix_ncol(block) != p || matrix_nrow(cross) != q ||
                matrix_ncol(cross) != p)
            error("Invalid Schur block dimensions");

        SEXP chol = PROTECT(duplicate(block));
        F77_CALL(dpotrf)("U", &p, REAL(chol), &p, &info FCONE);
        check_info(info, "dpotrf");
        for (int column = 0; column < p; ++column)
            for (int row = column + 1; row < p; ++row)
                REAL(chol)[row + p * column] = 0.0;

        SEXP inverse_cross = PROTECT(allocMatrix(REALSXP, p, q));
        double *inverse = REAL(inverse_cross);
        double *between = REAL(cross);
        for (int column = 0; column < q; ++column)
            for (int row = 0; row < p; ++row)
                inverse[row + p * column] = between[column + q * row];
        F77_CALL(dpotrs)("U", &p, &q, REAL(chol), &p, inverse, &p,
                         &info FCONE);
        check_info(info, "dpotrs");

        F77_CALL(dgemm)("N", "N", &q, &q, &p, &minus_one,
                        between, &q, inverse, &p, &one,
                        REAL(schur), &q FCONE FCONE);
        SET_VECTOR_ELT(cholesky, g, chol);
        UNPROTECT(2);
    }

    SEXP core_cholesky = PROTECT(duplicate(schur));
    F77_CALL(dpotrf)("U", &q, REAL(core_cholesky), &q, &info FCONE);
    check_info(info, "Schur dpotrf");
    for (int column = 0; column < q; ++column)
        for (int row = column + 1; row < q; ++row)
            REAL(core_cholesky)[row + q * column] = 0.0;

    SEXP output = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(output, 0, core_cholesky);
    SET_VECTOR_ELT(output, 1, cholesky);
    UNPROTECT(4);
    return output;
}

SEXP cdrgam_schur_factor_sparse(SEXP system, SEXP core_indices,
                                 SEXP block_indices) {
    SEXP dimensions = R_do_slot(system, install("Dim"));
    SEXP column_pointers = R_do_slot(system, install("p"));
    SEXP row_indices = R_do_slot(system, install("i"));
    SEXP values = R_do_slot(system, install("x"));
    int n = INTEGER(dimensions)[0];
    int q = length(core_indices);
    int group_count = length(block_indices);

    int *core_position = (int *) R_alloc((size_t) n, sizeof(int));
    int *block_number = (int *) R_alloc((size_t) n, sizeof(int));
    int *block_position = (int *) R_alloc((size_t) n, sizeof(int));
    for (int j = 0; j < n; ++j) {
        core_position[j] = -1;
        block_number[j] = -1;
        block_position[j] = -1;
    }
    for (int j = 0; j < q; ++j)
        core_position[INTEGER(core_indices)[j] - 1] = j;
    for (int g = 0; g < group_count; ++g) {
        SEXP indices = VECTOR_ELT(block_indices, g);
        for (int j = 0; j < length(indices); ++j) {
            int global = INTEGER(indices)[j] - 1;
            block_number[global] = g;
            block_position[global] = j;
        }
    }

    SEXP core = PROTECT(allocMatrix(REALSXP, q, q));
    memset(REAL(core), 0, (size_t) q * q * sizeof(double));
    SEXP blocks = PROTECT(allocVector(VECSXP, group_count));
    SEXP crosses = PROTECT(allocVector(VECSXP, group_count));
    for (int g = 0; g < group_count; ++g) {
        int p = length(VECTOR_ELT(block_indices, g));
        SEXP block = PROTECT(allocMatrix(REALSXP, p, p));
        SEXP cross = PROTECT(allocMatrix(REALSXP, q, p));
        memset(REAL(block), 0, (size_t) p * p * sizeof(double));
        memset(REAL(cross), 0, (size_t) q * p * sizeof(double));
        SET_VECTOR_ELT(blocks, g, block);
        SET_VECTOR_ELT(crosses, g, cross);
        UNPROTECT(2);
    }

    int *pointers = INTEGER(column_pointers);
    int *rows = INTEGER(row_indices);
    double *x = REAL(values);
    for (int column = 0; column < n; ++column) {
        for (int entry = pointers[column]; entry < pointers[column + 1]; ++entry) {
            int row = rows[entry];
            double value = x[entry];
            int core_row = core_position[row];
            int core_column = core_position[column];
            int row_block = block_number[row];
            int column_block = block_number[column];
            if (core_row >= 0 && core_column >= 0) {
                REAL(core)[core_row + q * core_column] = value;
                REAL(core)[core_column + q * core_row] = value;
            } else if (row_block >= 0 && row_block == column_block) {
                SEXP block = VECTOR_ELT(blocks, row_block);
                int p = matrix_nrow(block);
                int local_row = block_position[row];
                int local_column = block_position[column];
                REAL(block)[local_row + p * local_column] = value;
                REAL(block)[local_column + p * local_row] = value;
            } else if (core_row >= 0 && column_block >= 0) {
                SEXP cross = VECTOR_ELT(crosses, column_block);
                REAL(cross)[core_row + q * block_position[column]] = value;
            } else if (core_column >= 0 && row_block >= 0) {
                SEXP cross = VECTOR_ELT(crosses, row_block);
                REAL(cross)[core_column + q * block_position[row]] = value;
            } else if (row != column) {
                error("Sparse Schur layout contains coupling between blocks");
            }
        }
    }

    SEXP dense_factor = PROTECT(cdrgam_schur_factor(core, blocks, crosses));
    SEXP output = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(output, 0, VECTOR_ELT(dense_factor, 0));
    SET_VECTOR_ELT(output, 1, VECTOR_ELT(dense_factor, 1));
    SET_VECTOR_ELT(output, 2, crosses);
    UNPROTECT(5);
    return output;
}

SEXP cdrgam_schur_solve(SEXP core_indices, SEXP block_indices,
                         SEXP crosses, SEXP block_cholesky,
                         SEXP core_cholesky, SEXP rhs) {
    int n = matrix_nrow(rhs);
    int nrhs = matrix_ncol(rhs);
    int q = length(core_indices);
    int group_count = length(block_indices);
    int info;
    const double one = 1.0;
    const double zero = 0.0;
    const double minus_one = -1.0;

    SEXP output = PROTECT(allocMatrix(REALSXP, n, nrhs));
    memset(REAL(output), 0, (size_t) n * nrhs * sizeof(double));
    SEXP adjusted_sexp = PROTECT(allocMatrix(REALSXP, q, nrhs));
    double *adjusted = REAL(adjusted_sexp);
    double *input = REAL(rhs);
    for (int column = 0; column < nrhs; ++column)
        for (int row = 0; row < q; ++row)
            adjusted[row + q * column] =
                input[(INTEGER(core_indices)[row] - 1) + n * column];

    SEXP block_solutions = PROTECT(allocVector(VECSXP, group_count));
    for (int g = 0; g < group_count; ++g) {
        SEXP indices = VECTOR_ELT(block_indices, g);
        SEXP cross = VECTOR_ELT(crosses, g);
        SEXP chol = VECTOR_ELT(block_cholesky, g);
        int p = length(indices);
        SEXP solved = PROTECT(allocMatrix(REALSXP, p, nrhs));
        for (int column = 0; column < nrhs; ++column)
            for (int row = 0; row < p; ++row)
                REAL(solved)[row + p * column] =
                    input[(INTEGER(indices)[row] - 1) + n * column];
        F77_CALL(dpotrs)("U", &p, &nrhs, REAL(chol), &p,
                         REAL(solved), &p, &info FCONE);
        check_info(info, "block dpotrs");
        F77_CALL(dgemm)("N", "N", &q, &nrhs, &p, &minus_one,
                        REAL(cross), &q, REAL(solved), &p, &one,
                        adjusted, &q FCONE FCONE);
        SET_VECTOR_ELT(block_solutions, g, solved);
        UNPROTECT(1);
    }

    F77_CALL(dpotrs)("U", &q, &nrhs, REAL(core_cholesky), &q,
                     adjusted, &q, &info FCONE);
    check_info(info, "core dpotrs");
    for (int column = 0; column < nrhs; ++column)
        for (int row = 0; row < q; ++row)
            REAL(output)[(INTEGER(core_indices)[row] - 1) + n * column] =
                adjusted[row + q * column];

    for (int g = 0; g < group_count; ++g) {
        SEXP indices = VECTOR_ELT(block_indices, g);
        SEXP cross = VECTOR_ELT(crosses, g);
        SEXP chol = VECTOR_ELT(block_cholesky, g);
        SEXP solved = VECTOR_ELT(block_solutions, g);
        int p = length(indices);
        SEXP correction = PROTECT(allocMatrix(REALSXP, p, nrhs));
        F77_CALL(dgemm)("T", "N", &p, &nrhs, &q, &one,
                        REAL(cross), &q, adjusted, &q, &zero,
                        REAL(correction), &p FCONE FCONE);
        F77_CALL(dpotrs)("U", &p, &nrhs, REAL(chol), &p,
                         REAL(correction), &p, &info FCONE);
        check_info(info, "correction dpotrs");
        for (int column = 0; column < nrhs; ++column)
            for (int row = 0; row < p; ++row)
                REAL(output)[(INTEGER(indices)[row] - 1) + n * column] =
                    REAL(solved)[row + p * column] -
                    REAL(correction)[row + p * column];
        UNPROTECT(1);
    }

    UNPROTECT(3);
    return output;
}
