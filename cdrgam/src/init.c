#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP cdrgam_schur_factor(SEXP, SEXP, SEXP);
extern SEXP cdrgam_schur_factor_sparse(SEXP, SEXP, SEXP);
extern SEXP cdrgam_schur_solve(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cdrgam_schur_solve_batched(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef call_methods[] = {
    {"cdrgam_schur_factor", (DL_FUNC) &cdrgam_schur_factor, 3},
    {"cdrgam_schur_factor_sparse", (DL_FUNC) &cdrgam_schur_factor_sparse, 3},
    {"cdrgam_schur_solve", (DL_FUNC) &cdrgam_schur_solve, 6},
    {"cdrgam_schur_solve_batched", (DL_FUNC) &cdrgam_schur_solve_batched, 6},
    {NULL, NULL, 0}
};

void attribute_visible R_init_cdrgam(DllInfo *dll) {
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
