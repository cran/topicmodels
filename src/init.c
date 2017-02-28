#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP rGibbslda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
		      SEXP, SEXP, SEXP, SEXP);
extern SEXP rctm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rlda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"rGibbslda", (DL_FUNC) &rGibbslda, 13},
    {"rctm",      (DL_FUNC) &rctm,       9},
    {"rlda",      (DL_FUNC) &rlda,       9},
    {NULL, NULL, 0}
};

void R_init_topicmodels(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
