#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP _RGraphM_run_graph_match(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_RGraphM_run_graph_match", (DL_FUNC) &_RGraphM_run_graph_match, 3},
    {NULL, NULL, 0}
};

void R_init_RGraphM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
