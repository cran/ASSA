#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(dbar_)(void *, void *, void *, void *);
extern void F77_NAME(autocov_)(void *, void *, void *);
extern void F77_NAME(yi_)(void *, void *, void *,void *,void *, void *, void *);
extern void F77_NAME(yyi_)(void *, void *, void *,void *,void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"dbar_", (DL_FUNC) &F77_NAME(dbar_), 4},
    {"autocov_", (DL_FUNC) &F77_NAME(autocov_), 3},
    {"yi_", (DL_FUNC) &F77_NAME(yi_), 7},
    {"yyi_", (DL_FUNC) &F77_NAME(yyi_), 7},
    {NULL, NULL, 0}
};

void R_init_ASSA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
