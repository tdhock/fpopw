#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void colibri_op_R_c(void *, void *, void *, void *, void *, void *, void *);
extern void colibri_op_weights_R_c(void *, void *, void *, void *, void *, void *, void *, void *);
extern void colibri_sn_R_c(void *, void *, void *, void *, void *, void *, void *, void *);
extern void colibri_sn_weights_nomemory_R_c(void *, void *, void *, void *, void *, void *, void *);
extern void colibri_sn_weights_R_c(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"colibri_op_R_c",                  (DL_FUNC) &colibri_op_R_c,                  7},
    {"colibri_op_weights_R_c",          (DL_FUNC) &colibri_op_weights_R_c,          8},
    {"colibri_sn_R_c",                  (DL_FUNC) &colibri_sn_R_c,                  8},
    {"colibri_sn_weights_nomemory_R_c", (DL_FUNC) &colibri_sn_weights_nomemory_R_c, 7},
    {"colibri_sn_weights_R_c",          (DL_FUNC) &colibri_sn_weights_R_c,          9},
    {NULL, NULL, 0}
};

void R_init_fpopw(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
