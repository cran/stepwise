#include <R.h>
#include <Rdefines.h>
#include <string.h>
SEXP Rmaxchi(SEXP argv, SEXP breaks, SEXP RWinHalfWidth, SEXP RpermReps);
SEXP Rphylpro(SEXP argv, SEXP breaks, SEXP RWinHalfWidth, SEXP RpermReps);


SEXP stepwise(SEXP argv, SEXP breaks, SEXP RWinHalfWidth, SEXP RpermReps) {
    char pargv[9];
    SEXP val;
    PROTECT(argv = AS_CHARACTER(argv));
    strcpy(pargv, CHAR(STRING_ELT(argv, 0)));

    if (strcmp("Rmaxchi", pargv) == 0) {
	val = (SEXP)Rmaxchi(argv, breaks, RWinHalfWidth, RpermReps);
    }
    else if (strcmp("Rphylpro", pargv) == 0) {
	val = (SEXP)Rphylpro(argv, breaks, RWinHalfWidth, RpermReps);
    }
    UNPROTECT(1);
    return(val);
}
