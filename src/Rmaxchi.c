#include <stdlib.h>
#include "init.h"
/* #include "getMaxChi.h" */
#include <R.h>
#include <Rdefines.h>

SEXP Rmaxchi(SEXP argv, SEXP breaks, SEXP RWinHalfWidth, SEXP RpermReps)
{
    int i, *pairmem1, *pairmem2;
    double *quants, *chisqs; 
    int numsig=0, *winlocs, *polyposn;
    char *pargv[2];
    int *pbreaks;
    int val_len = 6;
    SEXP Rpolyposn, Rchisqs, Rwinlocs, Rpairmem1, Rpairmem2, Rquants, val, val_nm;
    double *pRchisqs;
    int *pRpolyposn, *pRwinlocs;

    PROTECT(argv = AS_CHARACTER(argv));
    pargv[0] = R_alloc(strlen(CHAR(STRING_ELT(argv, 0))), sizeof(char));
    pargv[1] = R_alloc(strlen(CHAR(STRING_ELT(argv, 1))), sizeof(char));
    strcpy(pargv[0], CHAR(STRING_ELT(argv, 0)));
    strcpy(pargv[1], CHAR(STRING_ELT(argv, 1)));

    PROTECT(breaks = AS_INTEGER(breaks));
    pbreaks = INTEGER_POINTER(breaks);
    PROTECT(RWinHalfWidth = AS_INTEGER(RWinHalfWidth));
    winHalfWidth = INTEGER_POINTER(RWinHalfWidth)[0];
    PROTECT(RpermReps = AS_INTEGER(RpermReps));
    permReps = INTEGER_POINTER(RpermReps)[0];

    numsig =  main_maxchi(2, pargv, pbreaks, &polyposn, &winlocs, &chisqs, 
&pairmem1, &pairmem2, &quants);

    UNPROTECT(4);
    if (numsig == 0) return(R_NilValue);

    PROTECT(Rchisqs = NEW_NUMERIC(numsig));
    PROTECT(Rwinlocs = NEW_INTEGER(numsig));

    PROTECT(Rpairmem1 = allocVector(STRSXP,numsig));
    PROTECT(Rpairmem2 = allocVector(STRSXP,numsig));
    pRchisqs = NUMERIC_POINTER(Rchisqs);
    pRwinlocs = INTEGER_POINTER(Rwinlocs);
    for(i = 0; i < numsig; i++) {
	pRchisqs[i] = chisqs[i];
	pRwinlocs[i] = winlocs[i];
        SET_STRING_ELT(Rpairmem1, i, mkChar(sequenceLabels[pairmem1[i]]));
        SET_STRING_ELT(Rpairmem2, i, mkChar(sequenceLabels[pairmem2[i]]));
    }   
    PROTECT(Rpolyposn = NEW_INTEGER(numBases));
    pRpolyposn = INTEGER_POINTER(Rpolyposn);
    for(i = 0; i < numBases; i++) 
	pRpolyposn[i] = polyposn[i] +1;
    PROTECT(Rquants = NEW_NUMERIC(2));
    NUMERIC_POINTER(Rquants)[0] = quants[0];
    NUMERIC_POINTER(Rquants)[1] = quants[1];

    PROTECT(val_nm = allocVector(STRSXP,6));
    SET_STRING_ELT(val_nm, 0, mkChar("polyposn"));
    SET_STRING_ELT(val_nm, 1, mkChar("chisqs"));
    SET_STRING_ELT(val_nm, 2, mkChar("winlocs"));
    SET_STRING_ELT(val_nm, 3, mkChar("pairmem1"));
    SET_STRING_ELT(val_nm, 4, mkChar("pairmem2"));
    SET_STRING_ELT(val_nm, 5, mkChar("quants"));

    PROTECT(val = allocVector(VECSXP, val_len));
    SET_VECTOR_ELT(val, 0, Rpolyposn);
    SET_VECTOR_ELT(val, 1, Rchisqs);
    SET_VECTOR_ELT(val, 2, Rwinlocs);
    SET_VECTOR_ELT(val, 3, Rpairmem1);
    SET_VECTOR_ELT(val, 4, Rpairmem2);
    SET_VECTOR_ELT(val, 5, Rquants);
    setAttrib(val, R_NamesSymbol, val_nm);

    UNPROTECT(8);
    return(val);
}
