#include <stdlib.h>
#include "init.h"
#include "getPhylpro.h"
#include <R.h>
#include <Rdefines.h>

SEXP Rphylpro(SEXP argv, SEXP breaks, SEXP RWinHalfWidth, SEXP RpermReps) {

    int numsig=0, val_len=5, i, *winlocs, *polyposn, *pbreaks, *targetseqs, *pRpolyposn, *pRwinlocs;
    char *pargv[2], pargvSize[2][15]; 
    double *quants, *corrs, *pRcorrs, *pRquants;
    SEXP Rpolyposn, Rcorrs, Rwinlocs, Rtargetseqs, Rquants, val_nm, val=R_NilValue;
    *pargv = &pargvSize[0][0];

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

    numsig =  mainContinue(2, pargv, pbreaks, &polyposn, &winlocs, &corrs, &targetseqs, &quants);

    UNPROTECT(4);

    if (numsig == 0) return(R_NilValue);

    PROTECT(Rcorrs = NEW_NUMERIC(numsig));
    PROTECT(Rwinlocs = NEW_INTEGER(numsig));
    PROTECT(Rtargetseqs = allocVector(STRSXP, numsig));
    pRcorrs = NUMERIC_POINTER(Rcorrs);

    pRwinlocs = INTEGER_POINTER(Rwinlocs);

    for(i = 0; i < numsig; i++) {
        pRcorrs[i] = corrs[i];
	pRwinlocs[i] = winlocs[i];
        SET_STRING_ELT(Rtargetseqs, i, mkChar(sequenceLabels[targetseqs[i]]));
    }   
    PROTECT(Rpolyposn = NEW_INTEGER(numBases));
    pRpolyposn = INTEGER_POINTER(Rpolyposn);
    for(i = 0; i < numBases; i++) 
        pRpolyposn[i] = polyposn[i]+1;
    PROTECT(Rquants = NEW_NUMERIC(2));
    NUMERIC_POINTER(Rquants)[0] = quants[0];
    NUMERIC_POINTER(Rquants)[1] = quants[1];

    PROTECT(val_nm = allocVector(STRSXP, val_len)); 
    SET_STRING_ELT(val_nm, 0, mkChar("polyposn"));
    SET_STRING_ELT(val_nm, 1, mkChar("corrs"));
    SET_STRING_ELT(val_nm, 2, mkChar("winlocs"));
    SET_STRING_ELT(val_nm, 3, mkChar("target.seqs"));
    SET_STRING_ELT(val_nm, 4, mkChar("quants"));

    PROTECT(val = allocVector(VECSXP, val_len));
    SET_VECTOR_ELT(val, 0, Rpolyposn);
    SET_VECTOR_ELT(val, 1, Rcorrs);
    SET_VECTOR_ELT(val, 2, Rwinlocs);
    SET_VECTOR_ELT(val, 3, Rtargetseqs);
    SET_VECTOR_ELT(val, 4, Rquants);
    setAttrib(val, R_NamesSymbol, val_nm);

    UNPROTECT(7);
    return(val);
}
