#ifndef _INIT_H_
#define _INIT_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define LABEL_LENGTH 10 /* length of sequence labels -- same as in phylip */
#define STRLEN 100
#define FALSE 0
#define TRUE 1
#define NA -1  /*Our usual missing value code*/
#define HUGE_INT 1000000
#define BIG_NUMBER 1.0e+30
#define SMALL_NUMBER 0.000000001
#define MIN_LENGTH 5 /*min seq length to do distance calcs on */
#define PERMREPS 1000
/* macro definitions */
#define min(a,b) (a<b ? a : b)
#define max(a,b) (a>b ? a : b)
#define TRUE 1
#define FALSE 0

/* global variables */
extern int numBases, nseqs;
extern int winHalfWidth, permReps;
extern char **sequenceLabels;
extern char **sequences;

/*prototypes*/
void SetSeed (int seed);
double rndu (void);
void readPhylipData(int argc, char **argv);
int anyGapsAtSite(int j);
int isVariable(int j);
void printSequences();
int *readOtherData(int *pnumBreaks);
char **makeCharArray(int nrows, int ncols);
double **make_double_mat(int nrows, int ncols);
int **make_int_mat(int nrows, int ncols);
int *make_initial_index(void);

void reduceToUniqueSeqs(void);
int *findUniqueSeqs(void);
int newSeq(char *sequence, int *uniqueSeqInd, int uniqueSeqNum);

int *copyIntVec(int *vec, int len);
int *findEndpoints(int *index, int *breaks, int numBreaks);
void free_double_mat(double **mat, int nrows);
void free_int_mat(int **mat, int nrows);
void freeCharArray(char **cArray, int nrows);

int eof(FILE *f); /* from Phylip */
int eoln(FILE *f); /* from Phylip */

int **dynamicArray(int nrows, int ncolumns); 

#endif /* _INIT_H_ */
