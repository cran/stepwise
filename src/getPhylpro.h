#include <math.h>

extern int numBases, nseqs;
extern char **sequenceLabels;
extern char **sequences;

/***********************************************************************
                      Begin Function prototypes                       */

double cor(int *X, int *Y, int dim);
double getMinCor(int *polyposn, int nsites, int winHalfWidth);
void doPhylpro(int *polyposn, int nsites, int winHalfWidth, double nullquant,
		     int maxnum, int *numsigOut, double *corsOut, 
		     		     int *winlocsOut, int *targetOut);

/*                     End function prototypes
***********************************************************************/


