#include <math.h>

extern int numBases, nseqs;
extern char **sequenceLabels;
extern char **sequences;

/***********************************************************************
                      Begin Function prototypes                       */

double mychisq(int c,int d, int WH);
int diff(int seq1, int seq2, int *range, int *polyposn);
double getMaxChi(int *polyposn, int nsites, int winHalfWidth); 
void doChi(int *polyposn, int nsites, int winHalfWidth, double nullquant,
		      int maxnum, int *numsigOut, double *chisqsOut, 
                    int *winlocsOut, int *pairmem1Out, int *pairmem2Out);

/*                     End function prototypes
***********************************************************************/


