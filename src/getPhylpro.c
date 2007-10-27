#include "getPhylpro.h"
#include "init.h"
#include <stdlib.h>

static int diff(int seq1, int seq2, int *range, int *polyposn);

/***********************************************************************/
double getMinCor(int *polyposn, int nsites, int winHalfWidth) 
{
/*******************************************************************
  Compute the negative of the min cor stat over all possible target sequences 
  and over all possible positions for the centre of the sliding window with 
  half-width winHalfWidth.

  Input: *polyposn    - vector of positions of polysites in segment of alignment 
                        being considered
         nsites       - number of polysites in segment of alignment being considered
   	 winHalfWidth - window half width
  Output: minCorOut   - the min correlation over all window placements and all 
                        target sequences       
******************************************************************/

  /******Begin variable declarations******/
  double minCor=1.0; /* initialize minimum phylpro correlation to be 1 */
  double myCor;
  
  /* Create vectors to hold indices of first and last sites in left
     half of sliding window (lSites) and right half (rSites)*/
  int lSites[2], rSites[2];
  int i,j,k; /* indices for looping */
  int **lprdiffs, **rprdiffs, *lDist, *rDist;

  lprdiffs=(int **)dynamicArray(nseqs, nseqs); /*matrix of left pairdiffs*/ 
  rprdiffs=(int **)dynamicArray(nseqs, nseqs); /*matrix of right pairdiffs*/ 

  lDist=(int *)malloc(sizeof(int)*(nseqs-1));
  rDist=(int *)malloc(sizeof(int)*(nseqs-1));

/*  int lDist[nseqs-1]; *//*left distance vector for a target seq*/
/*  int rDist[nseqs-1]; *//*right distance vector for a target seq*/
  /*******End variable declarations******/

  /* Loop over all positions for sliding window. 
     Let window be centred between adjacent sites and let i be the site 
     just *after* the window centre so that i belongs to the right 
     fragment and i-1 to the left fragment. Thus window consists of sites 
     123 456 i-winHalfWidth to i+winHalfWidth-1. E.G. window half width
      x  ix
     of 2 and i=4 means the window is from 4-2=2 to 4+2-1=5
     First placement of window centre is at i=winHalfWidth+1 and last
     placement is at nsites-winHalfWidth+1
     *BUT* we're in C now so placement of first window is at
     i=winHalfWidth and last is at nsites-winHalfWidth.  */

  for(i=winHalfWidth; i<=(nsites-winHalfWidth); i++){ /*window loop*/
    lSites[0]=(i-winHalfWidth); /*index of posn of 1st polysite in left fragment*/
    lSites[1]=(i-1); /*index of posn of last polysite in left fragment*/
    rSites[0]=i; /*index of posn of first polysite in right fragment*/
    rSites[1]=(i+winHalfWidth-1); /*index of posn last polysite in right fragment*/
    for(j=0; j<(nseqs-1); j++)  /*pairs loop part I*/
       for(k=j+1; k<nseqs; k++){ /*pairs loop part II*/
           lprdiffs[j][k]=diff(j,k,lSites,polyposn);
           rprdiffs[j][k]=diff(j,k,rSites,polyposn);
        } /*end k-loop*/
    /*end j-loop and complete the left and right difference matrices*/
    for(j=0; j<nseqs; j++) { /*loop over "target" seqs to compute corrs*/
       /*Pairwise difference matrices are upper triangular so the 
         distance vectors we want are based on "L"'s in the difference
         matrices with the vertical part of the L being els [0:(j-1),j] 
         and the horizontal part being els [j,(j+1):(nseqs-1)]. 
         Copy diffs into lDist and rDist vectors, then get correlation. */
       for(k=0;k<j;k++) { /*do the vertical part of the L*/
         lDist[k]=lprdiffs[k][j];
         rDist[k]=rprdiffs[k][j];
       } /*end for k*/
       for(k=(j+1); k<nseqs; k++) { /*do the horizontal part of the L*/
         lDist[k-1]=lprdiffs[j][k];
         rDist[k-1]=rprdiffs[j][k];
       } /*end for k*/
       myCor=cor(lDist,rDist,(nseqs-1));
       if(myCor<minCor) minCor=myCor;
    } /*finish loop over target seqs*/
  } /*finish the loop through possible window locations*/

  freedynamicArray(lprdiffs);
  freedynamicArray(rprdiffs);
  free(lDist);
  free(rDist);
     
  return(minCor);
}

/***********************************************************************/
double cor(int *X, int *Y, int dim)
{
  int i;
  double SXX, SXY, SYY, Xbar, Ybar, corr;

  SXX=SXY=SYY=Xbar=Ybar=0.0;
  for(i=0;i<dim;i++) {
      SXX += X[i]*X[i];
      SXY += X[i]*Y[i];
      SYY += Y[i]*Y[i];
      Xbar += X[i];
      Ybar += Y[i];
  }        
  Xbar /= ((double) dim);
  Ybar /= ((double) dim);
  SXX = (SXX-dim*Xbar*Xbar)/(dim-1.0); /*now this is var(X)*/
  SYY = (SYY-dim*Ybar*Ybar)/(dim-1.0);
  if(SXX>0 && SYY>0) {
    corr = (SXY - dim*Xbar*Ybar)/( sqrt(SXX)*sqrt(SYY)*(dim-1.0) );
  } else {
    corr = HUGE_VAL; /*not defined so set to infty*/
  }
  return(corr);
}
/***********************************************************************/
static int diff(int seq1, int seq2, int *range, int *polyposn)
{
	  int i; /*for looping*/
	    int diffCount=0; /*counter for differences*/

	      for(i=range[0]; i<=range[1]; i++) 
	          if(sequences[seq1][polyposn[i]]!= sequences[seq2][polyposn[i]]) diffCount++;

		    return(diffCount);
}
/***********************************************************************/
void doPhylpro(int *polyposn, int nsites, int winHalfWidth, double nullquant,
               int maxnum, int *numsigOut, double *corsOut, int *winlocsOut,
               int *targetOut)

{
/*******************************************************************
  Compute the min cor stats over all possible target
  sequences and over all possible positions for the centre of the sliding 
  window with half-width winHalfWidth. Returns numsig, the number of 
  significant mincor statistics, and up to maxnum number of values of 
  these statistics, the associated target sequences and the window locations.
  The suffix "Out" for variables in the function call means they're outputs.

 Input: *polyposn - vector of positions of polysites in segment of alignment 
         nsites - number of polysites in segment of alignment 
         winHalfWidth - window half width
   	 nullquant - the quantile from the null distribution
 	 maxnum - maximum number of significant corr stats to return

 Output: numsigOut - number of significant min cor stats found
	 corsOut - vector of cor stats found (up to maxnumIn of them)
	 winlocsOut - vector of window locations found (up to maxnumIn)
	 targetOut - vector of target seqs with significant min cor stat
******************************************************************/

  /******Begin variable declarations******/
  double myCor;

 int numsig= *(numsigOut); /* working copy of number of signif corrs found*/ 

  /* Create vectors to hold indices of first and last sites in left
     half of sliding window (lSites) and right half (rSites)*/
  int lSites[2], rSites[2];
  int i,j,k; /* indices for looping */
  int *lDist; /*left distance vector for a target seq*/
  int *rDist; /*right distance vector for a target seq*/

  int **lprdiffs=(int **)dynamicArray(nseqs, nseqs); /*matrix of left pairdiffs*/ 
  int **rprdiffs=(int **)dynamicArray(nseqs, nseqs); /*matrix of right pairdiffs*/ 

  lDist=(int *)malloc(sizeof(int)*(nseqs-1)); 
  rDist=(int *)malloc(sizeof(int)*(nseqs-1));

  /*******End variable declarations******/

  /* Loop over all positions for sliding window. 
     Let window be centred between adjacent sites and let i be the site 
     just *after* the window centre so that i belongs to the right 
     fragment and i-1 to the left fragment. Thus window consists of sites 
     123 456 i-winHalfWidth to i+winHalfWidth-1. E.G. window half width
      x  ix
     of 2 and i=4 means the window is from 4-2=2 to 4+2-1=5
     First placement of window centre is at i=winHalfWidth+1 and last
     placement is at nsites-winHalfWidth+1
     *BUT* we're in C now so placement of first window is at
     i=winHalfWidth and last is at nsites-winHalfWidth.  */

  for(i=winHalfWidth; i<=(nsites-winHalfWidth); i++){ /*window loop*/
    lSites[0]=(i-winHalfWidth); /*first site in left fragment*/
    lSites[1]=(i-1); /*last site in left fragment*/
    rSites[0]=i; /*first site in right fragment*/
    rSites[1]=(i+winHalfWidth-1); /*last site in right fragment*/
    for(j=0; j<(nseqs-1); j++)  /*pairs loop part I*/
       for(k=j+1; k<nseqs; k++){ /*pairs loop part II*/
           lprdiffs[j][k]=diff(j,k,lSites,polyposn);
           rprdiffs[j][k]=diff(j,k,rSites,polyposn);
        } /*end k-loop*/
    /*end j-loop and complete the left and right difference matrices*/
    for(j=0; j<nseqs; j++) { /*loop over "target" seqs to compute corrs*/
       /*Pairwise difference matrices are upper triangular so the 
         distance vectors we want are based on "L"'s in the difference
         matrices with the vertical part of the L being els [0:(j-1),j] 
         and the horizontal part being els [j,(j+1):(nseqs-1)]. 
         Copy diffs into lDist and rDist vectors, then get correlation. */
       for(k=0;k<j;k++) { /*do the vertical part of the L*/
         lDist[k]=lprdiffs[k][j];
         rDist[k]=rprdiffs[k][j];
       } /*end for k*/
       for(k=(j+1); k<nseqs; k++) { /*do the horizontal part of the L*/
         lDist[k-1]=lprdiffs[j][k];
         rDist[k-1]=rprdiffs[j][k];
       } /*end for k*/
       myCor=cor(lDist,rDist,(nseqs-1));
       if(myCor<nullquant) {
	   /*numsig starts as 0 so is appropriate to index els
             of the output vectors. After adding to the output
             vectors increment numsig to reflect the sig stat 
	     found and for use in the next iteration. */
           corsOut[numsig]=myCor; 
	   /*indices are in C so they're 1 less than what would
	     make sense in S -- add 1 to each*/
	   winlocsOut[numsig]=polyposn[i-1]+1; 
	   targetOut[numsig]=(j+1);
	   numsig++; 
       } /*end if*/
    } /*finish loop over target seqs*/
  } /*finish the loop through possible window locations*/
  freedynamicArray(lprdiffs);
  freedynamicArray(rprdiffs);
  free(lDist);
  free(rDist);
     
  numsigOut[0]=numsig;
}
