#include "getMaxChi.h"
#include "init.h"

/* NB: We return the site just BEFORE a breakpoint (as in TREEVOLVE). 
See code for C-function doChi. */

/***********************************************************************/
double getMaxChi(int *polyposn, int nsites, int winHalfWidth) 
{
/*
  Compute the maximum chi-squared stat over all possible pairs of
  sequences and over all possible positions for the centre of the sliding 
  window with half-width winHalfWidth. 

  Input: *polyposn - vector of positions of polysites in segment of alignment being considered
	 nsites - number of polysites in segment of alignment being considered
	 winHalfWidth - window half width
  Output: maxChiOut - the max chisquare stat over all window placements
                      and all pairs of sequences                     */

  /******Begin variable declarations******/
  double maxChi=0.0; /* initialize maximum chisqare stat to be 0*/
  double myChi;
  
  /* Create vectors to hold indices of first and last sites in left
     half of sliding window (lSites) and right half (rSites)*/
  int lSites[2], rSites[2];
  int i,j,k; /* indices for looping */

  int **lprdiffs, **rprdiffs;
  lprdiffs=(int **)dynamicArray(nseqs, nseqs); /*matrix of left pairdiffs*/
  rprdiffs=(int **)dynamicArray(nseqs, nseqs); /*matrix of right pairdiffs*/

  /*******End variable declarations******/

  /* Loop over all positions for sliding window. 
     Let window be centred between adjacent sites and let i be the site 
     just *after* the window centre so that i belongs to the right 
     fragment and i-1 to the left fragment. Thus window consists of sites 
     123c456 i-winHalfWidth to i+winHalfWidth-1. E.G. window half width
      x  ix
     of 2 and i=4 means the window is from 4-2=2 to 4+2-1=5
     First placement of window centre is at i=winHalfWidth+1 and last
     placement is at nsites-winHalfWidth+1
     *BUT* we're in C now so placement of first window is at
     i=winHalfWidth and last is at nsites-winHalfWidth.  */

  for(i=winHalfWidth; i<=(nsites-winHalfWidth); i++){ /*window loop*/
    lSites[0]=(i-winHalfWidth); /*index of posn of first polysite in left fragment*/
    lSites[1]=(i-1); /*index of posn of last polysite in left fragment*/
    rSites[0]=i; /*index of posn of first polysite in right fragment*/
    rSites[1]=(i+winHalfWidth-1); /*index of posn of last polysite in right fragment*/
    for(j=0; j<(nseqs-1); j++)  /*pairs loop part I*/
       for(k=j+1; k<nseqs; k++){ /*pairs loop part II*/
           lprdiffs[j][k]=diff(j,k,lSites,polyposn);
           rprdiffs[j][k]=diff(j,k,rSites,polyposn);
           myChi=mychisq(lprdiffs[j][k],rprdiffs[j][k],winHalfWidth);
           if(myChi>maxChi)  maxChi=myChi;
        } /*end k-loop*/
    /*end j-loop and complete the left and right difference matrices*/
  } /*finish the loop through possible window locations*/
     
  return(maxChi);
}

/***********************************************************************/
double mychisq(int c,int d, int WH){
   int a,b,num,den;

   a=WH-c; b=WH-d;
   num=(a*d - b*c)*(a*d - b*c);
   if(num==0) return(0.0);
   else {
   den=WH*WH*(a+b)*(c+d); 
   return( 2*WH*num/((double) den) );
   } 
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
void doChi(int *polyposn, int nsites, int winHalfWidth, double nullquant,
		     int maxnum, int *numsigOut,
		     double *chisqsOut, int *winlocsOut,
		     int *pairmem1Out, int *pairmem2Out)
{
/* Compute the chi-squared stats over all possible pairs of
  sequences and over all possible positions for the centre of the sliding 
  window with half-width winHalfWidth. 
  Returns numsig, the number of significant chi-squared statistics,
  and up to maxnum number of values of these statistics, the associated 
  pair of sequences and the window locations.
  The suffix "Out" means the variable is an output.

 Input: *polyposn - vector of positions of polysites in segment of alignment 
         nsites - number of polysites in segment of alignment 
         winHalfWidth - window half width
	 nullquant - the quantile from the null distribution
	 maxnum - maximum number of significant chisq stats to return
 Output: numsigOut - number of significant chisq stats found
	 chisqsOut - vector of chisq stats found (up to maxnumIn of them)
	 winlocsOut - vector of window locations found (up to maxnumIn)
	              NB: winloc is site to **right** of break -- see below
	 pairmem1Out - vector of first members of pairs of sequences with
  	               significant chisq stat
         pairmem2Out - vector of second pair members */

  /******Begin variable declarations******/
  double myChi;

  int numsig= *(numsigOut); /* working copy of number of signif chisqs found*/
  
  /* Create vectors to hold indices of first and last sites in left
     half of sliding window (lSites) and right half (rSites)*/
  int lSites[2], rSites[2];
  int i,j,k; /* indices for looping */
  int **lprdiffs, **rprdiffs;
  lprdiffs=(int **)dynamicArray(nseqs, nseqs);
  rprdiffs=(int **)dynamicArray(nseqs, nseqs);

/*  int lprdiffs[nseqs][nseqs]; *//* matrix of left pairdiffs*/
/*  int rprdiffs[nseqs][nseqs]; *//*matrix of right pairdiffs*/
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
           myChi=mychisq(lprdiffs[j][k],rprdiffs[j][k],winHalfWidth);
           if(myChi>nullquant) {
	      /*numsig starts as 0 so is appropriate to index els
	        of the output vectors. After adding to the output
		vectors increment numsig to reflect the sig stat 
		found and for use in the next iteration. */
              if(numsig<maxnum){
	         chisqsOut[numsig]=myChi; 
	         /*indices are in C so they're 1 less than what would
	           make sense in S -- add 1 to each*/
	         winlocsOut[numsig]=polyposn[i-1]+1; 
		     /*NB +1 is to change C-style index to natural indexing.
		       This returns site to *left* of the bp (as in TREEVOLVE);
		       change to polyposn[i]+1 if site to *right* is wanted */
                 pairmem1Out[numsig]=j;
                 pairmem2Out[numsig]=k;
              }
	      numsig++; 
	   } /*end if*/
        } /*end k-loop*/
    /*end j-loop and complete the left and right difference matrices*/
  } /*finish the loop through possible window locations*/
     
  *(numsigOut)=numsig; /* update the output variable */
}
