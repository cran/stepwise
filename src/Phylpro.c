#include "getPhylpro.h"
#include "init.h"
#include <stdlib.h>
/*-----------------------------------------------------------*/
/*                                                           */
/* Jinko Graham                     Brad McNeney             */
/* jgraham@stat.sfu.ca       and    mcneney@stat.sfu.ca      */
/*               Department of Statistics                    */
/*               Simon Fraser University                     */
/*               Burnaby, BC Canada V5A 1S6                  */
/*                          and                              */
/* Francoise Seillier-Moiseiwitsch                           */
/* seillier@math.umbc.edu                                    */
/* Bioinformatics Research Center                            */
/* Department of Mathematics and Statistics                  */
/* University of Maryland-Baltimore County                   */
/* Baltimore MD 21250, U.S.A.                                */
/*                                                           */
/* NB: Code for reading in data in Phylip format is taken    */
/* from Phylip version 3.56c.                                */
/* (c) Copyright 1993 by Joseph Felsenstein.                 */
/*-----------------------------------------------------------*/

/*prototypes*/

int main_phylpro(int argc, char **argv, int *breaks, int **polyposn, int 
**winlocs, double **corrs, int **targetseqs, double **quants);
static int checkSegs(int *endpoints);
static double getOverallPhylpro(int *polyposn, int *endpoints);
static void overalldoPhylpro(int *polyposn, int *endpoints, double nullquant, int maxnum,
                  int *numsigOut, double *corrsOut, int *winlocsOut,
                  int *targetseqsOut);
static double *getNullquant(int *polyposn, int *endpoints);
static double *getNulldist(int *polyposn, int *endpoints);
static double *quantile(double *permDistn, double *props, int num);
static void sort(double *vec);
static int compareDbl(const void *key1, const void *key2); /* Needed in sort */
static void permutePosn(int *polyposn, int *endpoints);
static void permute(int *polyposn, int start, int stop);
static void printResults(double *quants, int numsig,double *corrs,
                  int *winlocs, int *targetseqs);
static int siteSpecificSummary(int numsig, double *corrs, int *winlocs, 
		int *targetseqs, double *siteCorrs, 
		int *siteWinlocs, char **targnames);


/* -------------- global variables --------------- */
int numBases, numBreaks, nseqs;
int winHalfWidth, permReps;
char **sequenceLabels;
char **sequences;
static char programName[] = "Stepwise Phylpro";
static char version[] = "0.1";


/*=================================================================*/

int main_phylpro(int argc, char **argv, int* breaks, int **Rpolyposn, int 
**Rwinlocs, double **Rcorrs, int **Rtargetseqs, double **Rquants)
{
    int *polyposn;
    int *winlocs;
    double *corrs;
    int *targetseqs;
    double *quants;

    int i, orignseqs, counter=0, noValidSegs=0, maxnum, numsig=0;
  int *endpoints; /*vector of indices for polyposns just before bps plus ends of sequence*/

  /*----------------------------------------------------------------*/
  /* 1. Read in sequence data to a global variable called sequences */
  SetSeed(time(NULL));
  readPhylipData(argc, argv); /* Reads in sequence data in Phylip format */
  orignseqs=nseqs; /*backup copy of number of sequences read in*/

  /* 3. More preparation: find unique sequences and polymorphic sites */
  reduceToUniqueSeqs(); /*strip out duplicate sequences: 
		alters global var nseqs and sequences matrix*/
  polyposn = make_initial_index(); /*Indices of variable, non-gap sites:
				  modifies the global variable numBases*/
  if (!strcmp("phylpro", argv[0]) || !strcmp("./phylpro", argv[0])) {
  fprintf(stdout,"\nThere are %d unique sequences in the %d provided.\n\n",nseqs,orignseqs);
  fprintf(stdout,"There are %d ungapped polymorphic sites:\n",
	  numBases);
  for(i=0; i<numBases; i++) { /* print a list of polymorphic sites */
    fprintf(stdout,"%d ", polyposn[i]+1);
    counter++;
    if(counter>19) { 
       fprintf(stdout,"\n");
       counter=0;
    }
  }
  fprintf(stdout,"\n");
  }
  /* Find the endpoints of segments defined by previously found breaks, */
  /* where endpoints are indices for the polyposn just before a breakpoint. */
  endpoints = findEndpoints(polyposn,breaks,numBreaks); 
  noValidSegs = checkSegs(endpoints); /* returns 1 if all segments are 
					    smaller than window width */
  if(noValidSegs == 1) {
    printf("\nNo further steps are possible: specified window half-width larger than number\n");
    printf("of polymorphic sites in all segments defined by previously-declared breaks\n\n");
  } else {
    /*----------------------------------------------------------------*/
    /* 4. Get the 90th and 95th percentiles of the permutation distribution 
          for the extreme value over the alignment. */

    quants = getNullquant(polyposn, endpoints);
    /*----------------------------------------------------------------*/
    /* 5. Find the observed values that exceed the quantile.          */
    /* First set up space for the output -- there are nseqs target
       sequences, and there are numBases-2*winHalfWidth+1 possible
       breakpoints (e.g. if there are 5 bases and winHalfWidth is
       2 there are two possible breakpoints). */
    maxnum = nseqs*(numBases-2*winHalfWidth+1);
    corrs = (double *)malloc(maxnum*sizeof(double));
    winlocs = (int *)malloc(maxnum*sizeof(int));
    targetseqs = (int *)malloc(maxnum*sizeof(int));
    overalldoPhylpro(polyposn, endpoints, quants[0], maxnum, &numsig, 
                     corrs, winlocs, targetseqs);
    /*----------------------------------------------------------------*/
    /* 6. Print out the results.                                      */
    if (!strcmp("phylpro", argv[0]) || !strcmp("./phylpro", argv[0])) {
	printResults(quants,numsig,corrs,winlocs,targetseqs);
    }
  } /* end else */

  *Rpolyposn = polyposn;
  *Rwinlocs = winlocs;
  *Rcorrs = corrs;
  *Rtargetseqs = targetseqs;
  *Rquants = quants;
  return(numsig);
}

int main(int argc, char **argv)
{
  int *breaks;  /*vector of break points*/
  int *polyposn;  /*vector of posns of polysites in original alignment */
  double *quants; 
  int *targetseqs, *winlocs;
  double *corrs;

  /* 0. Print out description of the program */
  printf("-----------------------------------------------\n");
  printf("%s, version %s output\n",programName,version);
  printf("-----------------------------------------------\n");

  /*----------------------------------------------------------------*/
  /* 2. Read in other params such as window half-widths and previous breaks */
  breaks = readOtherData(&numBreaks);
  /*----------------------------------------------------------------*/
  main_phylpro(argc, argv, breaks, &polyposn, &winlocs, &corrs, &targetseqs, &quants);
  return(1);
} /* end main */
  

/*===================================================================*/
/* Functions used in the main program */

static int checkSegs(int *endpoints)
{
  int i,k;
  int noValidSegs=1; /* init to 1=true */

  /* The endpoints vector contains indices of polymorphic sites just 
   * before breakpoints *except* first and last elements which are
   * 0 and numBases-1, respectively, indicating the ends of the
   * alignment. The first segment (call it seg_0) from endpoints[0] 
   * to endpoints[1] includes endpoints[0] as an element. Subsequent
   * segments seg_i, i>0, are from endpoints[i]+1 to endpoints[i+1]
   * inclusive. */

  /* test the first segment */
  k=endpoints[1]-endpoints[0]+1; /* number of indices in seg_0 */
  if(k>= 2*winHalfWidth) {
    noValidSegs=0;
  }
  /* test remaining segments, if any */
  for(i=1; i<numBreaks+1; i++) {
    k=endpoints[i+1]-endpoints[i]; /* number of indices in seg_i */
    if(k>= 2*winHalfWidth) {
      noValidSegs=0;
    }
  }
  return(noValidSegs);

}

static double getOverallPhylpro(int *polyposn, int *endpoints)
{
  /* Run Phylrpo over all segments of the alignment                 */

  int i, j, k, *myvec;
  double newcor;
  double mincor=1.0;

  myvec = (int *)malloc(numBases*sizeof(int));

  /* As noted above in checkSegs, the first segment is special and needs
   * to be treated separately */
  k=0;
  for(j=endpoints[0]; j<=endpoints[1]; j++) {
    myvec[k]=polyposn[j];
    k++;
  }
  if(k>= 2*winHalfWidth) {
    newcor = getMinCor(myvec,k,winHalfWidth);
    if(newcor<mincor)
       mincor=newcor;
  }
  /* Now do the remaining segments, if any */
  for(i=1; i<numBreaks+1; i++) {
    k=0;
    for(j=(endpoints[i]+1); j<=endpoints[i+1]; j++) {
      myvec[k]=polyposn[j];
      k++;
    }
    if(k>= 2*winHalfWidth) {
      newcor = getMinCor(myvec,k,winHalfWidth);
      if(newcor<mincor)
	  mincor=newcor;
    }
  }
  return(mincor);
}

static void overalldoPhylpro(int *polyposn, int *endpoints, double nullquant, 
                  int maxnum, int *numsigOut, double *corrsOut, 
		  int *winlocsOut, int *targetseqsOut)
{
  /* Run doPhylpro over all segments of the alignment                */

  int i, j, k, *myvec;

  myvec = (int *)malloc(numBases*sizeof(int));
  /* As noted above in checkSegs, the first segment is special and 
     needs to be treated separately */
   k=0;
   for(j=endpoints[0]; j<=endpoints[1]; j++) {
       myvec[k]=polyposn[j];
       k++;
   }
    if(k>= 2*winHalfWidth) {
      /* numsigOut is a running counter of the number of significant
         correlations found so far. This will be modified in doPhylpro.
         In doPhylpro, all indexing of the vectors corrsOut, winlocsOut, 
         targetseqsOut is done via numsigOut */
      doPhylpro(myvec, k, winHalfWidth, nullquant, maxnum, numsigOut,
            corrsOut, winlocsOut, targetseqsOut);
    }

  for(i=1; i<numBreaks+1; i++) { /* loop over remaining segments */
    k=0;
    for(j=endpoints[i]+1; j<=endpoints[i+1]; j++) {
      myvec[k]=polyposn[j];
      k++;
    } 
    if(k>= 2*winHalfWidth) {
      /* numsigOut is a running counter of the number of significant
         correlations found so far. This will be modified in doPhylpro.
         In doPhylpro, all indexing of the vectors corrsOut, winlocsOut, 
         targetseqsOut is done via numsigOut */
      doPhylpro(myvec, k, winHalfWidth, nullquant, maxnum, numsigOut,
            corrsOut, winlocsOut, targetseqsOut);
    }
  }
}

/**********************************************************************/
/* getNullquant and related functions                                 */

static double *getNullquant(int *polyposn, int *endpoints)
{
  /* Find the 90th and 95th percentiles of the permutation null distribution
     of maxima. Calls getNulldist and quantile */
  double *permDistn, *nullquants, *props;

  props = (double *)malloc(sizeof(double)*2);
  props[0]=0.10; /* for phylpro so minimum correlation */
  props[1]=0.05; 
  permDistn = getNulldist(polyposn, endpoints);
  nullquants = quantile(permDistn, props, 2); /* ask for the 2 quantiles */
  return(nullquants);
}

static double *getNulldist(int *polyposn, int *endpoints)
{
  /* Return permutation null distribution of maxima */
  /* Make use of the global variable permReps */
  double *permDistn;
  int i;
  int *mypolyposn; /* to be filled with a permutation */
  mypolyposn = copyIntVec(polyposn,numBases); 


  permDistn = (double *)malloc(permReps*sizeof(double));
  for(i=0; i<permReps; i++) {
    permutePosn(mypolyposn, endpoints); /* permute sites within segments */
    permDistn[i]=getOverallPhylpro(mypolyposn, endpoints);
  }
  return(permDistn);
}

static void permutePosn(int *polyposn, int *endpoints)
{
  /* Permute polymorphic sites (given by polyposn) within 
     segments (given by endpoints). Calls permute.  */

  int j;

  /*first segment is a special case since it starts at endpoints[0]
     rather than endpoints[0]+1. Recall that endpoints[0] is the index of
     the first polymorphic site in the vector polyposn of positions
     in the original alignment; similarly endpoints[1] is the index
     of the second polysite in the vector polyposn of positions in
     the original alignment.   */
  permute(polyposn,endpoints[0],endpoints[1]); 

  /*there are numBreaks+2 endpoints*/
  for(j=1;j<(numBreaks+1);j++) 
      permute(polyposn,endpoints[j]+1,endpoints[j+1]);
}


static void permute(int *polyposn, int start, int stop)
{
/* permute by drawing numbers w/o replacement from an ``urn'' and    */
/* storing in the appropriate places in the vector ``polyposn''      */
  int i,j, urnElement, urnLength=(stop-start+1);
  int *urn=(int *)malloc(sizeof(int)*(stop-start+1));

  /*Copy the appropriate elements of polyposn into the urn*/
  for(i=start;i<=stop;i++)
    urn[i-start]=polyposn[i];

  for(i=start;i<=stop;i++) {
    urnElement = floor(urnLength*rndu());  /* randomly select an element */
    polyposn[i]=urn[urnElement];
    for(j=urnElement+1;j<urnLength;j++)
       urn[j-1]=urn[j];
    urnLength--;
  }
}

static double *quantile(double *permDistn, double *props, int num)
{
  /* Sort the permutation distribution and return the quantiles */
  double *myquants;
  int i, indx;

  myquants = (double *)malloc(sizeof(double)*num);
  sort(permDistn);
  for(i=0;i<num;i++) {
    indx=floor(props[i]*permReps)-1; /* index of our quantile; -1 for C-indexing. */
    myquants[i]=permDistn[indx]; 
  }
  return(myquants);
}

static void sort(double *vec)
{

  /* Just a wrapper for the standard C library qsort function */
  qsort(vec, permReps, sizeof(double), compareDbl); 

}

static int compareDbl(const void *key1, const void *key2)
{

  if(*(double *)key1 < *(double *)key2) 
    return -1;
  else if(*(double *)key1 > *(double *)key2) 
    return 1;
  else 
    return 0;
}

/**********************************************************************/

static void printResults(double *quants, int numsig, double *corrs,
                  int *winlocs, int *targetseqs)
{
  /* quants is a vector with 10th and 5th percentiles of perm'tn distn
   * respectively. There are numsig significant results in corrs, winlocs
   * and target
   */
  int i, sigSites;
  double *siteCorrs;
  int *siteWinlocs;
  char **targnames; 

siteCorrs = (double *)malloc(sizeof(double)*numsig);
siteWinlocs = (int *)malloc(sizeof(int)*numsig);
targnames = (char **)malloc(sizeof(char *)*numsig);

  for(i=0;i<numsig;i++) 
	  targnames[i]=(char *)malloc(sizeof(char)*500);
  if(numsig>0) {
    sigSites = siteSpecificSummary(numsig, corrs, winlocs, targetseqs,
        siteCorrs, siteWinlocs, targnames);
    printf("-----------------------------------------------\n");
    printf("There were %d site-specific minimum correlation statistics significant at the\n", sigSites);
    printf("10 percent level (10th percentile = %5.3f, 5th percentile = %5.3f):\n\n", quants[0],quants[1]);
    printf("Number Location  MinCor   targets\n");
    for(i=0; i<sigSites; i++) {
      if(siteCorrs[i]<quants[1]) {
	 printf("%6d   %6d  %5.3f*  %s\n",
           (i+1), siteWinlocs[i], siteCorrs[i], targnames[i]);
      } else {
        printf("%6d   %6d  %5.3f   %s\n",
           (i+1), siteWinlocs[i], siteCorrs[i],  targnames[i]);
      }
    }
    printf("------------------------------------------------\n");
    printf("Notes - \"Location\" is the polymorphic site just before the proposed breakpoint.\n");
    printf("      - MinCor statistics significant at the 5 percent level indicated by a * \n\n");
  } else {
    printf("  No significant minimum correlation statistics found:\n\n");
  }
}
  
static int siteSpecificSummary(int numsig, double *corrs, int *winlocs, 
		int *targetseqs, double *siteCorrs, 
		int *siteWinlocs, char **targnames)
{
  /* Find most significant signal at each site. Report all target seq names that
   * tie for most significant signal. */

  int i; 
  int sigSites=0;

  siteCorrs[0]=corrs[0]; 
  siteWinlocs[0]=winlocs[0];
  sprintf(targnames[sigSites],"%s",sequenceLabels[targetseqs[0]]);
  for(i=1;i<numsig;i++) {
    if(winlocs[i] != siteWinlocs[sigSites]) { /* found a new site */
      sigSites++;
      siteCorrs[sigSites]=corrs[i]; 
      siteWinlocs[sigSites]=winlocs[i];
      sprintf(targnames[sigSites],"%s",sequenceLabels[targetseqs[i]]);
    } else if(corrs[i] < siteCorrs[sigSites]) {
      siteCorrs[sigSites]=corrs[i];
      sprintf(targnames[sigSites],"%s",sequenceLabels[targetseqs[i]]);
    } else if(corrs[i] == siteCorrs[sigSites]) {
      /* concatenate this pair on end of current list of pairs */
      sprintf(targnames[sigSites],"%s,%s",targnames[sigSites],
		                    sequenceLabels[targetseqs[i]]);
    }
  }
  return(sigSites+1);
}

/*-------------------------------------------------------------------*/
