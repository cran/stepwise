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

#include "init.h"
#include "getMaxChi.h"

/*prototypes*/

int main_maxchi(int argc, char **argv, int *breaks, int **polyposn, int 
**winlocs, double **chisqs, int **pairmem1, int **pairmem2, double 
**quants);
static int checkSegs(int *endpoints);
static double getOverallMaxChi(int *polyposn, int *endpoints);
static void overalldoChi(int *polyposn, int *endpoints, double nullquant, int maxnum,
                  int *numsigOut, double *chisqsOut, int *winlocsOut,
                  int *pairmem1Out, int *pairmem2Out);
static double *getNullquant(int *polyposn, int *endpoints);
static double *getNulldist(int *polyposn, int *endpoints);
static double *quantile(double *permDistn, double *props, int num);
static void sort(double *vec);
static int compareDbl(const void *key1, const void *key2); /* Needed in sort */
static void permutePosn(int *polyposn, int *endpoints);
static void permute(int *polyposn, int start, int stop);
static void printResults(double *quants, int numsig,double *chisqs,
                  int *winlocs, int *pairmem1, int *pairmem2);
static int siteSpecificSummary(int numsig, double *chisqs, int *winlocs, 
		int *pairmem1, int *pairmem2, double *siteChisqs, 
		int *siteWinlocs, char **pairs);


/* -------------- global variables --------------- */
extern int numBases, numBreaks, nseqs;
extern int winHalfWidth, permReps;
extern char **sequenceLabels;                          
extern char **sequences;                               

#ifdef C_PROGRAM
static char programName[] = "Stepwise MaxChi";
static char version[] = "0.1-1";
#endif

/*=========================================================================*/

int main_maxchi(int argc, char **argv, int* breaks, int **Rpolyposn, int 
**Rwinlocs, double **Rchisqs, int **Rpairmem1, int **Rpairmem2, double 
**Rquants) 
{

  int *polyposn, *winlocs;
  int *pairmem1, *pairmem2;
  double *chisqs, *quants;
  int i, orignseqs, counter=0, noValidSegs=0, numsig=0;
  int *endpoints; /*vector of indices for polyposns just before bps plus ends of sequence*/
  int maxnum;

  /*----------------------------------------------------------------*/
  /* 1. Read in sequence data to a global variable called sequences */
  SetSeed(time(NULL));
  readPhylipData(argc, argv); /* Reads in sequence data in Phylip format */
  orignseqs=nseqs; /*backup copy of number of sequences read in*/


  /*----------------------------------------------------------------*/
  /* 3. More preparation: find unique sequences and polymorphic sites */
  reduceToUniqueSeqs(); /*strip out duplicate sequences: 
		alters global var nseqs and sequences matrix*/

  polyposn = make_initial_index(); /*Indices of variable, non-gap sites:
				     modifies the global variable numBases*/

  if (!strcmp("maxchi", argv[0]) || !strcmp("./maxchi", argv[0])) {
    fprintf(stdout,"\nThere are %d unique sequences in the %d provided.\n\n",nseqs, orignseqs);
    fprintf(stdout,"There are %d ungapped polymorphic sites:\n",
	    numBases);
    for(i=0; i<numBases; i++) {/* print out a list of polymorphic sites */
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
    /*-------------------------------------------------------------*/
    /* 4. Get the 90th and 95th percentiles of the permutation distribution 
       for the extreme value over the alignment. */

    quants = getNullquant(polyposn, endpoints);
    /*-------------------------------------------------------------*/
    /* 5. Find the observed values that exceed the quantile.       */
    /* First set up space for the output -- there are nseqs choose 2
       pairs of sequences, and there are numBases-2*winHalfWidth+1 
       possible breakpoints (e.g. if there are 5 bases and winHalfWidth is 2 there are two possible bps). */
    maxnum = (nseqs*(nseqs-1)/2)*(numBases-2*winHalfWidth+1);
    chisqs = (double *)malloc(maxnum*sizeof(double));
    winlocs = (int *)malloc(maxnum*sizeof(int));
    pairmem1 = (int *)malloc(maxnum*sizeof(int));
    pairmem2 = (int *)malloc(maxnum*sizeof(int));

    overalldoChi(polyposn, endpoints, quants[0], maxnum, &numsig, chisqs, winlocs, pairmem1, pairmem2);
    /*-------------------------------------------------------------*/
    if (!strcmp("maxchi", argv[0]) || !strcmp("./maxchi", argv[0])) {
       /* 6. Print out the results.  */
          printResults(quants,numsig,chisqs,winlocs,pairmem1,pairmem2);  
    }
  } /* end else */

  *Rpolyposn =  polyposn;
  *Rwinlocs  =  winlocs;
  *Rpairmem1 =  pairmem1;
  *Rpairmem2 =  pairmem2;
  *Rchisqs   =  chisqs;
  *Rquants   =  quants;
  return(numsig);
}

#ifdef C_PROGRAM
int main(int argc, char **argv)
{
  int *polyposn;  /*vector of posns of polysites in original alignment */
  double *quants; 
  int *pairmem1, *pairmem2, *winlocs;
  double *chisqs; 
  int *breaks;  /*vector of break points*/

  /* 0. Print out description of the program */
  printf("-----------------------------------------------\n");
  printf("%s, version %s output\n",programName,version);
  printf("-----------------------------------------------\n");
  /*----------------------------------------------------------------*/
  /* 2. Read in other params such as window half-widths and previous breaks */
  breaks = readOtherData(&numBreaks);

  main_maxchi(argc, argv, breaks, &polyposn, &winlocs, &chisqs, &pairmem1, 
&pairmem2, &quants);

  return(0);
  } /* end main */  

#endif

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

/*===================================================================*/
static double getOverallMaxChi(int *polyposn, int *endpoints)
{
  /* Run MaxChi over all segments of the alignment                 */

  int i, j, k, *myvec;
  double newchi;
  double maxchi=0.0;

  myvec = (int *)malloc(numBases*sizeof(int));

  /* As noted above in checkSegs, the first segment is special and needs
   * to be treated separately */
  k=0;
  for(j=endpoints[0]; j<=endpoints[1]; j++) {
    myvec[k]=polyposn[j];
    k++;
  }
  if(k>= 2*winHalfWidth) {
    newchi = getMaxChi(myvec,k,winHalfWidth);
    if(newchi>maxchi)
       maxchi=newchi;
  }
  /* Now do the remaining segments, if any */
  for(i=1; i<numBreaks+1; i++) {
    k=0;
    for(j=(endpoints[i]+1); j<=endpoints[i+1]; j++) {
      myvec[k]=polyposn[j];
      k++;
    }
    if(k>= 2*winHalfWidth) {
      newchi = getMaxChi(myvec,k,winHalfWidth);
      if(newchi>maxchi)
	  maxchi=newchi;
    }
  }
  return(maxchi);
}

/*===================================================================*/
static void overalldoChi(int *polyposn, int *endpoints, double nullquant, 
                  int maxnum,
                  int *numsigOut, double *chisqsOut, int *winlocsOut,
                  int *pairmem1Out, int *pairmem2Out)

{
  /* Run doChi over all segments of the alignment                 */

  int i, j, k, *myvec;

  myvec = (int *)malloc(numBases*sizeof(int));
  /* As noted above in checkSegs, the first segment is special and needs
   * to be treated separately */
  k=0;
  for(j=endpoints[0]; j<=endpoints[1]; j++) {
    myvec[k]=polyposn[j];
    k++;
  }
  if(k>= 2*winHalfWidth) {
    doChi(myvec, k, winHalfWidth, nullquant, maxnum, numsigOut,
            chisqsOut, winlocsOut, pairmem1Out, pairmem2Out);
  }
  /* Now do the remaining segments, if any */
  for(i=1; i<numBreaks+1; i++) {
    k=0;
    for(j=endpoints[i]+1; j<=endpoints[i+1]; j++) {
      myvec[k]=polyposn[j];
      k++;
    } 
    if(k>= 2*winHalfWidth) {
      /* numsigOut is a running counter of the number of significant
         chisquares found so far. This will be modified in doChi.
         In doChi, all indexing of the vectors chisqsOut, winlocsOut, 
         pairmem1Out, pairmem2Out is done via numsigOut */
      doChi(myvec, k, winHalfWidth, nullquant, maxnum, numsigOut,
            chisqsOut, winlocsOut, pairmem1Out, pairmem2Out);
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
  props[0]=0.90;
  props[1]=0.95;
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
    permDistn[i]=getOverallMaxChi(mypolyposn, endpoints);
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

  int *urn = (int *)malloc((stop-start+1) * sizeof(int));

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

static void printResults(double *quants, int numsig, double *chisqs,
                  int *winlocs, int *pairmem1, int *pairmem2)
{
  /* quants is a vector with 90th and 95th percentiles of perm'tn distn
   * respectively. There are numsig significant results in chisqs, winlocs
   * pairmem1 and pairmem2
   */
  int i, sigSites;

double *siteChisqs = (double *)malloc(numsig * sizeof(double));
int *siteWinlocs = (int *)malloc(numsig * sizeof(int));
char **pairs = (char **)malloc(numsig * sizeof(char *));

  for(i=0;i<numsig;i++) 
	  pairs[i]=(char *)malloc(sizeof(char)*500);
  if(numsig>0) {
    sigSites = siteSpecificSummary(numsig, chisqs, winlocs, pairmem1,
	       pairmem2, siteChisqs, siteWinlocs, pairs);
    printf("-----------------------------------------------\n");
    printf("There were %d site-specific MaxChi statistics significant at the\n", sigSites);
    printf("10 percent level (90th percentile = %5.3f, 95th percentile = %5.3f):\n\n", quants[0],quants[1]);
    printf("Number Location  MaxChi   pairs\n");
    for(i=0; i<sigSites; i++) {
      if(siteChisqs[i]>quants[1]) {
	printf("%6d   %6d  %5.3f*  %s\n",
	       (i+1), siteWinlocs[i], siteChisqs[i], pairs[i]);
      } else {
	printf("%6d   %6d  %5.3f   %s\n",
	       (i+1), siteWinlocs[i], siteChisqs[i],  pairs[i]);
      }
    }

    printf("------------------------------------------------\n");
    printf("Notes - \"Location\" is the polymorphic site just before the proposed breakpoint.\n");
    printf("      - MaxChi statistics significant at the 5 percent level indicated by a * \n\n");
  } else {
    printf("  No significant MaxChi statistics found.\n\n");
  }
  return;
}
  
static int siteSpecificSummary(int numsig, double *chisqs, int *winlocs, 
		int *pairmem1, int *pairmem2, double *siteChisqs, 
		int *siteWinlocs, char **pairs)
{
  /* Find most significant signal at each site. Report all pairs that
   * tie for most significant signal in pairs. */

  int i; 
  int sigSites=0;

  siteChisqs[0]=chisqs[0]; 
  siteWinlocs[0]=winlocs[0];
  sprintf(pairs[sigSites],"(%s:%s)",sequenceLabels[pairmem1[0]],
                                    sequenceLabels[pairmem2[0]]);
  for(i=1;i<numsig;i++) {
    if(winlocs[i] != siteWinlocs[sigSites]) { /* found a new site */
      sigSites++;
      siteChisqs[sigSites]=chisqs[i]; 
      siteWinlocs[sigSites]=winlocs[i];
      sprintf(pairs[sigSites],"(%s:%s)",sequenceLabels[pairmem1[i]],
                                    sequenceLabels[pairmem2[i]]);
    } else if(chisqs[i] > siteChisqs[sigSites]) {
      siteChisqs[sigSites]=chisqs[i];
      sprintf(pairs[sigSites],"(%s:%s)",sequenceLabels[pairmem1[i]],
                                    sequenceLabels[pairmem2[i]]);
    } else if(chisqs[i] == siteChisqs[sigSites]) {
      /* concatenate this pair on end of current list of pairs */
      sprintf(pairs[sigSites],"%s\n\t\t\t  (%s:%s)",pairs[sigSites],
		                    sequenceLabels[pairmem1[i]],
                                    sequenceLabels[pairmem2[i]]);
    }
  }

  return(sigSites+1);
}

/*-------------------------------------------------------------------*/
