#include "init.h"
#include <time.h>
#include <ctype.h>

static int z_rndu=137;
static unsigned w_rndu=13757;

void SetSeed (int seed)
{
   z_rndu = 170*(seed%178) + 137;
   w_rndu=seed;
}

double rndu (void)
{
   static int x_rndu=11, y_rndu=23;
   double r;

   x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
   y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
   z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
   if (x_rndu<0) x_rndu+=30269;
   if (y_rndu<0) y_rndu+=30307;
   if (z_rndu<0) z_rndu+=30323;
   r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
   return (r-(int)r);
}
            

/* Functions to make/initialize data structures and read in data */
/*-----------------------------------------------------------------------*/
/* Read Phylip-style parameter file, either from inputfile given on the
 * command-line, or from the default file "infile" */
/*-----------------------------------------------------------------------*/
void readPhylipData(int argc, char **argv)
{
  int  i,j;
  char ch;
  char infileName[STRLEN];
  FILE *infile;
  int done;
	
  if(argc>1) 
    strcpy(infileName,argv[1]);
  else
    strcpy(infileName,"infile");
    
  if( (infile=fopen(infileName,"r")) == NULL ) {
    fprintf(stderr,"%s: Error in opening input file %s\n",infileName,
	    argv[1]);
    exit(0);
  }
  if(feof(infile)){
    fprintf(stderr, "%s: Unable to read input from %s\n",argv[0],infileName);
    exit(0);
  }
  fscanf(infile, "%d %d", &nseqs, &numBases);
  fscanf(infile, "%*[^\n]"); 
  ch=getc(infile); /* get past line-break */
  sequenceLabels = makeCharArray(nseqs,LABEL_LENGTH);
  sequences = makeCharArray(nseqs,numBases);

  /* most of the rest of this function is taken/adaped from Phylip */
  for(i=0;i<nseqs;i++) {
    for(j=0;j<LABEL_LENGTH;j++) {
      if (eof(infile) || eoln(infile)){
	fprintf(stderr,"ERROR: END-OF-LINE OR END-OF-FILE IN");
        fprintf(stderr," THE MIDDLE OF A SPECIES NAME\n");
        exit(0);
      }
      sequenceLabels[i][j] = getc(infile);
    }

    j = done = 0;
    while (!done && !eof(infile)) {
      while (j < numBases && !(eoln(infile) || eof(infile))) {
	ch = getc(infile);
	if (ch == ' ' || (ch >= '0' && ch <= '9'))
	  continue;
	ch=toupper(ch);
	if (strchr("ABCDGHKMNRSTUVWXY?O-.",ch) == NULL){
	  printf("ERROR: BAD BASE:%c AT POSITION%5d OF SPECIES %3d\n",
		 ch, j, i);
	  exit(-1);
	}
	j++;
	if (ch == '.')
	  ch = sequences[0][j - 1];
          sequences[i][j - 1] = ch;
        }
        if (j < numBases) {
          fscanf(infile, "%*[^\n]");
          getc(infile);
	} else if (j == numBases)
          done = 1;
    }
    fscanf(infile, "%*[^\n]");
    getc(infile);
  }
  if (!strcmp("maxchi", argv[0]) || !strcmp("./maxchi", argv[0])) 
    fprintf(stdout, "\nRead data from %s: sample size %d with %d bases\n",
	    infileName,nseqs,numBases);
}


void printSequences()
{
  int i,j;

  for(i=0;i<nseqs;i++) {
    for(j=0;j<numBases;j++)
      fputc(sequences[i][j],stdout);
    printf("\n");
  }
}

  


/*-------------------------------------------------------------------------*/
/* Function to read numBreaks and bnewVec and return bnewVec               */
/*-------------------------------------------------------------------------*/
int *readOtherData(int *pnumBreaks)
{
  int i, *bnewVec;

  fprintf(stderr,"\nEnter number of breaks: ");
  scanf("%d", pnumBreaks);
  fprintf(stdout,"\n%d prior breaks\n ",*pnumBreaks);

  /*Allocate memory for bnewVec */
  bnewVec = (int *)malloc((*pnumBreaks)*sizeof(int));
  
  if(*pnumBreaks > 0){
    fprintf(stderr,"\nEnter the %d ordered site(s) just before the break(s): ",
	    (*pnumBreaks));
    for(i=0;i<(*pnumBreaks);i++) {
      scanf("%d", &(bnewVec[i]));
      bnewVec[i]--; /*convert to C-style indexing*/
    }
    fprintf(stdout,"Breaks entered: ");
    for(i=0;i<(*pnumBreaks);i++) {
      fprintf(stdout,"%d ",bnewVec[i]+1);
    }
    fprintf(stdout,"\n");
  }

  fprintf(stderr,"\nEnter window half width to use: ");
  scanf("%d", &winHalfWidth);
  fprintf(stdout,"\nwindow half-width=%d sites.\n ",winHalfWidth);
  fprintf(stderr,"\nEnter number of MC reps to use for the permutation distribution: ");
  scanf("%d", &permReps);
  fprintf(stdout,"\n%d MC reps for the permutation distribution\n ",permReps);
  return(bnewVec);
}



char **makeCharArray(int nrows, int ncols) 
{
  char **newArray;
  int i;

  if((newArray = (char **)malloc(nrows*sizeof(char*)))==NULL){
    fprintf(stderr,"Out of memory allocating character array\n");
    exit(0);
  }
  for(i=0;i<nrows;i++) {
    if((newArray[i] = (char *)malloc((ncols+1)*sizeof(char)))==NULL){
      fprintf(stderr,"Out of memory allocating character array\n");
      exit(0);
    }
    newArray[i][ncols]='\0'; /* to terminate with null char so 
				newArray[i] behaves like a string */
  }
  return(newArray);
}



double ** make_double_mat(int nrows, int ncols)
{
  int i;
  double **mat;

  if((mat = (double **)malloc(nrows*sizeof(double *)))==NULL){
    fprintf(stderr,"Out of memory allocating double-precision matrix\n");
    exit(0);
  }
  for(i=0;i<nrows;i++)
    if((mat[i] = (double *)malloc(ncols*sizeof(double)))==NULL){
      fprintf(stderr,"Out of memory allocating double-precision matrix\n");
      exit(0);
    }
  return(mat);
}

int ** make_int_mat(int nrows, int ncols)
{
  int i;
  int **mat;

  if((mat = (int **)malloc(nrows*sizeof(int *)))==NULL){
    fprintf(stderr,"Out of memory allocating integer matrix\n");
    exit(0);
  }
  for(i=0;i<nrows;i++)
    if((mat[i] = (int *)malloc(ncols*sizeof(int)))==NULL){
      fprintf(stderr,"Out of memory allocating integer matrix\n");
      exit(0);
    }
  return(mat);
}

int *copyIntVec(int *vec, int len)
{
  int i;
  int *newvec = (int *)malloc(len*sizeof(int));

  for(i=0;i<len;i++)
    newvec[i]=vec[i];

  return(newvec);
}
    

int *make_initial_index( )
{
  int j, currBase=0;
  int *index=(int*)(malloc(numBases*sizeof(int)));

  for(j=0;j<numBases; j++) {
    if(!anyGapsAtSite(j) && isVariable(j)) {
      index[currBase]=j;
      currBase++;
    }
  }
  numBases=currBase; /*number of bases after removing gaps*/

  return(index);
}


int anyGapsAtSite(int j)
{
  int i;

  for(i=0;i<nseqs;i++)
    if(sequences[i][j]=='-')
      return(TRUE);
  return(FALSE); /*If we get here, no gaps found*/
}

int isVariable(int j)
{
  int i;

  for(i=1;i<nseqs;i++)
    if(sequences[i][j]!=sequences[0][j])
      return(TRUE); /*we've found variation*/
  return(FALSE); /*If we get here, all bases at the site are the same*/
}

/*******************************************************************/
/* code to reduce the original alignment to unique sequences only  */
/*******************************************************************/

void reduceToUniqueSeqs()
{
  int *uniqueSeqInd, i;

  /*First find the indeces of the unique sequences in the alignment*/
  uniqueSeqInd=findUniqueSeqs(); /*alters global variable nseqs to be the
                                   number of unique sequences*/

  /*Now cycle thru alignment and ``remove'' duplicates by copying 
    the unique sequences into the first nseqs rows of the sequences 
    matrix that holds the alignment.*/
  for(i=0;i<nseqs;i++){
    if(uniqueSeqInd[i]!=i){ /*then we need to do some copying*/
      strcpy(sequences[i],sequences[uniqueSeqInd[i]]); 
          /*copies sequences[uniqueSeqInd[i]] into sequences[i]*/
       strcpy(sequenceLabels[i],sequenceLabels[uniqueSeqInd[i]]); 
    }
  }

  free(uniqueSeqInd);
}

/* Return a vector with the indeces of the unique sequences in the alignment*/
int *findUniqueSeqs()
{
  int i, uniqueSeqNum;
  int *uniqueSeqInd = (int*)malloc(nseqs*sizeof(int));

  uniqueSeqInd[0]=0;  /*first seq. in alignment is first in list of uniques*/
  uniqueSeqNum=1; /*1 unique so far*/

  for(i=1;i<nseqs;i++) {
    if(newSeq(sequences[i],uniqueSeqInd,uniqueSeqNum)) {
      uniqueSeqInd[uniqueSeqNum]=i;
      uniqueSeqNum++;
    }
  }
  nseqs=uniqueSeqNum; /*number of unique sequences*/

  return(uniqueSeqInd);
}

/* Test if sequence is a new one. Search list of unique sequences found
   so far. If found return FALSE, if not return TRUE. */
int newSeq(char *sequence, int *uniqueSeqInd, int uniqueSeqNum)
{
  int i;

  for(i=0;i<uniqueSeqNum;i++)
    if(strcmp(sequences[uniqueSeqInd[i]],sequence)==0) /*then they match*/
      return(FALSE);

  return(TRUE); /*sequence not found so must be a new one*/
}

int *findEndpoints(int *polyposn, int *breaks, int numBreaks)
{
  int i, counter=0, *endpoints;

  endpoints=(int*) malloc((numBreaks+2)*sizeof(int));
  endpoints[0]=0; 
  if(numBreaks>0){
    for(i=1;i<(numBases-1);i++) /*no need to check ends of sequence*/
      if(polyposn[i]<=breaks[counter] && polyposn[i+1] > breaks[counter]) {
        endpoints[counter+1]=i;
        counter++;
      }
  }
  endpoints[numBreaks+1]=numBases-1; 

  return(endpoints);
}


void free_double_mat(double **mat, int nrows)
{
  int i;

  for(i=0;i<nrows;i++)
    free(mat[i]);
  free(mat);
}

void free_int_mat(int **mat, int nrows)
{
  int i;

  for(i=0;i<nrows;i++)
    free(mat[i]);
  free(mat);
}

void freeCharArray(char **cArray, int nrows) 
{
  int i;

  for(i=0;i<nrows;i++)
    free(cArray[i]);
  free(cArray);
}



/* Next two functions from phylip */
int eof(FILE *f)
{
    register int ch;

    if (feof(f))
        return 1;
    if (f == stdin)
        return 0;
    ch = getc(f);
    if (ch == EOF)
        return 1;
    ungetc(ch, f);
    return 0;
}

int eoln(FILE *f)
{
    register int ch;

    ch = getc(f);
    if (ch == EOF)
        return 1;
    ungetc(ch, f);
    return (ch == '\n');
}

/***********************************************************************/ 
/*  Allocate two-demensional int array dynamically                     */
int **dynamicArray(int nrows, int ncolumns) {
     int i;
     int **array = (int **)malloc(nrows * sizeof(int *));
     if(!array) {
        printf("Out of memory\n");
        exit(1);
     }
     array[0] = (int *)malloc(nrows * ncolumns * sizeof(int));
     if (!array[0]){
        printf("Out of memory\n");
        exit(1);
     }
     for(i = 1; i < nrows; i++)
        array[i] = array[0] + i * ncolumns;
     return array;
}

