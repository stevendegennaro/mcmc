#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "evolve.h"
#include "structures.h"
#include "loadModels.h"

#define  COL_MAX        1000                // max number of cluster stars 
#define  FS_NUM          0.0		            // mcmc outputs negative mass for field star 
#define  CLUS_READ         9
#define  CLUS_STAT_WRITE  10
#define  MASS1_FILE       11
#define  MASS2_FILE       12
#define  STARS_STAT_WRITE 13
#define  CMD_FILE         16

// Used by evolve.c
double ltau[2], wdLogTeff[2];
int aFilt=-1;

// Used by a bunch of different functions.
int    verbose, needMassNow=1, useFilt[FILTS];//, numFilts;

void openOutputFiles(FILE **filePtr, char *filename, int fileType);
double intlFinalMassReln(double zamsMass, int IFMR);

int main(void)

{
  
  /////////////////////////////////
  /////// Declare Variables ///////
  /////////////////////////////////
  
  int    j, nr, p, m, row, filt, isClusterMember[2][COL_MAX], cm[COL_MAX],
  stuckMass[2][COL_MAX], stuckParam[NPARAMS], tempCols, useParam[NPARAMS],
  numStepSizeMass[2][COL_MAX],iMag;
  double mass[2][COL_MAX], prevMass[2][COL_MAX], mean = 0.0,minMag,maxMag,
  sumMass[2][COL_MAX], sumSquaresMass[2][COL_MAX], minMass[2][COL_MAX], maxMass[2][COL_MAX], 
  thisStep, meanStepSizeMass[2][COL_MAX], meanMass[2][COL_MAX], variance=0.0, sigma, N, csProb;
  double param[NPARAMS], meanParam[NPARAMS],varParam[NPARAMS],minParam[NPARAMS],maxParam[NPARAMS],
  prevParam[NPARAMS], sumStepSizeParam[NPARAMS];
  double sumPhot[COL_MAX][FILTS],sumSquaresPhot[COL_MAX][FILTS];
  char   line[240], filename[100], **starName;
  
  const char *paramNames[] = {"logAge    ",
    "Y         ",
    "[Fe/H]    ",
    "modulus   ",
    "absorption"};
  
  FILE   *rMassPtr[2], *rClusterPtr, *rDataPtr;
  FILE   *wClusterStatPtr, *wStarsStatPtr, *wCmdPtr;
  
  struct cluster theCluster;
  //struct star *stars;
  
  initCluster(&theCluster);
  
  
  ////////////////////////////////////////////
  /////// Open files to read and write ///////
  ////////////////////////////////////////////
  
  printf("Enter file name containing color/magnitude data:\n> ");
  scanf("%s",filename);
  if((rDataPtr = fopen(filename,"r")) == NULL) {
    printf("***Error: file %s was not found.***\n",filename);
    printf("[Exiting...]\n");
    exit(1);
  }
  
  printf("Enter minimum and maximum magnitude of MS to use and band (0,1, or 2):\n> ");
  scanf("%lf %lf %d", &minMag, &maxMag,&iMag);
  if(iMag<0||iMag>2){
    printf("***Error: %d not a valid magnitude index.  Choose 0,1,or 2.***\n",iMag);
    printf("[Exiting...]\n");
    exit(1);
  }   
  
  printf("\n Enter mcmc file name : ");
  scanf("%s",filename);
  openOutputFiles(&rClusterPtr, filename, CLUS_READ);
  openOutputFiles(&wClusterStatPtr, filename, CLUS_STAT_WRITE);
  openOutputFiles(&wStarsStatPtr, filename, STARS_STAT_WRITE);
  openOutputFiles(&wCmdPtr, filename, CMD_FILE);
  
  
  ///////////////////////////////////////
  //////// Determine how many stars /////
  ////////// and allocate memory ////////
  ///////////////////////////////////////
  
  for(m = 0 ; m < 2 ; m++){
    openOutputFiles(&rMassPtr[m], filename, MASS1_FILE + m);
    fscanf(rMassPtr[m],"%*s %*s %*s %*s %*s %*s %d",((m) ? &(theCluster.nStars) : &tempCols));           // Read in the number of stars 
  }
  
  if(theCluster.nStars < 1) {
    printf("\n need at least a one column file - exiting.\n");
    exit(1);
  }
  if(theCluster.nStars > COL_MAX) {
    printf("\n Exceeded the limit of the number of cluster stars in a file - exiting.\n");
    exit(1);
  }
  
  if(theCluster.nStars != tempCols){
    printf("\n Different numbers of stars in .mass1 and .mass2 files");
    exit(1);
  }
  
  if((theCluster.stars = (struct star *) malloc(theCluster.nStars * sizeof(struct star))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");
  
  if((starName = (char **) calloc(theCluster.nStars, sizeof(char*))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");
  for(j = 0; j < theCluster.nStars; j++) {
    if((starName[j] = (char *) calloc(100 , sizeof(char))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
  }
  
  for(j = 0; j < theCluster.nStars; j++) initStar(&(theCluster.stars[j]));
  
  ////////////////////////////////////////
  //// Determine which parameters are ////
  //// included in the .cluster file /////
  ////////////////////////////////////////
  
  // Skip first line of .cluster file
  fgets(line,240,rClusterPtr);
  
  //Read in second line
  fscanf(rClusterPtr,"%*s ");
  for(p=0;p<NPARAMS;p++){
    fscanf(rClusterPtr,"%d %lf ",&(useParam[p]),&(param[p]));
  }
  fgets(line,240,rClusterPtr);
  
  
  ////////////////////////////////////////////
  ////// Determine which filters to use //////
  ////// read header from .scatter file //////
  ////////////////////////////////////////////
  
  /*   
   fgets(line, 300, rDataPtr);                               // skip first header line
   
   for(filt=0;filt<FILTS;filt++){
   fscanf(rDataPtr,"%d ",&useFilt[filt]);
   if(useFilt[filt]){
   theCluster.evoModels.numFilts++;
   //printf("** %d **\n",theCluster.evoModels.numFilts);
   aFilt = filt;                          // Sets this to a band we know we are using (for evolve) 
   }
   }
   fgets(line,240,rDataPtr);				// and next header line 
   */
  
  char *pch, sig[]="sig\0", comp[]="   \0";
  
  //Parse the header of the file to determine which filters are being used
  fgets(line, 300, rDataPtr); // Read in the header line
  //printf("%s",line);
  
  // split the string on these delimiters into "tokens"
  pch = strtok (line," "); // Ignore the first token, which is "id", and move 
  pch = strtok (NULL, " "); // to the next (which should be the first filter name)
  printf("filterSet = %d\n",findFilterSet(pch));   // figure out which filter set we're using
  filt=0;
  
  while (pch != NULL)
  {                              
    strncpy(comp,pch,3);                  // copy the first three letters into the dummy string 'comp'
    if(strcmp(comp,sig)==0) break;        // and check to see if they are 'sig'.  If they are, there are no more filters
    
    for(filt=0;filt<FILTS;filt++){                 // Otherwise check to see what this filter's name is
      //printf("%d * %s * %s *\n", filt, pch, getFilterName(filt));
      if(strcmp(pch,getFilterName(filt))==0){
        useFilt[filt] = 1;
        theCluster.evoModels.numFilts++;
        if (aFilt < 0) aFilt = filt;                         // Sets this to a band we know we are using (for evolve)
        break;
      }
    }
    pch = strtok (NULL, " ");
  }
  
  //for(filt=0;filt<FILTS;filt++) printf("%d ",useFilt[filt]);
  //printf("\n");
  
  ///////////////////////////////////
  ///// Read in other variables /////
  ///////// and load models /////////
  ///////////////////////////////////
  
  printf("\n Enter M_wd_up : ");
  scanf("%lf", &theCluster.M_wd_up);
  
  printf("\n Run in verbose mode (0=no, 1=yes, 2=YES) ?");
  scanf("%d",&verbose);
  if(verbose < 0 || verbose > 2) verbose = 1;		// give standard feedback if incorrectly specified 
  
  chooseModels(&theCluster);
  loadModels(&theCluster);                                      // read in stellar evol & WD models 
  
  //////////////////////////////
  /////// Output headers ///////
  //////////////////////////////
  
  fprintf(wStarsStatPtr, "     Star    ");
  for(m = 1 ; m <= 2 ; m++){
    fprintf(wStarsStatPtr, "meanmass%d  sigmass%d   minmass%d   maxmass%d stuck%d meanStep%d ",m,m,m,m,m,m);
    if(m==1) fprintf(wStarsStatPtr, "meanWDMass%d sigWDmass%d ",m,m);
    fprintf(wStarsStatPtr, "CSprob%d Niter%d   ",m,m);
  }
  for(filt = 0; filt < FILTS; filt++){
    if(useFilt[filt]){
      char tmpString[100] = "mean\0";
      strcat(tmpString,getFilterName(filt));
      fprintf(wStarsStatPtr,"%8s ",tmpString);
      tmpString[0] ='\0';
      strcat(tmpString,"var");
      strcat(tmpString,getFilterName(filt));
      fprintf(wStarsStatPtr,"  %8s   ",tmpString);
    }
  }
  fprintf(wStarsStatPtr," pop\n");
  
  
  ////////////////////////////////
  ///// Initialize variables /////
  ////////////////////////////////
  
  m=0;
  nr=0;
  row=0;
  for(p = 0; p < NPARAMS; p++){
    minParam[p] = HUGE_VAL;
    maxParam[p] = -HUGE_VAL;
	  stuckParam[p] = 0;
	  meanParam[p] = 0;
	  varParam[p] = 0;
	  sumStepSizeParam[p] = 0;
  }
  for(j = 0 ; j < theCluster.nStars ; j++){
    for(filt = 0 ; filt < FILTS ; filt++){
      sumPhot[j][filt] = 0.0;
      sumSquaresPhot[j][filt]= 0.0;
    }
  }
  
  for(m = 0 ; m < 2 ; m++){
    for(j = 0 ; j < COL_MAX ; j++) {                                     // initialize variables 
      sumMass[m][j]          = 0.0;
      sumSquaresMass[m][j]   = 0.0;
      isClusterMember[m][j]  = 0;
      minMass[m][j]          = 100;                                 // impossibly high min mcmc value for each mass 
      maxMass[m][j]          = -100;                                // impossibly low max mcmc value for each mass 
      prevMass[m][j]         = -100.;                                 // impossible value, should never be used 
      stuckMass[m][j]        = 0;
      meanStepSizeMass[m][j] = 0.;
      numStepSizeMass[m][j]  = 0;
    }
  }
  
  /////////////////////////////////
  ///// Determine which stars /////
  /////   were used by mcmc   /////
  /////////////////////////////////
  for(j = 0; j < theCluster.nStars; j++){
    do{
      if((nr = fscanf(rDataPtr,"%s ", starName[j])) == EOF) break;
      for(p=0;p<theCluster.evoModels.numFilts;p++) fscanf(rDataPtr, "%lf ",&(theCluster.stars[j].obsPhot[p]));
      for(p=0;p<theCluster.evoModels.numFilts;p++) fscanf(rDataPtr, "%*f ");
      fscanf(rDataPtr, "%*f %*f %d %*d %*d",&theCluster.stars[j].status[0]);
      fgets(line,240,rDataPtr);
      //printf("%s %f %f %f\n",starName[j], theCluster.stars[j].obsPhot[iMag], minMag, maxMag);
    } while(theCluster.stars[j].status[0]==MSRG && (theCluster.stars[j].obsPhot[iMag] < minMag || theCluster.stars[j].obsPhot[iMag] > maxMag));
    //printf("%s \n",starName[j]);
  }
  //exit(1);
  
  
  ////////////////////////////
  ///////// Main loop ////////
  ////////////////////////////

  
  while(1){
    // Read a line in from both mass files
    for(j = 0; j < theCluster.nStars; j++){
      if((nr = fscanf(rMassPtr[0],"%lf",&mass[0][j])) == EOF) break;
      fscanf(rMassPtr[1],"%lf",&mass[1][j]);
      //printf("%d %f %f\n",j,mass[0][j],mass[1][j]);
    }    
    if(nr == EOF) break;
    
    fscanf(rClusterPtr,"%*d "); //Skip the first entry in each line of .cluster file (which is just iteration #)
    // Store stats for each cluster parameter
    for(p = 0; p < NPARAMS; p++){ 
      if(useParam[p]) fscanf(rClusterPtr,"%lf ",&param[p]);
      meanParam[p] = (row*meanParam[p] + param[p])/(row + 1);        // Calculate mean and variance recursively 
      if(row > 0) varParam[p]  = row*varParam[p]/(row + 1) + pow((param[p] - meanParam[p]),2)/row;
      if(param[p] > maxParam[p]) maxParam[p] = param[p];
      if(param[p] < minParam[p]) minParam[p] = param[p];
      
      thisStep = fabs(param[p] - prevParam[p]);                      // accumulate correlation and sticking stats 
      if(thisStep < EPS) stuckParam[p]++;
      sumStepSizeParam[p] += thisStep;
      prevParam[p] = param[p];
      theCluster.parameter[p] = param[p];
    }
    
    
    // Store stats for each mass for each star
    for(j = 0 ; j < theCluster.nStars ; j++) {
      if(mass[0][j] > EPS) cm[j] = 1;
      else cm[j] = 0;
      if(cm[j]){                                                      // if it is a cluster star 
        for(m = 0 ; m < 2 ; m++){                                     // For each component 
          sumMass[m][j]        += mass[m][j];                         // for mean, sigma, N 
          sumSquaresMass[m][j] += mass[m][j] * mass[m][j];
          
          if(minMass[m][j] > mass[m][j]) minMass[m][j] = mass[m][j];  // min & max 
          if(maxMass[m][j] < mass[m][j]) maxMass[m][j] = mass[m][j];
          
          if(prevMass[m][j]>0){
            thisStep = fabs(mass[m][j] - prevMass[m][j]);               // accumulate correlation and sticking stats 
            if(thisStep < EPS) stuckMass[m][j]++;
            meanStepSizeMass[m][j] += thisStep;
            numStepSizeMass[m][j]++;
          }
          
          prevMass[m][j] = mass[m][j];
          isClusterMember[m][j]++;                                    // Keep track of how many iterations we're averaging for each star 
          
          theCluster.stars[j].U = mass[0][j];
          theCluster.stars[j].massRatio = mass[1][j]/mass[0][j];
        }        
      }
    }
    evolve(&theCluster,-1);
    for(j = 0 ; j < theCluster.nStars ; j++) {
      if(cm[j]){                                                      // if it is a cluster star 
        for(filt = 0 ; filt < theCluster.evoModels.numFilts ; filt++){
          sumPhot[j][filt] += theCluster.stars[j].photometry[filt];
          sumSquaresPhot[j][filt] +=  theCluster.stars[j].photometry[filt] * theCluster.stars[j].photometry[filt];
        }
      }
    }
    row++;
  }
  
  /////////////////////////////////////////////
  //////// Calculate stats and output /////////
  /////////////////////////////////////////////
  
  // Calculate mass stats for each star, for each component
  for(j = 0 ; j < theCluster.nStars ; j++) {                                  // calc mean, sigma values for each star 
    
    fprintf(wStarsStatPtr,"%9s  ",starName[j]);
    for(m = 0 ; m < 2 ; m++){
      N = (double) isClusterMember[m][j];
      if(N > 1.0) {
        variance = (N*sumSquaresMass[m][j] - sumMass[m][j]*sumMass[m][j]) / (N * (N-1));
        sigma    = sqrt(fabs(variance));
      }
      else sigma = 0.0;
      if(N >= 0.5) {
        meanMass[m][j]            = sumMass[m][j] / N;                      // else mean = 0, anyway 
        if(numStepSizeMass[m][j] > 0) meanStepSizeMass[m][j]   /= numStepSizeMass[m][j];
      }
      
      csProb = N / row;                                                      // row = mcmc iterations, N = # of times it was a cluster member 
      
      // Output mass stats
      fprintf(wStarsStatPtr,"%11.5f %9.5f %10.5f %10.5f %6d  %8.6f  ",
              meanMass[m][j], sigma, minMass[m][j], maxMass[m][j], stuckMass[m][j], meanStepSizeMass[m][j]);
      if(m==0){
        if(theCluster.stars[j].status[0]==WD) 
          fprintf(wStarsStatPtr,"%10.5f  %9.3e  ", intlFinalMassReln(meanMass[m][j], theCluster.evoModels.IFMR), 
                  intlFinalMassReln(meanMass[m][j], theCluster.evoModels.IFMR) -
                  intlFinalMassReln(meanMass[m][j] - sigma, theCluster.evoModels.IFMR));
        else
          fprintf(wStarsStatPtr,"%10.5f  %9.3e  ",0.0,0.0);
      }
      fprintf(wStarsStatPtr,"%6.4f %6.0f ", csProb, N);
      
      
      // Calculate photometry stats
      if(m==1){
        //printf("*** %d ***\n",theCluster.evoModels.numFilts);
        for(filt=0;filt<theCluster.evoModels.numFilts;filt++){
          if(N > 1.0) {
            variance = (N*sumSquaresPhot[j][filt] - sumPhot[j][filt]*sumPhot[j][filt]) / (N * (N-1));
          }
          if(N >= 0.5) {
            mean           = sumPhot[j][filt] / N;
          }
          else mean = 99.999;
          //printf("***%10f %10.4e***",mean, variance);
          fprintf(wStarsStatPtr,"%10f %10.4e ",mean < 99 ? mean : 99.999, mean < 99 ? variance : 0.0);
        }
        fprintf(wStarsStatPtr," 1 \n");
      }
    }
  }
  
  // Output stats for cluster parameters 
  fprintf(wClusterStatPtr, "Parameter     mean     sigma       min      max    stuck  meanStep\n");
  for(p = 0 ; p < NPARAMS ; p++) fprintf(wClusterStatPtr,"%s %8.5f %9.5f %10.5f %8.5f %6d  %8.5f\n",paramNames[p],meanParam[p],
                                         sqrt(varParam[p]),minParam[p],maxParam[p],stuckParam[p], sumStepSizeParam[p]/(row - 1));
  
  
  ////////////////////////////////////
  //// Create a simulated cluster ////
  //// based on the cluster stats ////
  //// and output for comparison  ////
  ////////////////////////////////////
  
  //theCluster.nStars = (int)(theCluster.M_wd_up * 100) + 1;
  for(p=0;p<NPARAMS;p++) theCluster.parameter[p] = meanParam[p];
  
  //void *tempAlloc;		// temporary for allocation
  //if((tempAlloc = (void *) realloc(theCluster.stars, theCluster.nStars * sizeof(struct star))) == NULL)
   // perror("MEMORY ALLOCATION ERROR \n");
  //else
  //  theCluster.stars = (struct star *) tempAlloc;
  
  //for(j=0;j<theCluster.nStars;j++){
  //  theCluster.stars[j].U = j*.01;
   // theCluster.stars[j].massRatio = 0.0;
 // }
  
  returnIso(&theCluster);
  
  fprintf(wCmdPtr," mass1 stage1");
  for(filt = 0 ; filt < FILTS ; filt++) if(useFilt[filt]) fprintf(wCmdPtr,"          %s",getFilterName(filt));
  fprintf(wCmdPtr,"\n");
  
  double prevV=100;
  for(j=0;j<theCluster.nStars;j++){
          fprintf(wCmdPtr,"%6.2f %6d ",theCluster.stars[j].U,theCluster.stars[j].status[0]);
          for(filt = 0 ; filt < theCluster.evoModels.numFilts ; filt++) fprintf(wCmdPtr,"%10f ",theCluster.stars[j].photometry[filt]);
          fprintf(wCmdPtr,"\n");
  }
  
  /*
  double prevV=100;
  for(j=0;j<theCluster.nStars;j++){
    if(theCluster.stars[j].photometry[0]<90){
      if(theCluster.stars[j].status[0]==MSRG){
        if(prevV > theCluster.stars[j].photometry[0]){
          fprintf(wCmdPtr,"%6.2f %6d ",theCluster.stars[j].U,theCluster.stars[j].status[0]);
          for(filt = 0 ; filt < theCluster.evoModels.numFilts ; filt++) fprintf(wCmdPtr,"%10f ",theCluster.stars[j].photometry[filt]);
          fprintf(wCmdPtr,"\n");
          prevV = theCluster.stars[j].photometry[0];
        }
        else{
          prevV = 0.0;
          continue;
        }
      }
      else{
        fprintf(wCmdPtr,"%5.2f %3d ",theCluster.stars[j].U,theCluster.stars[j].status[0]);
        for(filt = 0 ; filt < theCluster.evoModels.numFilts ; filt++) fprintf(wCmdPtr,"%10f ",theCluster.stars[j].photometry[filt]);
        fprintf(wCmdPtr,"\n");
      }
    }
  }
  */
  
  fclose(wClusterStatPtr);
  fclose(rClusterPtr);
  fclose(rMassPtr[0]);
  fclose(rMassPtr[1]);
  fclose(wStarsStatPtr);
  return(0);
}


void openOutputFiles(FILE **filePtr, char *filename, int fileType)
{
  char tmpfile[100];
  const char *mode = "r";
  
  strcpy(tmpfile,filename);
  switch (fileType) {
      // output files differentiated by file name extension
    case CLUS_READ        : strcat(tmpfile,".cluster");
      break;
    case MASS1_FILE       : strcat(tmpfile,".mass1");
      break;
    case MASS2_FILE       : strcat(tmpfile,".mass2");
      break;
    case CLUS_STAT_WRITE  : strcat(tmpfile,".cluster.stat");
      mode = "w";
      break;
    case STARS_STAT_WRITE : strcat(tmpfile,".stars.stat");
      mode = "w";
      break;
    case CMD_FILE         : strcat(tmpfile,".cmd");
      mode = "w";
      break;
    default               : printf("***Error: Bad file choice in openOutputFiles().***\n");
      printf("[Exiting...]\n");
      exit(0);
      
  }
  
  printf("Reading file  : %s (%s)\n",tmpfile,mode);
  if((*filePtr = fopen(tmpfile,mode)) == NULL) {
    printf("***Error: File %s was not available. %s ***\n",tmpfile,mode);
    printf("[Exiting...]\n");
    exit(1);
  }
}
