#include <stdio.h> 
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "mcmc.h"
#include "evolve.h"
#include "loadModels.h"
#include "densities.h"
#include "decide.h"
#include "samplers.h"
#include "leastSquares.h"
#include "mt19937ar.h"
#include "linInterp.h"

////////////////////////////////////////////
//////// declare global variables //////////
////////////////////////////////////////////

// TO BE MADE TO INPUT LATER
double filterPriorMin[FILTS];
double filterPriorMax[FILTS];

// Used by other methods in mcmc.c (can be made static?)
int    B0, B2, T, inloopback=0, inagebeta=0;
char   w_file[100];

// Used by evolve.c
double ltau[2];
int aFilt=-1;

// Used in densities.c.
//double priorMean[NPARAMS];//, priorVar[NPARAMS];
extern double ageLimit[2]; // Defined in evolve.c, set in the appropriate model during loadModels.

// Used by a bunch of different functions.
int    verbose = 0, needMassNow = 0, useFilt[FILTS];//, numFilts = 0;

// For random number generator (mt19937ar.c)
unsigned long  mt[N], seed;
int            mti=N+1;


////////////////////////////////////////////
////////////////////////////////////////////
// MAIN FUNCTION
////////////////////////////////////////////
////////////////////////////////////////////
int main()
{
  printf("\n ***You are running BASE08/MCMC version %.1f.***\n",VERSION);
  printf("\n\nMCMC: Markov chain Monte Carlo sampler\n");
  printf("      for Bayesian analysis of star cluster data\n\n");
  
  struct chain mc;
  struct mcmcControl ctrl;
  
  initCluster(&mc.clust);
  initStepSizes(&mc.clust);
  initMcmcControl(&mc, &ctrl);
   
  struct block initialWanderBlock      = {INITIAL_WANDER_BLOCK, runWander, ctrl.B0};
  struct block massRatioBlock          = {MASS_RATIO_BLOCK, runMassVsMassRatio, ctrl.T};
  struct block ageBlock                = {AGE_BLOCK, runMassVsAge, ctrl.T};
  struct block modulusBlock            = {MODULUS_BLOCK, runModulus, ctrl.T};
  struct block fehBlock                = {FEH_BLOCK, runFeH, ctrl.T};
  struct block yBlock                  = {Y_BLOCK, runY, ctrl.T};
  struct block turnFSSamplingOffBlock  = {TURN_FS_SAMPLING_OFF_BLOCK, runFSSamplingOff, 0};
  struct block fieldStarPriorBlock     = {FIELD_STAR_PRIOR_BLOCK, runFieldStarPriorCalc, ctrl.T};
  struct block turnFSSamplingOnBlock   = {TURN_FS_SAMPLING_ON_BLOCK, runFSSamplingOn, 0};
  struct block finalWanderBlock        = {FINAL_WANDER_BLOCK, runWander, ctrl.B2};
  struct block endBurnInBlock          = {END_BURN_IN_BLOCK, runEndBurnIn, 0};
  struct block postBurnInBlock         = {POST_BURN_IN_BLOCK, runPostBurnIn, ctrl.nIter}; 

  // Initialize array of chain blocks
  // This defines the sequence of the burn-in 
  struct block chainBlocks[] = {
    initialWanderBlock,
    massRatioBlock,
    ageBlock,
    ageBlock,
    modulusBlock,
    modulusBlock,
    fehBlock,
    yBlock,
    turnFSSamplingOffBlock,
    fieldStarPriorBlock,
    turnFSSamplingOnBlock,
    initialWanderBlock,
    massRatioBlock,
    ageBlock,
    modulusBlock,
    modulusBlock,
    fehBlock,
    yBlock,
    turnFSSamplingOffBlock,
    fieldStarPriorBlock,
    turnFSSamplingOnBlock,
    finalWanderBlock,
    endBurnInBlock,
    postBurnInBlock
  };
  if (!ctrl.runBurnIn) {
    if (ctrl.rBetaInFile == NULL) {
      printf("***Error:  Can't skip burnin unless you have a .beta.in file\n");
      exit(0);
    }
    else {
      // only run the postBurnIn
      chainBlocks[0] = postBurnInBlock;
    }
  }
  
  
  setBurnIter(&ctrl, chainBlocks);
  readCmdData(&mc, &ctrl);
  allocateChainMemory(&mc);
  initChain(&mc, &ctrl);
  printOutputFileHeaders(&mc, &ctrl);
  printf("[Running MCMC...]\n");
  
  //int end = 0;
  
  /********* MAIN LOOP ********/
  int iter = 0;
  int iBlock = 0;
  int moreBlocks = 1;
  while (moreBlocks) {
    ctrl.currentBlock = chainBlocks[iBlock].type;
    chainBlocks[iBlock].run(&iter, chainBlocks[iBlock].nIter, &mc, &ctrl);
    if (chainBlocks[iBlock].type == POST_BURN_IN_BLOCK){
      //if(end)
        moreBlocks = 0;
      //else end = 1;
    }
    else
      iBlock++;
  }
  /********* END MAIN LOOP ********/
	  
  closePostBurnInFiles(&ctrl);
  freeChainArrays(&mc);
  return(0);
}
////////////////////////////////////////////
////////////////////////////////////////////
/////////// END MAIN FUNCTION //////////////
////////////////////////////////////////////
////////////////////////////////////////////



/////////////////////////////////////
// Close post burn-in output files //
/////////////////////////////////////
void closePostBurnInFiles(struct mcmcControl *ctrl) {
  fclose(ctrl->rData);
  fclose(ctrl->wMass1File[0]);
  fclose(ctrl->wMass2File[0]);
  fclose(ctrl->wClusterFile[0]);
  //fclose(ctrl->wDebugFile[0]);
  if(ctrl->rBetaInFile != NULL) fclose(ctrl->rBetaInFile);
}


///////////////////////
// Free chain arrays //
///////////////////////
void freeChainArrays(struct chain *mc) {
  free(mc->clust.stars);
  free(mc->acceptMass);
  free(mc->rejectMass);
  free(mc->acceptMassRatio);
  free(mc->rejectMassRatio);
  free(mc->isFieldStarNow);
  free(mc->isClusterStarNow);
}


/////////////////////////
////// Sample Step //////
/////////////////////////
void sampleStep(int iter, struct chain *mc, struct mcmcControl *ctrl, int SAMPLE_TYPE){
  int fsOnState = ctrl->fsSamplingOn;
  
  if(mc->clust.priorVar[SAMPLE_TYPE] > EPSILON) {
    // If calculating the correlation between age and mass (ie in AGE_BLOCK), then
    // don't sample YYY, FEH, MOD, or ABS
    if (ctrl->currentBlock != AGE_BLOCK || SAMPLE_TYPE <= AGE) {
      //if(ctrl->currentBlock != AGE_MODULUS_BLOCK || SAMPLE_TYPE != MOD){
      mc->clust = decideClust(mc, fsOnState, SAMPLE_TYPE);
      //}
    }
  }
  else (mc->acceptClust[SAMPLE_TYPE])++;


} // sampleStep


/////////////////////////
///// Open Outputs //////
/////////////////////////
void openOutputFiles(FILE **file_ptr, FILE **file_bi_ptr, char *filename, int length, int WHICH_FILE, int outputBurnIn)
{
  char tmpfile[100];
  char tmpfile2[100];
  if(length > 80) {
    printf("***Error: File name too long.***\n");
    printf("[Exiting...]\n");
    exit(1);
  }
  strcpy(tmpfile,filename);
  switch (WHICH_FILE) {
      // output files differentiated by file name extension
      case CLUS_FILE   : strcat(tmpfile,".cluster");
                         break;
      case MASS1_FILE  : strcat(tmpfile,".mass1");
                         break;
      case MASS2_FILE  : strcat(tmpfile,".mass2");
                         break;
      case DEBUG_FILE  : strcat(tmpfile,".debug");
                         break;
      case BETA_FILE   : strcat(tmpfile,".betas");
                         break;
      default          : printf("***Error: Bad file choice in openOutputFiles().***\n");
                         printf("[Exiting...]\n");
                         exit(0);

  }
  if(WHICH_FILE == BETA_FILE){
    strcpy(tmpfile2,tmpfile);
    strcat(tmpfile2,".in");
    if((*file_ptr = fopen(tmpfile2,"r")) == NULL) {
      printf("There is no beta.in file...'%s'\n",tmpfile2);
    }
    strcat(tmpfile,".out");
    if((*file_bi_ptr = fopen(tmpfile,"w")) == NULL) {
      printf("***Error: File %s was not available for writing.***\n",tmpfile);
      printf("[Exiting...]\n");
      exit(1);
    }
    return;
  }
  else{
    if((*file_ptr = fopen(tmpfile,"w")) == NULL) {
      printf("***Error: File %s was not available for writing.***\n",tmpfile);
      printf("[Exiting...]\n");
      exit(1);
    }
    // write burn-in iterations to *.burnin file
    if(outputBurnIn){
      strcat(tmpfile,".burnin");
      if((*file_bi_ptr = fopen(tmpfile,"w")) == NULL) {
        printf("***Error: File %s was not available for writing.***\n",tmpfile);
        printf("[Exiting...]\n");
        exit(1);
      }
    }
  }
} // openOutputFiles


//////////////////////////////////////////////////////////
// Set starting values for cluster parameter step sizes //
//////////////////////////////////////////////////////////
void initStepSizes(struct cluster *pClust)
{
   pClust->stepSize[AGE]    = 0.05;
   pClust->stepSize[FEH]    = 0.05;
   pClust->stepSize[MOD]    = 0.03;
   pClust->stepSize[ABS]    = 0.02;
   pClust->stepSize[YYY]    = 0.002;
} // initStepSizes


///////////////////////////
//// Adjust Step Sizes ////
///////////////////////////
double adjustStepSizes(int *accept, int *reject){

  int i, iupper = 4, ilower = 7;
  double fracAccept;
  double fa[7] = {0.9, 0.7, 0.5, 0.3, 0.05, 0.15, 0.2};
  double w[7] = {2.0, 1.8, 1.5, 1.2, 0.5, 1/1.8, 1/1.5};

  if(*accept + *reject < 200) return 1.0;

  fracAccept = *accept / (double) (*accept + *reject);
  *accept = 0;
  *reject = 0;

  for(i=0;i<iupper;i++){
    if(fracAccept > fa[i]) return w[i];
  }
  for(i = 4; i < ilower; i++){
    if(fracAccept < fa[i]) return w[i];
  }
  return 1.0;
} // adjustStepSizes


///////////////////////////////////////////////////////////////////
// Store masses and parameter values in masses and params arrays //
///////////////////////////////////////////////////////////////////
void storeBetaMass(int index, double **masses, double *params, struct chain *mc, struct mcmcControl *ctrl, int TYPE){
  int j;
  if((index % ctrl->increment) == 0) {
    params[index/ctrl->increment] = getParameter(&mc->clust,TYPE);
    for(j = 0; j < mc->clust.nStars; j++) {
      if(TYPE == AGE){
        masses[j][index/ctrl->increment] = getMass1(&mc->clust.stars[j], &mc->clust);
      }
      else{
        masses[j][index/ctrl->increment] = mc->clust.stars[j].U + 
          mc->clust.stars[j].beta[TYPE][0] * (mc->clust.parameter[TYPE] - mc->clust.mean[TYPE]);
      }
    }
  }
} // storeBetaMass


///////////////////////////////////
// Decorrelate masses and params //
/////////////////////////////////// 
void findBetaMass(double **masses, double *params, double tempMean, struct chain *mc, 
                  int TYPE, double lower, double upper)
{
  double massTemp[mc->clust.nStars];
  int j;

  for(j = 0; j < mc->clust.nStars; j++) {
    massTemp[j] = getMass1(&mc->clust.stars[j],&mc->clust);
    mc->clust.stars[j].beta[TYPE][0] = leastSquaresBeta(params, masses[j], N_LSQ_VALUES);
    if((mc->clust.stars[j].beta[TYPE][0] > upper) || (mc->clust.stars[j].beta[TYPE][0] < lower)) 
      mc->clust.stars[j].beta[TYPE][0] = 0.0;
    mc->acceptMass[j] = 0;
    mc->rejectMass[j] = 0;
  }
  mc->clust.mean[TYPE] = tempMean;
  for(j = 0; j < mc->clust.nStars; j++)
    setMass1(&mc->clust.stars[j],&mc->clust,massTemp[j]);
  
} // findBetaMass


/////////////////////////////////////////
// Calculate good mass starting values //
///////////////////////////////////////// 
void getStartMass(struct cluster *pCluster, int iMag){

  int j;
  struct cluster tempCluster = cloneCluster(pCluster);
  int m;
  
  returnIso(&tempCluster);

  // Find starting primary mass using magnitude
  for(j=0;j<pCluster->nStars;j++){
    if(pCluster->stars[j].useDuringBurnIn==0) pCluster->stars[j].U = 0.0;
    m = 0;
    while(pCluster->stars[j].obsPhot[iMag] < tempCluster.stars[m].photometry[iMag] 
          && m < tempCluster.nStars - 2 && tempCluster.stars[m+1].status[0] == 1) m++;
      m--;
    if(pCluster->stars[j].status[0] != WD && pCluster->stars[j].status[0] != BD){      
      // If the input mass is equal to zero, find a good starting mass
      if(fabs(pCluster->stars[j].U) < EPS){
        pCluster->stars[j].U = linInterpExtrap(tempCluster.stars[m].photometry[iMag], tempCluster.stars[m+1].photometry[iMag], 
                                               tempCluster.stars[m].U,                tempCluster.stars[m+1].U, 
                                               pCluster->stars[j].obsPhot[iMag]);
      }
      if(m==0) m++;
      //Set the initial mcmc step size
      if(fabs(tempCluster.stars[m+1].U - tempCluster.stars[m-1].U)/2 < pCluster->stars[j].UStepSize)
        pCluster->stars[j].UStepSize = fabs(tempCluster.stars[m+1].U - tempCluster.stars[m-1].U)/2;
    }
  }
  
  //for(j=1;j<tempCluster.nStars-2;j++) printf("%d %f\n",j,fabs(tempCluster.stars[j+1].U - tempCluster.stars[j-1].U)/2);
  
  evolve(pCluster,-1);
    
  freeCluster(&tempCluster);
} // getStartMass


////////////////////////////////////////////////////////////
// "Wandering" period: update chain and adjust step sizes //
////////////////////////////////////////////////////////////
void runWander(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl) {  
  int i, j;
  int startIter = *iter;
  
  for (; *iter < startIter + nIter; (*iter)++) {
     if (ctrl->outputBurnIn) printChainValues(*iter, mc, ctrl);
    for(j = 0; j < mc->clust.nStars; j++) {
      if(mc->clust.stars[j].isFieldStar) mc->isFieldStarNow[j]++;
      else mc->isClusterStarNow[j]++;
    }

    i = *iter - startIter;
        
    updateChain(*iter, mc, ctrl);

    if(i % 200 == 0) {
      adjustChainStepSizes(mc, ctrl);
    } 
  }
  
  //for(j=0;j<mc->clust.nStars;j++) printf("%d %f\n",j,mc->clust.stars[j].UStepSize);exit(1);
  
} // runWander


/////////////////////////////////////
// Decorrelate mass and mass ratio //
/////////////////////////////////////  
void runMassVsMassRatio(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl) {
  
  double tempMean[NPARAMS] ={0};
  
  int j;
  double **masses, **massRatios;
  //double *meanMassRatioTemp;
  if((masses = (double **) calloc(mc->clust.nStars, sizeof(double*))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");
  if((massRatios = (double **) calloc(mc->clust.nStars, sizeof(double*))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");
  for(j = 0; j < mc->clust.nStars; j++) {
    if((masses[j] = (double *) calloc(N_LSQ_VALUES , sizeof(double))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
    if((massRatios[j] = (double *) calloc(N_LSQ_VALUES , sizeof(double))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
  }
  //if((meanMassRatioTemp = (double *) calloc(mc->clust.nStars, sizeof(double))) == NULL)
  //  perror("MEMORY ALLOCATION ERROR \n");
  
    
  //////////////////////////////////////
  // Main loop for runMassVsMassRatio //
  //////////////////////////////////////
  int i;
  int startIter = *iter;
  for (; *iter < startIter + nIter; (*iter)++) {
    if (ctrl->outputBurnIn) printChainValues(*iter, mc, ctrl);
    for(j = 0; j < mc->clust.nStars; j++) {
      if(mc->clust.stars[j].isFieldStar) mc->isFieldStarNow[j]++;
      else mc->isClusterStarNow[j]++;
    }
    
    //Calculates temporary mean values so it can hold these factors constant in runMassVsAge
    tempMean[MOD] += getParameter(&mc->clust,MOD);
    tempMean[FEH] += getParameter(&mc->clust,FEH);
    tempMean[YYY] += getParameter(&mc->clust,YYY);
    tempMean[ABS] += getParameter(&mc->clust,ABS);
        
    i = *iter - startIter;
    for(j = 0; j < mc->clust.nStars; j++) {
      if(i % ctrl->increment == 0) {
        masses[j][i/ctrl->increment] = getMass1(&mc->clust.stars[j], &mc->clust);
        massRatios[j][i/ctrl->increment] = mc->clust.stars[j].massRatio;
      }
     // meanMassRatioTemp[j] += mc->clust.stars[j].massRatio;
    }
        
    updateChain(*iter, mc, ctrl);
    
    if(i % 200 == 0) {
      adjustChainStepSizes(mc, ctrl);
    } 
  }
  //////////////////////////////////////////
  // End main loop for runMassVsMassRatio //
  ////////////////////////////////////////// 
    
  tempMean[MOD] /= ctrl->T;
  tempMean[FEH] /= ctrl->T;
  tempMean[YYY] /= ctrl->T;
  tempMean[ABS] /= ctrl->T;

  double massTemp[mc->clust.nStars];
  for(j = 0; j < mc->clust.nStars; j++) {

    /* power law fit */
    //mc->clust.stars[j].meanMassRatio = meanMassRatioTemp[j] / ctrl->T;
    massTemp[j] = getMass1(&mc->clust.stars[j],&mc->clust);
    mc->clust.stars[j].betaMassRatio[1] = powerLaw(massRatios[j],masses[j],N_LSQ_VALUES,&(mc->clust.stars[j].betaMassRatio[0]));
  
    /* linear fit */
    // meanMassRatioTemp[j] /= ctrl->T;
    // massTemp = mc->clust.stars[j].U + mc->clust.stars[j].betaMassRatio[0] * 
    //   (mc->clust.stars[j].massRatio - mc->clust.stars[j].meanMassRatio);
    // 
    // mc->clust.stars[j].betaMassRatio[0] = leastSquaresBeta(massRatios[j],masses[j], N_LSQ_VALUES);
    // if(mc->clust.stars[j].betaMassRatio[0] > 0 || mc->clust.stars[j].betaMassRatio[0] < -100) mc->clust.stars[j].betaMassRatio[0] = 0.0;
    // mc->clust.stars[j].U = massTemp - mc->clust.stars[j].betaMassRatio[0] * 
    //   (mc->clust.stars[j].massRatio - meanMassRatioTemp[j]);
    // mc->clust.stars[j].meanMassRatio     = meanMassRatioTemp[j];
  }
  
  double tempAbs = getParameter(&mc->clust,ABS);
  
  mc->clust.mean[FEH] = tempMean[FEH];
  mc->clust.mean[MOD] = tempMean[MOD];
  mc->clust.mean[YYY] = tempMean[YYY];
  mc->clust.mean[ABS] = tempMean[ABS];
  
  setParameter(&mc->clust,ABS, tempAbs);
  
  for(j = 0; j < mc->clust.nStars; j++) {
    //meanMassRatioTemp[j]       = 0.0;
    mc->acceptMassRatio[j]     = 0;
    mc->rejectMassRatio[j]     = 0;
    mc->clust.stars[j].UStepSize    *= 2;
    setMass1(&mc->clust.stars[j],&mc->clust, massTemp[j]);
  }
  
  //free(meanMassRatioTemp);
  for(i = 0; i < mc->clust.nStars; i++) {
    free(masses[i]);
    free(massRatios[i]);
  }
  free(masses);
  free(massRatios);
} // runMassVsMassRatio


//////////////////////////////
// Decorrelate mass and age //
//////////////////////////////
void runMassVsAge(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl) {
  
  if (mc->clust.priorVar[YYY] > EPSILON) setParameter(&mc->clust, YYY, mc->clust.mean[YYY]);
  if (mc->clust.priorVar[FEH] > EPSILON) setParameter(&mc->clust, FEH, mc->clust.mean[FEH]);
  if (mc->clust.priorVar[MOD] > EPSILON) setParameter(&mc->clust, MOD, mc->clust.mean[MOD]);
  if (mc->clust.priorVar[ABS] > EPSILON) setParameter(&mc->clust, ABS, mc->clust.mean[ABS]);
  evolve(&mc->clust, -1);

  int j;
  double logAges[N_LSQ_VALUES] = {0};
  double tempMeanAge = 0.0;
  double **masses;
  if((masses = (double **) calloc(mc->clust.nStars, sizeof(double*))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");
  for(j = 0; j < mc->clust.nStars; j++) {
    if((masses[j] = (double *) calloc(N_LSQ_VALUES , sizeof(double))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
  }
  
  ////////////////////////////////
  // Main loop for runMassVsAge //
  ////////////////////////////////  
  int i;
  int startIter = *iter;
  for (; *iter < startIter + nIter; (*iter)++) {

    if (ctrl->outputBurnIn) printChainValues(*iter, mc, ctrl);
    for(j = 0; j < mc->clust.nStars; j++) {
      if(mc->clust.stars[j].isFieldStar) mc->isFieldStarNow[j]++;
      else mc->isClusterStarNow[j]++;
    }

    i = *iter - startIter;
    tempMeanAge += getParameter(&mc->clust,AGE);
    storeBetaMass(i, masses, logAges, mc, ctrl, AGE);
      
    updateChain(*iter, mc, ctrl);
    
    if(i % 200 == 0) {
      adjustChainStepSizes(mc, ctrl);
    } 
  }
  ////////////////////////////////////
  // End main loop for runMassVsAge //
  ////////////////////////////////////  
  
  tempMeanAge /= ctrl->T;
  findBetaMass(masses, logAges, tempMeanAge, mc, AGE, -HUGE_VAL, 0.0);
  mc->acceptClust[AGE] = 0;
  mc->rejectClust[AGE] = 0;
  
  for(j = 0; j < mc->clust.nStars; j++) {
    free(masses[j]);
  }
  free(masses);
} // runMassVsAge


//////////////////////////////////////////////////////////////
// Decorrelate mass and modulus, and absorption and modulus //
//////////////////////////////////////////////////////////////
void runModulus(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl) {  
  
  if(mc->clust.priorVar[MOD] > EPSILON) { 
    
    int j;
    double tempMeanMod = 0.0;
    double tempMeanAbs = 0.0;
    double moduli[N_LSQ_VALUES] = {0};
    double absorptions[N_LSQ_VALUES] = {0};
    double **masses;
    if((masses = (double **) calloc(mc->clust.nStars, sizeof(double*))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
    for(j = 0; j < mc->clust.nStars; j++) {
      if((masses[j] = (double *) calloc(N_LSQ_VALUES , sizeof(double))) == NULL)
        perror("MEMORY ALLOCATION ERROR \n");
    }
    
    //////////////////////////////
    // Main loop for runModulus //
    //////////////////////////////  
    int i;
    int startIter = *iter;
    for (; *iter < startIter + nIter; (*iter)++) {

      if (ctrl->outputBurnIn) printChainValues(*iter, mc, ctrl);
      for(j = 0; j < mc->clust.nStars; j++) {
        if(mc->clust.stars[j].isFieldStar) mc->isFieldStarNow[j]++;
        else mc->isClusterStarNow[j]++;
      }

      i = *iter - startIter;
      
      // store values
      tempMeanMod += getParameter(&mc->clust,MOD);
      storeBetaMass(i, masses, moduli, mc, ctrl, MOD);
      if(mc->clust.priorVar[ABS] > EPSILON) {
        if(i % ctrl->increment == 0) {
          absorptions[i/ctrl->increment] =  getParameter(&mc->clust,ABS); 
        }
        tempMeanAbs += getParameter(&mc->clust,ABS);
      }

      updateChain(*iter, mc, ctrl);

      if(i % 200 == 0) {
        adjustChainStepSizes(mc, ctrl);
      } 
    }
    //////////////////////////////////
    // End main loop for runModulus //
    //////////////////////////////////
        
    // Calculate beta[MOD][0] (correlation between modulus and each stellar mass) 
    // and increase step size (since it should be sampling better)
    tempMeanMod /= ctrl->T;    
    findBetaMass(masses, moduli, tempMeanMod, mc, MOD, -HUGE_VAL, HUGE_VAL);
    for(j = 0; j < mc->clust.nStars; j++){
      if(mc->clust.stars[j].status[0]!=WD && mc->clust.stars[j].beta[MOD][0] < 0){
        double tempMass = getMass1(&mc->clust.stars[j],&mc->clust);
        mc->clust.stars[j].beta[MOD][0] = 0.0;
        setMass1(&mc->clust.stars[j],&mc->clust,tempMass);
      }
      mc->clust.stars[j].UStepSize *= 2;
    }
    // if(!inloopback) mc->clust.stepSize[MOD] *= 2;
    if(!ctrl->fsSamplingOn) mc->clust.stepSize[MOD] *= 2;
    mc->acceptClust[MOD] = 0;
    mc->rejectClust[MOD] = 0;

    if(mc->clust.priorVar[ABS] > EPSILON){
      // Calculate betamabs (correlation between modulus and absorption)
      double absTemp = getParameter(&mc->clust,ABS);
      mc->clust.betamabs = leastSquaresBeta(moduli, absorptions, N_LSQ_VALUES);
      setParameter(&mc->clust,ABS,absTemp);
      
      // if(!inloopback) mc->clust.stepSize[ABS] *= 3;
      if(!ctrl->fsSamplingOn) mc->clust.stepSize[ABS] *= 3;
      mc->clust.mean[ABS] = tempMeanAbs / ctrl->T;
      mc->acceptClust[ABS] = 0;
      mc->rejectClust[ABS] = 0;
    }
      
    for(j = 0; j < mc->clust.nStars; j++) {
      free(masses[j]);
    }
    free(masses);
  }
  else *iter += nIter;
} // runModulus


/////////////////////////////////////////////////
// Decorrelate FeH and masses, and FeH and abs //
/////////////////////////////////////////////////
void runFeH(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl) {
  
  if(mc->clust.priorVar[FEH] > EPSILON) {
    
    int j;
    double tempMeanAbs = 0.0;
    double FeHs[N_LSQ_VALUES] = {0};
    double absorptions[N_LSQ_VALUES] = {0};
    double **masses;
    if((masses = (double **) calloc(mc->clust.nStars, sizeof(double*))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
    for(j = 0; j < mc->clust.nStars; j++) {
      if((masses[j] = (double *) calloc(N_LSQ_VALUES , sizeof(double))) == NULL)
        perror("MEMORY ALLOCATION ERROR \n");
    }
    
    //////////////////////////
    // Main loop for runFeH //
    //////////////////////////  
    int i;
    int startIter = *iter;
    for (; *iter < startIter + nIter; (*iter)++) {

      if (ctrl->outputBurnIn) printChainValues(*iter, mc, ctrl);
      for(j = 0; j < mc->clust.nStars; j++) {
        if(mc->clust.stars[j].isFieldStar) mc->isFieldStarNow[j]++;
        else mc->isClusterStarNow[j]++;
      }

      i = *iter - startIter;
      
      // store values
      storeBetaMass(i, masses, FeHs, mc, ctrl, FEH);
      if(mc->clust.priorVar[ABS] > EPSILON) {
        if(i % ctrl->increment == 0) {
          absorptions[i/ctrl->increment] =  getParameter(&mc->clust,ABS); 
        }
        tempMeanAbs += getParameter(&mc->clust,ABS);
      }

      updateChain(*iter, mc, ctrl);

      if(i % 200 == 0) {
        adjustChainStepSizes(mc, ctrl);
      } 
    }
    //////////////////////////////
    // End main loop for runFeH //
    //////////////////////////////
    
    // Calculate beta[FEH][0] -- the correlation between FeH and mass
    // And betaFabs -- the correlation between FeH and absorption
    if(mc->clust.priorVar[ABS] > EPSILON) {
      double absTemp = getParameter(&mc->clust,ABS);
      mc->clust.betaFabs = leastSquaresBeta(FeHs, absorptions, N_LSQ_VALUES);
      setParameter(&mc->clust,ABS,absTemp);

      mc->acceptClust[ABS] = 0;
      mc->rejectClust[ABS] = 0;
      mc->clust.stepSize[ABS] *= 5;
    }

    double tempMeanFeH = mc->clust.mean[FEH];
    findBetaMass(masses, FeHs, tempMeanFeH, mc, FEH, 0.0, 1.0);
    mc->acceptClust[FEH] = 0;
    mc->rejectClust[FEH] = 0;

    evolve(&mc->clust, -1);
    
    for(j = 0; j < mc->clust.nStars; j++) {
      free(masses[j]);
    }
    free(masses);
  }
  else *iter += nIter;
} // runFeH


///////////////////////////
// Decorrelate Y and FeH //
///////////////////////////
void runY(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl) {
  // evolve(&mc->clust, -1);
  
  if(mc->clust.priorVar[YYY] > EPSILON && mc->clust.priorVar[FEH] > EPSILON) {
    double FeHs[N_LSQ_VALUES] = {0};
    double Ys[N_LSQ_VALUES] = {0};
    
    /////////////////////////////////////
    // Main loop for runY
    /////////////////////////////////////  
    int i, j;
    int startIter = *iter;
    for (; *iter < startIter + nIter; (*iter)++) {

      if (ctrl->outputBurnIn) printChainValues(*iter, mc, ctrl);
      for(j = 0; j < mc->clust.nStars; j++) {
        if(mc->clust.stars[j].isFieldStar) mc->isFieldStarNow[j]++;
        else mc->isClusterStarNow[j]++;
      }

      i = *iter - startIter;
      
      // store values
      if(i % ctrl->increment == 0) {
        Ys[i/ctrl->increment] = getParameter(&mc->clust,YYY);//fullY(&mc->clust); 
        FeHs[i/ctrl->increment] = getParameter(&mc->clust,FEH);
      }

      updateChain(*iter, mc, ctrl);

      if(i % 200 == 0) {
        adjustChainStepSizes(mc, ctrl);
      } 
    }
    ////////////////////////////
    // End main loop for runY //
    ////////////////////////////
  
    double yTemp = getParameter(&mc->clust,YYY);
    mc->clust.betaFY = leastSquaresBeta(FeHs,Ys, N_LSQ_VALUES);
    setParameter(&mc->clust,YYY,yTemp);
    mc->acceptClust[YYY] = 0;
    mc->rejectClust[YYY] = 0;
  }
  else *iter += nIter;
} // runY


/*
/////////////////////////////////
// Decorrelate Age and Modulus //
/////////////////////////////////
void runAgeVsModulus(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl) {

  
  int increment = nIter/N_LSQ_VALUES;
  
  // For the moment, we only need this to work on single field WDs, not clusters
  if(mc->clust.nStars > 1) return;
  
  if(mc->clust.priorVar[MOD] > EPSILON) {
    double ages[N_LSQ_VALUES] = {0};
    double moduli[N_LSQ_VALUES] = {0};
 
    ///////////////////////////////////
    // Main loop for runAgeVsModulus //
    ///////////////////////////////////  
    int i, j;
    int startIter = *iter;
    for (; *iter < startIter + nIter; (*iter)++) {
      
      if (ctrl->outputBurnIn) printChainValues(*iter, mc, ctrl);
      for(j = 0; j < mc->clust.nStars; j++) {
        if(mc->clust.stars[j].isFieldStar) mc->isFieldStarNow[j]++;
        else mc->isClusterStarNow[j]++;
      }
      
      i = *iter - startIter;
      
      // store values
      if(i % increment == 0) {
        ages[i/increment] = getParameter(&mc->clust,AGE); 
        moduli[i/increment] = getParameter(&mc->clust,MOD);
        //if(mc->clust.parameter[MOD] < 6.0) mc->clust.parameter[MOD] += .01;
        //else if(mc->clust.parameter[MOD] > 2.0) mc->clust.parameter[MOD] -= .01;
        //printf("%d %f\n",i,mc->clust.parameter[MOD]);
      }
      
      updateChain(*iter, mc, ctrl);
      
      if(i % 200 == 0) {
        adjustChainStepSizes(mc, ctrl);
      } 
    }
    ///////////////////////////////////////
    // End main loop for runAgeVsModulus //
    ///////////////////////////////////////
    
    double ageTemp = getParameter(&mc->clust, AGE);
    mc->clust.betaAgeMod[0] = brokenBeta(moduli, ages, N_LSQ_VALUES, &mc->clust.betaAgeMod[1]);
    setParameter(&mc->clust, AGE, ageTemp);
    mc->acceptClust[AGE] = 0;
    mc->rejectClust[AGE] = 0;
    
    evolve(&mc->clust, -1);
    
  }
  else {
    // if the prior variance of MOD < EPSILON, run the sampler anyway
    runWander(iter, nIter, mc, ctrl);
  }
  

  double ageTemp = getParameter(&mc->clust, AGE);
  mc->clust.betaAgeMod[0] = 4.574000000;
  mc->clust.betaAgeMod[1] = -0.008737569;
  mc->clust.betaAgeMod[2] = 0.369287485;
  setParameter(&mc->clust, AGE, ageTemp);
  double massTemp = getMass1(&mc->clust.stars[0],&mc->clust);
  mc->clust.stars[0].beta[MOD][0] = 0;
  setMass1(&mc->clust.stars[0],&mc->clust,massTemp);

  printf("Done\n");fflush(stdout);
} // runAgeVsModulus

*/


//////////////////////////////////
// Turn field star sampling off //
//////////////////////////////////
void runFSSamplingOff(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl) {
  // don't allow field star sampling while saving values for the field star prior distributions
  int j;
  if (ctrl->fsSamplingOn) {
    ctrl->fsSamplingOn = 0;
    for (j = 0; j < mc->clust.nStars; j++) {
      if (mc->isClusterStarNow[j] / (double) (mc->isFieldStarNow[j] + mc->isClusterStarNow[j]) > 0.5) 
        mc->clust.stars[j].useDuringBurnIn = 1;
      else 
        mc->clust.stars[j].useDuringBurnIn = 0;
    }
  }  
} // runFSSamplingOff


/////////////////////////////////
// Turn field star sampling on //
/////////////////////////////////
void runFSSamplingOn(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl) {
  
  int j;
  for(j=0;j<mc->clust.nStars;j++) {
    mc->isFieldStarNow[j] = 0;
    mc->isClusterStarNow[j] = 0;
    if(!mc->clust.stars[j].useDuringBurnIn) mc->clust.stars[j].isFieldStar = 1;
  }
  ctrl->fsSamplingOn = 1;
} // runFSSamplingOn


/////////////////////////////////////////////////////////////////////
// Calculate means and variances of field star prior distributions //
/////////////////////////////////////////////////////////////////////
void runFieldStarPriorCalc(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl) {
  // Memory allocation
  double *meanMassRatioTemp, *meanUTemp, *varMassRatioTemp, *varUTemp;
  if((varUTemp = (double *) calloc(mc->clust.nStars, sizeof(double))) == NULL)
     perror("MEMORY ALLOCATION ERROR \n");
  if((varMassRatioTemp = (double *) calloc(mc->clust.nStars, sizeof(double))) == NULL)
     perror("MEMORY ALLOCATION ERROR \n");
  if((meanUTemp = (double *) calloc(mc->clust.nStars, sizeof(double))) == NULL)
     perror("MEMORY ALLOCATION ERROR \n");
  if((meanMassRatioTemp = (double *) calloc(mc->clust.nStars, sizeof(double))) == NULL)
     perror("MEMORY ALLOCATION ERROR \n");  
  
  /////////////////////////////////////////
  // Main loop for runFieldStarPriorCalc //
  /////////////////////////////////////////
  int i, j;
  int startIter = *iter;
  double m;  
  for (; *iter < startIter + nIter; (*iter)++) {

    if (ctrl->outputBurnIn) printChainValues(*iter, mc, ctrl);

    i = *iter - startIter;
  
    // store values
    if(i % ctrl->increment == 0) {
      for (j = 0; j < mc->clust.nStars; j++) {
        m = i/ctrl->increment + 1;
        meanUTemp[j] = ((m - 1) * meanUTemp[j] + mc->clust.stars[j].U)/m;
        if(m > 1) {
           varUTemp[j] = ((m - 1) * varUTemp[j] / m) + 
             (mc->clust.stars[j].U - meanUTemp[j]) * (mc->clust.stars[j].U - meanUTemp[j]) / (m - 1);
        }
        meanMassRatioTemp[j] = ((m - 1) * meanMassRatioTemp[j] + mc->clust.stars[j].massRatio) / m;
        if(m > 1) {
          varMassRatioTemp[j] = ((m - 1) * varMassRatioTemp[j] / m) + 
           (mc->clust.stars[j].massRatio - meanMassRatioTemp[j]) * (mc->clust.stars[j].massRatio - meanMassRatioTemp[j]) / (m - 1);
        }
      }
    }

    updateChain(*iter, mc, ctrl);

    if(i % 200 == 0) {
      adjustChainStepSizes(mc, ctrl);
    } 
  }
  /////////////////////////////////////////////
  // End main loop for runFieldStarPriorCalc //
  /////////////////////////////////////////////
  
  ///////////////////////////////////////////////
  // Set the mean and var of the distributions //
  ///////////////////////////////////////////////
  for(j = 0; j < mc->clust.nStars; j++) {
    mc->clust.stars[j].meanMassRatio = meanMassRatioTemp[j];
    mc->clust.stars[j].meanU = meanUTemp[j];
    mc->clust.stars[j].varU = varUTemp[j];
    mc->clust.stars[j].varMassRatio = varMassRatioTemp[j];
  }

  free(meanMassRatioTemp);
  free(meanUTemp);
  free(varMassRatioTemp);
  free(varUTemp);
}  // runFieldStarPriorCalc


////////////////////////////////////////////////////////////
// Write slopes and close files at the end of the burn-in //
////////////////////////////////////////////////////////////
void runEndBurnIn(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl) {
  int p;
  for(p=0;p<NPARAMS;p++) 
    fprintf(ctrl->wBetaOutFile,"    param[%d]     mean[%d]     step[%d]",p,p,p);

  fprintf(ctrl->wBetaOutFile," betamabs betaFabs betaFY\n");
  writeClust(ctrl->wBetaOutFile, &mc->clust); // function defined in structures.h

  fprintf(ctrl->wBetaOutFile,"        U massRatio isWD isFS CMprior CMprop ");
  for(p=0;p<NPARAMS;p++) 
    fprintf(ctrl->wBetaOutFile," beta[%d]",p);
  fprintf(ctrl->wBetaOutFile,"  betaMR[0]  betaMR[1]    meanMR     varMR     Ustep MRstep       meanU       varU   use isCSNow isFSnow\n");
  int j;
  for(j=0;j<mc->clust.nStars;j++){
    writeStar(ctrl->wBetaOutFile, &mc->clust.stars[j]); // function defined in structures.h
    fprintf(ctrl->wBetaOutFile,"%7d %7d\n",mc->isClusterStarNow[j],mc->isFieldStarNow[j]);
  }
  fclose(ctrl->wBetaOutFile);
  if(ctrl->outputBurnIn){
    fclose(ctrl->wMass1File[1]);
    fclose(ctrl->wMass2File[1]);
    fclose(ctrl->wClusterFile[1]);
    //fclose(ctrl->wDebugFile[1]);
  }
  printf("Burnin complete...\n");
} // runEndBurnIn


////////////////////////////////////////
// Update the chain after the burn-in //
////////////////////////////////////////
void runPostBurnIn(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl) {
  int i, j;
  int startIter = *iter;
  for (; *iter < startIter + nIter * ctrl->thin; (*iter)++) {
    i = *iter - startIter;
    if (i % ctrl->thin == 0) {
      printChainValues(*iter, mc, ctrl);
      for(j = 0; j < mc->clust.nStars; j++) {
        if(mc->clust.stars[j].isFieldStar) mc->isFieldStarNow[j]++;
        else mc->isClusterStarNow[j]++;
      }
    }      
    updateChain(*iter, mc, ctrl);
  }
} // runPostBurnIn


///////////////////////////////////////////////
// read control parameters from input stream //
///////////////////////////////////////////////
void initMcmcControl(struct chain *mc, struct mcmcControl *ctrl) {

  double priorSigma;
  
  ctrl->verbose = 0;
  ctrl->sampleVarScale = 0;
  mc->clust.evoModels.numFilts = 0;
    
  ////////////////////////////////////////////////////////////////////////////////////////
  /////////// Query user for number of steps, burn-in details, random seed ///////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  
  printf("How many steps after the burn-in?\n> ");
  scanf("%d",&ctrl->nIter);

  printf("How many steps for the initial wandering period (B0)?\n> ");
  scanf("%d",&ctrl->B0);
  if(ctrl->B0 <= 0) {
    printf("***Warning: B0 must be greater than 0.***\n");
    printf("[Setting B0 to default value of 1000...]\n");
    ctrl->B0 = 1000;
  }
  printf("How many steps for the final step size calculation period (B2)?\n> ");
  scanf("%d", &ctrl->B2);
  if(ctrl->B2 <= 0) {
    printf("***Warning: B2 must be greater than 0.***\n");
    printf("[Setting B2 to default value of 1000...]\n");
    ctrl->B2 = 1000;
  }
  printf("How many steps for calculating correlations (T)?\n> ");
  scanf("%d",&ctrl->T);
  if(ctrl->T % N_LSQ_VALUES || ctrl->T <= 0) {
    printf("***Warning: T must be a (positive) multiple of %d.***\n", N_LSQ_VALUES);
    // because INCREMENT (ie, how often values are saved when calculating
    // correlations) is defined in mcmc.h as T/N_LSQ_VALUES
    printf("[Setting T to default value of %d...]\n", 2*N_LSQ_VALUES);
    ctrl->T = 2 * N_LSQ_VALUES;
  }
  
  ctrl->increment = ctrl->T / N_LSQ_VALUES;

  printf("Run the entire burnin?:\n> ");
  scanf("%d",&ctrl->runBurnIn);

  printf("Output the entire burnin?:\n> ");
  scanf("%d",&ctrl->outputBurnIn);

  if (!ctrl->runBurnIn) ctrl->outputBurnIn = 0;

  printf("Save every ? iteration? (minimum is 1)\n> ");
  scanf("%d", &ctrl->thin);
  if(ctrl->thin < 1) {
    ctrl->thin = 1;
  }

  printf("Enter an integer seed:\n> ");
  scanf("%ld",&seed);
  init_genrand(seed);
  
  /////////////////////////////////////////////////////////////
  //////////// Query user for priors and sigmas ///////////////
  /////////////////////////////////////////////////////////////  
  
  printf("Enter mean and sigma of metallicity prior (minimum is 0.0):\n> ");
  scanf("%lf %lf",&mc->clust.priorMean[FEH],&priorSigma);
  mc->clust.priorVar[FEH] = priorSigma*priorSigma;
  //if(mc->clust.priorVar[FEH] < 0.0) mc->clust.priorVar[FEH] = 0.0;

  printf("Enter mean and sigma of modulus prior (minimum is 0.0):\n> ");
  scanf("%lf %lf",&mc->clust.priorMean[MOD],&priorSigma);
  mc->clust.priorVar[MOD] = priorSigma*priorSigma;
  //if(mc->clust.priorVar[MOD] < 0.0) mc->clust.priorVar[MOD] = 0.0;

  printf("Enter mean and sigma of absorption prior (minimum is 0.0):\n> ");
  scanf("%lf %lf",&mc->clust.priorMean[ABS],&priorSigma);
  mc->clust.priorVar[ABS] = priorSigma*priorSigma;
  //if(mc->clust.priorVar[ABS] < 0.0) mc->clust.priorVar[ABS] = 0.0;

  printf("Enter a starting value for the log(age) of the cluster:\n> ");
  scanf("%lf",&ctrl->initialAge);

  // This line is important even though the age prior is not gaussian because
  // the code checks in several places for priorVar[XXX] > 0
  mc->clust.priorVar[AGE] = 1.0;

  //////////////////////////////////////////////////////////////////////////////
  ////////////// open model file, choose model set, and load models ////////////
  //////////////////////////////////////////////////////////////////////////////

  chooseModels(&(mc->clust));
  loadModels(&(mc->clust));

  if(mc->clust.evoModels.mainSequenceEvol == CHABHELIUM) {
    printf("Enter mean and sigma of He prior (minimum is 0):\n> ");
    scanf("%lf %lf", &mc->clust.priorMean[YYY], &priorSigma);
    mc->clust.priorVar[YYY] = priorSigma*priorSigma;
//    if(mc->clust.priorVar[YYY] < 0.0) mc->clust.priorVar[YYY] = 0.0;
  }
  else {
    mc->clust.priorMean[YYY] = 0.0;
    mc->clust.priorVar[YYY] = 0.0;
  }

  //////////////////////////////////////////////////////////////////////////////
  /////////// open files for reading (data) and writing (mcmc) /////////////////
  //////////////////////////////////////////////////////////////////////////////

  // open the data file 
  char filename[100];
  printf("Enter file name containing color/magnitude data:\n> ");
  scanf("%s", filename);
  if((ctrl->rData = fopen(filename,"r")) == NULL) {
    printf("***Error: file %s was not found.***\n",filename);
    printf("[Exiting...]\n");
    exit(1);
  }

  printf("Enter minimum and maximum magnitude of MS to use and band:\n> ");
  scanf("%lf %lf %d", &ctrl->minMag, &ctrl->maxMag, &ctrl->iMag);
  /*
  if(ctrl->iMag < 0 || ctrl->iMag >= mc->clust.evoModels.numFilts){
    printf("***Error: %d not a valid magnitude index.\n", ctrl->iMag);
    printf("Choose a number 0<= N <= %d less than.***\n",mc->clust.evoModels.numFilts);
    printf("[Exiting...]\n");
    exit(1);
  }
*/
  // open output files 
  printf("Enter output file name:\n> ");
  scanf("%s", filename);
  openOutputFiles(&ctrl->wClusterFile[0], &ctrl->wClusterFile[1], filename, strlen(filename)-1, CLUS_FILE,  ctrl->outputBurnIn);
  openOutputFiles(&ctrl->wMass1File[0],   &ctrl->wMass1File[1],   filename, strlen(filename)-1, MASS1_FILE, ctrl->outputBurnIn);
  openOutputFiles(&ctrl->wMass2File[0],   &ctrl->wMass2File[1],   filename, strlen(filename)-1, MASS2_FILE, ctrl->outputBurnIn);
  //openOutputFiles(&ctrl->wDebugFile[0],   &ctrl->wDebugFile[1],   filename, strlen(filename)-1, DEBUG_FILE, ctrl->outputBurnIn);
  openOutputFiles(&ctrl->rBetaInFile,     &ctrl->wBetaOutFile,    filename, strlen(filename)-1, BETA_FILE,  ctrl->outputBurnIn);

  /////////////////////////////////////////////////////////////////////////
  ////////// Decide whether to skip the burnin period or not //////////////
  /////////////////////////////////////////////////////////////////////////

  // if (ctrl->runBurnIn) ctrl->iStart = 0;
  ctrl->iStart = 0;
  ctrl->fsSamplingOn = 0;
  if (!ctrl->runBurnIn) {
    if (ctrl->rBetaInFile == NULL) {
      printf("***Error:  Can't skip burnin unless you have a .beta.in file\n");
      exit(0);
    }
    else {
      // inloopback = 1;
      ctrl->fsSamplingOn = 1;
      // ctrl->iStart = TOTAL_BURNIN;
    }
  }
  
  /////////////////////////////////////////////////////////////
  ////////// Initialize filter prior mins and maxes ///////////
  /////////////////////////////////////////////////////////////

  int j;
  for(j = 0; j < FILTS; j++) {
    ctrl->filterPriorMin[j] = 1000;
    ctrl->filterPriorMax[j] = -1000;
  }
   
} // initMcmcControl


//////////////////////////////////////
// Set number of burn-in iterations //
//////////////////////////////////////
void setBurnIter(struct mcmcControl *ctrl, struct block *chainBlocks) {
  int burnIter = 0;
  int iBlock = 0;
  int moreBlocks = 1;
  
  while (moreBlocks) {
    if (chainBlocks[iBlock].type != POST_BURN_IN_BLOCK) {
      burnIter += chainBlocks[iBlock].nIter;
      iBlock++;
    }
    else moreBlocks = 0;
  }
  
  ctrl->burnIter = burnIter;

} // setBurnIter


///////////////////
// Read CMD data //
///////////////////
void readCmdData(struct chain *mc, struct mcmcControl *ctrl) {
  char line[300];
  double tempSigma;
  int filt;
  char *pch, sig[]="sig\0", comp[]="   \0";
  
  //if((mc->stars = (struct star *) calloc(ALLOC_CHUNK, sizeof(struct star))) == NULL)
  //  perror("MEMORY ALLOCATION ERROR \n");
  
  if((mc->clust.stars = (struct star *) calloc(ALLOC_CHUNK, sizeof(struct star))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");

  
  //printf("%p %p\n",mc->clust.stars, mc->stars);
  //mc->clust.stars = mc->stars;
  //printf("%p %p\n",mc->clust.stars, mc->stars);
  
  //printf("%d %d\n",mc->clust.stars[0].status[0],mc->stars[0].status[0]);
  //exit(1);
  
  //Parse the header of the file to determine which filters are being used
  fgets(line, 300, ctrl->rData); // Read in the header line
  
  pch = strtok (line," "); // split the string on these delimiters into "tokens"
  filt=0;
  while (pch != NULL)
  {
    pch = strtok (NULL, " "); // Ignore the first token (which is "id") and move 
                              // to the next (which should be the first filter name)
    strncpy(comp,pch,3);                  // copy the first three letters into the dummy string 'comp'
    if(strcmp(comp,sig)==0) break;        // and check to see if they are 'sig'.  If they are, there are no more filters
    
    for(filt=0;filt<FILTS;filt++){                 // Otherwise check to see what this filter's name is
      if(strcmp(pch,getFilterName(filt))==0){
        ctrl->useFilt[filt] = 1;
        mc->clust.evoModels.numFilts++;
        if (aFilt < 0) aFilt = filt;                         // Sets this to a band we know we are using (for evolve)
        break;
      }
    }    
  }
  
  //printf("nFilts = %d afilt = %d\n",mc->clust.evoModels.numFilts, aFilt);
  //for(filt=0;filt<FILTS;filt++) printf("%8s %d\n",getFilterName(filt),ctrl->useFilt[filt]);
  
  // This loop reads in photometry data
  // It also reads a best guess for the mass
  int nr, j = 0,i;
  int moreStars = 1;           // true
  void *tempAlloc;             // temporary for allocation
  while (moreStars) {
    if((j%ALLOC_CHUNK) == 0) {
      if((tempAlloc = (void *) realloc(mc->clust.stars, (j+ALLOC_CHUNK) * sizeof(struct star))) == NULL)
        perror("MEMORY ALLOCATION ERROR \n");
      else
        mc->clust.stars = (struct star *) tempAlloc;
    }
    nr = fscanf(ctrl->rData,"%*s");
    if(nr == EOF) break;
    for(i = 0; i < mc->clust.evoModels.numFilts; i++) {
      nr = fscanf(ctrl->rData, "%lf", &(mc->clust.stars[j].obsPhot[i]));
      if(mc->clust.stars[j].obsPhot[i] < ctrl->filterPriorMin[i]) ctrl->filterPriorMin[i] = mc->clust.stars[j].obsPhot[i];
      if(mc->clust.stars[j].obsPhot[i] > ctrl->filterPriorMax[i]) ctrl->filterPriorMax[i] = mc->clust.stars[j].obsPhot[i];
    }
    // copy to global variables
    for(i = 0; i < mc->clust.evoModels.numFilts; i++) {
      filterPriorMin[i] = ctrl->filterPriorMin[i];
      filterPriorMax[i] = ctrl->filterPriorMax[i];
    }
    for(i = 0; i < mc->clust.evoModels.numFilts; i++) {
      nr = fscanf(ctrl->rData, "%lf",&tempSigma);
      mc->clust.stars[j].variance[i]=tempSigma*fabs(tempSigma); 
      // The fabs() keeps the sign of the variance the same as that input by the user for sigma
      // Negative sigma (variance) is used to signal "don't count this band for this star"
    }
    nr = fscanf(ctrl->rData,"%lf %lf %d %lf %d %*d", 
          &(mc->clust.stars[j].U), &(mc->clust.stars[j].massRatio),&(mc->clust.stars[j].status[0]),
          &(mc->clust.stars[j].clustStarPriorDens),&(mc->clust.stars[j].useDuringBurnIn));
    if(mc->clust.stars[j].status[0]==WD || mc->clust.stars[j].status[0]==BD ||
       (mc->clust.stars[j].obsPhot[ctrl->iMag] >= ctrl->minMag && mc->clust.stars[j].obsPhot[ctrl->iMag] <= ctrl->maxMag)){
      j++;
    }
  }
  mc->clust.nStars = j;
  
  // copy to global values
  for (i = 0; i < FILTS; i++) useFilt[i] = ctrl->useFilt[i];

  /*
  //Temporary kludge to scatter brown dwarfs
  for(j=0;j<mc->clust.nStars;j++){
    
     for(i=0;i<mc->clust.evoModels.numFilts;i++)
       if(mc->clust.stars[j].variance[i]> 1e-9)
         mc->clust.stars[j].obsPhot[i]  = gen_norm(mc->clust.stars[j].obsPhot[i],sqrt(mc->clust.stars[j].variance[i]));
    printf("%d %f ",j, getMass1(&mc->clust.stars[j],&mc->clust));
    for(i=0;i<mc->clust.evoModels.numFilts;i++) printf("%f %f ",mc->clust.stars[j].obsPhot[i],sqrt(mc->clust.stars[j].variance[i]));
    printf("\n");
  }
  */
//exit(1);  


} // readCmdData


///////////////////////
// Memory allocation //
///////////////////////
void allocateChainMemory(struct chain *mc) {
   if((mc->acceptMass = (int *) calloc(mc->clust.nStars, sizeof(int))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
     
   if((mc->rejectMass = (int *) calloc(mc->clust.nStars, sizeof(int))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
   
   if((mc->acceptMassRatio = (int *) calloc(mc->clust.nStars, sizeof(int))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
   
   if((mc->rejectMassRatio = (int *) calloc(mc->clust.nStars, sizeof(int))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
   
   if((mc->isFieldStarNow = (int *) calloc(mc->clust.nStars, sizeof(int))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
   
   if((mc->isClusterStarNow = (int *) calloc(mc->clust.nStars, sizeof(int))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
} // allocateChainMemory


//////////////////////
// Initialize chain //
//////////////////////
void initChain(struct chain *mc, const struct mcmcControl *ctrl) {  
  int p;
  for (p = 0; p < NPARAMS; p++) {
    mc->acceptClust[p] = mc->rejectClust[p] = 0;
  }
  
  // If there is no beta in file, initialize everything to prior means
  if(ctrl->rBetaInFile == NULL){
    setParameter(&mc->clust, FEH, mc->clust.priorMean[FEH]);
    setParameter(&mc->clust, MOD, mc->clust.priorMean[MOD]);
    setParameter(&mc->clust, ABS, mc->clust.priorMean[ABS]);
    setParameter(&mc->clust, YYY, mc->clust.priorMean[YYY]);
    setParameter(&mc->clust, AGE, ctrl->initialAge);

    mc->clust.mean[AGE]         = getParameter(&mc->clust, AGE);
    mc->clust.mean[YYY]         = getParameter(&mc->clust, YYY);
    mc->clust.mean[MOD]         = getParameter(&mc->clust, MOD);
    mc->clust.mean[FEH]         = getParameter(&mc->clust, FEH);
    mc->clust.mean[ABS]         = getParameter(&mc->clust, ABS);
    mc->clust.betamabs          = 0.0;
    mc->clust.betaFabs          = 0.0;
  }
  // If there is a beta.in file, use it to initialize
  else {
    char line[300];
    fgets(line, 300, ctrl->rBetaInFile);
    readClust(ctrl->rBetaInFile, &mc->clust);
    fgets(line, 300, ctrl->rBetaInFile);
    fgets(line, 300, ctrl->rBetaInFile);
  }

  int i, j;
  for (j = 0; j < mc->clust.nStars; j++) {
    mc->clust.stars[j].meanMassRatio         = 0.0; 
    mc->clust.stars[j].isFieldStar           = 0;
    mc->clust.stars[j].clustStarProposalDens = mc->clust.stars[j].clustStarPriorDens; // Use prior prob of being clus star 
    mc->clust.stars[j].UStepSize             = 0.001;       // within factor of ~2 for most main sequence stars 
    mc->clust.stars[j].massRatioStepSize     = 0.001;
    for(i = 0 ; i < NPARAMS ; i++){ 
      mc->clust.stars[j].beta[i][0]          = 0.0;
      mc->clust.stars[j].beta[i][1]          = 0.0;
    }
    mc->clust.stars[j].betaMassRatio[0]      = 0.0;
    mc->clust.stars[j].betaMassRatio[1]      = 0.0;
    mc->clust.stars[j].meanU                 = 0.0;
    mc->clust.stars[j].varU                  = 0.0;
    for(i = 0 ; i < 2 ; i++) mc->clust.stars[j].wdType[i] = DA;
    for(i = 0; i < mc->clust.evoModels.numFilts; i++) {
      mc->clust.stars[j].photometry[i]       = 0.0;
      //mc->clust.stars[j].variance[i]     = 0.0;
    }
    if(ctrl->rBetaInFile != NULL){
      readStar(ctrl->rBetaInFile,&(mc->clust.stars[j]));
      if(ctrl->runBurnIn) mc->clust.stars[j].isFieldStar = 0;
    }
    // find photometry for initial values of currentClust and mc->clust.stars
    evolve(&mc->clust, j);
    if(mc->clust.stars[j].status[0] == WD && ctrl->rBetaInFile == NULL) {
       mc->clust.stars[j].UStepSize = 0.05; // use larger initial step size for white dwarfs
       mc->clust.stars[j].massRatio = 0.0;
    }
  }

  //j=140;
  //printf("%f ** %f\n",mc->clust.stars[j].U,mc->clust.stars[j].UStepSize);
  //mc->clust.stars[j].U = 0.0;
  // Calculate good mass starting values
  getStartMass(&mc->clust, ctrl->iMag);
    //printf("%f ** %f\n",mc->clust.stars[j].U,mc->clust.stars[j].UStepSize);

} // initChain


///////////////////////////////
// Print output file headers //
///////////////////////////////
void printOutputFileHeaders(struct chain *mc, const struct mcmcControl *ctrl) {
  const char *paramNames[] = {"    logAge",
                              "         Y",
                              "    [Fe/H]",
                              "   modulus",
                              "absorption"};
  
   fprintf(ctrl->wClusterFile[0],    "# Burnin = %d\n", ctrl->burnIter);
   fprintf(ctrl->wMass1File[0],      "# Burnin = %d numStars = %d\n", ctrl->burnIter, mc->clust.nStars);
   fprintf(ctrl->wMass2File[0],      "# Burnin = %d numStars = %d\n", ctrl->burnIter, mc->clust.nStars);
   if(ctrl->outputBurnIn){
     fprintf(ctrl->wClusterFile[1],    "# Burnin = %d\n", ctrl->burnIter);
     fprintf(ctrl->wMass1File[1],      "# Burnin = %d numStars = %d\n", ctrl->burnIter, mc->clust.nStars);
     fprintf(ctrl->wMass2File[1],      "# Burnin = %d numStars = %d\n", ctrl->burnIter, mc->clust.nStars);
   }

   // Output info so makeCMD knows what parameters are being sampled on
   fprintf(ctrl->wClusterFile[0],    "#");
   int p;
   for (p = 0; p < NPARAMS; p++){
     if(mc->clust.priorVar[p]> EPSILON) 
       fprintf(ctrl->wClusterFile[0],"%d %f ", 1, getParameter(&mc->clust,p));
     else
       fprintf(ctrl->wClusterFile[0],"%d %f ", 0, getParameter(&mc->clust,p));
   }
   fprintf(ctrl->wClusterFile[0],"\n    i ");
   for (p = 0; p < NPARAMS; p++) {
     //printf("%d %s %f\n",p,paramNames[p],mc->clust.priorVar[p]);
     if(mc->clust.priorVar[p]> EPSILON)
       fprintf(ctrl->wClusterFile[0], "%s " , paramNames[p]);
   }
   if (ctrl->sampleVarScale) fprintf(ctrl->wClusterFile[0], "varScale ");
   fprintf(ctrl->wClusterFile[0],"\n");
   
   if(ctrl->outputBurnIn){
     fprintf(ctrl->wClusterFile[1],    "#");
     for (p = 0; p < NPARAMS; p++){
       if(mc->clust.priorVar[p]> EPSILON) 
         fprintf(ctrl->wClusterFile[1],"%d %f ",1,getParameter(&mc->clust,p));
       else
         fprintf(ctrl->wClusterFile[1],"%d %f ",0,getParameter(&mc->clust,p));
     }
     fprintf(ctrl->wClusterFile[1],"\n    i ");
     for (p = 0; p < NPARAMS; p++) {
       if(mc->clust.priorVar[p]> EPSILON)
         fprintf(ctrl->wClusterFile[1],"%s ",paramNames[p]);
     }
     if (ctrl->sampleVarScale) fprintf(ctrl->wClusterFile[1], "varScale ");    
     fprintf(ctrl->wClusterFile[1],"\n");
   }
} // printOutputFileHeaders


/////////////////////////////////////////////
// Output current values of all parameters //
/////////////////////////////////////////////
void printChainValues(int iter, struct chain *mc, const struct mcmcControl *ctrl) {
  int inBurnIn = iter < ctrl->burnIter;
  //printf("%d %d\n",iter, ctrl->burnIter);fflush(stdout);
  fprintf(ctrl->wClusterFile[inBurnIn], "%5d ", iter);
  //  printf("line 2\n");fflush(stdout);
  
  int p;
  for(p=0;p<NPARAMS;p++){
    if(mc->clust.priorVar[p] > EPSILON){
    fprintf(ctrl->wClusterFile[inBurnIn], "%10.6f ", getParameter(&mc->clust, p));
    }
  }
  if (ctrl->sampleVarScale) 
    fprintf(ctrl->wClusterFile[inBurnIn], "%10.6f ", mc->clust.varScale);
  fprintf(ctrl->wClusterFile[inBurnIn], " \n");
  
  int j;
  for(j = 0; j < mc->clust.nStars; j++) {
    if(mc->clust.stars[j].isFieldStar) {
      fprintf(ctrl->wMass1File[inBurnIn], "-%8.6f ", getMass1(&mc->clust.stars[j],&mc->clust));
      fprintf(ctrl->wMass2File[inBurnIn], "-%8.6f ", getMass2(&mc->clust.stars[j],&mc->clust));
      //fprintf(ctrl->wDebugFile[inBurnIn], "%9.6f ",mc->clust.stars[j].photometry[0]);
      //fprintf(ctrl->wDebugFile[inBurnIn], "%9.6f ",mc->clust.stars[j].U);
    }
    else {
      fprintf(ctrl->wMass1File[inBurnIn], "%9.6f ", getMass1(&mc->clust.stars[j],&mc->clust));
      fprintf(ctrl->wMass2File[inBurnIn], "%9.6f ", getMass2(&mc->clust.stars[j],&mc->clust));
      //fprintf(ctrl->wDebugFile[inBurnIn], "%9.6f ",mc->clust.stars[j].photometry[0]);
      //fprintf(ctrl->wDebugFile[inBurnIn], "%9.6f ",mc->clust.stars[j].U);
    }
  }

  fprintf(ctrl->wMass1File[inBurnIn], "\n");
  fprintf(ctrl->wMass2File[inBurnIn], "\n");
  //fprintf(ctrl->wDebugFile[inBurnIn],"\n");
} // printChainValues


////////////////////////////////////////////
// Update parameter values for entire chain
////////////////////////////////////////////
void updateChain(int iter, struct chain *mc, struct mcmcControl *ctrl) {
  // Individual star sampling
  if(ctrl->fsSamplingOn) decideFieldStar(&mc->clust);//, ctrl->wDebugFile[1]);
  decideMass(mc);
  decideMassRatio(mc);

  // cluster parameter sampling
  sampleStep(iter, mc, ctrl, AGE);
  sampleStep(iter, mc, ctrl, FEH);
  sampleStep(iter, mc, ctrl, MOD);
  sampleStep(iter, mc, ctrl, ABS);
  if(mc->clust.evoModels.mainSequenceEvol == CHABHELIUM) sampleStep(iter, mc, ctrl, YYY);

  // variance scale sampling
  if(ctrl->sampleVarScale && ctrl->fsSamplingOn) updateVarScale(&mc->clust);
} // updateChain


///////////////////////////////////
// Adjust step sizes for a chain //
///////////////////////////////////
void adjustChainStepSizes(struct chain *mc, const struct mcmcControl *ctrl){ 
  mc->clust.stepSize[AGE] *= adjustStepSizes(&mc->acceptClust[AGE],&mc->rejectClust[AGE]);

  if(ctrl->currentBlock != MASS_RATIO_BLOCK) {
    if(mc->clust.priorVar[MOD] > EPSILON)
      mc->clust.stepSize[MOD] *= adjustStepSizes(&mc->acceptClust[MOD],&mc->rejectClust[MOD]);
    if(mc->clust.priorVar[FEH] > EPSILON)
      mc->clust.stepSize[FEH] *= adjustStepSizes(&mc->acceptClust[FEH],&mc->rejectClust[FEH]);
    if(mc->clust.priorVar[ABS] > EPSILON)
      mc->clust.stepSize[ABS] *= adjustStepSizes(&mc->acceptClust[ABS],&mc->rejectClust[ABS]);
    if(mc->clust.evoModels.mainSequenceEvol == CHABHELIUM)
      mc->clust.stepSize[YYY] *= adjustStepSizes(&mc->acceptClust[YYY],&mc->rejectClust[YYY]);
  }
  
  int j;
  for(j = 0; j < mc->clust.nStars; j++) {
    mc->clust.stars[j].UStepSize         *= adjustStepSizes(&mc->acceptMass[j], &mc->rejectMass[j]);
    mc->clust.stars[j].massRatioStepSize *= adjustStepSizes(&mc->acceptMassRatio[j], &mc->rejectMassRatio[j]);
    if(mc->clust.stars[j].massRatioStepSize > 1.0) mc->clust.stars[j].massRatioStepSize = 1.0;
  }
} // adjustChainStepSizes
