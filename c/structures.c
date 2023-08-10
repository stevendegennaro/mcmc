#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "evolve.h"
#include "structures.h"

static double getAge(struct cluster *pCluster);
static double getY  (struct cluster *pCluster);
static double getFeH(struct cluster *pCluster);
static double getMod(struct cluster *pCluster);
static double getAbs(struct cluster *pCluster);
static void   setAge(struct cluster *pCluster, double newAge);
static void   setY  (struct cluster *pCluster, double newY);
static void   setFeH(struct cluster *pCluster, double newFeH);
static void   setMod(struct cluster *pCluster, double newMod);
static void   setAbs(struct cluster *pCluster, double newAbs);

static double (*getThisParameter[NPARAMS])(struct cluster*) =
{
	&getAge,
	&getY,
	&getFeH,
	&getMod,
  &getAbs
};

static void (*setThisParameter[NPARAMS])(struct cluster*, double) =
 {
   &setAge,
   &setY,
   &setFeH,
   &setMod,
   &setAbs
 };

void initCluster(struct cluster *pCluster){
  
  int p;
  pCluster->nStars = 0;
  pCluster->M_wd_up = 8.0;
  pCluster->betamabs = 0.0;
  pCluster->betaFabs = 0.0;
  pCluster->betaFY = 0.0;
  pCluster->AGBt_zmass = 0.0;
  for(p=0;p<3;p++) pCluster->betaAgeMod[p] = 0.0;
  for(p=0;p<NPARAMS;p++){
    pCluster->stepSize[p] = 0.0;
    pCluster->mean[p] = 0.0;
    pCluster->parameter[p] = 0.0;
    pCluster->priorVar[p] = 0.0;
    pCluster->priorMean[p] = 0.0;
  }
  pCluster->varScale = 1.0;
  initModels(&pCluster->evoModels);
}

void initModels(struct model *models){
	models->evoModel = 0;
  models->brownDwarfEvol = 0;
	models->mainSequenceEvol = 0;
	models->IFMR = 0;
	models->WDcooling = 0;
	models->WDatm = 0;
	models->filterSet = 0;
  models->numFilts = 0;
  models->needFS = 0;
  models->minMass = 0.0;
};

void initStar(struct star *pStar) {

  int i;
  pStar->id = 0;
  for(i = 0;i< FILTS; i++){
    pStar->obsPhot[i]=0.0;
    pStar->photometry[i]=0.0;
    pStar->variance[i]=0.0;
  }
  for(i = 0; i < 2; i++){
    pStar->status[i]=0;
    pStar->massNow[i]= 0.0;
    pStar->wdLogTeff[i]=0.0;
    pStar->wdType[i]=0;
    pStar->betaMassRatio[i] = 0.0;
  }
  for(i = 0; i < NPARAMS; i++){
    pStar->beta[i][0] = 0.0;
    pStar->beta[i][1] = 0.0;
  }
  pStar->U = 0.0;
  pStar->massRatio=0.0;
  pStar->isFieldStar=0;
  pStar->useDuringBurnIn=0;
  pStar->clustStarPriorDens=0.0;
  pStar->clustStarProposalDens=0.0;
  pStar->meanMassRatio = 0.0;
  pStar->varMassRatio = 0.0;
  pStar->UStepSize = 0.0;
  pStar->massRatioStepSize = 0.0;
  pStar->meanU = 0.0;
  pStar->varU = 0.0;
  pStar->boundsFlag = 0;
}

double getMass1(struct star *pStar, struct cluster *pCluster){

  double mass;

  mass = pStar->U + pStar->beta[AGE][0]*(getParameter(pCluster,AGE) - pCluster->mean[AGE]) + 
    pStar->beta[MOD][0]*(getParameter(pCluster,MOD) - pCluster->mean[MOD]) + 
    pStar->beta[FEH][0]*(getParameter(pCluster,FEH) - pCluster->mean[FEH]) + 
    pStar->beta[YYY][0]*(getParameter(pCluster,YYY) - pCluster->mean[YYY]) + 
    pStar->betaMassRatio[0]*pow(pStar->massRatio,pStar->betaMassRatio[1]);

  return mass;

}

double getMass2(struct star *pStar, struct cluster *pCluster){

  double mass2;
  mass2 = getMass1(pStar,pCluster)*pStar->massRatio;
  return mass2;

}

void setMass1(struct star *pStar, struct cluster *pCluster, double newMass){
    
  pStar->U = newMass - (pStar->beta[AGE][0]*(getParameter(pCluster,AGE) - pCluster->mean[AGE]) + 
  pStar->beta[MOD][0]*(getParameter(pCluster,MOD) - pCluster->mean[MOD]) + 
  pStar->beta[FEH][0]*(getParameter(pCluster,FEH) - pCluster->mean[FEH]) + 
  pStar->beta[YYY][0]*(getParameter(pCluster,YYY) - pCluster->mean[YYY]) + 
  pStar->betaMassRatio[0]*pow(pStar->massRatio,pStar->betaMassRatio[1]));
  
}

void setMass2(struct star *pStar,struct cluster *pCluster, double newMass){
  pStar->massRatio = newMass/getMass1(pStar, pCluster);
  
}


void writeStar(FILE *pFile, struct star *pStar){

  int p;

  fprintf(pFile,"%9.6f   %7.5f    %d    %d  %6.3f %6.3f  ",
          pStar->U,pStar->massRatio, pStar->status[0], pStar->isFieldStar,
          pStar->clustStarPriorDens,pStar->clustStarProposalDens);
  for(p=0;p<NPARAMS;p++) fprintf(pFile,"%7.3f ",pStar->beta[p][0]);
  fprintf(pFile," %6.3f %6.3f %6.3e %9.3e ",pStar->betaMassRatio[0], pStar->betaMassRatio[1], pStar->meanMassRatio,pStar->varMassRatio);
  fprintf(pFile,"%9.3e %6.3f   %9.6f  %6.3e     %d",
          pStar->UStepSize, pStar->massRatioStepSize, pStar->meanU, pStar->varU, pStar->useDuringBurnIn);

}

void readStar(FILE *pFile, struct star *pStar){
  
  int p;

  fscanf(pFile,"%lf %lf %d %d %lf %lf",
         &(pStar->U),&(pStar->massRatio), &(pStar->status[0]), &(pStar->isFieldStar),
         &(pStar->clustStarPriorDens),&(pStar->clustStarProposalDens));
  for(p=0;p<NPARAMS;p++) fscanf(pFile,"%lf ",&(pStar->beta[p][0]));
  fscanf(pFile," %lf %lf %lf %lf ",
         &(pStar->betaMassRatio[0]), &(pStar->betaMassRatio[1]), &(pStar->meanMassRatio),&(pStar->varMassRatio));
  fscanf(pFile,"%lf %lf %lf %lf %d %*d %*d",
         &(pStar->UStepSize), &(pStar->massRatioStepSize),
         &(pStar->meanU),&(pStar->varU), &(pStar->useDuringBurnIn));
}


void readClust(FILE *pFile, struct cluster *pCluster){

  int p;

  for(p=0;p<NPARAMS;p++) fscanf(pFile,"%lg %lg %lg", &(pCluster->parameter[p]), &(pCluster->mean[p]), &(pCluster->stepSize[p]));
  fscanf(pFile,"%lf %lf %lf",&(pCluster->betamabs),&(pCluster->betaFabs),&(pCluster->betaFY));

}



void writeClust(FILE *pFile, struct cluster *pCluster){

  int p;

  for(p=0;p<NPARAMS;p++) fprintf(pFile," %11.3e %11.3e %11.3e",pCluster->parameter[p],pCluster->mean[p],pCluster->stepSize[p]);
  fprintf(pFile," %8.3f %8.3f %8.3f\n",pCluster->betamabs,pCluster->betaFabs,pCluster->betaFY);

}
/*
// Copies the relevant portions of pStarFrom into pStarTo 
// after a call to evolve and an accept by mcmc
void quickCopy(struct star *pStarFrom, struct star *pStarTo){

  int i = 0;

  for(i=0;i<FILTS;i++) pStarTo->photometry[i] =  pStarFrom->photometry[i];
  for(i=0;i<2;i++) pStarTo->status[i] = pStarFrom->status[i];

}
*/

// Creates a duplicate of a cluster structure, including allocating memory for the stars array
struct cluster cloneCluster(struct cluster *pCluster){
  
  struct cluster newCluster = (*pCluster);
  struct star *newStars;
  
  if((newStars = (struct star *) malloc(newCluster.nStars * sizeof(struct star))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");
  
  newCluster.stars = newStars;
  return newCluster;
  
}

void freeCluster(struct cluster *pCluster){
  free(pCluster->stars);
  free(pCluster);
}

double getParameter(struct cluster *pCluster, int TYPE)
{
  return getThisParameter[TYPE](pCluster);
}

void setParameter(struct cluster *pCluster, int TYPE, double newValue)
{
  return setThisParameter[TYPE](pCluster,newValue);
}

char filterNames[FILTS][10];

void setFilterNames(int filterSet){
	
	int f;
	for(f=0;f<FILTS;f++) strcpy(filterNames[f],"\0");
	if(filterSet == UBVRIJHK){
		strcat(filterNames[0],"U\0");
		strcat(filterNames[1],"B\0");
		strcat(filterNames[2],"V\0");
		strcat(filterNames[3],"R\0");
		strcat(filterNames[4],"I\0");
		strcat(filterNames[5],"J\0");
		strcat(filterNames[6],"H\0");
		strcat(filterNames[7],"K\0");
	}
	else if(filterSet == ACS){
		strcat(filterNames[0],"F435W\0");
		strcat(filterNames[1],"F475W\0");
		strcat(filterNames[2],"F550M\0");
		strcat(filterNames[3],"F555W\0");
		strcat(filterNames[4],"F606W\0");
		strcat(filterNames[5],"F625W\0");
		strcat(filterNames[6],"F775W\0");
		strcat(filterNames[7],"F814W\0");
	}
	else if(filterSet == SDSS){
		strcat(filterNames[0],"u\0");
		strcat(filterNames[1],"g\0");
		strcat(filterNames[2],"r\0");
		strcat(filterNames[3],"i\0");
		strcat(filterNames[4],"z\0");
		strcat(filterNames[5],"J\0");
		strcat(filterNames[6],"H\0");
		strcat(filterNames[7],"K\0");
	}
	else{
		printf("\nfilterSet %d not available.  Exiting.\n",filterSet);
		exit(1);
	}
  
  strcat(filterNames[8], "IrB\0");
  strcat(filterNames[9], "IrR\0");
  strcat(filterNames[10],"SpB1\0");
  strcat(filterNames[11],"SpB2\0");
  strcat(filterNames[12],"SpB3\0");
  strcat(filterNames[13],"SpB4\0");
}

//setFilterNames or findFilterSet must be called first
char *getFilterName(int index){
  return filterNames[index];
}

// This assumes that the filter names have no overlap between sets.
// If at some point in the future, that turns out to no longer be true,
// this will have to be modified accordingly --SD
int findFilterSet(char *aFilter){

  int i,filt;

  for(i=0;i<3;i++){
    setFilterNames(i);
    for(filt=0;filt<8;filt++){
      if(strcmp(aFilter,getFilterName(filt))) return i;
    }
  }
  
  printf("Filterset not found.  Exiting.\n");
  exit(1);
}


static double getAge(struct cluster *pCluster){
  /*
  if(getParameter(pCluster,MOD) < pCluster->betaAgeMod[0])
    return pCluster->betaAgeMod[1]*(getParameter(pCluster,MOD)-pCluster->betaAgeMod[0]) + pCluster->parameter[AGE];
  else
    return pCluster->betaAgeMod[2]*(getParameter(pCluster,MOD)-pCluster->betaAgeMod[0]) + pCluster->parameter[AGE];
   */
  return pCluster->parameter[AGE];
  
}

static double getY(struct cluster *pCluster){
  double Y;
  Y =  pCluster->parameter[YYY] + pCluster->betaFY*(getParameter(pCluster,FEH)-pCluster->mean[FEH]);
  return Y;
}

static double getFeH(struct cluster *pCluster){
  return pCluster->parameter[FEH];
}

static double getMod(struct cluster *pCluster){
  return pCluster->parameter[MOD];
}

static double getAbs(struct cluster *pCluster){
  double abs;
  abs =  pCluster->parameter[ABS] 
    + pCluster->betaFabs*(getParameter(pCluster,FEH)-pCluster->mean[FEH])
    + pCluster->betamabs*(getParameter(pCluster,MOD)-pCluster->mean[MOD]);
  return abs;
}

static void   setAge(struct cluster *pCluster, double newAge){
  pCluster->parameter[AGE] = newAge;

  //if(getParameter(pCluster,MOD) < pCluster->betaAgeMod[0])
  //  pCluster->parameter[AGE] = newAge - pCluster->betaAgeMod[1]*(getParameter(pCluster,MOD)-pCluster->betaAgeMod[0]);
  //else
  //  pCluster->parameter[AGE] = newAge - pCluster->betaAgeMod[2]*(getParameter(pCluster,MOD)-pCluster->betaAgeMod[0]);
}
static void   setY(struct cluster *pCluster, double newY){
  pCluster->parameter[YYY] = newY - pCluster->betaFY * (getParameter(pCluster,FEH) - pCluster->mean[FEH]);
}
static void   setFeH(struct cluster *pCluster, double newFeH){
  pCluster->parameter[FEH] = newFeH;
}
static void   setMod(struct cluster *pCluster, double newMod){
  pCluster->parameter[MOD] = newMod;
}
static void   setAbs(struct cluster *pCluster, double newAbs){
  pCluster->parameter[ABS] = newAbs 
    - pCluster->betamabs * (getParameter(pCluster,MOD) - pCluster->mean[MOD])
    - pCluster->betaFabs * (getParameter(pCluster,FEH) - pCluster->mean[FEH]);
}


//NOTE: nEntries and nFilts need to be set in the calling program before this will work properly.
void allocateGlobalIso(struct globalIso *newIso) {

  if(newIso->nEntries == 0 || newIso->nFilts == 0){
    printf("Cannot allocate memory for global isochrone with %d entries and/or %d filts\n",newIso->nEntries, newIso->nFilts);
    exit(1);
  }
  
  if(newIso->mass != NULL){
    printf("globalIso memory already allocated\n");
    return;
  }
  
  if((newIso->mass = (double *) calloc(newIso->nEntries, sizeof(double))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");
  
  if((newIso->massNow = (double *) calloc(newIso->nEntries, sizeof(double))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");

  if((newIso->eep = (int *) calloc(newIso->nEntries, sizeof(int))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");

  int j=0;
  if((newIso->mag = (double **) calloc(newIso->nEntries, sizeof(double*))) == NULL)
    perror("MEMORY ALLOCATION ERROR \n");
  for(j = 0; j < newIso->nEntries; j++) {
    if((newIso->mag[j] = (double *) calloc(newIso->nFilts, sizeof(double))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
  }
}

// Swaps two mass entries in a global isochrone (n and n-1)
void swapGlobalEntries(struct globalIso *thisIso, int n, int useFilt[FILTS]){
  
  int filt, tempEep;
  double tempMass, tempMag[thisIso->nFilts];
  
  tempMass = thisIso->mass[n];
  tempEep  = thisIso->eep[n];
  for(filt=0;filt<thisIso->nFilts;filt++) if(useFilt[filt]) tempMag[filt] = thisIso->mag[n][filt];
  
  thisIso->mass[n] = thisIso->mass[n-1];
  thisIso->eep[n]  = thisIso->eep[n-1];
  for(filt=0;filt<thisIso->nFilts;filt++) if(useFilt[filt]) thisIso->mag[n][filt] = thisIso->mag[n-1][filt];
  
  thisIso->mass[n-1] = tempMass;
  thisIso->eep[n-1]  = tempEep;
  for(filt=0;filt<thisIso->nFilts;filt++) if(useFilt[filt]) thisIso->mag[n-1][filt] = tempMag[filt];
  
}

void outputGlobalIso(struct globalIso *thisIso, FILE *wPtr){
  int m,filt;
  fprintf(wPtr,"#age = %.4f logAge = %.4f FeH = %.4f z = %.4f Y = %.4f AGBt = %.5f nEntries = %d nFilts = %d\n",
          thisIso->age,thisIso->logAge,
          thisIso->FeH, thisIso->z, thisIso->Y, 
          thisIso->AgbTurnoffMass,
          thisIso->nEntries, thisIso->nFilts);

  //Add filternames and fix filts output
  fprintf(wPtr,"eep mass massNow ");
  for(filt=0;filt<thisIso->nFilts;filt++) fprintf(wPtr,"%s ",getFilterName(filt));
  fprintf(wPtr,"\n");
  
  for(m=0;m<thisIso->nEntries;m++){
    fprintf(wPtr,"%d %f %f ",thisIso->eep[m],thisIso->mass[m], thisIso->massNow[m]);
    for(filt=0;filt<thisIso->nFilts;filt++) fprintf(wPtr,"%f ",thisIso->mag[m][filt]);
    fprintf(wPtr,"\n");
  }
  fflush(wPtr);
}


//Functions for interfacing with Python
int getNPARAMS(){
    return NPARAMS;
}

int getFILTS(){
    return FILTS;
}


struct star * getStarPtr(struct cluster *pCluster, int i)
{
    if(i < pCluster->nStars) return &pCluster->stars[i];
    else return NULL;
}

struct cluster * allocCluster(int nStars){
    struct cluster *pCluster = (struct cluster *) malloc(sizeof(struct cluster));
    initCluster(pCluster);
    pCluster->nStars = nStars;

    if((pCluster->stars = (struct star *) calloc(pCluster->nStars, sizeof(struct star))) == NULL)
        perror("MEMORY ALLOCATION ERROR \n");
    for(int i = 0; i < pCluster->nStars; i ++){ 
        initStar(&pCluster->stars[i]);
    }

    return pCluster;
}

void reallocStarsArray(struct cluster *pCluster){
  /// Does not copy whatever stars are in the old array. Should only be called when
  /// generating an entirely new population of stars.
    if((pCluster->stars = (struct star *) realloc(pCluster->stars, pCluster->nStars * sizeof(struct star))) == NULL)
        perror("MEMORY ALLOCATION ERROR \n");
    for(int i = 0; i < pCluster->nStars; i ++){ 
        initStar(&pCluster->stars[i]);
    }
    return;
}
