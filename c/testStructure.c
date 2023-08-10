#include <stdio.h>
#include <stdlib.h>

#define NPARAMS              5
#define FILTS                14

struct model{
    int evoModel;
    int brownDwarfEvol;
    int mainSequenceEvol;
    int IFMR;
    int WDcooling;
    int WDatm;
    int filterSet;
    int numFilts;
    int needFS;
};

struct star{
  double  obsPhot[FILTS]; 
  double  photometry[FILTS]; 
  double  variance[FILTS];
  double  useFilt[FILTS];
  double  U;
  double  massRatio;            // massRatio = secondary mass / primary mass (between 0 and 1)
  int     status[2];
  int     wdType[2];

  int     isFieldStar;
  int     useDuringBurnIn;              // switch whether to use star to burn in cluster parameters
  double  clustStarPriorDens;       // prior probability that the star is a cluster star
  double  clustStarProposalDens;    // proposal density for steps to the cluster star model

  double  beta[NPARAMS][2];
  double  betaMassRatio[2];
  double  meanU;
  double  varU;
  double  meanMassRatio;
  double  varMassRatio;

  double  UStepSize;
  double  massRatioStepSize;
  int     boundsFlag;
  double  wdLogTeff[2];
  double  massNow[2];                   // Actual current masses of each component (i.e. not zams_mass)
};

struct cluster {
  int    nStars;
  int    photometrySet;
  double M_wd_up;
  double parameter[NPARAMS];
  double stepSize[NPARAMS];
  double mean[NPARAMS];
  double priorVar[NPARAMS];
  double priorMean[NPARAMS];
  double betamabs;
  double betaFabs;
  double betaFY;
  double betaAgeMod[3];
  double AGBt_zmass;
  double varScale;
  struct model evoModels;
  struct star *stars;
};

void initStar(struct star *pStar);
void initCluster(struct cluster *pCluster);
void initModels(struct model *models);

void initStar(struct star *pStar) {

  int i;

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

void initCluster(struct cluster *pCluster){
  int p;
  pCluster->nStars = 0;
  pCluster->photometrySet = -1;
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
};



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

struct cluster * createCluster(int nStars){
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