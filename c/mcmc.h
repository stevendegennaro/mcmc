/*** Last update: 19jun06 ***/

#ifdef MCMC_H
  /* the file has been included already */
#else
#define MCMC_H

#define CLUS_FILE   10
#define MASS1_FILE  11
#define MASS2_FILE  12
#define DEBUG_FILE  14
#define BETA_FILE   15

#define EPSILON		1.e-15

#define N_LSQ_VALUES  1000
#define	INCREMENT     (T/N_LSQ_VALUES)

#define ALLOC_CHUNK   5

#include "evolve.h"
#include "structures.h"

void openOutputFiles(FILE **file_ptr, FILE **file_bi_ptr, char *filename, int length, int WHICH_FILE, int outputBurnIn);
void initCurrentStars(int j, int numFilts);
void initStepSizes(struct cluster *pClust);
double adjustStepSizes(int *accept, int *reject);
void storeBetaMass(int index, double **masses, double *params, struct chain *mc, struct mcmcControl *ctrl, int TYPE);
void findBetaMass(double **masses, double *params, double tempMean, struct chain *mc, 
                  int TYPE, double lower, double upper);
void getStartMass(struct cluster *pCluster, int index);
void initMcmcControl(struct chain *mc, struct mcmcControl *ctrl);
void setBurnIter(struct mcmcControl *ctrl, struct block *chainBlocks);
void readCmdData(struct chain *mc, struct mcmcControl *ctrl);
void allocateChainMemory(struct chain *mc);
void initChain(struct chain *mc, const struct mcmcControl *ctrl);
void printOutputFileHeaders(struct chain *mc, const struct mcmcControl *ctrl);
void printChainValues(int iter, struct chain *mc, const struct mcmcControl *ctrl);
void adjustChainStepSizes(struct chain *mc, const struct mcmcControl *ctrl);
void sampleStep(int iter, struct chain *mc, struct mcmcControl *ctrl, int SAMPLE_TYPE);
void updateChain(int iter, struct chain *mc, struct mcmcControl *ctrl);
void closePostBurnInFiles(struct mcmcControl *ctrl);
void freeChainArrays(struct chain *mc);

// updating functions for different blocks of the chain
void runWander(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);
void runMassVsMassRatio(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);
void runMassVsAge(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);
//void runAgeVsModulus(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);
void runModulus(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);
void runFeH(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);
void runY(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);
void runFSSamplingOff(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);
void runFieldStarPriorCalc(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);
void runFSSamplingOn(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);
void runEndBurnIn(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);
void runPostBurnIn(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);

#endif
