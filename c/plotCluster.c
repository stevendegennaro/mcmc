// last update, 26jun08 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include "mt19937ar.h"
#include "evolve.h"
#include "structures.h"
#include "loadModels.h"

// Used by evolve.c
double ltau[2];
int    aFilt=0;

// Used by a bunch of different functions.
int   verbose, needMassNow=1, useFilt[FILTS];

// For random # generator (mt19937ar.c)
unsigned long mt[N], seed=0;
int           mti=N+1;

int    i, filt, nStars, cmpnt, nPops,fault;
int    filt, fault;
char   w_file[100];
FILE   *w_ptr;

double drawFromIMF(void);
double genrand_res53(void);
void   updateCount(struct star *pStar, int cmpnt);

int setGlobals(struct cluster *pCluster, int newSeed, int newVerbose)

{
  
  for(filt=0; filt<8; filt++) useFilt[filt] = 1; 	// calculate all of U-K
  for(filt=8; filt<FILTS; filt++) useFilt[filt] = 0; 	// but not the other crap

  printf("\n ***You are running BASE08/simCluster version %.1f.***\n",VERSION);

  seed = newSeed;
  printf("Seed: %ld\n",seed);
  init_genrand(seed);
  verbose = newVerbose;
  printf("Verbose: %i\n",verbose);

  if(verbose < 0 || verbose > 2) verbose = 1;		// give standard feedback if incorrectly specified 
  fault = chooseModels(pCluster);
  if(fault){
    return 1;
  }
  loadModels(pCluster);
  return 0;
}

int getUseFilt(int filt){
  return useFilt[filt];
}

void setUseFilt(int filt, int state){
  useFilt[filt] = state;
}

double getLTau(int i){
  return ltau[i];
}

