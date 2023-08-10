#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "evolve.h"
#include "loadModels.h"
#include "msRgbEvol.h"
#include "gBergMag.h"
#include "wdCooling.h"
#include "gBaraffeMag.h"

// Declared in parent program (mcmc, simCluster, makeCMD)
extern int    verbose, needMassNow;

int chooseModels(struct cluster *theCluster)

// **************************************************************************************
// last update: 29may12
// 
// Query the user for model set #s.
// **************************************************************************************

{
  //printf("useFilt address loadmodels: %p\n",(void *)&useFilt);
  // printf("needMassNow address loadmodels: %p\n",(void *)&needMassNow);
  // printf("verbose address loadmodels: %p\n",(void *)&verbose);
  //char path[100] = "../models/\0";
  printf("Main Sequence Model: %i\n",theCluster->evoModels.mainSequenceEvol);
  if(theCluster->evoModels.mainSequenceEvol < 0 || theCluster->evoModels.mainSequenceEvol > 3) {
      printf("***Error: No models found for main sequence evolution model %d.***\n",theCluster->evoModels.mainSequenceEvol);
      printf("[Exiting...]\n");
      return 1;
  }
  printf("Filter Set: %i\n",theCluster->evoModels.filterSet);
  if(theCluster->evoModels.filterSet < 0 || theCluster->evoModels.filterSet > 2) {
      printf("***Error: No models found for filter set %d.***\n",theCluster->evoModels.filterSet);
      printf("[Exiting...]\n");
      return 1;
  }

  setFilterNames(theCluster->evoModels.filterSet);
  printf("IFMR: %i\n",theCluster->evoModels.IFMR);
  if(theCluster->evoModels.IFMR < 0 || theCluster->evoModels.IFMR > 3) {
      printf("***Error: No models found for IFMR %d.***\n",theCluster->evoModels.IFMR);
      printf("[Exiting...]\n");
      return 1;
  }
  printf("Brown Dwarf Model: %i\n",theCluster->evoModels.brownDwarfEvol);
  if(theCluster->evoModels.IFMR < 0 || theCluster->evoModels.IFMR > 1) {
      printf("***Error: No models found for brown dwarfs %d.***\n",theCluster->evoModels.brownDwarfEvol);
      printf("[Exiting...]\n");
      return 1;
  }

  
  theCluster->evoModels.WDcooling = WOOD;
  theCluster->evoModels.WDatm = BERGERON;
  
  return 0;
}
 
// Call appropriate function to load that module's models
void loadModels(struct cluster *theCluster){
   
  printf("\nReading models...\n");

  //setModels(theCluster, modelSet);
  char path[100] = "../models/\0";

  loadMSRgbModels(theCluster, path);
  loadWDCool(path, theCluster->evoModels.WDcooling);
  loadBergeron(path,  theCluster->evoModels.filterSet);    
  
  if(theCluster->evoModels.brownDwarfEvol == BARAFFE)
    loadBaraffe(path);

  if(theCluster->evoModels.mainSequenceEvol == YALE) theCluster->evoModels.minMass = 0.4;
  else if(theCluster->evoModels.mainSequenceEvol == DSED) theCluster->evoModels.minMass =  0.25;
  else theCluster->evoModels.minMass = 0.15;

  printf("Models read.\n");

}
