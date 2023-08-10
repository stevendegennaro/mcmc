#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "evolve.h"
#include "linInterp.h"
#include "wdCooling.h"

static int    nIso;
static double *wdMasses;
static struct wdCoolingCurve *wdCurves;

void loadWDCool(char *path, int modelSet){

  int i, j;
  FILE *pCoolingModels;
  char tempFile[100]="\0",line[240];
  double tempAge,tempTeff,tempMass,tempRadius,lastMass=0.0;
  void *tempAlloc;		// temporary for allocation

  // Allocate memory for the cooling curves dynamically so that
  // any number of curves can be read in (up to limits of memory)
  if((wdCurves = (struct wdCoolingCurve *) calloc(1, sizeof(struct wdCoolingCurve))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");
  if((wdMasses = (double *) calloc(1, sizeof(double))) == NULL)
      perror("MEMORY ALLOCATION ERROR \n");

  strcat(tempFile, path);
  if(modelSet == WOOD)
	strcat(tempFile,"xb.comb\0");
  else{
	  printf("\nCooling models do not exist.  Exiting...\n");
	  exit(1);
  }
		 
  if((pCoolingModels = fopen(tempFile,"r")) == NULL) {
    printf("\n\n file %s was not found - exiting\n",tempFile);
    exit(1);
  }

  fgets(line,240,pCoolingModels);                                       /* after header line, read in Wood model file */
  i=0;
  j=0;
  while(fgets(line,240,pCoolingModels) != NULL) {
    sscanf(line,"%*d %lf %*f %*f %*f %lf %lf %*f %*f %*f %lf",
           &tempAge,&tempRadius,&tempTeff,&tempMass);
    if(i==0 && j==0) lastMass = tempMass;
    // If the mass for this entry isn't the same as the mass
    // of the last entry, it's a new cooling curve.  Re-allocate
    // memory for the wdCurves array and start storing in the
    // next entry.
    if(tempMass != lastMass){
      wdCurves[i].length = j;
      i++;
      j=0;
      if((tempAlloc = (void *) realloc(wdCurves, (i+1) * sizeof(struct wdCoolingCurve))) == NULL)
        perror("MEMORY ALLOCATION ERROR \n");
      else
        wdCurves = (struct wdCoolingCurve *) tempAlloc;
      if((tempAlloc = (void *) realloc(wdMasses, (i+1) * sizeof(double))) == NULL)
        perror("MEMORY ALLOCATION ERROR \n");
      else
        wdMasses = (double *) tempAlloc;
    }
    wdCurves[i].logTeff[j] = tempTeff;
    wdCurves[i].mass = wdMasses[i] = tempMass;
    wdCurves[i].logRadius[j] = tempRadius;
    wdCurves[i].logAge[j] = log10(tempAge);
    lastMass = tempMass;
    j++;
  }
  nIso = i + 1;
  wdCurves[nIso-1].length = j;
}



double wdMassToTeffAndRadius(double logAge, double wdPrecLogAge, double wdMass, double *thisWDLogRadius){

  int binarySearch (double *searchArray, int size, double searchItem);

  double wdCoolLogAge = 0.0,newTeff=0.0;
  int massIndex = -1,m,ageIndex;
  double ageTeff[2]={0,0};
  double ageRadius[2]={0,0};
  
  if(logAge > wdPrecLogAge) {		// age limit check: otherwise wdCoolLogAge
    wdCoolLogAge = log10(pow(10., logAge) - pow(10., wdPrecLogAge));
  }
  else { 			      	// mcmc.c can cause this by adjusting masses and ages 
	  (*thisWDLogRadius) = 0.0;
	  return 0.0;		               	// no need to calculate anything, return to evolve.c here 
  }

  massIndex = binarySearch(wdMasses,nIso,wdMass);
  if(massIndex < 0){
    printf("Error in binary search on mass (wdCooling.c)\n");
    exit(1);
  }

  //For each mass entry, interpolate in age
  for(m=massIndex;m<massIndex+2;m++){
    ageIndex = -1;
    ageIndex = binarySearch(wdCurves[m].logAge,wdCurves[m].length,wdCoolLogAge);
    if(ageIndex<0){
      printf("Error in binary search on age (wdCooling.c)\n");
      exit(1);
    }
    
    ageTeff[m-massIndex] = linInterpExtrap(wdCurves[m].logAge[ageIndex],  wdCurves[m].logAge[ageIndex+1],
                                           wdCurves[m].logTeff[ageIndex], wdCurves[m].logTeff[ageIndex+1],
                                           wdCoolLogAge);
	
	ageRadius[m-massIndex] = linInterpExtrap(wdCurves[m].logAge[ageIndex],  wdCurves[m].logAge[ageIndex+1],
                                           wdCurves[m].logRadius[ageIndex], wdCurves[m].logRadius[ageIndex+1],
                                           wdCoolLogAge);
	
  }

  // Now interpolate in mass
  newTeff=linInterpExtrap(wdMasses[massIndex],wdMasses[massIndex+1],ageTeff[0],ageTeff[1],wdMass);
  (*thisWDLogRadius)=linInterpExtrap(wdMasses[massIndex],wdMasses[massIndex+1],ageRadius[0],ageRadius[1],wdMass);

  return newTeff;

}
