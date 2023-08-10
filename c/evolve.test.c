#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "evolve.h"
#include "msRgbEvol.h"
#include "gBaraffeMag.h"

extern int    verbose, useFilt[FILTS], aFilt, needMassNow;
extern double ltau[2];

// Used by sub-methods of msRgbEvol (gGirMag, gChabMag, etc...) and wdEvol (gBergMag)
double globalMags[FILTS];
double ageLimit[2];

struct globalIso isochrone;

void evolve(struct cluster *pCluster, struct star *stars, int index)

{

  int    filt, i;

  //Allocate memory to global isochrone(if it hasn't been allocated already)
  if(isochrone.mass == NULL){
    isochrone.nEntries = 370;
    isochrone.nFilts = FILTS;
    allocateGlobalIso(&isochrone);
  }
    
  deriveAgbTipMass(pCluster);                          // determine AGBt ZAMS mass, to find evol state    
  
  for(i=0;i<isochrone.nEntries;i++){
    stars[i].U = isochrone.mass[i];
    stars[i].boundsFlag = isochrone.eep[i];
    for(filt=0;filt<isochrone.nFilts;filt++) stars[i].photometry[filt] = isochrone.mag[i][filt];
  }
  
  pCluster->nStars = isochrone.nEntries;
}