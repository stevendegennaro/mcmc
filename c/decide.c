/*** Decides whether to accept proposed Metropolis-Hastings steps ***/
/*** last update:    29jun06     ***/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mcmc.h"
#include "samplers.h"
#include "densities.h"
#include "decide.h"
#include "mt19937ar.h"
#include "evolve.h"
#include "structures.h"

//extern int    numFilts;

/*** Decides whether to accept a proposed jump between field star and cluster star models ***/
void decideFieldStar(struct cluster *pCluster)
{
  int j;
  double u, alpha, post1, post2;
  struct cluster clust2 = (*pCluster);
  struct star star2;
  
  clust2.stars = &star2;
  
  for(j=0;j < pCluster->nStars;j++){                   // For each star,
    clust2.stars[0] = pCluster->stars[j];              // copy stars to new array,
    propFieldStar(&clust2.stars[0]);                   // propose a new field star status,
    
    post1 = logPost1Star(pCluster,j);
    post2 = logPost1Star(&clust2,0);
    
    if(fabs(post2 + HUGE_VAL) < EPS) continue;         // Proposed star no good.  Leave pCluster->stars alone.

    // Don't bother with this if the prior and proposal densities will cancel out
    if(fabs(pCluster->stars[j].clustStarPriorDens - pCluster->stars[j].clustStarProposalDens) > EPS){
      // Multiply by prior probability of star's FS status
      if(pCluster->stars[j].isFieldStar) post1 += log(1 - pCluster->stars[j].clustStarPriorDens);
      else post1 += log(pCluster->stars[j].clustStarPriorDens);
      if(clust2.stars[0].isFieldStar) post2 += log(1 - clust2.stars[0].clustStarPriorDens);
      else post2 += log(clust2.stars[0].clustStarPriorDens);
      // Divide by proposal densities
      if(pCluster->stars[j].isFieldStar) post1 -= log(1-pCluster->stars[j].clustStarProposalDens);
      else post1 -= log(pCluster->stars[j].clustStarProposalDens);
      if(clust2.stars[0].isFieldStar) post2 -= log(1-clust2.stars[0].clustStarProposalDens);
      else post2 -= log(clust2.stars[0].clustStarProposalDens);
    }

    alpha = post2 - post1;
    u = genrand_res53();
    if(u < 1.e-15) u = 1.e-15;
    u = log(u);
  
    // Accept proposed star (copy back to the stars array in he cluster)
    if(u < alpha) pCluster->stars[j].isFieldStar = clust2.stars[0].isFieldStar; 
  }
}

/*** Decides whether to accept a proposed mass ***/
void decideMass(struct chain *mc)
{
  int j;
  double u, alpha, post1, post2;
  struct star stars2[mc->clust.nStars];
  struct star *stars1 = mc->clust.stars;
  
  for(j = 0; j < mc->clust.nStars; j++){               // For each star,
    stars2[j] = mc->clust.stars[j];                    // copy stars to new array,
    propMass(&stars2[j]);                              // propose a new mass,
    stars2[j].boundsFlag = 0;                          // and set the boundsFlag to zero
  }
  
  //printf("%p\n",mc->clust.stars);
  
  mc->clust.stars = stars2;
  evolve(&mc->clust, -1);                              // Evolve all the (proposed) stars at once
  
  for(j = 0; j < mc->clust.nStars; j++){               // Accept or reject each star individually 
    mc->clust.stars = stars2;
    if(mc->clust.stars[j].boundsFlag || getMass1(&mc->clust.stars[j], &mc->clust)<EPS) mc->rejectMass[j]++;      // Proposed star no good.  Leave stars1 alone.
    else{
      post2 = logPost1Star(&mc->clust,j);
      if(fabs(post2 + HUGE_VAL) < EPS) mc->rejectMass[j]++; // Proposed star no good.  Leave stars1 alone.
      else{
        mc->clust.stars = stars1;
        post1 = logPost1Star(&mc->clust,j);
        alpha = post2;
        //        alpha += log(getMass1(&stars2[j],pCluster)); 
        alpha -= post1;
        //        alpha -= log(getMass1(&stars1[j],pCluster));
        
        u = genrand_res53();
        if(u < 1.e-15) u = 1.e-15;
        u = log(u);
        
        if(u < alpha) {
          mc->acceptMass[j]++;                             // Accept proposed star
          mc->clust.stars[j] = stars2[j];                  // And copy back to the stars1 array
        }
        else mc->rejectMass[j]++;                          // Proposed star no good.  Leave stars1 alone.
      }
    }
  }
  mc->clust.stars = stars1;
}

/*** Decides whether to accept a proposed mass ratio ***/
void decideMassRatio(struct chain *mc)
{
  double u, alpha, post1, post2;
  int j;
  
  struct star stars2[mc->clust.nStars];
  struct star *stars1 = mc->clust.stars;

  for(j = 0; j < mc->clust.nStars; j++){               // For each star,
    stars2[j] = mc->clust.stars[j];                    // copy stars to new array,
    propMassRatio(&stars2[j]);                         // propose a new mass ratio,
    stars2[j].boundsFlag = 0;                          // and set the boundsFlag to zero
  }

  mc->clust.stars = stars2;
  evolve(&mc->clust, -1);                              // Evolve all the (proposed) stars at once
  
  for(j = 0; j < mc->clust.nStars; j++){               // Accept or reject each star individually
    mc->clust.stars = stars2;
    //if(mc->clust.stars[j].status[0] == WD || mc->clust.stars[j].status[0] != BD) mc->rejectMassRatio[j]++;
    if(stars2[j].boundsFlag || getMass1(&mc->clust.stars[j], &mc->clust) < EPS) mc->rejectMassRatio[j]++;       // Proposed star no good.  Leave stars1 alone.
    else if(mc->clust.stars[j].massRatio > 1.0 || mc->clust.stars[j].massRatio < 0.0) mc->rejectMassRatio[j]++; // Proposed star no good.  Leave stars1 alone.
    else{
      post2 = logPost1Star(&mc->clust, j);
      if(fabs(post2 + HUGE_VAL) < EPS) mc->rejectMassRatio[j]++;                            // Proposed star no good.  Leave stars1 alone.
      else{
        mc->clust.stars = stars1;
        post1 = logPost1Star(&mc->clust, j);
        alpha = post2;
        //        alpha += log(getMass1(&stars2[j],pCluster));
        alpha -= post1;
        //        alpha -= log(getMass1(&stars1[j],pCluster));
        
        u = genrand_res53();
        if(u < 1.e-15) u = 1.e-15;
        u = log(u);
        
        if(u < alpha) {
          mc->acceptMassRatio[j]++;                             // Accept proposed star
          mc->clust.stars[j] = stars2[j];                       // And copy back to the stars1 array 
        }
        else mc->rejectMassRatio[j]++;                          // Proposed star no good.  Leave stars1 alone.
      }
    }
  }
  mc->clust.stars = stars1;
}

// Decides whether to accept a proposed cluster property
struct cluster decideClust(struct chain *mc, const int FS_ON_STATE, const int SAMPLE_TYPE)
{
  int j;
   double u, alpha, post1=0.0, post2=0.0;
   
   struct star stars2[mc->clust.nStars];
   struct cluster clust2;

   clust2 = mc->clust;
   clust2.stars = stars2;
   propClustParam(&clust2,SAMPLE_TYPE);		// propose a new value

   post1 = logPriorClust(&mc->clust);
   post2 = logPriorClust(&clust2);

   if(fabs(post2 + HUGE_VAL) < EPS) {
     mc->rejectClust[SAMPLE_TYPE]++;
     return mc->clust;
   }

   for(j = 0; j < mc->clust.nStars; j++) {
     clust2.stars[j] = mc->clust.stars[j];
     if(getMass1(&mc->clust.stars[j],&clust2)<EPS){
       mc->rejectClust[SAMPLE_TYPE]++;
       return mc->clust;
     }
     clust2.stars[j].boundsFlag = 0;
   }

   evolve(&clust2, -1);

   for(j = 0; j < mc->clust.nStars; j++){
     if(clust2.stars[j].boundsFlag){
       mc->rejectClust[SAMPLE_TYPE]++;
       return mc->clust;
     }
     if(FS_ON_STATE || mc->clust.stars[j].useDuringBurnIn) {

       post1 += logPost1Star(&mc->clust, j);
       post2 += logPost1Star(&clust2, j);

       if(fabs(post2 + HUGE_VAL) < EPS) {
         mc->rejectClust[SAMPLE_TYPE]++;
         return mc->clust;
       }
     }
   }
   alpha = (post2 - post1);
   
   u = genrand_res53();
   if(u < 1.e-15) u = 1.e-15;
   u = log(u);
   
   if(u < alpha) {
     for(j = 0 ; j < mc->clust.nStars ; j++) mc->clust.stars[j] = clust2.stars[j];
     clust2.stars = mc->clust.stars;
     mc->acceptClust[SAMPLE_TYPE]++;
     return clust2;
   }
   else {
     mc->rejectClust[SAMPLE_TYPE]++;
     return mc->clust;
   }
}

// Draw a new varScale from a scaled Inv-gamma distribution.
// This assumes that the prior distribution for the varScale parameter is Inv-chisq(prior_df)
void updateVarScale(struct cluster *pCluster)
{
  int nClustStars = 0, i, j;
  double scale = 0.0;
  double prior_df = 3.0;
  for (i = 0; i < pCluster->nStars; i++) {
    if (!pCluster->stars[i].isFieldStar) {
      nClustStars++;
      for (j = 0; j < pCluster->evoModels.numFilts; j++) {
        if(pCluster->stars[i].variance[j] > 0) scale += 0.5 * ((pCluster->stars[i].photometry[j] - pCluster->stars[i].obsPhot[j]) * (pCluster->stars[i].photometry[j] - pCluster->stars[i].obsPhot[j]) / pCluster->stars[i].variance[j]);
      }
    }
  }
  scale += 0.5;
  double g = gamdev((pCluster->evoModels.numFilts * nClustStars + prior_df) / 2.0);

  pCluster->varScale = scale / g;
}
