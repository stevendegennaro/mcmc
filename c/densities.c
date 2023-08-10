/******* Density Routines *******/
/******* last update:  29jun06      *******/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mcmc.h"
#include "evolve.h"
#include "structures.h"
#include "densities.h"

#define  SQR(z) ((z)*(z))

static double logMassNorm = 0.0;
static int calcMassNorm = 0;

extern double filterPriorMin[FILTS], filterPriorMax[FILTS];
extern double ageLimit[2];

double logPriorMass(struct cluster *pCluster, int index)
// Compute log prior density
{
  const double mf_sigma = 0.67729, mf_mu = -1.02;
  double mass1, log_m1, mass2, logPrior=0.0; 
  double loglog10 = log(log(10));

  if(pCluster->stars[index].status[0] == BD) return 0.0;
  
  // Calculate the mass normalization factor once so that we don't have to
  // calculate it every time we want the mass prior 
  if(!calcMassNorm){
    double p,q,c;
    double tup,tlow;
    
    p = mf_mu + mf_sigma*mf_sigma*log(10);
    tup = (log10(pCluster->M_wd_up)-p)/(mf_sigma);
    tlow = (-1-p)/mf_sigma;
    q = exp(-(mf_mu*mf_mu-p*p)/(2*mf_sigma*mf_sigma));
    c = 1/(q*mf_sigma*sqrt(2*M_PI)*(Phi(tup)-Phi(tlow)));
     
    logMassNorm = log(c);
    calcMassNorm = 1;
  }

  mass1 = getMass1(&pCluster->stars[index],pCluster);
  mass2 = getMass2(&pCluster->stars[index],pCluster);
  
  if(mass1 > 0.1 && mass1 <= pCluster->M_wd_up) {
    if(pCluster->stars[index].isFieldStar) {
      logPrior = logTDens(pCluster->stars[index].U, pCluster->stars[index].meanU, pCluster->stars[index].varU, DOF);
      if(pCluster->stars[index].status[0] != WD) logPrior += logTDens(pCluster->stars[index].massRatio, pCluster->stars[index].meanMassRatio, pCluster->stars[index].varMassRatio, DOF);
      return logPrior;
    }
    else {
      log_m1 = log10(mass1);
      logPrior = logMassNorm + -0.5*SQR(log_m1-mf_mu)/(SQR(mf_sigma)) - log(mass1)-loglog10;
      return logPrior;
    }
  }
  else return -HUGE_VAL;
}

// Compute log prior density for cluster properties
double logPriorClust(struct cluster *pCluster)
{   
  if(getParameter(pCluster,AGE) < ageLimit[0]) return -HUGE_VAL;		 // these are possible, we just don't have models for them YET
  else if (getParameter(pCluster,AGE) > ageLimit[1]) return -HUGE_VAL;	 // appropriate for the MS/RGB models but not the WDs
  double prior = 0.0;
  if (getParameter(pCluster,ABS) < 0.0) return -HUGE_VAL;
  if (pCluster->priorVar[FEH] > EPSILON) prior += (-0.5) * SQR(getParameter(pCluster,FEH) - pCluster->priorMean[FEH]) / pCluster->priorVar[FEH];
  if (pCluster->priorVar[MOD] > EPSILON) prior += (-0.5) * SQR(getParameter(pCluster,MOD) - pCluster->priorMean[MOD]) / pCluster->priorVar[MOD];
  if (pCluster->priorVar[ABS] > EPSILON) prior += (-0.5) * SQR(getParameter(pCluster,ABS) - pCluster->priorMean[ABS]) / pCluster->priorVar[ABS];
  if (pCluster->priorVar[YYY] > EPSILON) prior += (-0.5) * SQR(getParameter(pCluster,ABS) - pCluster->priorMean[YYY]) / pCluster->priorVar[YYY];

  return prior;
}

double logLikelihood(struct star *pStar, int numFilts)
// Computes log likelihood 
{
  int    i;
  double likelihood = 0.0;
  for(i = 0; i < numFilts; i++) {
    if(pStar->isFieldStar) {
      if(filterPriorMin[i] <= pStar->obsPhot[i] && pStar->obsPhot[i] <= filterPriorMax[i])
        likelihood -= log(filterPriorMax[i] - filterPriorMin[i]);
      else {
        likelihood = -HUGE_VAL;
        return likelihood;
      }
    }
    else{
      if(pStar->variance[i] > 1e-9)
        likelihood -= 0.5 * (log(2 * M_PI * pStar->variance[i]) 
                             + (SQR(pStar->photometry[i] - pStar->obsPhot[i]) / pStar->variance[i]));
    }
  }
  return likelihood;
}

double scaledLogLike(struct star *pStar, int numFilts, double varScale)
// Computes log likelihood
{
  int    i;
  double likelihood = 0.0;
  for(i = 0; i < numFilts; i++) {
    if(pStar->isFieldStar) {
      if(filterPriorMin[i] <= pStar->obsPhot[i] && pStar->obsPhot[i] <= filterPriorMax[i])
        likelihood -= log(filterPriorMax[i] - filterPriorMin[i]);
      else {
        likelihood = -HUGE_VAL;
        return likelihood;
      }
    }
    else{
      if(pStar->variance[i] > 1e-9){
        likelihood -= 0.5 * (log(2 * M_PI * varScale * pStar->variance[i]) 
                              + (SQR(pStar->photometry[i] - pStar->obsPhot[i]) / (varScale * pStar->variance[i])));
        //printf("%f %d %f %f %f\n",pStar->U,i,pStar->photometry[i],pStar->obsPhot[i],likelihood);
      }

    }
  }
  return likelihood;
}


double logPost1Star(struct cluster *pCluster, int index)
// Compute posterior density for 1 star:
{
  double likelihood = 0.0, logPrior = 0.0;
   
  logPrior = logPriorMass(pCluster, index);
  if(fabs(logPrior + HUGE_VAL) < EPSILON) return(logPrior);
  likelihood = scaledLogLike(&pCluster->stars[index], pCluster->evoModels.numFilts, pCluster->varScale);
  if(fabs(likelihood + HUGE_VAL) < EPSILON) return(likelihood);
  return (logPrior + likelihood);
}


// computes normal distribution Phi(x) (integral from -Inf to x of normal density)
// taken from: http://www.jstatsoft.org/v11/i04/v11i04.pdf 
double Phi(double x)
{
  long double s = x, t = 0, b = x, q = x * x, i = 1;
  while(s != t) s = (t = s) + (b *= q / (i += 2));
  return 0.5 + s * exp(-0.5 * q - 0.91893853320467274178L);
}

// do not call this routine with nu = 2, which wouldn't make much sense anyway
double logTDens(double x, double mean, double var, double nu)
{
  double logp = 0;
  double s;

  //if(fabs(nu-6.0)<EPSILON){
    s = sqrt(DOF / (var * (DOF - 2)));
    logp = log(s) + GAMMA6 - 3.5 * log(1 + pow(s * (x - mean), 2) / DOF);
    //}
    /*
  else{
    s = sqrt(nu/(var*(nu-2)));
    logp = log(s) + gammln((nu+1)/2) - gammln(nu/2) - 0.5*log(nu*M_PI) - ((nu+1)/2)*log(1+pow(s*(x-mean),2)/nu);
  }
    */
  return logp;
}
