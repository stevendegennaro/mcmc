/* densities.h */

#ifdef DENSITIES_H
  /* the file has been included already */
#else
#define DENSITIES_H

#define DOF                  6.0
#define GAMMA6              -2.0590305444197083417635   /* GAMMA for DOF=6 */

double logPriorMass(struct cluster *pCluster, int index);
double logPriorClust(struct cluster *p_Clust);
double logLikelihood(struct star *pStar, int numFilts);
double scaledLogLike(struct star *pStar, int numFilts, double varScale);
double logPost1Star(struct cluster *pClust, int index);
double Phi(double x);
double logTDens(double x, double mean, double var, double nu);
#endif
