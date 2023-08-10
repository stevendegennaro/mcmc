/* decide.h */

#ifdef DECIDE_H
  /* the file has been included already */
#else

#define DECIDE_H

void decideFieldStar(struct cluster *pCluster);
void decideMass(struct chain *mc);
void decideMassRatio(struct chain *mc);
struct cluster decideClust(struct chain *mc, const int FS_ON_STATE, const int SAMPLE_TYPE);
void updateVarScale(struct cluster *pCluster);
         
#endif
