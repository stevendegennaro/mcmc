#if defined( LOAD_MODELS_H )
  /* the file has been included already */
#else
#define LOAD_MODELS_H

void loadModels(struct cluster *theCluster);
int chooseModels(struct cluster *theCluster);

#endif
