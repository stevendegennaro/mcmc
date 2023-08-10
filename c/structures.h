#ifdef STRUCT_H
  /* the file has been included already */
#else
#define STRUCT_H

//Be careful adding new sample types.  Their order is important.
//There are a few places in the code that test for SAMPLE_TYPE > or <.
#define MASS                -2
#define AGE_DURING_WANDER   -1          // age sampling during AGE_WANDER (defined in samplers.h)
#define AGE	                 0          // age sampling
#define YYY		               1	        // helium sampling
#define FEH	                 2          // metallicity sampling
#define MOD	                 3          // modulus sampling
#define ABS                  4           // absorption sampling;

#define NPARAMS              5

/*** Define a structure model that houses information about the evolution model ***/
struct model{
	int evoModel;
  int brownDwarfEvol;
	int mainSequenceEvol;
	int IFMR;
	int WDcooling;
	int WDatm;
	int filterSet;
  int numFilts;
  int needFS;
  double minMass;
};

struct cluster {
  int    nStars;
  double M_wd_up;
  double parameter[NPARAMS];
  double stepSize[NPARAMS];
  double mean[NPARAMS];
  double priorVar[NPARAMS];
  double priorMean[NPARAMS];
  double betamabs;
  double betaFabs;
  double betaFY;
  double betaAgeMod[3];
  double AGBt_zmass;
  double varScale;
  struct model evoModels;
  struct star *stars;
};

// Define a structure star that houses all star properties
struct star{
  int     id;
  double  obsPhot[FILTS]; 
  double  photometry[FILTS]; 
  double  variance[FILTS];
  double  useFilt[FILTS];
  double  U;
  double  massRatio;			// massRatio = secondary mass / primary mass (between 0 and 1)
  int     status[2];
  int     wdType[2];

  int     isFieldStar;
  int     useDuringBurnIn;              // switch whether to use star to burn in cluster parameters
  double  clustStarPriorDens;		// prior probability that the star is a cluster star
  double  clustStarProposalDens;	// proposal density for steps to the cluster star model

  double  beta[NPARAMS][2];
  double  betaMassRatio[2];
  double  meanU;
  double  varU;
  double  meanMassRatio;
  double  varMassRatio;

  double  UStepSize;
  double  massRatioStepSize;
  int     boundsFlag;
  double  wdLogTeff[2];
  double  massNow[2];                   // Actual current masses of each component (i.e. not zams_mass)
};


struct globalIso {
  double FeH;
  double logAge;
  double age;
  double Y;
  double z;
  double *mass;
  double *massNow;
  double **mag;  //Dimensions are [NMASSES][NFILTS];
  double AgbTurnoffMass;
  int nEntries;
  int nFilts;
  int *eep;
};

struct chain {
  struct cluster  clust;
  double temperature;
  int acceptClust[NPARAMS];
  int rejectClust[NPARAMS];
  int *acceptMass;
  int *rejectMass;
  int *acceptMassRatio;
  int *rejectMassRatio;
  int *isFieldStarNow;
  int *isClusterStarNow;
};

enum blockType {
  INITIAL_WANDER_BLOCK,
  MASS_RATIO_BLOCK,
  AGE_BLOCK,
  MODULUS_BLOCK,
  AGE_MODULUS_BLOCK,
  FEH_BLOCK,
  Y_BLOCK,
  TURN_FS_SAMPLING_OFF_BLOCK,
  FIELD_STAR_PRIOR_BLOCK,
  TURN_FS_SAMPLING_ON_BLOCK,
  FINAL_WANDER_BLOCK,
  END_BURN_IN_BLOCK,
  POST_BURN_IN_BLOCK
};

struct mcmcControl {
  FILE *rData;
  FILE *wMass1File[2];
  FILE *wMass2File[2];
  FILE *wClusterFile[2];
  FILE *wDebugFile[2]; 
  FILE *rBetaInFile;
  FILE *wBetaOutFile;
  int fsSamplingOn;
  int sampleVarScale;
  int B0;                   // length of first "wandering" period
  int B2;                   // length of final "wandering" period
  int T;                    // length of decorrelation periods
  int increment;            // = ctrl.T / N_LSQ_VALUES
  int nIter;                // number of post burn-in iterations
  int burnIter;             // total number of burn-in iterations
  int runBurnIn;
  int outputBurnIn;
  int thin;
  //double priorMean[NPARAMS];
  //double priorVar[NPARAMS];
  double initialAge;
  double minMag;
  double maxMag;
  int iMag;
  int iStart;
  double filterPriorMin[FILTS];
  double filterPriorMax[FILTS];
  int verbose;
  int useFilt[FILTS];
  //int numFilts;
  enum blockType currentBlock;
};

struct block {
  enum blockType type;
  void (*run)(int *iter, int nIter, struct chain *mc, struct mcmcControl *ctrl);
  int nIter;
};


// Helper functions for star and cluster structures.
void initStar(struct star *pStar);
void initCluster(struct cluster *pCluster);
void initModels(struct model *models);
void readStar(FILE *pFile, struct star *pStar);
void writeStar(FILE *pFile, struct star *pStar);
void readClust(FILE *pFile, struct cluster *pCluster);
void writeClust(FILE *pFile, struct cluster *pCluster);
struct cluster cloneCluster(struct cluster *pCluster);
void freeCluster(struct cluster *pCluster);
double getMass1(struct star *pStar, struct cluster *pCluster);
double getMass2(struct star *pStar,struct cluster *pCluster);
void setMass1(struct star *pStar, struct cluster *pCluster, double newMass);
void setMass2(struct star *pStar,struct cluster *pCluster, double newMass);
//void quickCopy(struct star *pStarFrom, struct star *pStarTo);
double getParameter (struct cluster *pCluster, int TYPE);
void setParameter (struct cluster *pCluster, int TYPE, double newValue);
void setFilterNames(int filterSet);
char *getFilterName(int index);
int findFilterSet(char *aFilter);
void allocateGlobalIso(struct globalIso *newIso);
void swapGlobalEntries(struct globalIso *thisIso, int n, int useFilt[FILTS]);
void outputGlobalIso(struct globalIso *thisIso, FILE *wPtr);


//Helper functions for Python
int getNPARAMS();
int getFILTS();
struct star * getStarPtr(struct cluster *pCluster, int i);
struct cluster * allocCluster(int nStars);
void reallocStarsArray(struct cluster *pCluster);
//struct cluster * freeCluster(struct cluster *pCluster);

#endif


