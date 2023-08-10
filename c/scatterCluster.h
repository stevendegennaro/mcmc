#if defined( SCATTER_H )
  // the file has been included already 
#else

double signalToNoise(double mag, double exptime, int filt);
int readLine(FILE* filePtr);
int magCutoff(int firstFilt, double brightLimit, double faintLimit);
int stageCutoff(int isFS);
int scatterPhot(double limitSigToNoise);
double gen_norm(double mean, double std_dev);  // Returns a normal rv
int outputScatter(FILE* w_ptr, int isFS, double clusterMemberPrior, int popID);
#endif