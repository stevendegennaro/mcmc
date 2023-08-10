#if defined( WDCOOL_H )
  /* the file has been included already */
#else
#define WDCOOL_H

#define MAX_WD_MODEL        300		       	/* Biggest WD model file (to date) has ~<100 entries */

struct wdCoolingCurve {
  double mass;
  int length;
  double logRadius[MAX_WD_MODEL];
  double logAge[MAX_WD_MODEL];
  double logTeff[MAX_WD_MODEL];
};

void loadWDCool(char *path, int modelSet);
double wdMassToTeffAndRadius(double logAge, double wdPrecLogAge, double wdMass, double *thisWDLogRadius);
#endif
