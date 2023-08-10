#if defined( GDSEDMAG_H )

#else
#define GDSEDMAG_H

#define N_DSED_Z           9                    // number of metallicities in DSED isochrones
#define N_DSED_AGES        53                   // number of ages in DSED isochonres (+1, because the 1Gyr isocrhone is repeated
#define N_DSED_FILTS       8
#define MAX_DSED_ENTRIES   370
#define LN_10              2.3025850929940459


void   loadDsed(char *path, int filterSet);
double deriveDsedAgbTip(double newFeH, double newY, double newAge);
double getDsedMags(double zamsMass);
double wdPrecLogAgeDsed(double thisFeH, double thisY, double zamsMass);
#endif
