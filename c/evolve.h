#if defined( EVOLVE_H )
  // the file has been included already 
#else
#define EVOLVE_H

#define EPS                0.000001
#define FILTS              14		       	// UBVRIJHK, respectively 
#define MAG_LIST           26		       	// list length for mag vs. S/N 
#define MODEL_LIST         118		       	// age/Mv/V-I entries in WD model array 

#define MSRG               1                    // Main-sequence or red giant 
#define WD                 3                    // White dwarf 
#define NSBH               4                    // Neutron star or black hole 
#define BD                 5                    // Brown dwarf/planet
#define DNE                9                    // Does not exist (zams mass > 100 Msun or < 0.0001) 

// Numbers for different populations
#define MAXPOPS           10
#define FSPOP             (MAXPOPS - 1)
#define BDPOP             (MAXPOPS - 2)

// Filter sets
#define UBVRIJHK           0
#define ACS                1
#define SDSS               2

// Define model sets
//Main sequence evolution
#define GIR				         0
#define CHABHELIUM		     1
#define YALE			         2
#define DSED			         3

//IFMR
#define WEIDEMANN		       0
#define WILLIAMS		       1
#define SALARISLIN         2
#define SALARISPW          3

//WD Cooling
#define WOOD			         0

//WD Atmosphere
#define BERGERON		       0

//Brown Dwarf/Planet
#define NONE               0
#define BARAFFE            1

#define ZAMS               0
#define NOW                1

#define LOW                0
#define HIGH               1

#define DA                 0
#define DB                 1


#include "structures.h"

#define VERSION            1.0				// last update: 23may11

void evolve(struct cluster *pCluster, int index);
void returnIso(struct cluster *pCluster);
double getAgeLimit(int index);
double getFeHLimit(int index);

#endif
