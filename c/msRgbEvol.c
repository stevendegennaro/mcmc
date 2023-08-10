#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "msRgbEvol.h"

/****************************************************************************************
last update: 20jul10

Pulled all MS model related functions into one file to facilitate adding MSRGB models.

To do so:
(1) import the new header file in msRgbEvol.h
(2) all MS models must delcare four methods:
	1. loadXXX(), which takes a char array specifying the path where the models are stored and an int indicating
		what filter set to use and returns void
	2. deriveAgbTipMassXXX(), which takes three double arguments (FeH, Y, logAge) and returns a double (AGB tip zams mass)
	3. getXXXMags(), which takes one double argmuent (zams mass) and returns a double (mass now--most models don't use this 
                                                                                     and just return the zams mass)
	4. wdPrecLogAgeXXX(), which takes three double arguments (FeH, Y, zamsMass) and returns a double (precursor age)
(3) Add each method to the end of the appropriate function array defined below
(4) Add new model set indicator to those defined in evolve.h (has to be the next number for the function arrays to 
															  work right)

****************************************************************************************/

//Static array of pointers to the individual functions
static void (*loadMSRgbModelFunctions[])(char*, int) =
	{
    &loadGirardi,
    &loadChaboyer,
    &loadYale,
    &loadDsed
  };

static double (*deriveAgbTipMassFunctions[])(double, double, double) = 
	{
    &deriveGirAgbTip,
    &deriveChabAgbTip,
    &deriveYYAgbTip,
    &deriveDsedAgbTip
  };

static double (*getMagsFunctions[])(double) = 
	{
    &getGirardiMags,
    &getChaboyerMags,
    &getYaleMags,
    &getDsedMags
  };

static double (*wdPrecLogAgeFunctions[])(double, double, double) = 
	{
    &wdPrecLogAgeGir,
    &wdPrecLogAgeChaboyer,
    &wdPrecLogAgeYY,
    &wdPrecLogAgeDsed
  };

double msRgbEvol(struct cluster *pCluster, double zamsMass)

/****************************************************************************************
last update: 20jul10

Perform interpolation via  calls to getGirardiMags() or similar.   
The former does 3-D interpolation of the Girardi isochrones.

deriveAgbTipMass() needs to be called first
****************************************************************************************/

{
  double massNow=0.0;
  massNow = (*getMagsFunctions[(pCluster->evoModels).mainSequenceEvol])(zamsMass);
  return massNow;
}


void deriveAgbTipMass(struct cluster *pCluster)

/****************************************************************************************
last update: 20jul07

Derive AGBt mass (actually the ZAMS mass for the appropriate AGBt star) for a given 
metallicity, age, and helium abundance (if necessary).  Uses an array of pointers to functions.
Functions must return a double and have three double arguments (newFeH, newY, newLogAge)

Array indices are defined in evolve.h
****************************************************************************************/

{	
		//printf("%f %f %f\n",getParameter(pCluster,FEH),  getParameter(pCluster,YYY), getParameter(pCluster,AGE));
		pCluster->AGBt_zmass = (*deriveAgbTipMassFunctions[(pCluster->evoModels).mainSequenceEvol])(getParameter(pCluster,FEH),  getParameter(pCluster,YYY), getParameter(pCluster,AGE));
  	return;
}



double wdPrecLogAge(struct cluster *pCluster, double zamsMass)

/****************************************************************************************
last update: 20jul10

Derive WD precursor age for a given metallicity, calling in turn wd_prec_g_lage to 
interpolate among Girardi or wd_prec_c_lage to interpolate among Chaboyer isochrones in 
mass and age.

Distributed most of the code to the respective subroutines, leaving only those to be 
modified for different model sets.
****************************************************************************************/

{
	double wdPrecLogAge=0.0;
	wdPrecLogAge = (*wdPrecLogAgeFunctions[(pCluster->evoModels).mainSequenceEvol])(getParameter(pCluster,FEH),  getParameter(pCluster,YYY), zamsMass);
	return wdPrecLogAge;
}


void loadMSRgbModels(struct cluster *pCluster, char *path)
{

		// The Yale models have issues creating field stars in simCluster, so if you need field stars 
	// and you are using the Yale models, also load the DSED models to create the field stars
	//if(pCluster->evoModels.mainSequenceEvol == YALE && pCluster->evoModels.needFS) 
	(*loadMSRgbModelFunctions[DSED])(path, pCluster->evoModels.filterSet);
	
	//printf("******************************\n");
	//fflush(stdout);
	//printf("Loading modelset %d from path %s \n",pCluster->evoModels.mainSequenceEvol,path);
	(*loadMSRgbModelFunctions[pCluster->evoModels.mainSequenceEvol])(path, pCluster->evoModels.filterSet);
}
