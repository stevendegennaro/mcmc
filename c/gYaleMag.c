#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "evolve.h"
#include "linInterp.h"
#include "binSearch.h"
#include "gYaleMag.h"
#include "binSearch.h"

// Fortran functions found in YYinterp.f
//extern void loadyyf_();
//extern void deriveiso_(double *FeH, double *age, double masses[],
 //                       int *numM, double yy[][N_YY_FILTS], int *flag);

extern int useFilt[FILTS];
extern double globalMags[FILTS];

// Defined in parent program (mcmc, simCluster,makeCMD)
extern int    verbose;

// Defined in evolve.c
extern double ageLimit[2];
extern double FeHLimit[2];
extern struct globalIso isochrone;
static int  loaded = 0;

//global variables
static int    iZ, iAge;
static double yyFeH[N_YY_Z],yyZ[N_YY_Z];
static double yyLogAge[N_YY_Z][N_YY_AGES],yyAge[N_YY_Z][N_YY_AGES];
static struct yyIsochrone yyIso[N_YY_Z][N_YY_AGES];
static double yyAGBt[N_YY_Z][N_YY_AGES];
static char   tempFile[100];
static struct globalIso tempIso[2];
static double coeff[6][8];

//Static funtions
static void   initIso(struct yyIsochrone *newIso);
static void   getFileName(char *path, int z);
static void   eepset(double x[], double y[], int nstep, int igrd, double *slope);
static int    ToffM(double x[4], double yy[4], double *xp, double *ybar, int iorder);
static void   convertColorsToMags(struct yyIsochrone *iso, double param[MAX_YY_ENTRIES][N_YY_PARAMS]);
static double feh2z(double FeH);
static void   intpolZ(int iZ, int iAge, double newZ);
static void   intpolAge(int iAge, double newAge);
//static void   outputYYIso(struct yyIsochrone *iso, FILE *ptr, int byMag);

void loadYale(char *path, int filterSet){

  //loadyyf_();
  //exit(1);
  
  FILE *pYY;// = NULL;
	int z, a,n,p;
	char line[500], temp[50], *lineCopy = line;
  int m=0,im=0,ii=0;
  
  int meep=0,kipnm=0;
  double eep[MAX_YY_ENTRIES],xeep[MAX_YY_ENTRIES],tlkip=0.0,blkip=0.0,slope=0.0;
  double temar[MAX_YY_ENTRIES][N_YY_PARAMS],xxeep,xim;
  double param[MAX_YY_ENTRIES][N_YY_PARAMS];
    
  if(loaded && verbose) {printf("Yale models already loaded.\n"); return;}
  else loaded = 1;

  if(filterSet != UBVRIJHK){
    printf("\nFilter set %d not available on YY models.  Exiting...\n",filterSet);
    exit(1);
  }  
  
  for(z=0 ; z < N_YY_Z ; z++) {                                    // foreach Dsed metallicity/isochrone file 
		yyFeH[z] = 0.0;

		for(a=0 ; a < N_YY_AGES ; a++) {                               // initialize age/boundary pointers 
			yyLogAge[z][a] = 0.0;
			yyAge[z][a]    = 0.0;
			initIso(&(yyIso[z][a]));                                        // initialize array of model parameters 
    }

    getFileName(path, z);							 // work on one Dsed model at a time 
    if((pYY = fopen(tempFile,"r")) == NULL) {                    // open file 
      printf("\n\n file %s was not found - exiting\n",tempFile);
      exit(1);
    }
    
    // Read header line
    fgets(line,500,pYY);
    strncpy(temp,&line[2],8);
    temp[8]='\0';
    yyZ[z] = atof(temp);
    strncpy(temp,&line[51],9);
    temp[9]='\0';
    yyFeH[z] = atof(temp);
    for(a=0 ; a < N_YY_AGES ; a++){
      yyIso[z][a].FeH = yyFeH[z];
      yyIso[z][a].z   = yyZ[z];
    }
    
    a=-1;
    while(fgets(line,500,pYY) != NULL) {  // YY model for all ages 
      if(line[0] == 'a'){
        a++;
        strncpy(temp,&line[9],6);
        temp[6]='\0';
        yyIso[z][a].age = yyAge[z][a] = atof(temp);
        yyIso[z][a].logAge = yyLogAge[z][a] = log10(yyIso[z][a].age*1e9);
        strncpy(temp,&line[16],3);
        temp[3]='\0';
        yyIso[z][a].nEntries = atoi(temp);
        m = 0;
      }
      else if(line[1] != '\n') {
        //       m,t,l,g,mv,ub,bv,vr,vi,vj,vh,vk,vl,vm,n1,n135,n3
        //  M/Msun     logT  logL/Ls   logg    Mv     U-B    B-V    V-R    V-I    V-J    V-H    V-K    V-L    V-M     #(x=-1)    #(x=1.35)    #(x=3)
        p=0;
        lineCopy = line;
        while (sscanf (lineCopy, "%s%n",temp, &n) == 1 ) {
          param[m][p] = atof(temp);
          lineCopy += n;
          p++;
        }
        
        // >>> YCK modified for the eep
        //--------------yi----------------
        if(m <= 0){
          tlkip=param[m][1];
          blkip=param[m][2];
          eep[m]=1.0;
        }
        else{
          eep[m]=eep[m-1]
          + sqrt(100*SQR(param[m][1]-tlkip)
                   + SQR(param[m][2]-blkip));  
          tlkip=param[m][1];
          blkip=param[m][2];
        }
        xeep[m]=param[m][1];
        if(yyIso[z][a].nEntries < MAX_YY_ENTRIES) kipnm = yyIso[z][a].nEntries;
        m++;
      }
      else{ //i.e., at the end of each age entry
        
        // >>> YCK modified for the eep
        slope=3.0e-2;
        // to find the turnoff mass to utilize it as an anchor point
        if(yyIso[z][a].nEntries != kipnm){
            eepset(eep,xeep,yyIso[z][a].nEntries,kipnm,&slope);
        }
        
        for(im=0;im<yyIso[z][a].nEntries;im++){
          xim=im*slope;
          xeep[im] = atan(xim) *
            (eep[yyIso[z][a].nEntries-1]-eep[0])/(atan(slope*(yyIso[z][a].nEntries-1))) +eep[0];
        }
        
        meep=0;
        for(im=0;im<yyIso[z][a].nEntries;im++){
          xxeep=xeep[im];
          meep = binarySearch(eep,yyIso[z][a].nEntries,xxeep);
          
          if(meep >= (yyIso[z][a].nEntries-2)) meep = yyIso[z][a].nEntries-2;
          if(meep <= 0) meep=0;
          for(ii=0;ii<N_YY_PARAMS;ii++){
            temar[im][ii]= POLLIN(eep[meep], param[meep][ii],
                                  eep[meep+1],param[meep+1][ii],xxeep);
          }
        }
        
        // normalized
        for(im=0;im<yyIso[z][a].nEntries;im++){
          for(ii=0;ii<N_YY_PARAMS;ii++){
            param[im][ii]=temar[im][ii];
          }
        }
        
        convertColorsToMags(&(yyIso[z][a]),param);

        yyAGBt[z][a] = yyIso[z][a].AgbTurnoffMass = yyIso[z][a].mass[yyIso[z][a].nEntries-1];
      }
    }
  }
  
  /////////////////////////////////////////////////////////////////////
  // Read in the coefficients needed to calculate the precurser ages //
  /////////////////////////////////////////////////////////////////////  
  
  // Open coeff file for reading
  strcpy(tempFile, path);
  strcat(tempFile,"YYiso/yyAGBtcoeff.dat\0");
  
  //fscanf(pModelList,"%s",tempFile);
  if((pYY = fopen(tempFile,"r")) == NULL) {
    printf("\n\n file %s was not found - exiting\n",tempFile);
    exit(1);
  }
  
  int i,j;
  // Read in the coefficients needed to calculate the precurser ages
  for(i=0;i<6;i++){
    for(j=0;j<8;j++){
      fscanf(pYY,"%lf ",&coeff[i][j]);
    }
  }
  
  //Set the min and max age in this model set (for use in densities.c)
  ageLimit[0] = yyLogAge[0][0];
  ageLimit[1] = yyLogAge[0][N_YY_AGES - 1];
  FeHLimit[0] = yyFeH[0];
  FeHLimit[1] = yyFeH[N_YY_Z - 1];
  
  for(i=0;i<2;i++){
    tempIso[i].nEntries = MAX_YY_ENTRIES;
    tempIso[i].nFilts   = N_YY_FILTS;
    allocateGlobalIso(&tempIso[i]); 
  }
   
  fclose(pYY);

  
}

/**********************************************************************************
// Last updated 11apr11--SD
// Interpolates between isochrones for two ages using linear interpolation
// Must run loadYale() first for this to work.
// Currently ignores newY
**********************************************************************************/

double deriveYYAgbTip(double newFeH, double newY, double newLogAge){
  
  double newAge=pow(10,newLogAge)/1e9;
  double newZ;
  
  newZ = feh2z(newFeH);
  
  iAge=-1;
  iZ=-1;
  
  if(newAge < yyAge[0][0]) {
    if(verbose) printf("\n Requested age (%.3f) too young. (gYaleMag.c)",newLogAge);
    return 0.0;
  }
  if(newAge > yyAge[N_YY_Z-1][N_YY_AGES-1]) {
    if(verbose) printf("\n Requested age (%.3f) too old. (gYaleMag.c)",newLogAge);
    return 0.0;
  }
  if(newZ < yyZ[0]){
    if(verbose) printf("\n Requested Z (%.3f, FeH = %.3f) too low. (gYaleMag.c)",newZ,newFeH);
    return 0.0;
  }
  if(newZ > yyZ[N_YY_Z - 1]){
    if(verbose) printf("\n Requested Z (%.3f, FeH = %.3f) too high. (gYaleMag.c)",newZ,newFeH);
    return 0.0;
  }
  
  
  // Find the values for each parameter that we will be interpolating 
  // between and calculate the interpolation coefficients.
  iZ = binarySearch(yyZ,N_YY_Z,newZ);
  iAge = binarySearch(yyAge[iZ],N_YY_AGES,newAge);

  iZ--;
  if(iZ >= N_YY_Z-4) iZ=N_YY_Z-4;
  if(iZ <= 0) iZ=0;
//    newY = dydz*(newZ-zp)+yp;
  
  intpolZ(iZ, iAge, newZ);    
  intpolAge(iAge,newAge);  
  
  isochrone.age    = newAge;
  isochrone.logAge = log10(newAge*1e9);
  isochrone.FeH    = newFeH;
  isochrone.z      = newZ;
  isochrone.AgbTurnoffMass = isochrone.mass[isochrone.nEntries - 1];

  return isochrone.AgbTurnoffMass;
}


static void intpolZ(int iZ, int iAge, double newZ){
  
  int i, p, m;

  //Will eventually interpolate in age between these two isochrones
  for(i=0;i<2;i++){
    for(m = 0 ; m < yyIso[iZ+1][i+iAge].nEntries ; m++){
      tempIso[i].mass[m] = 
      CUBEINT(log(yyZ[iZ]),
              yyIso[iZ][i+iAge].mass[m],
              log(yyZ[iZ+1]),
              yyIso[iZ+1][i+iAge].mass[m],
              log(yyZ[iZ+2]),
              yyIso[iZ+2][i+iAge].mass[m],
              log(yyZ[iZ+3]),
              yyIso[iZ+3][i+iAge].mass[m],
              log(newZ));
      for(p=0 ; p<N_YY_FILTS ; p++){
        tempIso[i].mag[m][p] = 
        CUBEINT(log(yyZ[iZ]),
                yyIso[iZ][i+iAge].mag[m][p],
                log(yyZ[iZ+1]),
                yyIso[iZ+1][i+iAge].mag[m][p],
                log(yyZ[iZ+2]),
                yyIso[iZ+2][i+iAge].mag[m][p],
                log(yyZ[iZ+3]),
                yyIso[iZ+3][i+iAge].mag[m][p],
                log(newZ));
      }
    }
    
    tempIso[i].nEntries = yyIso[iZ+1][iAge+i].nEntries;
    tempIso[i].age = yyIso[iZ][iAge+i].age;
    tempIso[i].logAge = yyIso[iZ][iAge+i].logAge;
  }
  
  return;
}


void intpolAge(int iAge, double newAge){
  
  int iYYm, m, p;
  
  //SD -- Isochrones for the lower ages have fewer entries (~35 instead of 140)
  //If this age is right on the border, just use the two higher age entries
  //that have all 140 entries (and extrapolate, technically)
  iYYm = tempIso[0].nEntries;//yyIso[iZ][iAge+1].nEntries;
  if(yyIso[iZ][iAge].nEntries < iYYm) iAge++;

  for(m = 0 ; m < iYYm ; m++){
    isochrone.mass[m] = POLLIN(tempIso[0].age,
                                tempIso[0].mass[m],
                                tempIso[1].age,
                                tempIso[1].mass[m],
                                newAge);    
    
    for(p = 0 ; p < N_YY_FILTS ; p++){
      isochrone.mag[m][p] = POLLIN(tempIso[0].age,
                                    tempIso[0].mag[m][p],
                                    tempIso[1].age,
                                    tempIso[1].mag[m][p],
                                    newAge);
    }
  }
  
  isochrone.nEntries = iYYm;  
  return;
}

static void convertColorsToMags(struct yyIsochrone *iso, double param[MAX_YY_ENTRIES][N_YY_PARAMS]){
  
  int m,p;
  
  for(m=0;m<iso->nEntries;m++){
    iso->mass[m]   = param[m][0]; // Mass
    iso->mag[m][2] = param[m][4]; // V = Mv
    iso->mag[m][1] = param[m][6] + iso->mag[m][2]; // B = (B-V) + V
    iso->mag[m][0] = param[m][5] + iso->mag[m][1]; // U = (U-B) + B
    for(p=7;p<14;p++){ // R = V - (V-R), etc
      iso->mag[m][p-4] = iso->mag[m][2] - param[m][p];
    }
  }
  iso->AgbTurnoffMass = iso->mass[iso->nEntries - 1];
}

static void eepset(double x[], double y[], int nstep, int igrd, double *slope){
  double xp=0.0,ybar=0.0,ymax=0.0;
  double tomass, totemp,anchor,xhi,xlo,eeptag;
  int ikip=1, i=0;
  int rv = -1;
  
  //To find the turnoff mass
  for(i=0;i<nstep;i++){
    if(y[i] >= ymax){
      ikip=i;
      ymax=y[i];
    }
  }
  
  if(ToffM(&(x[ikip-1]),&(y[ikip-1]),&xp,&ybar,0)){
    xp = x[ikip];
  }
  else{
    rv=0;
    if(xp < x[ikip]){
      rv = ToffM(&(x[ikip-2]),&(y[ikip-2]),&xp,&ybar,1);
    }
    if(rv){
      xp = x[ikip];
    }
    else{
      if(xp < x[ikip] && xp > x[ikip+1]) printf("possibly incorrect\n");
    }
  }
  
  tomass=xp;
  totemp=ybar;
  
  // Got the turnoff mass
  // To find the slope
  anchor=(tomass-x[0])/(x[nstep-1]-x[0]);
  xhi=5.0;
  xlo=0.005;
  (*slope)=0.5*(xhi+xlo);
  
  while(1){
    eeptag=atan((*slope)*(igrd-1))/atan((*slope)*(nstep-1));
    if(fabs(anchor-eeptag) <= 0.5e-7) break;
    if(eeptag >= anchor) xhi=(*slope);
    else xlo=(*slope);
    (*slope)=0.5*(xhi+xlo);
    if(fabs((xhi-xlo)/(*slope)) <= 0.5e-12) break;
  }
  
  eeptag = atan((*slope)*(igrd-1))/atan((*slope)*(nstep-1));
  
  return;
}

//SD not to self, makesure iorder is fed in as one less than in Fortran code
static int ToffM(double x[4], double yy[4], double *xp, double *ybar, int iorder){
  int n=4,k,l,m;
  double y[n], a,b,c,s,xm,xbar;
  
  for(k=0;k<n;k++)
    y[k]=yy[k];
  
  for(k=0;k<n-1;k++){
    for(l=0;l<n-k-1;l++){
      y[l]=(y[l+1]-y[l])/(x[l+k+1]-x[l]);
    }
  }
  
  a=3.0*y[0];
  b=-2.0*y[0]*(x[3]+x[2]+x[1])+2.0*y[1];
  c=y[0]*(x[2]*x[3]+x[1]*x[3]+x[1]*x[2])-y[1]*(x[3]+x[2])+y[2];
  s=sqrt(b*b-4.0*a*c);
  
  (*xp)=(-b+s)/(2.0*a);
  xm=(-b-s)/(2.0*a);
  
  if((*xp) >= x[iorder] && (*xp) <= x[2])
    xbar=(*xp);
  else if(xm >= x[iorder] && xm <= x[2])
    xbar=xm;
  else{
    if (iorder >= 1)
      return 1;
    else{
      printf("Failed to find the turnoff mass (xp=%f, xm=%f).  Exiting.\n",(*xp),xm);
      exit(1);
    }
  }
  //    (*xp)=(-2.0d0*c) /(b+s)
  //  xm=(-2.0d0*c) /(b-s)
  (*ybar)=y[0];
  for(m=1;m<n;m++){
    (*ybar)=(*ybar)*(xbar-x[m])+y[m];
  }
  (*xp)=xbar;
  return 0;
}

static void initIso(struct yyIsochrone *newIso){
  
  int i,j;
  
  //newIso->nEntries = MAX_YY_ENTRIES;
  //newIso->nFilts   = N_YY_FILTS;
  //allocateGlobalIso(newIso);
  
  newIso->FeH=0.0;
  newIso->age=0.0;
  newIso->logAge=0.0;
  newIso->z=0.0;
  //newIso->y=0.0;
  newIso->nEntries=0;
  for(i=0;i<MAX_YY_ENTRIES;i++){
    newIso->mass[i]=0.0;
    for(j=0;j<N_YY_FILTS;j++)  newIso->mag[i][j]=0.0;
  }
}
/*
void outputYYIso(struct yyIsochrone *iso, FILE *ptr, int byMag){
  
  int m, filt;
  
  fprintf(ptr,"FeH = %f, z = %f, logAge = %f, age = %f, AgbTurnoffMass = %f, nEntries = %d\n",
                iso->FeH,iso->z,iso->logAge,iso->age,iso->AgbTurnoffMass,iso->nEntries);
  fprintf(ptr,"mass U B V R I J H K\n");
  for(m=0;m<iso->nEntries;m++){
    if(byMag){
      fprintf(ptr,"%6.3f ",iso->mass[m]);
      for(filt=0;filt<N_YY_FILTS;filt++)
        fprintf(ptr,"%6.3f ",iso->mag[m][filt]);
    }
    else{
      for(filt=0;filt<N_YY_PARAMS;filt++)
        fprintf(ptr,"%6.3f ",iso->param[m][filt]);
    }
    fprintf(ptr,"\n");
  }
}

*/
static void getFileName(char *path, int z){
  
  char fileNames[][19]={
    "76997z00001a0o2v2\0",
    "7697z0001a0o2v2\0",
    "7688z0004a0o2v2\0",
    "767z001a0o2v2\0",
    "758z004a0o2v2\0",
    "749z007a0o2v2\0",
    "74z01a0o2v2\0",
    "71z02a0o2v2\0",
    "65z04a0o2v2\0",
    "59z06a0o2v2\0",
    "53z08a0o2v2\0"
  };
  
  strcpy(tempFile,"\0");
	strcat(tempFile, path);
  strcat(tempFile,"YYiso/yy00l.x");
  strcat(tempFile,fileNames[z]);
  
}

static double feh2z(double FeH){
  
  double Z=0.0, FeH0=0.0, ZovX=0.0;
  double YYaf = 0.0;
  
  FeH0 = FeH - QUAD(0.,0.,0.3,FeHa2,0.6,FeHa4,YYaf);
  ZovX=pow(10.0,FeH0)*Zsun/Xsun;
  Z=ZovX*(1.0+dydz*zp-yp)/(1.0+ZovX*(1.0+dydz));
  
  return Z;
}

double   getYaleMags(double zamsMass){

/**************************************************************
// This is the subroutine that actually calculates the mags  //
// for a given logAge, FeH, and mass                         //
// Note: this should only be called by evolve after          //
//       deriveYYAgbTip() which does the bounds checking     //
//       and creates the isochrone                           //
**************************************************************/

  int    filt, m = 0;
  
  m = binarySearch(isochrone.mass, isochrone.nEntries, zamsMass);

  for(filt=0;filt < N_YY_FILTS;filt++) {
    if(!useFilt[filt]) continue;
    globalMags[filt] = linInterp(isochrone.mass[m],isochrone.mass[m+1],
                                 isochrone.mag[m][filt],isochrone.mag[m+1][filt],
                                 zamsMass);
  }
  
  return zamsMass;
}


double wdPrecLogAgeYY(double thisFeH, double thisY, double zamsMass)

/*************************************************************************************
last update: 12nov07

Determine WD precursor age by using coefficients derived from fitting to the current
Yale-Yonsei models.  If the models change, these will need to change as well.
Note that the appropriate AgbTurnoffMass mass and lifetime is not the ZAMS mass and lifetime of 
the star currently at the AgbTurnoffMass, but rather refers to the properties of the potentially
higher mass and younger AgbTurnoffMass star that was the WD precursor.
*************************************************************************************/

{

  int i,j;
  double  maCoeff[6], wdPrecLogAge;

  if(zamsMass < 1.0) zamsMass = 1.0;
  else if (zamsMass > 8.0) zamsMass = 8.0;
    wdPrecLogAge = 0.0;
    for(i = 0 ; i < 6 ; i++){
      maCoeff[i] = 0.0;
      for(j=0;j<8;j++) maCoeff[i] +=  coeff[i][j]*pow(thisFeH,j);
      wdPrecLogAge += maCoeff[i]*pow(zamsMass,i);
    }
    
  return wdPrecLogAge;

}

  


  
