#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "mt19937ar.h"
#include "evolve.h"
#include "structures.h"
#include "scatterCluster.h"

static double mass1, mass2, phot[FILTS], exptime[FILTS], sigma[FILTS];
static int stage1, stage2, starID, popID;


// int main()

// {
  
//   int    count[MAXPOPS], nr, nStars[MAXPOPS], wdCount, i, filt, np,
//   filterSet, firstFilt, nPops;
//   double limitSigToNoise, brightLimit, faintLimit, clusterMemberPrior;
//   char   filename[100], line[1000], aFilterName[10];
//   FILE   *r_ptr, *w_ptr;
  
  
//   printf("\n ***You are running BASE08/scatterCluster version %.1f.***\n",VERSION);
  
//   printf("\n Enter simulated cluster file name : ");
//   scanf("%s",filename);
//   if((r_ptr = fopen(filename,"r")) == NULL) {
//     printf("\n\n file %s was not found - exiting ",filename);
//     exit(1);
//   }
  
//   //Scan header line to figure out which photometry set is being used
//   //for(i=0;i<29;i++) fscanf(r_ptr,"%*s ");
//   //fscanf(r_ptr,"%s",aFilterName);
  
//   for(i=0;i<3;i++) fscanf(r_ptr,"%*s ");
//   fscanf(r_ptr,"%s",aFilterName);
//   i=0;
//   while(aFilterName[i]!='1') i++;
//   aFilterName[i] = '\0';
  
  
//   for(filterSet=0;filterSet<3;filterSet++){
//     setFilterNames(filterSet);
//     printf("%d %s %s\n",filterSet, aFilterName, getFilterName(0)); 
//     if(strcmp(aFilterName,getFilterName(0)) == 0) break;
//   }
//   printf("filterSet = %d\n",filterSet);
  
//   fgets(line,1000,r_ptr);		// remove rest of header line
  
//   printf("\n Enter hours of exposure for noise model for each of ");
//   for(filt=0;filt<FILTS;filt++) printf("%s ",getFilterName(filt));
//   //else            printf("\n Enter hours of exposure for noise model for each of band1 band2 ... band8");
//   printf("\n                       e.g., 2.3 1.0 1.2 0. 0. 0. 0. 0.");
//   printf("\n                       where 0. exposure time means unused band. ");
//   scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
//         &exptime[0],&exptime[1],&exptime[2],&exptime[3], &exptime[4],&exptime[5],&exptime[6],&exptime[7],
//         &exptime[8],&exptime[9],&exptime[10],&exptime[11],&exptime[12],&exptime[13]);
  
//   printf("\n Enter limiting signal-to-noise (e.g. 15) : ");
//   scanf("%lf",&limitSigToNoise);
    
//   printf("\n Enter an integer seed: ");
//   scanf("%ld",&seed);
  
//   printf("\n Enter number of populations (excluding brown dwarfs and field stars) :");
//   scanf("%d",&nPops);
  
//   if(nPops > MAXPOPS - 2){
//     printf("Maximum number of populations is %d.  Exiting.\n", MAXPOPS - 2);
//     exit(1);
//   }
//   printf("\n Enter the number of stars to keep for each population (e.g., 100 50 100) :");
//   for(np=0;np<MAXPOPS;np++){nStars[np]=0;}
//   for(np=0;np<nPops;np++){
//     scanf ("%d",&nStars[np]);
//     //printf("%d %d\n",np,nStars[np]);
//   }
  
//   printf("\n Enter number of brown dwarfs to include (e.g. 0): ");
//   scanf("%d",&nStars[BDPOP]);
//   if(nStars[BDPOP] < 0) nStars[BDPOP] = 0;
  
  
//   printf("\n Enter number of field stars to include (e.g. 0): ");
//   scanf("%d",&nStars[FSPOP]);
//   if(nStars[FSPOP] < 0) nStars[FSPOP] = 0;
    
//   printf("\n Enter bright and faint and cut-off mags, and their filter: ");
//   scanf("%lf %lf %d",&brightLimit, &faintLimit, &firstFilt);		// brightLimit primarily used to cut off RGB
  
//   printf("\n Enter output file name : ");
//   scanf("%s",filename);
//   if((w_ptr = fopen(filename,"w")) == NULL) {
//     printf("\n\n file %s not available for writing - exiting ",filename);
//     exit(1);
//   }
//   printf("\n");
  
//   clusterMemberPrior = 0;
//   for(np=0;np<nPops;np++) clusterMemberPrior += nStars[np];
//   clusterMemberPrior =  clusterMemberPrior / (clusterMemberPrior + nStars[MAXPOPS-1]);
//   if(clusterMemberPrior > 0.99) clusterMemberPrior = 0.99;
//   if(clusterMemberPrior < 0.01) clusterMemberPrior = 0.01;

//   init_genrand(seed);
  
//   fprintf(w_ptr,"    id ");
//   for(filt=0;filt<FILTS;filt++) if(exptime[filt] > EPS) fprintf(w_ptr,"%6s ",   getFilterName(filt));
//   for(filt=0;filt<FILTS;filt++) if(exptime[filt] > EPS) fprintf(w_ptr,"sig%-5s ",getFilterName(filt));
  
//   fprintf(w_ptr,"mass1 massRatio stage1 CMprior useDuringBurnIn pop\n");
  
//   for(np=0;np<MAXPOPS;np++) count[np] = 0;
//   wdCount = 0;
  
//   np=0;
  
//   while((nr = readLine(r_ptr)) != EOF) {
//     // If this is the first star of a new pop (or the first field star or brown dwarf
//     if(popID != np) np = popID;
    
//     // skip this star if it doesn't meet the various criteria to be included
//     if(magCutoff(firstFilt,brightLimit,faintLimit) == 0) continue;
//     if(stageCutoff(np==FSPOP) == 0) continue;
//     if(scatterPhot(limitSigToNoise) == 0) continue;
//     if(mass1 < 0.25 && mass2 < 0.25 && stage1 != BD) continue;	// limit on modelSet=4 is 0.4 Mo - what to do? 
//     wdCount += outputScatter(w_ptr, np==FSPOP, clusterMemberPrior, np);
//     count[np]++;
    
//     // If we have enough stars of this type
//     if(count[np] == nStars[np]) {
//       // If these are field stars, we're done
//       if(np == FSPOP) break;
//       // If not...
//       else {
//         // Skip lines in the file until you get to...
//         while((nr = fscanf(r_ptr,"%d ",&popID)) != EOF) {
//           //...the brown dwarfs...
//           if(popID==BDPOP){
//             np=BDPOP;
//             fgets(line,1000,r_ptr);
//             break;
//           }
//           //...or the field stars
//           else if(popID==FSPOP){
//             np=FSPOP;
//             fgets(line,1000,r_ptr);
//             break;
//           }
//           //or the next population
//           else if(popID==np+1){
//             np++;
//             fgets(line,1000,r_ptr);
//             break;
//           }
//           fgets(line,1000,r_ptr);
//         }
//         //fseek(r_ptr, -7, SEEK_CUR);
//       }
//     }
//   }
//   for(np=0;np<MAXPOPS;np++){
//     //printf("np = %d\n",np);
//     if(count[np] < nStars[np]) {
//       printf(" Warning - not enough (%d) stars in input file to keep desired number (%d) of stars ", count[np],nStars[np]);
//       if(np==BDPOP) printf("for brown dwarfs\n");
//       else if(np==FSPOP) printf("for field stars\n");
//       else printf("for population %d\n",np+1);
//     }
//   }
    
//   printf(" There %s %d WD%s in this scatter file, %s\n",wdCount == 1 ? "is" : "are",wdCount,wdCount ==1 ? "" : "s",filename);
  
//   fclose(r_ptr);
//   fclose(w_ptr);
  
//   return(0);
// }
  

static double s2nCoeffs[][2]=
{
  {9.33989,  0.3375778},	// U 
  {10.0478,  0.3462758},	// B 
  {10.48098, 0.368201 },	// V 
  {10.71151, 0.3837847},	// R 
  {10.61035, 0.3930941},	// I 
  {9.282385, 0.386258 },	// J 
  {9.197463, 0.3970419},	// H 
  {9.024068, 0.3985604},  // K
  {9.024068, 0.3985604},  // IRAC Blue
  {9.024068, 0.3985604},  // IRAC Red
  {9.024068, 0.3985604},  // Band 1
  {9.024068, 0.3985604},  // Band 2
  {9.024068, 0.3985604},  // Band 3
  {9.024068, 0.3985604}   // Band 4
};

double signalToNoise(double mag, double exptime, int filter)

/*
 This is an approximation to the results one would obtain in one hour with the KPNO
 4m + Mosaic (UBVRI) or Flamingos (JHK) per band, assuming dark time, seeing=1.1", 
 airmass=1.2, and then scaling from their by sqrt(exptime).  I further approximated
 the CCDTIME results (run on the NOAO webste) with linear fits of mag vs. log(S/N).
 */

{
  
  double s2n, logS2N;
  
  if(filter >= FILTS || filter < 0){
    printf("filter (%d) out of range - exiting\n",filter);
    exit(1);
  }
  
  // Scatter BD photometry at 5%
  if(filter > 7) return 1.0/0.05;
  
  logS2N  = s2nCoeffs[filter][0] - s2nCoeffs[filter][1] * mag;   
  s2n  = pow(10., logS2N);
  s2n *= sqrt(exptime);
  
  return(s2n);
  
}

int readLine(FILE* filePtr){
  
  int nr,filt;
  nr = fscanf(filePtr, "%d %d %lf ", &popID, &starID, &mass1);
  for(filt=0;filt<FILTS;filt++)  nr = fscanf(filePtr, "%*f ");
  nr = fscanf(filePtr, "%d %*f %*d %*f %*f %lf ",&stage1,&mass2);
  for(filt=0;filt<FILTS;filt++)  nr = fscanf(filePtr, "%*f ");
  nr = fscanf(filePtr, "%d %*f %*d %*f %*f ",&stage2);
  for(filt=0;filt<FILTS;filt++)  nr = fscanf(filePtr, "%lf ",&phot[filt]);
  
  //printf("%d %d %f %d %f %d\n",popID,starID,mass1,stage1,mass2,stage2);
  
  if(starID==EOF){
    printf("\nFatal error in readLine.  Exiting.\n");
    exit(1);
  }
  
  if(nr == EOF) return nr;
  return starID;
  
}

// Returns a 1 if the star is within the magnitude cutoff limits, 0 otherwise
int magCutoff(int firstFilt, double brightLimit, double faintLimit){
  if(phot[firstFilt] > 99. && stage1 != BD) return 0;         // check if real object.  if not, skip 
  if(phot[firstFilt] < brightLimit && stage1 != BD) return 0;	// typically used to remove red giants 
  if(phot[firstFilt] > faintLimit && stage1 != BD) return 0;	// typically used to remove WDs and lower MS 
  return 1;
}

// Returns a 1 if this star is a MS-MS binary or single WD, 0 otherwise (for cluster members)
// For field stars, it returns a 0 if the system contains a neutron star or black hole
int stageCutoff(int isFS){
  if(stage1 == NSBH || stage2 == NSBH) return 0;
  if(stage1 == WD && mass2 > 0.0 && !isFS) return 0;	// TEMPORARY KLUDGE -- ignore binaries of MS/RG + WDs and WD + WD
  if(stage2 == WD && mass1 > 0.0 && !isFS) return 0;
  return 1;
}

//Returns a 1 if this star has enough S/N, 0 otherwise
int scatterPhot(double limitSigToNoise){
  int filt;
  double sigToNoise;
  
  for(filt=(stage1==BD?8:0);filt<(stage1==BD?FILTS:8);filt++){
    if(exptime[filt] < EPS) continue;
    sigToNoise  = signalToNoise(phot[filt], exptime[filt], filt);
    if(sigToNoise < limitSigToNoise) {		// large photometric errors can lock mcmc during burnin 
      printf("Warning: star %4d, mass1=%.3f, stage1=%d, filter=%d, S/N = %.3f < %.1f (user limit) - skipping.\n",
             starID, mass1,stage1,filt,sigToNoise,limitSigToNoise);
      return 0;
    }
    
    sigma[filt] = 1./(sigToNoise);
    if(sigma[filt] < 0.005) sigma[filt] = 0.005;
    phot[filt] += gen_norm(0., sigma[filt]);
  }
  return 1;
}

//Returns a 1 for cluser member WDs, 0 otherwise (to keep track of the WD count for the cluster)
int outputScatter(FILE* w_ptr, int isFS, double clusterMemberPrior, int popID){
  int tempStage, filt;
  double tempMass, tempMassRatio;
  
  if(mass1 > mass2) {            // use to set starter mass for mcmc 
    tempMass      = mass1;       // (since higher mass star dominates the photometry) 
    tempMassRatio = mass2/mass1;
    tempStage     = stage1;
  }
  else {
    tempMass      = mass2;
    tempMassRatio = mass1/mass2;
    tempStage     = stage2;
  }
  
  fprintf(w_ptr,"%6d ",starID);
  for(filt=0;filt<FILTS;filt++) if(exptime[filt] > EPS) fprintf(w_ptr,"%6.3f ",phot[filt]);
  for(filt=0;filt<8;filt++){
    if(exptime[filt] > EPS){
      if(stage1 == BD) fprintf(w_ptr,"%8.5f ", -1.0);
      else fprintf(w_ptr,"%8.6f ",sigma[filt]);
    }
  }
  for(filt=8;filt<FILTS;filt++){
    if(exptime[filt] > EPS){
      if(stage1 == BD) fprintf(w_ptr,"%8.6f ",sigma[filt]);
      else fprintf(w_ptr,"%8.5f ",-1.0);
    }
  }
  
  fprintf(w_ptr,"%8.3f %6.3f %3d %6.3f   %d    %d\n",tempMass,(tempMassRatio>0.001?tempMassRatio:0.00),tempStage,clusterMemberPrior,!isFS,popID);
  
  if(tempStage == WD && !isFS) return 1;
  else return 0;
}



