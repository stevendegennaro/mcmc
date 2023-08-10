import structures as st
from parameters import Params
import pandas as pd

#   int    count[MAXPOPS], nr, nStars[MAXPOPS], wdCount, i, filt, np,
#   filterSet, firstFilt, nPops;
#   double limitSigToNoise, brightLimit, faintLimit, clusterMemberPrior;
#   char   filename[100], line[1000], aFilterName[10];
#   FILE   *r_ptr, *w_ptr;
	

params = Params()

print(f"***You are running BASE08/scatterCluster version {st.VERSION:.1f}.***\n")
	
filterSet = params.filterSet
print(f"filterSet = {filterSet}")

filename = params.outfile
with open(filename,"r") as f:
	print(f"Reading from file... {filename}")
	starsdf = pd.read_csv(f,sep=" ")
	pCluster = st.dataFrameToCluster(starsdf,filterSet)
	print(starsdf,len(starsdf.index))
	for j in range(pCluster.contents.nStars):
		print(st.c.getStarPtr(pCluster,j).contents.massRatio)
		
# 	printf("\n Enter hours of exposure for noise model for each of ");
# 	for(filt=0;filt<FILTS;filt++) printf("%s ",getFilterName(filt));
# 	//else            printf("\n Enter hours of exposure for noise model for each of band1 band2 ... band8");
# 	printf("\n                       e.g., 2.3 1.0 1.2 0. 0. 0. 0. 0.");
# 	printf("\n                       where 0. exposure time means unused band. ");
# 	scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
# 				&exptime[0],&exptime[1],&exptime[2],&exptime[3], &exptime[4],&exptime[5],&exptime[6],&exptime[7],
# 				&exptime[8],&exptime[9],&exptime[10],&exptime[11],&exptime[12],&exptime[13]);
	
# 	printf("\n Enter limiting signal-to-noise (e.g. 15) : ");
# 	scanf("%lf",&limitSigToNoise);
		
# 	printf("\n Enter an integer seed: ");
# 	scanf("%ld",&seed);
	
# 	printf("\n Enter number of populations (excluding brown dwarfs and field stars) :");
# 	scanf("%d",&nPops);
	
# 	if(nPops > MAXPOPS - 2){
# 		printf("Maximum number of populations is %d.  Exiting.\n", MAXPOPS - 2);
# 		exit(1);
# 	}
# 	printf("\n Enter the number of stars to keep for each population (e.g., 100 50 100) :");
# 	for(np=0;np<MAXPOPS;np++){nStars[np]=0;}
# 	for(np=0;np<nPops;np++){
# 		scanf ("%d",&nStars[np]);
# 		//printf("%d %d\n",np,nStars[np]);
# 	}
	
# 	printf("\n Enter number of brown dwarfs to include (e.g. 0): ");
# 	scanf("%d",&nStars[BDPOP]);
# 	if(nStars[BDPOP] < 0) nStars[BDPOP] = 0;
	
	
# 	printf("\n Enter number of field stars to include (e.g. 0): ");
# 	scanf("%d",&nStars[FSPOP]);
# 	if(nStars[FSPOP] < 0) nStars[FSPOP] = 0;
		
# 	printf("\n Enter bright and faint and cut-off mags, and their filter: ");
# 	scanf("%lf %lf %d",&brightLimit, &faintLimit, &firstFilt);		// brightLimit primarily used to cut off RGB
	
# 	printf("\n Enter output file name : ");
# 	scanf("%s",filename);
# 	if((w_ptr = fopen(filename,"w")) == NULL) {
# 		printf("\n\n file %s not available for writing - exiting ",filename);
# 		exit(1);
# 	}
# 	printf("\n");
	
# 	clusterMemberPrior = 0;
# 	for(np=0;np<nPops;np++) clusterMemberPrior += nStars[np];
# 	clusterMemberPrior =  clusterMemberPrior / (clusterMemberPrior + nStars[MAXPOPS-1]);
# 	if(clusterMemberPrior > 0.99) clusterMemberPrior = 0.99;
# 	if(clusterMemberPrior < 0.01) clusterMemberPrior = 0.01;

# 	init_genrand(seed);
	
# 	fprintf(w_ptr,"    id ");
# 	for(filt=0;filt<FILTS;filt++) if(exptime[filt] > EPS) fprintf(w_ptr,"%6s ",   getFilterName(filt));
# 	for(filt=0;filt<FILTS;filt++) if(exptime[filt] > EPS) fprintf(w_ptr,"sig%-5s ",getFilterName(filt));
	
# 	fprintf(w_ptr,"mass1 massRatio stage1 CMprior useDuringBurnIn pop\n");
	
# 	for(np=0;np<MAXPOPS;np++) count[np] = 0;
# 	wdCount = 0;
	
# 	np=0;
	
# 	while((nr = readLine(r_ptr)) != EOF) {
# 		// If this is the first star of a new pop (or the first field star or brown dwarf
# 		if(popID != np) np = popID;
		
# 		// skip this star if it doesn't meet the various criteria to be included
# 		if(magCutoff(firstFilt,brightLimit,faintLimit) == 0) continue;
# 		if(stageCutoff(np==FSPOP) == 0) continue;
# 		if(scatterPhot(limitSigToNoise) == 0) continue;
# 		if(mass1 < 0.25 && mass2 < 0.25 && stage1 != BD) continue;	// limit on modelSet=4 is 0.4 Mo - what to do? 
# 		wdCount += outputScatter(w_ptr, np==FSPOP, clusterMemberPrior, np);
# 		count[np]++;
		
# 		// If we have enough stars of this type
# 		if(count[np] == nStars[np]) {
# 			// If these are field stars, we're done
# 			if(np == FSPOP) break;
# 			// If not...
# 			else {
# 				// Skip lines in the file until you get to...
# 				while((nr = fscanf(r_ptr,"%d ",&popID)) != EOF) {
# 					//...the brown dwarfs...
# 					if(popID==BDPOP){
# 						np=BDPOP;
# 						fgets(line,1000,r_ptr);
# 						break;
# 					}
# 					//...or the field stars
# 					else if(popID==FSPOP){
# 						np=FSPOP;
# 						fgets(line,1000,r_ptr);
# 						break;
# 					}
# 					//or the next population
# 					else if(popID==np+1){
# 						np++;
# 						fgets(line,1000,r_ptr);
# 						break;
# 					}
# 					fgets(line,1000,r_ptr);
# 				}
# 				//fseek(r_ptr, -7, SEEK_CUR);
# 			}
# 		}
# 	}
# 	for(np=0;np<MAXPOPS;np++){
# 		//printf("np = %d\n",np);
# 		if(count[np] < nStars[np]) {
# 			printf(" Warning - not enough (%d) stars in input file to keep desired number (%d) of stars ", count[np],nStars[np]);
# 			if(np==BDPOP) printf("for brown dwarfs\n");
# 			else if(np==FSPOP) printf("for field stars\n");
# 			else printf("for population %d\n",np+1);
# 		}
# 	}
		
# 	printf(" There %s %d WD%s in this scatter file, %s\n",wdCount == 1 ? "is" : "are",wdCount,wdCount ==1 ? "" : "s",filename);
	
# 	fclose(r_ptr);
# 	fclose(w_ptr);
	
# 	return(0);
# }