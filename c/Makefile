CC=gcc
OFLAGS=-O3
CFLAGS=-O2 -g -pg -Wall
CFLAGS=-O3 -ffast-math -fstrict-aliasing -Wall
#CFLAGS=-I/Users/steven/anaconda3/envs/mcmc/include/python3.10 -I/Users/steven/anaconda3/envs/mcmc/in.10
LFLAGS=-g -pg
LFLAGS=

OBJS = binSearch.o gBaraffeMag.o gBergMag.o \
       gChabMag.o gDsedMag.o gGirMag.o gYaleMag.o evolve.o           \
       ifmr.o linInterp.o loadModels.o msRgbEvol.o structures.o       \
       wdCooling.o wdEvol.o

RAND_OBJS = gennorm.o mt19937ar.o

PLOT_OBJS = plotCluster.o

SCATTER_OBJS = scatterCluster.o

STAT_OBJS = densities.o samplers.o decide.o leastSquares.o poly.o solve.o mcmc.o

MAKE_CMD_OBJS = makeCMD.o 


all : $(OBJS) $(RAND_OBJS) $(SCATTER_OBJS) $(SIM_OBJS) $(STAT_OBJS) $(MAKE_CMD_OBJS)
	$(CC) $(LFLAGS) -o mcmc $(OBJS) $(RAND_OBJS) $(STAT_OBJS) 
	$(CC) $(LFLAGS) -o simCluster $(OBJS) $(RAND_OBJS) $(SIM_OBJS) 
	$(CC) $(LFLAGS) -o scatterCluster $(SCATTER_OBJS) $(RAND_OBJS) -lm
	$(CC) $(LFLAGS) -o makeCMD $(OBJS) $(MAKE_CMD_OBJS) 

mcmc : $(OBJS) $(STAT_OBJS) $(RAND_OBJS)
	$(CC) $(LFLAGS) -o mcmc $(OBJS) $(STAT_OBJS) $(RAND_OBJS)  

plotCluster : $(OBJS) $(PLOT_OBJS) $(RAND_OBJS)
	$(CC) -fPIC $(OBJS) $(PLOT_OBJS) $(RAND_OBJS) -shared -o plotCluster.so

#scatterCluster : $(SCATTER_OBJS) $(RAND_OBJS)
#	$(CC) -fPIC $(LFLAGS) -shared $(SCATTER_OBJS) $(RAND_OBJS) -o scatterCluster.so 

makeCMD : $(MAKE_CMD_OBJS) $(OBJS)
	$(CC) $(LFLAGS) -o makeCMD $(MAKE_CMD_OBJS) $(OBJS) 

# structures.so : structures.o
# 	gcc -fPIC structures.o -shared -o structures.so

binSearch.o : binSearch.c
	$(CC) -c $(CFLAGS) binSearch.c

decide.o : decide.c 
	$(CC) -c $(CFLAGS) decide.c

densities.o : densities.c 
	$(CC) -c $(CFLAGS) densities.c

drawFromIMF.o : drawFromIMF.c 
	$(CC) -c $(CFLAGS) drawFromIMF.c

evolve.o : evolve.c 
	$(CC) -c $(CFLAGS) evolve.c

gennorm.o : gennorm.c 
	$(CC) -c $(CFLAGS) gennorm.c

gBaraffeMag.o : gBaraffeMag.c
	$(CC) -c $(CFLAGS) gBaraffeMag.c

gBergMag.o : gBergMag.c
	$(CC) -c $(CFLAGS) gBergMag.c

gChabMag.o : gChabMag.c 
	$(CC) -c $(CFLAGS) gChabMag.c

gDsedMag.o : gDsedMag.c 
	$(CC) -c $(CFLAGS) gDsedMag.c

gGirmag.o : gGirMag.c 
	$(CC) -c $(CFLAGS) gGirMag.c

gYalemag.o : gYaleMag.c 
	$(CC) -c $(CFLAGS) gYaleMag.c

ifmr.o : ifmr.c 
	$(CC) -c $(CFLAGS) ifmr.c

leastSquares.o : leastSquares.c 
	$(CC) -c $(CFLAGS) leastSquares.c

linInterp.o :  linInterp.c 
	$(CC) -c $(CFLAGS) linInterp.c

loadModels.o : loadModels.c 
	$(CC) -c $(CFLAGS) loadModels.c

makeCMD.o : makeCMD.c 
	$(CC) -c $(CFLAGS) makeCMD.c

mcmc.o : mcmc.c
	$(CC) -c $(CFLAGS) mcmc.c

msRgbEvol.o : msRgbEvol.c
	$(CC) -c $(CFLAGS) msRgbEvol.c

mt19937ar.o : mt19937ar.c
	$(CC) -c $(CFLAGS) mt19937ar.c

poly.o : poly.c
	$(CC) -c $(CFLAGS) poly.c

samplers.o : samplers.c 
	$(CC) -c $(CFLAGS) samplers.c 

scatterCluster.o : scatterCluster.c 
	$(CC) -c $(CFLAGS) scatterCluster.c 

plotCluster.o : plotCluster.c 
	$(CC) -c $(CFLAGS) plotCluster.c 

solve.o : solve.c
	$(CC) -c $(CFLAGS) solve.c

structures.o : structures.c
	$(CC) -c $(CFLAGS) structures.c

wdEvol.o :  wdEvol.c
	$(CC) -c $(CFLAGS) wdEvol.c

wdCooling.o :  wdCooling.c 
	$(CC) -c $(CFLAGS) wdCooling.c

clean : 
	rm -f *.o *.so  mcmc *~ 

#gprof mcmc gmon.out > proc.out
