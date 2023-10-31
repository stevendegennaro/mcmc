import ctypes
from parameters import Params
import numpy as np
import pandas as pd
from drawFromIMF import drawFromIMF
from typing import TypeVar, List
from getCMacros import getCMacros

pdDataFrame = TypeVar('pd.core.frame.DataFrame')

# Read in the shared object file that contains all the variable,
# structures, and functions in C that will be needed by Python
c = ctypes.CDLL("../c/simCluster.so")

## Pull in and parse the #define macros in evolve.h and structures.h
macros = getCMacros("../c/structures.h")
macros.update(getCMacros("../c/evolve.h"))
for key in macros:
	exec(key + f" = {macros[key]}")

#### Use ctypes to create classes that mirror the structures
#### used by simCluster and MCMC for models, clusters, and stars

#### Model structure ####
class model(ctypes.Structure):
	_fields_ = [('evoModel',ctypes.c_int),
				('brownDwarfEvol',ctypes.c_int),
				('mainSequenceEvol',ctypes.c_int),
				('IFMR',ctypes.c_int),
				('WDcooling',ctypes.c_int),
				('WDatm',ctypes.c_int),
				('filterSet',ctypes.c_int),
				('numFilts',ctypes.c_int),
				('needFS',ctypes.c_int),
				('minMass',ctypes.c_double)]
	def __str__(self):
		string = "Model Parameters:\n"
		for i in range(len(self._fields_)):
			string += f"\t{self._fields_[i][0]:18}: "
			string += f" {getattr(self,self._fields_[i][0])}\n"
		return string

#### Model pointer ####
modelPtr = ctypes.POINTER(model)

#### Star Structure ####
class star(ctypes.Structure):	
	_fields_ = [('id',ctypes.c_int),
				('obsPhot',ctypes.c_double * FILTS), 
				('photometry',ctypes.c_double * FILTS), 
				('variance',ctypes.c_double * FILTS),
				('useFilt',ctypes.c_double * FILTS),
				('U',ctypes.c_double),
				('massRatio',ctypes.c_double),
				('status',ctypes.c_int * 2),
				('wdType',ctypes.c_int * 2),
				('isFieldStar',ctypes.c_int),
				('useDuringBurnIn',ctypes.c_int),
				('clustStarPriorDens',ctypes.c_double),
				('clustStarProposalDens',ctypes.c_double),
				('beta',(ctypes.c_double * 2) * NPARAMS),
				('betaMassRatio',ctypes.c_double * 2),
				('meanU',ctypes.c_double),
				('varU',ctypes.c_double),
				('meanMassRatio',ctypes.c_double),
				('varMassRatio',ctypes.c_double),
				('UStepSize',ctypes.c_double),
				('massRatioStepSize',ctypes.c_double),
				('boundsFlag',ctypes.c_int),
				('wdLogTeff',ctypes.c_double * 2),
				('massNow',ctypes.c_double * 2)]
	def __str__(self):
		string = f"Star Properites for star {self.id}:\n"
		for i in range(1,len(self._fields_)):
			string += f"\t{self._fields_[i][0]:22}: "
			if type(getattr(self,self._fields_[i][0])) == int or type(getattr(self,self._fields_[i][0])) == float:
				string += f" {getattr(self,self._fields_[i][0])}\n"
			else:
				for j in range(getattr(self,self._fields_[i][0])._length_):
					array = getattr(self,self._fields_[i][0])
					if j != 0: string += f"\t{'':22}  "
					string += f"{array[j]}\n"
				string += "\n"
		return string

#### Star pointer ####
starPtr = ctypes.POINTER(star)

#### Cluster Structure ####
class cluster(ctypes.Structure):
	_fields_ = [('nStars',ctypes.c_int),
				('M_wd_up',ctypes.c_double),
				('parameter',ctypes.c_double * NPARAMS),
				('stepSize',ctypes.c_double * NPARAMS),
				('mean',ctypes.c_double * NPARAMS),
				('priorVar',ctypes.c_double * NPARAMS),
				('priorMean',ctypes.c_double * NPARAMS),
				('betamabs',ctypes.c_double),
				('betaFabs',ctypes.c_double),
				('betaFY',ctypes.c_double),
				('betaAgeMod',ctypes.c_double * 3),
				('AGBt_zmass',ctypes.c_double),
				('varScale',ctypes.c_double),
				('evoModels',model),
				('stars',starPtr)]

	def __str__(self):
		string = "\nCluster Parameters:\n"
		for i in range(len(self._fields_) - 1):
			string += f"\t{self._fields_[i][0]:11}: "
			if type(getattr(self,self._fields_[i][0])) == int or type(getattr(self,self._fields_[i][0])) == float:
				string += f" {getattr(self,self._fields_[i][0])}\n"
			elif isinstance(getattr(self,self._fields_[i][0]),model):
				pass
			else:
				for j in range(getattr(self,self._fields_[i][0])._length_):
					array = getattr(self,self._fields_[i][0])
					string += f"{array[j]:7.4f} "
				string += "\n"
		string += "\n"
		string += str(self.evoModels)
		# for j in range(self.nStars):
		# 	string += str(c.getStarPtr(self,j).contents)
		return string

#### Cluster pointer ####
clusterPtr = ctypes.POINTER(cluster)

#### Argument and result types for functions defined in c ####
c.getStarPtr.restype = starPtr
c.getStarPtr.argtypes = [clusterPtr]

c.allocCluster.restype = clusterPtr
c.allocCluster.argtypes = [ctypes.c_int]

c.getMass1.restype = ctypes.c_double
c.getMass1.argtypes = [starPtr,clusterPtr]

c.getMass2.restype = ctypes.c_double
c.getMass2.argtypes = [starPtr,clusterPtr]

c.getFilterName.restype = ctypes.c_char_p
c.getFilterName.argtypes = [ctypes.c_int]

c.getLTau.restype = ctypes.c_double
c.getLTau.argtypes = [ctypes.c_int]

c.getFilterName.restype = ctypes.c_char_p

c.setMass1.restype = None
c.setMass1.argtypes = [starPtr,clusterPtr,ctypes.c_double]

c.setMass2.restype = None
c.setMass2.argtypes = [starPtr,clusterPtr,ctypes.c_double]

c.getAgeLimit.restype = ctypes.c_double
c.getAgeLimit.argtypes = [ctypes.c_int]

c.getFeHLimit.restype = ctypes.c_double
c.getFeHLimit.argtypes = [ctypes.c_int]

# Get the minimum mass for different model sets
def getMinMass(modelSet: int) -> float:
	if(modelSet == YALE): return  0.4
	if(modelSet == DSED): return 0.25
	return 0.15

# Get a list of the filter names
def getFilters():
	return [c.getFilterName(filt).decode() for filt in range(FILTS) if c.getUseFilt(filt)]

### Creates a cluster, allocates memory for it (in C), initializes
### memory for the cluster (which in turn initiailzes the memory for the
### array of stars. Returns a pointer to the cluster
### Does not evolve the cluster, just creates the structure in memory
def createCluster(params: Params) -> clusterPtr:
	pCluster = c.allocCluster(params.nSystems[0])
	pCluster.contents.evoModels.mainSequenceEvol = params.msRgbModels
	pCluster.contents.evoModels.filterSet = params.filterSet
	pCluster.contents.evoModels.IFMR = params.IFMR
	pCluster.contents.evoModels.brownDwarfEvol = params.BDmodel
	pCluster.contents.evoModels.minMass = getMinMass(pCluster.contents.evoModels.mainSequenceEvol)
	return pCluster

### Change to adifferent model set. Cluster needs to be re-evolved to update photometry
def switchModels(pCluster: clusterPtr, params: Params) -> None:
	pCluster.contents.evoModels.mainSequenceEvol = params.msRgbModels
	pCluster.contents.evoModels.filterSet = params.filterSet
	pCluster.contents.evoModels.IFMR = params.IFMR
	pCluster.contents.evoModels.brownDwarfEvol = params.BDmodel
	pCluster.contents.evoModels.minMass = getMinMass(pCluster.contents.evoModels.mainSequenceEvol)

### Set initial cluster parameters. Cluster needs to be evolved to get photometry
def setClusterParams(pCluster: clusterPtr, params: Params, pop: int):
	pCluster.contents.nStars = params.nSystems[pop]
	c.reallocStarsArray(pCluster)
	pCluster.contents.M_wd_up = params.WDMassUp[pop]
	pCluster.contents.parameter[AGE] = params.logClusAge[pop]
	pCluster.contents.parameter[YYY] = params.Y[pop]
	pCluster.contents.parameter[FEH] = params.Fe_H[pop]
	pCluster.contents.parameter[MOD] = params.distMod[pop]
	pCluster.contents.parameter[ABS] = params.Av[pop]

### Changes the values of cluster parameters without reallocating stars array
### Cluster needs to be re-evolved to get new photometry
def changeClusterParams(pCluster: clusterPtr, params: Params, pop: int):
	pCluster.contents.M_wd_up = params.WDMassUp[pop]
	pCluster.contents.parameter[AGE] = params.logClusAge[pop]
	pCluster.contents.parameter[YYY] = params.Y[pop]
	pCluster.contents.parameter[FEH] = params.Fe_H[pop]
	pCluster.contents.parameter[MOD] = params.distMod[pop]
	pCluster.contents.parameter[ABS] = params.Av[pop]

### Get a single isochrone directly from evolve.c
### The Girardi models don't store isocrhones, so we
### create one of our own for that model set (0)
def populateIso(pIso: clusterPtr) -> None:
	if pIso.contents.evoModels.mainSequenceEvol == 0:
		pIso.contents.nStars = int((pIso.contents.M_wd_up - pIso.contents.evoModels.minMass)/0.01)
		c.reallocStarsArray(pIso)
		c.deriveAgbTipMass(pIso)
		j = 0
		mass = pIso.contents.evoModels.minMass
		while mass < pIso.contents.AGBt_zmass:
			pIso.contents.stars[j].U = mass;
			pIso.contents.stars[j].massRatio = 0.0
			mass += 0.01
			j += 1
		mass = pIso.contents.AGBt_zmass + .005
		while mass < pIso.contents.M_wd_up:
			pIso.contents.stars[j].U = mass;
			pIso.contents.stars[j].massRatio = 0.0
			mass += 0.01
			j += 1
		c.evolve(pIso,-1)
	else:
		c.returnIso(pIso)


### Randomly draws from the IMF to create a cluster
def populateCluster(pCluster: clusterPtr, fractionBinary: float = 0.0) -> None:
	for j in range(pCluster.contents.nStars):	# for all systems in the cluster
		star = c.getStarPtr(pCluster,j)
		star.contents.id = j
		# create single stars in arbitrary mass order 
		while True:
			star.contents.U = drawFromIMF()
			if star.contents.U > pCluster.contents.evoModels.minMass:
				break
		## Need to evolve primaries to see if they are neutron stars or black holes
		c.evolve(pCluster, j)
		# create binaries among appropriate fraction 
		if star.contents.status[0] != NSBH \
			and star.contents.status[0] != WD:
			if np.random.random() < fractionBinary:
				massRatio = np.random.normal(1,0.5)
				while massRatio < 0 or massRatio > 1:
					massRatio = np.random.normal(1,0.5)
				star.contents.massRatio = massRatio
			else:
				star.contents.massRatio = 0.0
		else:
			star.contents.massRatio = 0.0

	#print(pCluster.contents)

	c.evolve(pCluster, -1)							# Evolve stars

### Create random field stars
def populateFieldStars(params: Params, minV: float, maxV: float) -> None:

	if params.nFieldStars > 0:
		minFeH = -2.5
		maxFeH = 0.56
		minAge = 8.4
		maxAge = 10.17

		pCluster = createCluster(params)
		setClusterParams(pCluster,params,0)
		pCluster.contents.evoModels.mainSequenceEvol = DSED
		minMass = pCluster.contents.evoModels.minMass

		pCluster.contents.nStars = 1
		c.reallocStarsArray(pCluster)
		starContents = c.getStarPtr(pCluster,0).contents
		starContents.status[0] = MSRG
		fieldStars = []

		for j in range(params.nFieldStars):
			while True:     
				# Draw a new age and metallicity
				pCluster.contents.parameter[AGE] = minAge + (maxAge-minAge) * np.random.random()
				pCluster.contents.parameter[FEH] = minFeH + (maxFeH-minFeH) * np.random.random()
				
				# Determine a new distance, weighted so 
				# there are more stars behind than in front
				pCluster.contents.parameter[MOD] = params.distMod[0] - 12.0 + np.log10(np.power(10,(np.power(np.power(26.0,3.0) * np.random.random(),1.0/3.0))))
				
				populateCluster(pCluster, 0.5)
							
				# print(j,minV,maxV,starContents.photometry[2])
				if (starContents.photometry[2] > minV and 
					starContents.photometry[2] < maxV and
					(starContents.photometry[1] 
						- starContents.photometry[2]) > -0.5 and
					(starContents.photometry[1] 
						- starContents.photometry[2]) < 1.7):
					break
			starContents.id = 9000 + j

			if j == 0:
				fieldstardf = clusterToDataFrame(pCluster)
			fieldStars.extend(clusterToList(pCluster))


		fieldstardf = pd.DataFrame(fieldStars,columns = fieldstardf.columns)
		columns = [c.getFilterName(filt).decode() for filt in range(FILTS) if c.getUseFilt(filt)]
		fieldstardf[fieldstardf[columns] > 99] = 99.999
		fieldstardf.insert(0,column = "pop",value=9)

		return fieldstardf

	else:
		return None

### Takes the photometry and other data from the cluster structure and puts it
### into a list in the same order as the old simCluster. Returns the list
def clusterToList(pCluster: clusterPtr) -> List[List]:

	data = []
	for j in range(pCluster.contents.nStars):
		star = c.getStarPtr(pCluster,j)
		line = [star.contents.id, 
				c.getMass1(star,pCluster),
				star.contents.status[0],
				star.contents.massNow[0] if star.contents.status[0] == 3 else 0.0,
				star.contents.wdLogTeff[0],
				0,
				c.getLTau(0),
				c.getMass2(star,pCluster),
				star.contents.status[1],
				star.contents.massNow[1] if star.contents.status[1] == 3 else 0.0,
				star.contents.wdLogTeff[1],
				0,
				c.getLTau(1)]
		line.extend([star.contents.photometry[filt] for filt in range(FILTS) if c.getUseFilt(filt)])
		data.append(line)
	return data

### Takes the photometry and other data from the cluster structure and puts it
### into a DataFrame in the same order as the old simCluster. Returns the DataFrame
def clusterToDataFrame(pCluster: clusterPtr) -> pdDataFrame:
	data = clusterToList(pCluster)
	columns =  ["id","mass1","stage1", "wdM1", "wdType1", "wdLogTeff1", "ltau1","mass2","stage2", "wdM2", "wdType2", "wdLogTeff2", "ltau2"]
	columns.extend([c.getFilterName(filt).decode() for filt in range(FILTS) if c.getUseFilt(filt)])
	df = pd.DataFrame(data, columns = columns)
	return df

### Takes a DataFrame in simCluster format and turns it into a cluster structure
### Puts the result in the cluster pointed to by pCluster
def dataFrameToCluster(pCluster: clusterPtr, data: pdDataFrame) -> None:
	c.initCluster(pCluster)
	pCluster.contents.nStars = len(data.index)
	c.reallocStarsArray(pCluster)
	data= data.reset_index()
	for j, row in data.iterrows():
		star = c.getStarPtr(pCluster,j)
		star.contents.id = int(row["id"])

		# star 1
		c.setMass1(star,pCluster,row["mass1"])
		star.contents.stage1 = int(row["stage1"])
		star.contents.massNow[0] = row["wdM1"]
		star.contents.wdType[0] = int(row["wdType1"])	
		star.contents.wdLogTeff[0] = row["wdLogTeff1"]

		# star 2
		c.setMass2(star,pCluster,row["mass2"])
		star.contents.stage1 = int(row["stage2"])
		star.contents.massNow[1] = row["wdM2"]
		star.contents.wdType[1] = int(row["wdType2"])	
		star.contents.wdLogTeff[1] = row["wdLogTeff2"]

		#photometry
		for filt in range(FILTS):
			if c.getUseFilt(filt):
				star.contents.photometry[filt] = row[c.getFilterName(filt).decode()]


### Drops every column in a cluster DataFrame that isn't a photometry column
def keepPhotometry(data: pdDataFrame) -> None:
	data.drop(columns = ["id","mass1","stage1", "wdM1", "wdType1", "wdLogTeff1", "ltau1","mass2","stage2", "wdM2", "wdType2", "wdLogTeff2", "ltau2"],inplace = True)


# The signal 2 noise coefficients for each band of photometry
# Used by scatterPhot to add noise to simulated cluster star photometry
# Originally defined in scatterCuster.c
s2nCoeffs = [
	[9.33989,  0.3375778],  # U 
	[10.0478,  0.3462758],  # B 
	[10.48098, 0.368201 ],  # V 
	[10.71151, 0.3837847],  # R 
	[10.61035, 0.3930941],  # I 
	[9.282385, 0.386258 ],  # J 
	[9.197463, 0.3970419],  # H 
	[9.024068, 0.3985604],  # K
	[9.024068, 0.3985604],  # IRAC Blue
	[9.024068, 0.3985604],  # IRAC Red
	[9.024068, 0.3985604],  # Band 1
	[9.024068, 0.3985604],  # Band 2
	[9.024068, 0.3985604],  # Band 3
	[9.024068, 0.3985604]   # Band 4
]

### Scatter photometry and return
def scatterPhot(simdf: pdDataFrame, params: Params) -> pdDataFrame:
	#For each filter that we are outputting
	scatterdf = simdf.copy()
	for filt in range(FILTS):
			if c.getUseFilt(filt):
				if params.exposures[filt] > 0.0:
					#Check if the star's photometry is above our desired limiting S/N
					logS2N  = s2nCoeffs[filt][0] - s2nCoeffs[filt][1] * scatterdf[c.getFilterName(filt).decode()]
					s2n = 10. ** logS2N
					s2n *= np.sqrt(params.exposures[filt])
					sigma = 1./(s2n)
					sigma[sigma < 0.005] = 0.005
					scatterdf["sig" + c.getFilterName(filt).decode()] = sigma
					scatterdf[c.getFilterName(filt).decode()] += np.random.normal(0.,sigma)
	return scatterdf

### Get rid of stars that don't meet the criteria we want for scatterCluster
def scatterClean(inputdf: pdDataFrame, params: Params) -> pdDataFrame:
	filter_col = [col for col in inputdf if col.startswith('sig')]
	s2n = 1/inputdf[filter_col]
	s2n["keep"] = s2n.apply(lambda row: all([(x > params.limitS2N) for x in row]), axis = 1)
	outputdf = inputdf[s2n["keep"]]
	brightLimit = params.brightLimit[0]
	faintLimit = params.brightLimit[1]
	filt = c.getFilterName(params.brightLimit[2]).decode()
	# print(params.limitS2N,params.brightLimit[0],params.brightLimit[1],params.brightLimit[2])
	return outputdf.loc[((outputdf[filt]>brightLimit) & (outputdf[filt]<faintLimit))| (outputdf["stage1"] == WD) | (outputdf["pop"] == 9)]

	# with pd.option_context('display.max_rows', None,):
	# 	print(outputdf.loc[((outputdf[filt]<brightLimit) &  (outputdf[filt]<faintLimit))| (outputdf["stage1"] == WD) | (outputdf["pop"] == 9)])
	# 