import ctypes
from parameters import Params
import numpy as np
import pandas as pd
import structures as st
import math
import sys
import mcmc_plot_functions as mc

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

def pow(x,y):
	return x ** y

def getSimHeader():
	header = ["pop", "id", "mass1"]
	for filt in range(st.FILTS):
		if st.c.getUseFilt(filt):
			header.append(f"{st.c.getFilterName(filt).decode()}1")
	header.extend(["stage1", "wdM1", "wdType1", "wdLogTeff1", "ltau1", "mass2"])
	for filt in range(st.FILTS):
		if st.c.getUseFilt(filt):
			header.append(f"{st.c.getFilterName(filt).decode()}2")
	header.extend(["stage2", "wdM2", "wdType2", "wdLogTeff2", "ltau2"])
	for filt in range(st.FILTS):
		if st.c.getUseFilt(filt):
			header.append(f"{st.c.getFilterName(filt).decode()}")
	return header

def getScatterHeader():
	header = ["id"]
	for filt in range(st.FILTS):
		if st.c.getUseFilt(filt):
			if params.exposures[filt] > 0.0:
				header.append(f"{st.c.getFilterName(filt).decode()}")
	for filt in range(st.FILTS):
		if st.c.getUseFilt(filt):
			if params.exposures[filt] > 0.0:
				header.append(f"sig{st.c.getFilterName(filt).decode()}")
	header.extend(["mass1", "massRatio", "stage1", "CMprior", "useDuringBurnIn", "pop"])
	return header

params = Params()
np.random.seed(params.seed)

#### Create the cluster ####
pCluster = st.createCluster(params)

#### Sets global variables in the simCluster C module
fault = st.c.setGlobals(pCluster,
						params.seed, 
						params.verbose)
if fault: sys.exit()
minMass = pCluster.contents.evoModels.minMass



with open(params.outfile,"w") as f:

	print(f"\nFilename: {params.outfile}")

	##### Create Stars #####
	for i in range(params.nPops):

		if(i >= st.MAXPOPS - 2):
			print(f"Maximum number of populations is {st.MAXPOPS - 2}.  Skipping to field stars")
			break


		# pCluster.contents.nStars = params.nSystems[i]
		# st.c.reallocStarsArray(pCluster)
		# pCluster.contents.M_wd_up = params.WDMassUp[i]
		# pCluster.contents.parameter[st.AGE] = params.logClusAge[i]
		# pCluster.contents.parameter[st.YYY] = params.Y[i]
		# pCluster.contents.parameter[st.FEH] = params.Fe_H[i]
		# pCluster.contents.parameter[st.MOD] = params.distMod[i]
		# pCluster.contents.parameter[st.ABS] = params.Av[i]

		st.setClusterParams(pCluster,params,i)
		# input as percentages, use as fractions 
		fractionBinary = params.percentBinary[i] / 100.
		fractionDB = params.percentDB[i] / 100.

		nStars = 0
		massTotal = 0.0 
		minV = 1000.0
		maxV = -1000.0
	
		################################
		##### Create cluster stars #####
		################################

		## Create the stars in the cluster, drawing masses from IMF
		## and mass ratios from a random distribution
		st.populateCluster(pCluster,fractionBinary)

		## Transfer the relevant values for all the stars into a dataframe
		## So we can output the final photometry at the end
		finaldf = st.clusterToDataFrame(pCluster)
		columns = [st.c.getFilterName(filt).decode() for filt in range(st.FILTS) if st.c.getUseFilt(filt)]
		finaldf[finaldf[columns] > 99] = 99.999

		## set mass ratios to zero
		for j in range(pCluster.contents.nStars):
			st.c.getStarPtr(pCluster,j).contents.massRatio = 0.0

		## evolve stars without secondaries and store results
		st.c.evolve(pCluster, -1)
		result = st.clusterToDataFrame(pCluster)
		newcolumns = [(column + "1") for column in columns]
		result.rename(columns = dict(zip(columns,newcolumns)), inplace = True)
		st.keepPhotometry(result)
		result[result > 99] = 99.999
		finaldf = finaldf.join(result)

		## set primary mass to current primary * mass ratio
		massRatios = finaldf['mass2']/finaldf['mass1']
		for j in range(pCluster.contents.nStars):
			st.c.getStarPtr(pCluster,j).contents.U = st.c.getStarPtr(pCluster,j).contents.U * massRatios[j]

		## evolve stars without primaries and store results
		st.c.evolve(pCluster, -1)
		result = st.clusterToDataFrame(pCluster)
		newcolumns = [(column + "2") for column in columns]
		result.rename(columns = dict(zip(columns,newcolumns)), inplace = True)
		st.keepPhotometry(result)
		result[result > 99] = 99.999
		finaldf = finaldf.join(result)

		#rearrange results and add pop column
		header = getSimHeader()
		finaldf.insert(0,column = "pop",value=i)

		minV = min(finaldf.loc[finaldf["V"] < 97,"V"])
		maxV = max(finaldf.loc[finaldf["V"] < 97,"V"])

		theCluster = pCluster.contents
		print(f"Properties for population {i+1}:")
		print(f"logClusAge     = {theCluster.parameter[st.AGE]:5.2f}")
		print(f"[Fe/H]         = {theCluster.parameter[st.FEH]:5.2f}")
		print(f"Y              = {theCluster.parameter[st.YYY]:5.2f}")
		print(f"modulus        = {theCluster.parameter[st.MOD]:5.2f}")
		print(f"Av             = {theCluster.parameter[st.ABS]:5.2f}")
		print(f"WDMassUp       = {theCluster.M_wd_up:4.1f}")
		print(f"fractionBinary = {fractionBinary:5.2f}")
			
		print("Totals:")
		print(f"nSystems       = {theCluster.nStars}")

		nStars = len(finaldf.loc[finaldf["stage2"] != st.DNE,"stage2"]) + len(finaldf["mass1"])
		print(f"nStars         = {nStars}")
		nMSRG  =  len(finaldf.loc[finaldf["stage2"] == st.MSRG,"stage2"]) + len(finaldf.loc[finaldf["stage1"] == 1,"stage1"])
		print(f"nMSRG          = {nMSRG}")
		nWD  =  len(finaldf.loc[finaldf["stage2"] == st.WD,"stage2"]) + len(finaldf.loc[finaldf["stage1"] == st.WD,"stage1"])
		print(f"nWD            = {nWD}")
		nNSBH = len(finaldf.loc[finaldf["stage2"] == st.NSBH,"stage2"]) + len(finaldf.loc[finaldf["stage1"] == st.NSBH,"stage1"])
		print(f"nNSBH          = {nNSBH}")
		massTotal  = sum(finaldf["mass1"]) + sum(finaldf["mass2"])
		print(f"massTotal      = {massTotal:5.2f}")
		MSRGMassTotal = sum(finaldf.loc[finaldf["stage1"] == st.MSRG,"mass1"]) + \
						sum(finaldf.loc[finaldf["stage2"] == st.MSRG,"mass2"])
		print(f"MSRGMassTotal  = {MSRGMassTotal:5.2f}")
		wdMassTotal = sum(finaldf.loc[finaldf["stage1"] == st.WD,"wdM1"]) + \
						sum(finaldf.loc[finaldf["stage2"] == st.WD,"wdM2"])
		print(f"wdMassTotal    = {wdMassTotal:5.2f}")
		print("")

	###########################
	### create field stars ####
	###########################

	minFeH = -2.5
	maxFeH = 0.56
	minAge = 8.4
	maxAge = 10.17

	print(f"nFieldStars: {params.nFieldStars}")
	pCluster.contents.nStars = 1
	tempMod = pCluster.contents.parameter[st.MOD]
	tempAbs = pCluster.contents.parameter[st.ABS]
	if(pCluster.contents.evoModels.mainSequenceEvol == st.YALE):
		pCluster.contents.evoModels.mainSequenceEvol = st.DSED
		st.c.loadModels(pCluster);
		minMass = pCluster.contents.evoModels.minMass
	st.c.reallocStarsArray(pCluster)
	starContents = st.c.getStarPtr(pCluster,0).contents
	starContents.status[0] = st.MSRG

	fieldStars = []
	for j in range(params.nFieldStars):
		while True:     
			# Draw a new age and metallicity
			pCluster.contents.parameter[st.AGE] = minAge + (maxAge-minAge) * np.random.random()
			pCluster.contents.parameter[st.FEH] = minFeH + (maxFeH-minFeH) * np.random.random()
			
			# Determine a new distance, weighted so 
			# there are more stars behind than in front
			pCluster.contents.parameter[st.MOD] = tempMod - 12.0 + math.log10(pow(10,(pow(pow(26.0,3.0) * np.random.random(),1.0/3.0))))
			
			st.populateCluster(pCluster, 0.5)
						
			if (starContents.photometry[2] > minV and 
				starContents.photometry[2] < maxV and
				(starContents.photometry[1] 
					- starContents.photometry[2]) > -0.5 and
				(starContents.photometry[1] 
					- starContents.photometry[2]) < 1.7):
				break
		starContents.id = 9000 + j

		if j == 0:
			fieldstardf = st.clusterToDataFrame(pCluster)
		fieldStars.extend(st.clusterToList(pCluster))

	fieldstardf = pd.DataFrame(fieldStars,columns = fieldstardf.columns)
	fieldstardf[fieldstardf[columns] > 99] = 99.999
	fieldstardf.insert(0,column = "pop",value=9)
	fieldstardf = fieldstardf.reindex(columns = header, fill_value=99.999)
	finaldf = pd.concat([finaldf,fieldstardf])

	finaldf = finaldf.reindex(columns = header)
	finaldf.reset_index(inplace=True,drop=True)
	finaldf.to_csv(f,sep=" ",index=False)

	#st.dataFrameToCluster(pCluster,finaldf)

# for filt in range(st.FILTS):
# 	for i in range(2):
# 		print(s2nCoeffs[filt][i])

###################################
######### Scatter Cluster #########
###################################

def scatterPhot(scatterdf = "scatterdf"):
	#For each filter that we are outputting
	for filt in range(st.FILTS):
			if st.c.getUseFilt(filt):
				if params.exposures[filt] > 0.0:
					#Check if the star's photometry is above our desired limiting S/N
					logS2N  = s2nCoeffs[filt][0] - s2nCoeffs[filt][1] * scatterdf[st.c.getFilterName(filt).decode()]
					s2n = 10. ** logS2N
					s2n *= math.sqrt(params.exposures[filt])
					scatterdf = scatterdf.loc[s2n > params.limitS2N]
					s2n = s2n.loc[s2n > params.limitS2N]
					sigma = 1./(s2n)
					sigma[sigma < 0.005] = 0.005
					scatterdf["sig" + st.c.getFilterName(filt).decode()] = sigma
					scatterdf[st.c.getFilterName(filt).decode()] += np.random.normal(0.,sigma) 
	return scatterdf


assert len(params.exposures) == st.FILTS, f"Number of exposure times {len(params.exposures)} not equal to number of filts {st.FILTS}"
limitSigToNoise = params.limitS2N
seed = params.seed
nPops = params.nPops
nStars = [int(params.nSystems[i]*params.keepRatio) for i in range(nPops)]
nStars.extend([int(params.nBrownDwarfs*params.keepRatio),int(params.nFieldStars*params.keepRatio)])
brightLimit, faintLimit, firstFilt = params.brightLimit

clusterMemberPrior  = sum(nStars[:-2]) / sum(nStars)
if clusterMemberPrior > 0.99: clusterMemberPrior = 0.99
if clusterMemberPrior < 0.01: clusterMemberPrior = 0.01

header = getScatterHeader()

with open(params.scatfile,"w") as f:
	#f.write(str(header))
	######  delete stars that don't meet the various criteria to be included ######
	### Main sequence stars that are above or below the brightness limits
	scatterdf = finaldf.loc[(finaldf["pop"] == 9) |
					((finaldf[st.c.getFilterName(firstFilt).decode()] < faintLimit) & 
					 (finaldf["stage1"] != st.WD)) &
					(finaldf[st.c.getFilterName(firstFilt).decode()] > brightLimit) | 
					(finaldf["stage1"] == st.WD)]

	### Remove neutron stars and black holes
	scatterdf = scatterdf.loc[~(scatterdf["stage1"] == st.NSBH) | (scatterdf["stage2"] == st.NSBH)]

	### Remove WD binaries (except for field stars)
	scatterdf = scatterdf.loc[~((scatterdf["stage1"] == st.WD) & 
								(scatterdf["mass2"] > 0.0) & 
								(scatterdf["pop"] != 9))]

	scatterdf = scatterPhot(scatterdf)
	scatterdf["massRatio"] = scatterdf["mass2"]/scatterdf["mass1"]
	scatterdf["CMprior"] = clusterMemberPrior
	scatterdf["useDuringBurnIn"] = 1
	scatterdf.loc[scatterdf["pop"] == 9,"useDuringBurnIn"] = 0
	scatterdf = scatterdf.reindex(columns=header)

	mc.HR(scatterdf,ebars=True)
	mc.HR(finaldf,overplot=True,color="red",type="l")
	mc.HR(scatterdf[scatterdf["useDuringBurnIn"]==0],overplot=True,color="blue")

	scatterdf.to_csv(f,sep=" ",index=False)


