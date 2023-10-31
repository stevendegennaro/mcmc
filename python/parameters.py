import numpy as np
from typing import List

class Params:
	def __init__(self):
		self.seed = 1
		self.verbose = 0
		self.nPops = 1
		self.nSystems = [200]		#create excess, as some will be RGB, HB, & AGB stars
		self.WDMassUp = [8.0]	
		self.percentBinary = [50.]
		self.percentDB = [0.]
		self.logClusAge = [9.1]
		self.Fe_H = [0.1]
		self.distMod = [10.0]
		self.Av = [0.1]
		self.Y = [0.29]

		#model set parameters
		self.msRgbModels = 0	# 0 = Girardi, 1 = Chaboyer-Dotter w/helium sampling, 2 = Yale-Yonsei, 3 = Dartmouth (DSED)
		self.filterSet = 0		# 0 = UBVRIJHK, 1 = ACS, 2 = SDSS + JHK
		self.IFMR =  0			# 0 = Weidemann, 1 = Williams, 2 = Salaris Linear, 3 = Salaris piecewise linear
		self.BDmodel = 0		# 0 = None, 1 = Baraffe

		#scatterCluster parameters
		self.brightLimit =  [13,22,2]		#brightLimit, faintLimit, firstFilt; clip RGB, faint WD/MS respectively; 0=bluest band
		self.limitS2N =  10
		self.exposures =  [.03,.02,.02,.02,0.1,0,0,0,0,0,0,0,0,0]
		self.nFieldStars = 50
		self.nBrownDwarfs = 0

		self.suffix = 'fe' + str(self.Fe_H[0])\
							+ '.a' + "{:.2f}".format(self.logClusAge[0]) \
							+ '.u' + str(self.WDMassUp[0])	\
							+ '.n' + str(self.nSystems[0]) \
							+ '.fs' + str(self.nFieldStars)	#filename setting more compact

		self.path = "../outputs/"
		self.updateFileStem()
		self.modelPath = "../models"

	def addPop(self,popParams: List):
		self.nPops += 1
		self.nSystems.append(popParams[0])
		self.WDMassUp.append(popParams[1])	
		self.percentBinary.append(popParams[2])
		self.percentDB.append(popParams[3])
		self.distMod.append(popParams[4])
		self.Av.append(popParams[5])
		self.logClusAge.append(popParams[6])
		self.Fe_H.append(popParams[7])
		self.Y.append(popParams[8])

	def setParams(self,popParams: List):
		print(popParams)
		self.nPops = 1
		self.nSystems = [popParams[0]]
		self.WDMassUp = [popParams[1]]	
		self.percentBinary = [popParams[2]]
		self.percentDB = [popParams[3]]
		self.distMod = [popParams[4]]
		self.Av = [popParams[5]]
		self.logClusAge = [popParams[6]]
		self.Fe_H = [popParams[7]]
		self.Y = [popParams[8]]

	def updateFileStem(self):
		self.suffix = 'fe' + "{:.2f}".format(self.Fe_H[0])\
					+ '.a' + "{:.2f}".format(self.logClusAge[0]) \
					+ '.u' + "{:.1f}".format(self.WDMassUp[0])	\
					+ '.n' + str(self.nSystems[0]) \
					+ '.fs' + str(self.nFieldStars)	#filename setting more compact

		self.file_stem  = self.path + 'sim.' + self.suffix
