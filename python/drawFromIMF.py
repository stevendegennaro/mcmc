import numpy as np
import ctypes

#*****************************************************************************************
# last update: 28sep05

# Draw a star randomly from the Miller-Scalo (1979, ApJS, 41, 513) IMF, using their gaussian 
# in log(M) form for age = 12 Gyr (c1 = 1.09, c2 = -1.02). Using a mass range of 0.1 to 100 
# Mo.  (Resetting and recalculating the constants each time this subroutine is called is 
# inefficient, but more calcuation time will be spent drawing the Gaussian random deviate 
# (typically done 2+ times per function call) and calculating the final inverse log.)
#*****************************************************************************************/
def drawFromIMF() -> float:
	mf_sigma = 0.67729
	mf_mu = -1.02;

	logMass = np.random.normal(mf_mu, mf_sigma)
	while (logMass < -4) or (logMass > 2.0):    # keep within mass (0.15 + EPS to 100 Mo) limits
		logMass = np.random.normal(mf_mu, mf_sigma)

	zamsMass = 10.0 ** logMass
	return zamsMass


# dfIMF = ctypes.CDLL("../c/drawFromIMF.so")

# x =[]
# y=[]
# for _ in range(1000000):
# 	x.append(drawFromIMF())
# 	y.append(dfIMF.drawFromIMF())

# import matplotlib.pyplot as plt

# # plt.hist(x,100,log=True)
# # plt.show()


