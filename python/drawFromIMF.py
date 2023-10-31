import numpy as np
import ctypes
import matplotlib.pyplot as plt

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
	mf_mu = -1.02

	logMass = np.random.normal(mf_mu, mf_sigma)
	while (logMass < -4) or (logMass > 2.0):    # keep within mass (0.15 + EPS to 100 Mo) limits
		logMass = np.random.normal(mf_mu, mf_sigma)

	zamsMass = 10.0 ** logMass
	return zamsMass

def plotIMF(log=False):
	
	def gaussian(x, mu, sig):
		return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)
	
	mf_sigma = 0.67729
	mf_mu = -1.02

	logmasses = np.arange(-4,2,.01)
	masses = np.power(10,logmasses)
	dens = gaussian(logmasses,mf_mu,mf_sigma) * 1000
	star_masses = np.array([drawFromIMF() for _ in range(10000)])

	if log:
		plt.plot(logmasses,dens,color="red")
		plt.hist(np.log10(star_masses),50,color="gray", edgecolor = "Black")
		plt.gca().set_xlabel("Log10(Mass)")
	else:
		# plt.plot(masses,dens,color="red")
		plt.hist(star_masses[star_masses<10],100,color="gray", edgecolor = "Black")
		plt.gca().set_xlabel("Mass")

	plt.gca().set_ylabel("N Stars")
	plt.title("Initial Mass Function")
	plt.show()


# dfIMF = ctypes.CDLL("../c/drawFromIMF.so")

# x =[]
# y=[]
# for _ in range(1000000):
# 	x.append(drawFromIMF())
# 	y.append(dfIMF.drawFromIMF())

# import matplotlib.pyplot as plt

# # plt.hist(x,100,log=True)
# # plt.show()


