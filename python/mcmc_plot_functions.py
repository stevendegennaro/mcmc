import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import TypeVar,Tuple
import math
PandasDataFrame = TypeVar('pd.core.frame.DataFrame')


#def set_pathname(path: str = "outputs/")->None:
#	global pathname
#	pathname = path


def import_from_sim_cluster(filename: str, pathname: str = "../outputs/") -> PandasDataFrame:
	fullfilename = pathname+filename

	with open(fullfilename,"r") as f:
		print(f"file {fullfilename} opened")
		data = pd.read_csv(f,delim_whitespace=True)
		data.set_index("id",inplace=True)

	return data

def HR(data:PandasDataFrame,
		band1: str="B",
		band2: str="V",
		marker: str=".",
		color="black",
		size=20,
		type: str="p",
		overplot: bool = False,
		ebars: bool = True):
	"""
	Takes a data frame that contains multiband stellar photometry 
	and plots a Hertzsprung-Russel diagram (aka color-magnitude diagram)
	with or without error bars.

	Returns the points and error bars as a tuple of points objects (so
	they can be hidden or removed in the parent program as needed)
	"""

	data=data[data[band2]<90]
	points = []

	if not overplot:
		plt.ion()
		plt.figure()
		plt.show()
		plt.gca().invert_yaxis()
		plt.gca().set_ylabel(band2)
		plt.gca().set_xlabel(band1 + " - " + band2)

	if type == "l":
		if 'stage2' in data.columns:
			data = data[data['stage2'] == 9]
		pops = data.value_counts(["pop"])
		groups = data.groupby("pop")
		popnames = [pop[0] for pop in pops.index if pop[0] != 9]
		for popname in popnames:
			group = groups.get_group(popname)
			group = group.sort_values(["mass1"])
			WDs = group[group["stage1"]==3]
			MS = group[group["stage1"]==1]
			points.append(plt.plot(WDs[band1]-WDs[band2],WDs[band2],c=color))
			points.append(plt.plot(MS[band1]-MS[band2],MS[band2],c=color))
	
	if type == "p":
		point = plt.scatter(data[band1]-data[band2],data[band2],marker=marker,c=color,s=size)
		points.append(point)
		if ebars:
			#first check to see if you have error bars in the file
			sig1 = "sig" + band1
			sig2 = "sig" + band2
			if sig1 in data.columns and sig2 in data.columns:
				sig_color = np.sqrt(data[sig1]**2 + data[sig2]**2)
				errors = plt.errorbar(data[band1]-data[band2],
										data[band2],
										yerr = data[sig2],
										xerr=sig_color,
										c=color,
										ls="none")
				points.append(errors)

	return points


def getHRlims(data:PandasDataFrame, band1: str, band2:str, ebars:bool = False) -> Tuple[Tuple]:
	miny = data.max(axis=0)[band2]
	maxy = data.min(axis=0)[band2]

	if ebars:
		#first check to see if you have error bars in the file
		sig1 = "sig" + band1
		sig2 = "sig" + band2
		if sig1 in data.columns and sig2 in data.columns:
			sig_color = np.sqrt(data[sig1]**2 + data[sig2]**2)
			minx = min(data[band1]-data[band2]-sig_color)
			maxx = max(data[band1]-data[band2]+sig_color)
		else:
			minx = min(data[band1]-data[band2])
			maxx = max(data[band1]-data[band2])
	else:
		minx = min(data[band1]-data[band2])
		maxx = max(data[band1]-data[band2])

	return (minx,maxx),(miny,maxy)

def magmag(data:PandasDataFrame,
		band1: str="B",
		band2: str="V",
		marker: str=".",
		color="black",
		size=20,
		type: str="p",
		overplot: bool = False,
		ebars: bool = True):
	"""
	Takes a data frame that contains multiband stellar photometry 
	and plots a Hertzsprung-Russel diagram (aka color-magnitude diagram)
	with or without error bars.

	Returns the points and error bars as a tuple of points objects (so
	they can be hidden or removed in the parent program as needed)
	"""

	data=data[data[band2]<90]
	points = []

	if not overplot:
		plt.ion()
		plt.figure()
		plt.gca().set_ylabel(band1)
		plt.gca().set_xlabel(band2)
		plt.gca().invert_yaxis()
		plt.gca().invert_xaxis()
		plt.show()


	if type == "l":
		if 'stage2' in data.columns:
			data = data[data['stage2'] == 9]
		pops = data.value_counts(["pop"])
		groups = data.groupby("pop")
		popnames = [pop[0] for pop in pops.index if pop[0] != 9]
		for popname in popnames:
			group = groups.get_group(popname)
			group = group.sort_values(["mass1"])
			WDs = group[group["stage1"]==3]
			MS = group[group["stage1"]==1]
			points.append(plt.plot(WDs[band1],WDs[band2],c=color))
			points.append(plt.plot(MS[band1],MS[band2],c=color))
	
	if type == "p":
		point = plt.scatter(data[band2],data[band1],marker=marker,c=color,s=size)
		points.append(point)
		if ebars:
			#first check to see if you have error bars in the file
			sig1 = "sig" + band1
			sig2 = "sig" + band2
			if sig1 in data.columns and sig2 in data.columns:
				errors = plt.errorbar(data[band1],
										data[band2],
										yerr = data[sig2],
										xerr = data[sig1],
										c=color,
										ls="none")
				points.append(errors)

	return points

