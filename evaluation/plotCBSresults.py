#!/usr/bin/env python

from __future__ import division
import sys
sys.path.append("..")
import os
import matplotlib
matplotlib.use('Agg')

import numpy as np
import subprocess
import multiprocessing as mp
import matplotlib.pyplot as plt
from scipy.stats import mode
import itertools
from hammlet_plotting import *
import scipy




def plotCBSresults(prefix, multiFile="", paletteFile=None, blocksFile="", ylim=None):
	filename=prefix+".csv"
	data = np.loadtxt(filename, skiprows=0, usecols=[1])
	means = np.loadtxt(prefix+"_MergeLevels_segmentation.csv", skiprows=1, usecols=[1], dtype=float)

	mergeLevelsStates = scipy.stats.rankdata(means)-1	
	mergeLevelsStates = -mergeLevelsStates.astype(int)-1
	# mimic the method='dense' keyword for older versions of scipy:
	uniqueStates = np.unique(mergeLevelsStates)
	nrStates = len(uniqueStates)
	print "Found %d states" % nrStates
	f=file("%s_nrStates.txt" % prefix, 'w')
	f.write(str(nrStates))
	f.close()
	for i in xrange(nrStates):
		mergeLevelsStates[mergeLevelsStates==uniqueStates[i]] = i

	
	calldir = os.path.dirname(os.path.realpath(__file__))
	if paletteFile==None:
		paletteFile=os.path.join(calldir, "..", "palette.txt")
		
	palette = np.loadtxt(paletteFile, comments="/", dtype=str)[::-1]
	
	
	
	if multiFile != "":
		multi=np.array([0]+list(np.loadtxt("all_bt474_multi.txt", usecols=[0], dtype=int))).cumsum()
		chrom=np.loadtxt("all_bt474_multi.txt", usecols=[1], dtype=str)

		for i in xrange(len(multi)-1):
			f = plt.figure(figsize=[8.0, 3.0])
			figname = "%s_CBS_%s_%dstates.png" % (prefix, chrom[i], nrStates)
			plotData(data, states=mergeLevelsStates, palette=palette, xlabel="Probe number", ylabel="Log2 ratio", zorder=2, start=multi[i], end=multi[i+1], ylim=ylim)
			print figname
			plt.savefig(figname, bbox_inches="tight")
			plt.close()
			
	else:
		f = plt.figure(figsize=[8.0, 3.0])
		figname = "%s_CBS_%dstates.png" % (prefix, nrStates)
		plotData(data, states=mergeLevelsStates, palette=palette, xlabel="Probe number", ylabel="Log2 ratio", zorder=2, ylim=ylim)
		if blocksFile != "":
			blockSizes = np.loadtxt(blocksFile)
			plotCoarsestSegmentation(blocks = blockSizes, zorder=1)
		print figname
		plt.savefig(figname, bbox_inches="tight")
		plt.close()

