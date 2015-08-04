#!/usr/bin/env python

from __future__ import division
import matplotlib
matplotlib.use('Agg')
import os
import sys

sys.path.append("..")
sys.path.append("../..")
hammlet = "python ../../../hammlet.py"	# this time, we go up one more level because files are in subdirectories

from Fmeasures import *
import numpy as np
import subprocess
import multiprocessing as mp
import matplotlib.pyplot as plt
from scipy.stats import mode
import itertools
from hammlet_plotting import *
import scipy

from plotCBSresults import *


nrThreads = 30
iterations=1000
burnin=900
states=5

hammlet_palette = np.loadtxt("../palette.txt", comments="/", dtype=str)



def bestFmeasure(data, inferredStates, trueStateSequence):
	"""trueStateSequence is assumed to be 0-1 vector, inferredStates ranked internally (it can be float or anythin like this, such as the means we get from MergeLevels)."""
	
	dist = 99999999999999999
	neutralInferredState = 0
	for s in np.unique(inferredStates):
		newdist = np.abs(np.mean(data[inferredStates==s]))
		if newdist < dist:
			neutralInferredState = s
			dist = newdist
	
	
	bestF = -1
	inferredStateSequence = np.zeros(len(trueStateSequence), dtype=int)
	sameAsNeutral = True	# check whether the best F-score comes from the state that is the closest to zero
	for s in np.unique(inferredStates):	# iterate over each potential neutral state
	
		inferredStateSequence[inferredStates!=s] = 1
		inferredStateSequence[inferredStates==s] = 0

		TP= np.logical_and(trueStateSequence, inferredStateSequence).sum()
		FP = np.logical_and(np.logical_not(trueStateSequence), inferredStateSequence).sum()
		FN = np.logical_and(trueStateSequence, np.logical_not(inferredStateSequence)).sum()
		
		if TP+FP==0:	# edge cases
			precision = 1	
		else:
			precision = TP/(TP+FP)
		if (TP+FN)==0:	# edge cases
			recall = 1	
		else:
			recall = TP/(TP+FN)
		
		currentF = 2* (precision * recall)/(precision+recall)
		if currentF > bestF:
			bestF = currentF
			sameAsNeutral = (s==neutralInferredState)
	
	return (bestF, sameAsNeutral)

	
	


def Fmeasures(prefices, outFilePrefix, burnin=0):
	HaMMLET_result = "# prefix\tF1\tneutral\n"
	CBS_result = "# prefix\tF1\tneutral\n"
	HaMMLET_F = np.zeros(len(prefices))
	CBS_F = np.zeros(len(prefices))
	fuzzyF = np.zeros(len(prefices))
	for i in xrange(len(prefices)):
		prefix = prefices[i]
		#print prefix
		
		dataFilename = prefix+".csv"
		data = np.loadtxt(dataFilename, ndmin=2)
		trueStateSequence = data[:,0].astype(int)
		emissions = data[:,1]
		# label neutral states with 0 and aberrant states with 1
		trueStateSequence[trueStateSequence!=2] = 0
		trueStateSequence[trueStateSequence==2] = 1
		trueStateSequence = 1-trueStateSequence


		### HaMMLET ###
		sampleFilename = prefix+"_statesequence.csv"
		
		samples = np.loadtxt(sampleFilename, skiprows=1+burnin)		
		inferredStates = mode(samples, axis=0)[0][0].astype(int)
		
		F, sameAsNeutral = bestFmeasure(emissions, inferredStates, trueStateSequence)
		HaMMLET_F[i] = F
		HaMMLET_result += "%s\t%f\t%s\n" % (prefix, F, sameAsNeutral)
		
		
		### CBS (= DNAcopy+MergeLevels) ###
		sampleFilename = prefix+"_MergeLevels_segmentation.csv"
		inferredStates = np.loadtxt(sampleFilename, skiprows=1, usecols=[1])	
		
		F, sameAsNeutral = bestFmeasure(emissions, inferredStates, trueStateSequence)
		CBS_F[i] = F
		CBS_result += "%s\t%f\t%s\n" % (prefix, F, sameAsNeutral)
		
		
	f = file(outFilePrefix+"_Fmeasures_HaMMLET.csv", 'w')
	f.write(HaMMLET_result)
	f.close()
	
	f = file(outFilePrefix+"_Fmeasures_CBS.csv", 'w')
	f.write(CBS_result)
	f.close()
		



def runCalls(prefix):
	try:
		print "\nStarting thread for %s" % prefix
		filename=prefix+".csv"
		blockFileName = prefix+"_blocks.csv"
		blocksFile = file(blockFileName, "w")
		count=1
		
		blocks = np.loadtxt(filename, usecols=[0]).astype(int)
		blockSizes=[]
		for i in xrange(1, len(blocks)):
			if blocks[i] != blocks[i-1]:
				blocksFile.write("%d\n" % count)
				blockSizes.append(count)
				count = 1
			else:
				count += 1
		blocksFile.write("%d\n" % count)
		blockSizes.append(count)
		blockSizes = np.asarray(blockSizes)
		blocksFile.close()
		
		
		data = np.loadtxt(prefix+".csv")[:,1]
		ylim = list(np.array([data.min(), data.max()])*1.5)
		
		# HaMMLET	
		os.system("%s --sequences   --iterations %d --states %d ../%dstates.txt '%s'" % (hammlet, iterations, states, states, filename))
		
		
		saveResultsFromFile(outFileName=prefix+"_HaMMLET.png", dataFiles=filename, sequenceFile=prefix+"_statesequence.csv",  nrStates=5, blocksFile=blockFileName, maxBlockSize=0, palette=hammlet_palette, burnin=burnin, allSequences=True, xlabel="Probe number", ylabel="Log2 ratio", ylim=ylim)


		
		### DNAcopy/CBS

		os.system("/usr/bin/time -f \"command\t%%C\nexit status\t%%x\nmax resident size [kB]\t%%M\navg resident size [kB]\t%%t\nuser CPU time [s]\t%%U\nprocess CPU time [s]\t%%S\nwall clock time [s]\t%%e\nmajor page faults\t%%F\nminor page faults\t%%R\n\" -o DNAcopy_benchmark.txt  Rscript %s %s" % ("../../DNAcopy_MergeLevels.R", prefix))		# "test.r"

		
		plotCBSresults(prefix, multiFile="",  blocksFile=prefix+"_blocks.csv", ylim=ylim)
		
	except:
		print "[ERROR] Failed on", prefix
		return 0
	return 0
	



root = os.path.abspath(os.getcwd())


prefix = "breakpoint_detection_and_merging"	
lowest=1	# 1
highest=500	#500
os.chdir(os.path.join(root, prefix))
prefices = ["sample%d" % s for s in xrange(lowest, highest+1)]
filenames = [p+".csv" for p in prefices]
pool = mp.Pool(nrThreads)	
pool.map(runCalls, prefices)	
Fmeasures(prefices, prefix, burnin)



prefix = "spatial_resolution_study"	
lowest = 1	# 1
highest = 100	# 100
os.chdir(os.path.join(root, prefix))
for testcase in ["below5", "5to10", "10to20", "above20"]:
	prefices = ["%s_sample_%d" % (testcase, x) for x in xrange(lowest, highest+1)]
	filenames = [p+".csv" for p in prefices]
	pool = mp.Pool(nrThreads)	
	pool.map(runCalls, prefices)
	Fmeasures(prefices, prefix+"_"+testcase, burnin)
	
prefix = "testing"	
dslowest = 1	# 1
dshighest = 500	# 500
slowest = 1	# 1
shighest = 20	# 20
os.chdir(os.path.join(root, prefix))
prefices = ["dataset_%d_sample_%d" % x for x in itertools.product(xrange(dslowest, dshighest+1), xrange(slowest, shighest+1))]
filenames = [p+".csv" for p in prefices]
pool = mp.Pool(nrThreads)	
pool.map(runCalls, prefices)

Fmeasures(prefices, prefix, burnin)	# calculate Fmeasures for all software



### create boxplot
os.chdir(root)
print "Creating boxplot..."
i=0 

def boxplot(filename):	# file contains a list of F-measures
	global i
	F = np.loadtxt(filename, usecols=[1])
	
	bp=plt.boxplot(F, positions = [i])
	plt.plot([i], F.mean(), 'rs')
	plt.setp(bp['boxes'], color='black')
	plt.setp(bp['whiskers'], color='black')
	plt.setp(bp['fliers'], color='black', marker=',')
	plt.setp(bp['medians'], color='black')
	ax = plt.gca()
	ax.set_ylim([0,1.05])
	i += 1
	return ax


files = [
	"breakpoint_detection_and_merging/breakpoint_detection_and_merging_Fmeasures_CBS.csv",
	"breakpoint_detection_and_merging/breakpoint_detection_and_merging_Fmeasures_HaMMLET.csv",
	"spatial_resolution_study/spatial_resolution_study_below5_Fmeasures_CBS.csv", "spatial_resolution_study/spatial_resolution_study_below5_Fmeasures_HaMMLET.csv", "spatial_resolution_study/spatial_resolution_study_5to10_Fmeasures_CBS.csv", "spatial_resolution_study/spatial_resolution_study_5to10_Fmeasures_HaMMLET.csv", "spatial_resolution_study/spatial_resolution_study_10to20_Fmeasures_CBS.csv", "spatial_resolution_study/spatial_resolution_study_10to20_Fmeasures_HaMMLET.csv", "spatial_resolution_study/spatial_resolution_study_above20_Fmeasures_CBS.csv", "spatial_resolution_study/spatial_resolution_study_above20_Fmeasures_HaMMLET.csv",
	"testing/testing_Fmeasures_CBS.csv",
	"testing/testing_Fmeasures_HaMMLET.csv"
	 ]

pos = [0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16]

fig = plt.figure()
bclr=["white", "grey"]
colors=[]
ax = plt.axes()
boxes=[]
means=[]
for i in xrange(len(files)):
	filename = files[i]
	F = np.loadtxt(filename, usecols=[1])	
	boxes.append(F)
	means.append(F.mean())
	plt.plot(pos[i], F.mean(), "ko", color=bclr[(i+1)%2],  zorder=len(files)+10)
	colors.append(bclr[i%2])
bp=plt.boxplot(boxes, positions=pos, patch_artist=True)
plt.setp(bp['whiskers'], color='black', linestyle='-')
plt.setp(bp['fliers'], color='black', marker='')
plt.setp(bp['medians'], color='black')
[bp["boxes"][i].set( color='black', facecolor = bclr[i%2] ) for i in xrange(len(files))]
[bp["medians"][i].set( color = bclr[(i+1)%2] ) for i in xrange(len(files))]


xlim = [-1,17]
ax = plt.gca()
ax.set_xlim(xlim)
ax.set_ylim([-0.05,1.05])
ax.set_xticks([3.5, 6.5, 9.5, 12.5])
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=18)
xtickNames = plt.setp(plt.gca(), xticklabels=["<5", "5-10", "10-20", ">20"])
plt.setp(xtickNames, fontsize=18)
plt.axvline(2, color='k')
plt.axvline(14, color='k')

ax2 = ax.twiny()
ax2.set_xlim(xlim)
ax2.set_xticks([0.5, 7.5, 15.5])
ax2.set_xticklabels(["BD&M", "SRS", "T"], fontsize=18)


plt.tight_layout()
plt.savefig('Fscores.pdf', bbox_inches="tight")


	


