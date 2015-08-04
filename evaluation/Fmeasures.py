#!/usr/bin/env python
# -*- coding:utf-8 -*-

## @package Fmeasures
# This package implements the micro- and macro-averaged F-measures by 
# Özgür, A., Özgür, L., & Güngör, T. (2005). Text Categorization with Class-Based and Corpus-Based Keyword Selection, 606–615. doi:10.1007/11569596_63 

from __future__ import division
import numpy as np
import itertools
import matplotlib.pyplot as plt
import copy
import os

from scipy import stats


# Edge cases are handled such that TP+FP=0 -> precision=1,  and TP+FN=0 -> recall=1.
def microFmeasure(TP, FP, FN):
	TP = np.asarray(TP)
	FP = np.asarray(FP)
	FN = np.asarray(FN)
	TPsum = TP.sum()
	TPFPsum = (TP+FP).sum()
	if TPFPsum == 0:	# edge cases
		precision = 1
	else:
		precision = TPsum/TPFPsum
	TPFNsum = (TP+FN).sum()	
	if TPFNsum == 0:	# edge cases
		recall = 1
	else:
		recall = TPsum/TPFNsum
	if precision+recall==0:
		return 0
	else:
		return 2*precision*recall/(precision+recall)


# Edge cases are handled such that TP+FP=0 -> precision=1,  and TP+FN=0 -> recall=1.
def macroFmeasure(TP, FP, FN):
	TP = np.asarray(TP)
	FP = np.asarray(FP)
	FN = np.asarray(FN)
	with np.errstate(invalid='ignore'):	# we expect invalid values for edge cases, but we take care of them
		precision = TP/(TP+FP)
		precision[(TP+FP==0)] = 1	# edge cases
		recall = TP/(TP+FN)
		recall[(TP+FN)==0] = 1	# edge cases
		Fi = 2*precision*recall/(precision+recall)
		Fi[(precision+recall)==0] = 0
		result = Fi.sum()/(len(Fi))
		return result




	
## Takes an array of labels 0<=label<nrLabels, and creates a matrix  nrLabels x len(data) that contains a 1 at [i,j] if the j-th label is i, and 0 otherwise.
def splitLabels(data, nrLabels):
	assert data.max() < nrLabels, data.max()
	classes = np.zeros((nrLabels, len(data)), dtype=bool)
	for c in xrange(nrLabels):
		classes[c,:] = data == c
	return classes



## Given an array of split labels, rearranges the rows by the ranks of  classValue.
def relabel(classes, classValue):
	return classes[np.argsort(classValue)]


	
## Returns an array of the true positives for each class (i-th entry is TP of class i).
def truePositives(inferredClasses, trueClasses):
	return np.logical_and(trueClasses, inferredClasses).sum(axis=1)



## Returns an array of the false positives for each class (i-th entry is FP of class i).
def falsePositives(inferredClasses, trueClasses):
	return np.logical_and(np.logical_not(trueClasses), inferredClasses).sum(axis=1)



## Returns an array of the false negatives for each class (i-th entry is FN of class i).
def falseNegatives(inferredClasses, trueClasses):
	return np.logical_and(trueClasses, np.logical_not(inferredClasses)).sum(axis=1)


def Fmeasure(TP, FP, FN, positiveLabel=1):
	TP = TP[positiveLabel]
	FP = FP[positiveLabel]
	FN = FN[positiveLabel]
	precision = TP/(TP+FP)
	recall = TP/(TP+FN)
	return 2* (precision*recall)/(precision+recall)



def majorityFmeasure(inferredClasses, trueClasses):
	"""Takes the most common label as negative, and all others as positive. This is used to calculate aberrations from the base case"""
	mostFrequentTrueClass = int(stats.mode(trueClasses)[0][0])
	trueAberrations = (trueClasses != mostFrequentTrueClass)
	
	mostFrequentInferredClass = int(stats.mode(inferredClasses)[0][0])
	inferredAberrations = (inferredClasses != mostFrequentInferredClass)
	
	TP = truePositives(inferredAberrations, trueAberrations)
	FP = falsePositives(inferredAberrations, trueAberrations)
	FN = falseNegatives(inferredAberrations, trueAberrations)
	return Fmeasure(TP, FP, FN, 1)

	

def majorityFmeasureFromSamples(simSampleFile, prefix, nrCNstates, iterations, burnin=0, plot=True, oneBased=False):
	# takes the true state sequence and state sequence samples and calculates the F-scores, assigning labels such that the sum of all F-scores is maximal
	trueStateSeq = np.loadtxt(simSampleFile, dtype="int", usecols=(0,))
	if oneBased:
		trueStateSeq-=1
	trueStateSeqByClass = splitLabels(trueStateSeq, nrCNstates)
	
		
	sampledStateSeqs = np.loadtxt(prefix+"_statesequence.csv", dtype="int", skiprows=1)
	
	nrSamples = sampledStateSeqs.shape[0]
	Fscores = np.zeros(nrSamples)
	FscoresTmp = np.zeros(nrSamples)
	
	# try all permutations of means and see which one is the largest
	
	for i in xrange(nrSamples):
		sampledStateSeqByClass = splitLabels(sampledStateSeqs[i,:], nrCNstates)
		TP = truePositives(sampledStateSeqByClass, trueStateSeqByClass)
		FP = falsePositives(sampledStateSeqByClass, trueStateSeqByClass)
		FN = falseNegatives(sampledStateSeqByClass, trueStateSeqByClass) 
		
		# calculate multi-class F-scores
		FscoresTmp[i] = majorityFmeasure(TP, FP, FN)
		
		if microFscoresTmp.sum()+macroFscoresTmp.sum() > microFscores.sum()+ macroFscores.sum():
			microFscores = copy.deepcopy(microFscoresTmp)
			macroFscores = copy.deepcopy(macroFscoresTmp)
		
	np.savetxt(prefix+"_majorityF.tsv", microFscores)
	if plot:
		plt.plot(range(burnin, iterations), microFscores, color="#d73027", label="micro-averaged")
		plt.plot(range(burnin, iterations), macroFscores, color="#1a9850", label="macro-averaged")
		plt.ylim([0,1.05])
		plt.xlabel("iteration")
		plt.ylabel("multi-class F-measures")
		plt.legend(loc="lower right")
		plt.tight_layout()
		plt.savefig(prefix+"_F-scores.png")
		plt.close()


def FmeasuresFromSamples(simSampleFile, prefix, nrCNstates, iterations, burnin=0, plot=True, oneBased=False):
	# takes the true state sequence and state sequence samples and calculates the F-scores, assigning labels such that the sum of all F-scores is maximal
	trueStateSeq = np.loadtxt(simSampleFile, dtype="int", usecols=(0,))
	if oneBased:
		trueStateSeq-=1
	trueStateSeqByClass = splitLabels(trueStateSeq, nrCNstates)
	
		
	sampledStateSeqs = np.loadtxt(prefix+"_statesequence.csv", dtype="int", skiprows=1)
	
	nrSamples = sampledStateSeqs.shape[0]
	microFscores = np.zeros(nrSamples)
	macroFscores = np.zeros(nrSamples)
	microFscoresTmp = np.zeros(nrSamples)
	macroFscoresTmp = np.zeros(nrSamples)
	
	# try all permutations of means and see which one is the largest
	for avg in itertools.permutations(range(nrCNstates)):
		for i in xrange(nrSamples):
			sampledStateSeqByClass = splitLabels(sampledStateSeqs[i,:], nrCNstates)
			
			sampledStateSeqByClass = relabel(sampledStateSeqByClass, avg)
			TP = truePositives(sampledStateSeqByClass, trueStateSeqByClass)
			FP = falsePositives(sampledStateSeqByClass, trueStateSeqByClass)
			FN = falseNegatives(sampledStateSeqByClass, trueStateSeqByClass) 

			# calculate multi-class F-scores
			microFscoresTmp[i] = microFmeasure(TP, FP, FN)
			macroFscoresTmp[i] = macroFmeasure(TP, FP, FN)
		if microFscoresTmp.sum()+macroFscoresTmp.sum() > microFscores.sum()+ macroFscores.sum():
			microFscores = copy.deepcopy(microFscoresTmp)
			macroFscores = copy.deepcopy(macroFscoresTmp)
		
	np.savetxt(prefix+"_microF.tsv", microFscores)
	np.savetxt(prefix+"_macroF.tsv", macroFscores)
	if plot:
		plt.plot(range(burnin, iterations), microFscores, color="#d73027", label="micro-averaged")
		plt.plot(range(burnin, iterations), macroFscores, color="#1a9850", label="macro-averaged")
		plt.ylim([0,1.05])
		plt.xlabel("iteration")
		plt.ylabel("multi-class F-measures")
		plt.legend(loc="lower right")
		plt.tight_layout()
		plt.savefig(prefix+"_F-scores.png")
		plt.close()


def FmeasureQuantilePlotsFromPrefices(prefixes,  iterations, hash="", xlabel="Iteration", ylabel="F-measure (quantiles)"):
	"""Takes a list of file prefixes and aggregates their micro- and macro-averaged F-measures into a quantile plot."""

	#boxplot has data in columns	
	Fmeasures = np.zeros((iterations, len(prefixes)))
	for i in xrange(len(prefixes)):
		prefix = prefixes[i]

		filename = prefix+".tsv"
		if os.path.exists(filename):
			fm=np.loadtxt(filename)
			assert len(fm)==iterations
			Fmeasures[:,i] = fm
		else:
			print "[ERROR] No such file", filename
	FmeasureQuantilePlots(Fmeasures,  iterations, hash, xlabel, ylabel)	# TODO FmeasureQuantilePlots takes filename now, not array


def FmeasureQuantiles(FmeasuresFile):
	quantilesFile = FmeasuresFile+".quantiles"
	if os.path.exists(quantilesFile):
		print quantilesFile, "already exists."
		return np.loadtxt(quantilesFile)
	else:
		Fmeasures = np.loadtxt(FmeasuresFile).T
		quantiles = range(5, 100, 5)
		percentiles = np.percentile(Fmeasures, quantiles, axis=1)
		np.savetxt(quantilesFile, percentiles)
		return percentiles

def FmeasureQuantilePlots(FmeasuresFile,  iterations, hash="", xlabel="Iteration", ylabel="F-measure (quantiles)", insetXlim=None, insetYlim=None, insetXticks=None, insetYticks=None, insetWidth="40%", insetHeight="40%", insetLoc=4):
	percentiles = FmeasureQuantiles(FmeasuresFile)
	quantiles = range(5, 100, 5)	# TODO move into quantiles file
	from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes
	inset = insetXlim != None
	
	ax = plt.gca()
	if inset:
		axins = inset_axes(ax, width=insetWidth, height=insetHeight, loc=4)
	
	for i in xrange(len(quantiles)):
		q = quantiles[i]
		if q==50:
			ax.plot(percentiles[i], color="black",zorder=len(quantiles), linewidth=2)
			if inset:
				axins.plot(percentiles[i], color="black",zorder=len(quantiles), linewidth=2)
		if q > 50:
			break
		
		mincolor = 0.1
		maxcolor=0.9
		ax.fill_between(range(iterations), percentiles[i], percentiles[-i-1], color=str(1-((q/100)/(maxcolor-mincolor)+mincolor)), linewidth=1, zorder=i)
		ax.plot(percentiles[i], color="black", linewidth=0.3, zorder=len(quantiles)+1)
		ax.plot(percentiles[-i-1], color="black", linewidth=0.3, zorder=len(quantiles)+1)
		
		if inset:
			axins.fill_between(range(iterations), percentiles[i], percentiles[-i-1], color=str(1-((q/100)/(maxcolor-mincolor)+mincolor)), linewidth=1, zorder=i)
			axins.plot(percentiles[i], color="black", linewidth=0.3, zorder=len(quantiles)+1)
			axins.plot(percentiles[-i-1], color="black", linewidth=0.3, zorder=len(quantiles)+1)
		
		
		
	ax.set_ylim([0,1.05])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	
	if inset:
		axins.set_xlim(insetXlim)
		axins.set_ylim(insetYlim)
		axins.xaxis.tick_top()
		axins.set_xticks(insetXticks)
		axins.set_yticks(insetYticks)
	plt.sca(ax)