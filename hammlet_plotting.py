#!/usr/bin/env python
# -*- coding:utf-8 -*-


###     Copyright 2014-2015 John Wiedenhoeft, Eric Brugel
###
###     This file is part of HaMMLET.
### 
###     HaMMLET is free software: you can redistribute it and/or modify
###     it under the terms of the GNU General Public License as published by
###     the Free Software Foundation, either version 3 of the License, or
###     (at your option) any later version.
### 
###     HaMMLET is distributed in the hope that it will be useful,
###     but WITHOUT ANY WARRANTY; without even the implied warranty of
###     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###     GNU General Public License for more details.
### 
###     You should have received a copy of the GNU General Public License
###     along with HaMMLET.  If not, see <http://www.gnu.org/licenses/>.


## @package Fmeasures
# This package implements the plotting routines for HaMMLET, including option for plotting known parameters from simulations.

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator 
import os




defaultColors = ["#a6cee3", "#1f78b4", "#fb9a99", "#e31a1c", "#b2df8a", "#33a02c", "#fdbf6f", "#ff7f00", "#cab2d6"]+['#ffed6f',
 '#ccebc5',
 '#bc80bd',
 '#d9d9d9',
 '#fccde5',
 '#b3de69',
 '#fdb462',
 '#80b1d3',
 '#fb8072',
 '#bebada',
 '#ffffb3',
 '#8dd3c7']

def lettercase(s):
	return s[0].upper()+s[1:]

def saveResultsFromFile(outFileName, dataFiles, sequenceFile,  nrStates, blocksFile=None, maxBlockSize=0, burnin=0, palette=defaultColors, allSequences=False, multipleOutputs=None, xlabel="Position", ylabel="Measurement"):
	
	if type(dataFiles) is str:
		dataFiles = [dataFiles]
		
	nrTracks = len(dataFiles)
	
	if allSequences:	# if we have all the state sequence samples
		seqs = np.loadtxt(sequenceFile, skiprows=1+burnin, dtype=int)
		nrValues = seqs.shape[1]
		marginals = np.zeros((nrStates, nrValues), dtype=int)
		for i in xrange(nrStates):
			marginals[i,:] = (seqs==i).sum(axis=0)	
	else:	# we already have the margins in seqFile
		marginals = np.loadtxt(sequenceFile, skiprows=1, dtype=int)
		nrValues = marginals.shape[1]

	blocks = None
	if blocksFile != None: 
		blocks = np.loadtxt(blocksFile)

	data = np.zeros((nrTracks, nrValues))
	ylim=[99999999999999999999999, -99999999999999999999999]
	for i in xrange(nrTracks):
		data[i, :] = np.loadtxt(dataFiles[i], dtype=float)[:,1]
		ylim[0] = min(ylim[0], data[i,:].min())
		ylim[1] = max(ylim[1], data[i,:].max())

	if multipleOutputs==None:
		saveResult(filename=outFileName, data=data, marginals=marginals, blocks=blocks, maxBlocksize=maxBlockSize, burnin=burnin, palette=palette, xlabel=xlabel, ylabel=ylabel)
	else:
		multipleOutputSizes = np.loadtxt(multipleOutputs, usecols=[0], dtype=int)
		multipleOutputSizes = np.asarray([0]+list(multipleOutputSizes)).cumsum()
		multipleOutputNames = np.loadtxt(multipleOutputs, usecols=[1], dtype=str)
		for i in xrange(len(multipleOutputSizes)-1):
			start = multipleOutputSizes[i]
			end = multipleOutputSizes[i+1]
			fn, fnext = os.path.splitext(outFileName)
			saveResult(filename="%s_%s%s" % (fn, multipleOutputNames[i], fnext), data=data, marginals=marginals, blocks=blocks, maxBlocksize=maxBlockSize, burnin=burnin, palette=palette, start=start, end=end, ylim=ylim, xlabel=xlabel, ylabel=ylabel)



def saveResult(filename, data, marginals, blocks=None, maxBlocksize=0, burnin=0, palette=defaultColors, start=0, end=None, ylim=None, xlabel="Position", ylabel="Measurement"):
	print filename
	if end==None:
		end = len(data[0])
	f = plt.Figure()
	f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)
	
	nrDataSets = data.shape[0]
	nrValues = data.shape[1]
	assert marginals.shape[1] == nrValues
	maxMargins = np.argmax(marginals, axis=0)
	
	# plot the data for each track in the data set
	for s in xrange(nrDataSets):
		plt.subplot(nrDataSets+1,1,s+1)
		if maxBlocksize>0:
			plotMaxBlocks(maxBlocksize, dataSize=nrValues, zorder=1)
		if blocks != None:
			plotCoarsestSegmentation(blocks, zorder=2)
		plotData(data[s, :], maxMargins, palette=palette, zorder=3, start=start, end=end, ylim=ylim, xlabel=xlabel, ylabel=ylabel)	
		plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower')) 

	# plot marginal distribution
	plt.subplot(nrDataSets+1,1,nrDataSets+1)
	stackedBarGraph(marginals, scale=True, palette=palette, offset=burnin, start=start, end=end, xlabel=xlabel, ylabel="marginal probabilities")
	f.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
	plt.savefig(filename, bbox_inches='tight')
	plt.close()


def plotData(data, states=None, palette=defaultColors, xlabel="position", ylabel="measurement", zorder=1, start=0, end=None, ylim=None): 
	if end==None:
		end = len(data)
	ax = plt.gca()
	if palette==None:
		color = ["black" for x in xrange(start, end)]
	else:
		color = [palette[i] for i in states[start:end]]
	ax.scatter(range(start, end), data[start:end], marker=".", linewidth=0, color=color, zorder=zorder)
	ax.set_xlim(start, end)
	if ylim != None:
		ax.set_ylim(ylim)
	ax.set_xlabel(lettercase(xlabel))
	ax.set_ylabel(lettercase(ylabel))
	return ax


def stackedBarGraph(M, x=None, scale=False, bins=False, palette=defaultColors, xlabel="position", ylabel="marginal counts", offset=0, start=0, end=None):	
	
	ax = plt.gca()
	nrSeries, nrPoints = M.shape	# M contains data in the rows, row 0 will be the lowest curve, row 1 stacked above and so forth
	if end==None:
		end = nrPoints
	nrPoints = end-start
	if x==None:
		x = np.array(range(start, end))
	elif bins:	
		x = np.cumsum(x)
	assert len(x) == nrPoints
	newM = np.zeros((nrSeries+1, 2+2*nrPoints))	# also add 0-series for fill_between reference
	newM[1:, 1:-2:2] = M[:,start:end]
	newM[1:, 2:-1:2] = M[:,start:end]
	if scale:
		newM /= np.sum(newM, axis=0).max()
	newM = np.cumsum(newM, axis=0)
	newX = np.zeros(2+2*nrPoints)
	if bins:	# if bins==True, x is interpreted as bin widths
		newX[2:-1:2]=x
		newX[3::2]=x
	else:			
		newX[0:-3:2]=x
		newX[1:-2:2]=x
		newX[-1] = max(x)+1
		newX[-2] = newX[-1]
		newX -= 0.5
	for i in xrange(nrSeries):
		ax.fill_between(newX, newM[i+1], newM[i], linewidth=0.0, color=palette[i])
	if scale:
		ax.set_ylim(0, 1)
	else:
		ax.set_ylim(0, max(newM[-1,:]))	
	ax.set_xlim(start, end)
	
	plt.xlabel(lettercase(xlabel))
	plt.ylabel(lettercase(ylabel))
	return ax



def plotCoarsestSegmentation(blocks, color="lightgrey", linewidth=1,  zorder=1):
	ax = plt.gca()
	blocks = blocks.cumsum()-0.5
	for i in xrange(len(blocks)):
		ax.axvline(blocks[i], color=color, zorder=zorder, linewidth=linewidth)
	return ax


def plotMaxBlocks(maxBlocksize, dataSize, color="lightgrey", linewidth=1, zorder=1):
	ax = plt.gca()
	if maxBlocksize>0:
		b = 0
		while b <=dataSize:
			ax.axvline(maxBlocksize*b, color=color, zorder=zorder, linewidth=linewidth)
			b+=1
	return ax



def saveParameters(filename, means, variances, trueMeans=None, trueVariances=None, burnin=0, palette=defaultColors):
	f = plt.Figure()
	plt.subplot(2,1,1)
	plotParameters(means, trueMeans, title="means", offset=burnin)
	plt.subplot(2,1,2)
	plotParameters(variances, trueVariances, title="variances", offset=burnin)
	plt.tight_layout()
	plt.savefig(filename, bbox_inches='tight')
	plt.close()

	

def plotParameters(theta, trueTheta=None, palette=defaultColors, title="", xlim=None, ylim=None, offset=0):
	ax = plt.gca()
	nrStates = theta.shape[0]
	nrIterations = theta.shape[1]
	if trueTheta != None:
		assert len(trueTheta)==nrStates
	for i in xrange(nrStates):
		ax.plot(range(-offset, nrIterations-offset), theta[i,:], linewidth=1, color=palette[i], zorder=2)
		if offset>0:
			ax.axvline(0, color="grey")
		if trueTheta != None:
			ax.axhline(trueTheta[i], color="grey", linewidth=1, zorder=1)
	ax.set_title(title)
	if xlim != None:
		ax.set_xlim = xlim
	if ylim != None:
		ax.set_ylim = ylim
	return ax


