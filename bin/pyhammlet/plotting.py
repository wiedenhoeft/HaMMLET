#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""Generic plotting functions used by HaMMLET and associated scripts and programs. """



from __future__ import print_function, division
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator 
import os
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes
from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm, Colormap








# like plt.imshow, but scaled down to a defined number of pixels 
def scaledImshow(matrix, cmap, norm,  maxNrPixels = 100000000, intType="uint16", *args, **kwargs):
	ax = plt.gca()
	height, width = matrix.shape
	nrPixels = height*width
	if nrPixels > maxNrPixels:
		try:
			from scipy.ndimage import zoom
		except:
			print("[ERROR] You must have SciPy installed to plot huge data!")
		scalingFactor = np.sqrt(maxNrPixels/nrPixels)		
		mat = zoom(matrix, scalingFactor, mode="nearest")
		assert mat.ndim==2
		ax.imshow(mat, cmap, norm,  *args, **kwargs)
	else:
		assert matrix.ndim==2
		ax.imshow(matrix, cmap, norm, *args, **kwargs)
	return ax









def lettercase(s):
	return s[0].upper()+s[1:]





### MATRIX PLOTTING, e.g. for marginals ###

def sortMatrix(
	matrix,
	axis=0,	# 1 for rows, 0 for columns
	reverse=False):
	return np.apply_along_axis(lambda x, r: np.array(sorted(x, reverse=r)), axis, matrix, not reverse)
		
	
def sortByFrequency(
	a,
	reverse=False):
	binCounts = np.bincount(a)
	idx=np.argsort(binCounts)
	if reverse:
		return np.repeat(np.array(range(max(a)+1))[idx], binCounts[idx])[::-1]
	else:
		return np.repeat(np.array(range(max(a)+1))[idx], binCounts[idx])


def sortMatrixByFrequency(
	m,
	axis = 0,
	reverse=False):
	"""Sorts matrix axes by the number of occurrences in each array, e.g. when sorting columns, each column will have its rarest element first and its most common entry last."""
	return np.apply_along_axis(sortByFrequency, axis, m, not reverse)




def plotMatrix(
	m,
	xlabel="Position along chromosome",
	ylabel="Marginal counts",
	xstretch=1,
	xmin=0,
	normalize=False,
	*args,
	**kwargs):
	
	ax = plt.gca()
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ymax, xmax = m.shape
	if normalize:
		ymax=1
	ext = [xmin, xmin+xmax*xstretch, 0, ymax]	# TODO include xstretch in xmin?
	scaledImshow(m, extent = ext, aspect="auto", origin="lower", interpolation="none",  *args, **kwargs)
	return ax







def matrixQuantilePlot(
	data, 
	quantiles = range(5, 100, 5),
	xlabel="Iteration", 
	ylabel="F-measure (quantiles)", 
	cmap="Blues",
	mincolor = 0.1,
	maxcolor=0.9,
	ylim = None,
	insetXlim=None, 
	insetYlim=None, 
	insetXticks=None, 
	insetYticks=None, 
	insetWidth="40%", 
	insetHeight="40%", 
	insetLoc=4):
	
	from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes
	inset = insetXlim != None
	percentiles = np.percentile(data, quantiles, axis=0)
	iterations = data.shape[1]
	
	
	ax = plt.gca()
	if inset:
		axins = inset_axes(ax, width=insetWidth, height=insetHeight, loc=4)
	
	colormap = plt.cm.get_cmap(cmap)
	
	for i in xrange(len(quantiles)):
		q = quantiles[i]
		if q == 50:
			ax.plot(percentiles[i], color="black",zorder=len(quantiles), linewidth=2)
			if inset:
				axins.plot(percentiles[i], color="black",zorder=len(quantiles), linewidth=2)
		if q > 50:
			break
		

		color = colormap(((q/100)/(maxcolor-mincolor)+mincolor))
		#getColor(value, cmap, minColor, maxColor)
		
		ax.fill_between(range(iterations), percentiles[i], percentiles[-i-1], color=color, linewidth=1, zorder=i)
		ax.plot(percentiles[i], color="black", linewidth=0.3, zorder=len(quantiles)+1)
		ax.plot(percentiles[-i-1], color="black", linewidth=0.3, zorder=len(quantiles)+1)
		
		if inset:
			axins.fill_between(range(iterations), percentiles[i], percentiles[-i-1], color=color, linewidth=1, zorder=i)
			axins.plot(percentiles[i], color="black", linewidth=0.3, zorder=len(quantiles)+1)
			axins.plot(percentiles[-i-1], color="black", linewidth=0.3, zorder=len(quantiles)+1)
		
	
	ymin=min([min(p) for p in percentiles])
	ymax = max([max(p) for p in percentiles])
	if ylim != None:
		ymin = min(ymin, ylim[0])
		ymax = min(ymax, ylim[1])
	
	margin = (ymax-ymin)*0.05
	ax.set_ylim([ymin-margin,ymax+margin])
	ax.set_xlim([-iterations/20, iterations+iterations/20 ])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	
	if inset:
		axins.set_xlim(insetXlim)
		axins.set_ylim(insetYlim)
		axins.xaxis.tick_top()
		axins.set_xticks(insetXticks)
		axins.set_yticks(insetYticks)
	#plt.sca(ax) 








# given data and an array of states of the same size, create a scatter plot for the start:end slice, colored according to the states
def plotData(
	data, 
	states=None,
	start=0, 
	end=None, 
	marker=".",
	linewidth=0,
	alpha=0.8,
	ylabel=None,
	cmap=None,
	norm=None,
	*args,
	**kwargs): 
	
	ax = plt.gca()
	
	if states is not None:
		assert cmap is not None, "Require colormap to plot data if a state sequence is provided!"
		assert norm is not None, "Require color normalization to plot data if a state sequence is provided!"
			
	#TODO chunksize
	#TODO multivariate
	if end is None:
		end = start+len(data)
	if states is None:
		states="k"
	else:
		states = states[start:end]

	ax.scatter(range(start, end), data[start:end], c=states,  marker=marker, linewidth=linewidth, alpha=alpha, cmap=cmap, norm=norm, *args)
	ax.set_ylabel(ylabel)
	ax.set_xlim([start, end])
	return ax



#even blocks are white, odd are dark
def plotBlockSizes(
	compressedBlocks, 
	start=0, 
	end=None,  
	chunkSize=1, 
	*args,
	**kwargs):
	ax = plt.gca()
	M=compressedBlocks[start:end].decompress().T
	# TODO more customizable way to plot boundaries (or not!)
	M[M==1] = 0
	plotMatrix(M,  xmin=start, xstretch=chunkSize,  *args, **kwargs)	# cmap="Blues", Spectral"
	return ax
	

# takes compressed marginal counts, extracts the [start:end] slice, expands the counts to a full matrix of iterations and plots it
def plotMarginals(
	compressedMarginalCounts, 
	cmap,
	norm,
	start=0, 
	end=None, 
	chunkSize=1, 
	*args,
	**kwargs):
	ax = plt.gca()
	counts = compressedMarginalCounts[start:end].decompress()
	
	# expand the counts to full array
	nrSegments, nrStates = counts.shape
	iterations = counts[0,:].sum()
	marginals = np.zeros((nrSegments, iterations), dtype=int)
	for i in xrange(nrSegments):
		marginals[i,:] = np.repeat(np.array(range(nrStates)), counts[i,:])
	plotMatrix(marginals.T,  xstretch=chunkSize, xmin=start, cmap=cmap, norm=norm, *args, **kwargs)
	del marginals
	del counts
	return ax



def plotSequences(
	sequences, # a list of RLE vectors
	cmap,
	norm,
	start=0, 
	end=None, 
	*args,
	**kwargs
	):
	
	ax = plt.gca()
	if end is None:
		end = sequences[0].size
	matrix = np.zeros([len(sequences), end-start], dtype=int)
	for i in xrange(len(sequences)):
		matrix[i] = sequences[i][start:end].decompress()
	#TODO vmin, vmax
	plotMatrix(matrix,  ylabel="Iteration", xmin=start, normalize=False, cmap=cmap, norm=norm, *args, **kwargs)
	return ax
	
	