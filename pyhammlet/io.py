#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""Basic I/O capabilities for HaMMLET, such as parsing compressed output etc."""

import numpy as np
import copy
from scipy import stats
import itertools
from bisect import bisect_right, bisect_left
from RLE import *


# TODO this reads marginals by iteration, not by class, which would be smaller	
def readMarginals(filename):
	M = np.loadtxt(filename, dtype=int, ndmin=2)
	segSizes = M[:,0]
	counts = M[:,1:]
	result = RunLengthArray(sizes=segSizes, array=counts)
	return result



#TODO start and end
def readCompressedStateSequences(seqFileName):	
	"""Reads compressed state sequence samples into a list of RLE'ed state sequences."""
	seqFile = file(seqFileName, "r")
	splitLines = [line.split() for line in filter(lambda x: not x.startswith("#"), seqFile.readlines())]
	nrIterations = len(splitLines)
	result = [[] for x in xrange(nrIterations)]
	for i in xrange(nrIterations):
		a, b = zip(*[x.split(":") for x  in splitLines[i]])
		assert len(a)==len(b)
		result[i] = RunLengthArray(sizes=np.array(a, dtype=int), array=np.array(b, dtype=int))
	return result
	


def readBlockSizes(filename):
	# read block sizes from file and create a compressed matrix in which each position corresponds to  (log of) the size of the block the value at this position is contained in. This yields a plot of the local compression.
	
	lines = [np.array(line.split(), int).cumsum() for line in file(filename).readlines()]
	for i in xrange(len(lines)-1):
		assert lines[i][-1] == lines[i+1][-1], "Block structure in input line %d does not match the previous ones in total size!" % (i+2)
	iterations = len(lines)
	
	# get all endpoints across all iterations
	ends = set()
	ends.update(*[set(L) for L in lines])
	ends = np.array(sorted(ends), dtype=int)
		
	# change end points to match size indices for run-length encoding
	for line in lines:
		e = 0
		for i in xrange(len(line)):
			while line[i] != ends[e]:
				e += 1
			e += 1
			line[i] = e
	for i in xrange(len(lines)-1):
		assert lines[i][-1] == lines[i+1][-1]
		
	# determine block labels
	blockedLabels = dict()	# position (start or end): set of forbidden labels
	validLabels=dict()	# valid label for (start, end)-tuples
	data = np.zeros((len(ends), iterations))
	iteration=0
	ends = subdiff(ends)
	prevLabel=0
	for line in lines:
		for i in xrange(len(line)):
			if i==0:
				start = 0
				end = line[0]
			else:
				start = line[i-1]
				end= line[i]
			label = np.log(sum(ends[start:end]))
			data[start:end, iteration] = label
		iteration += 1
	del blockedLabels

	return RunLengthArray(sizes=np.array(ends), array=data)





