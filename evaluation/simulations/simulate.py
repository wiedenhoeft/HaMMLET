#!/usr/bin/env python
import numpy as np
from numpy import sqrt
from numpy.random import randint, shuffle
import conf
from scipy.stats import norm
import os
import yaml

# PARAMETERS:
# nrCNstates: number of CN states
# nrSamples: length of the state sequence, for aCGH: how many probes along chromosome 
# nrGains
# nrLosses
# lenGains
# lenLosses
# nrDeNovoGains
# nrDeNovoLosses


def randomPartition(n, k, nonzero=True):
	"""Randomly split n into k summands. Nonzero means no summand can be zero."""
	if k <= 0:
		return np.array([n])
	if nonzero:
		assert n>=k
		splits = set()
		while len(splits) < k-1:	# to obtain k segments, we need k-1 splits
			splits.add(randint(1, n))
	else:
		assert n >=0
		splits = randint(0, n+1, size=k-1)
	splits = [0]+sorted(splits)+[n]
	result = np.array([splits[i+1]-splits[i] for i in xrange(k)])
	assert result.sum() == n, (result.sum(), n)
	return result

def embedFragments(totalFragLen, nrFrag, fragSym, totalSeqLen, seqSym):
	"""Create <nrFrag> random fragments of a sequence of <totalFragLen> copies of <fragSym>, and embed it randomly into a sequence of  <seqSym> s.t. the total length is <totalSeqLen>"""
	fragments = [[fragSym for i in xrange(x)] for x in randomPartition(totalFragLen, nrFrag)]
	bedding = [[seqSym for i in xrange(x)] for x in randomPartition(totalSeqLen-totalFragLen, nrFrag+1, False)]
	result = bedding[0]
	for i in xrange(len(fragments)):
		result += fragments[i]
		result += bedding[i+1]
	return result
	
	
def intPartitionKElem(n, k, lower=None, upper=None):
	if lower==None:
		lower = np.ones(k, dtype=int)
	if upper==None:
		upper = n*np.ones(k, dtype=int)
	"""Generates all integer solutions to sum_{i=1}^k x_i = n."""
	result = np.array(lower, dtype=int)
	result[-1] = n-k+1
	cumul = result.cumsum()
	result[-1] = n - cumul[-2]
	maxima = np.array(range(n+1), dtype=int)[-k:]	
	
	# make sure the values don't get so high that the sum to the right cannot reach n
	for i in reversed(xrange(1, k)):
		maxima[i-1] = maxima[i]-lower[i]
	maxima = np.minimum(maxima, upper)	# make sure individual values don't get too high
	if np.all(result<=maxima) and np.all(result>=lower) :
		yield result
		
	p = k-2	# pointer to the value we try to increase
	while True:
		while cumul[p] >= maxima[p]:	# check if we can increase result[p]
			p -= 1
			if p<0:
				return
		result[p] += 1
		cumul[p] += 1
		
		# fill everything to the right with the lowest possible values
		for i in xrange(p+1, k-1):
			result[i] = lower[i]
			cumul[i] = cumul[i-1]+lower[i]
		
		# set the last value as the difference to n, only yield if we are within bounds
		result[-1] = n-cumul[-2]
		if result[-1]>=lower[-1] and result[-1]<=upper[-1]:
			assert result.sum()==n, (result, cumul)
			assert np.all(result>=lower)
			assert np.all(result<=upper)
			yield result
		p=k-2	# the last value is set automatically, so the penultimate is the one we try to increase in the next step

		
		
			
if __name__ == "__main__":			
	conf = conf.Config()
	
	states = np.zeros(conf.nrSamples, dtype=int)
	gainsP = [np.ones(x, dtype=int)*3 for x in randomPartition(conf.lenGains, conf.nrGains)]
	lossesP = [np.ones(x, dtype=int)*1 for x in randomPartition(conf.lenLosses, conf.nrLosses)]
	neutralP = [np.ones(x, dtype=int)*2 for x in randomPartition(conf.nrSamples-conf.lenLosses-conf.lenGains, conf.nrLosses+conf.nrGains+1, False)]
	result = gainsP + lossesP + neutralP
	shuffle(result)
	states = np.concatenate(result)
	
	#restore numbers to allowed values
	states[states<1]=1
	states[states>conf.nrCNstates]=conf.nrCNstates
	
	data = np.zeros(conf.nrSamples)
	for c in xrange(conf.nrSamples):
		cn = states[c]
		data[c] = norm.rvs(loc=conf.means[cn-1], scale=sqrt(conf.variances[cn-1]))	# -1 because we use actual copy number (1-based), not index (0-based)
	data = np.dstack((states-1, data))[0,:,:]
	samplepath=os.path.abspath(conf.simSampleFile)
	np.savetxt(samplepath, data, fmt="%d\t%f")
	