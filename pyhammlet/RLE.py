#!/usr/bin/env python
# -*- coding:utf-8 -*-

import numpy as np
import copy
from scipy import stats
from bisect import bisect_right, bisect_left


def find_gt(a, x):
	'Find leftmost index if item greater than x'
	i = bisect_right(a, x)
	if i != len(a):
		return i
	raise ValueError


def find_ge(a, x):
    'Find leftmost index of item greater than or equal to x'
    i = bisect_left(a, x)
    if i != len(a):
        return i
    raise ValueError


# "subtractive difference" as inverse of "cumulative sum"
def subdiff(ends):
	result = copy.deepcopy(ends)
	result[1:] -= ends[:-1]
	return result







# A NumPy array, with a run-length encoding of its first axis.
class RunLengthArray:

	# self.ends
	# self.array
	
	def __init__(self, *args, **kwargs):
		if kwargs.has_key("array"):
			self.array = np.array(kwargs["array"])
			if kwargs.has_key("sizes"):
				self.ends = np.array(kwargs["sizes"], dtype=int)
			else:
				self.ends = np.ones(self.array.shape[0], dtype=int)
			self.ends = self.ends.cumsum();	# store the first position after the current block, for log-query
			self.size = self.ends[-1]
			
			assert len(self.ends) == len(self.array)

			assert len(self.ends) == self.array.shape[0], "Length of runlength array (%d) does not match the size of the data (%s)!" % (len(self.ends),  str(self.array.shape))
		else:
			self.size=0
			self.array=None
			self.ends=None
		
		
	def loadtxt(self, filename, dtype=int):
		array = np.loadtxt(filename, dtype=dtype)	
		self.__init__( sizes=array[:,0], array=array[:,1:])


	def compress(self):
		"""Compress the data as much as possible."""
		size = self.ends[-1]
		L = 0	# index for insertion
		for R in xrange(1,len(self.ends)):
			if not np.array_equal(self.array[L], self.array[R]):
				L += 1
				if L != R:
					self.array[L] = self.array[R]
			self.ends[L] = self.ends[R]
		L+=1
		self.ends = self.ends[0:L]
		self.array = self.array[0:L]
		assert len(self.ends)==len(self.array)
		assert self.ends[-1]==size

	# decompressed shape
	def shape(self):
		result = list(self.array.shape)
		result[0] = self.ends[-1]
		return tuple(result)
		
		
	def __str__(self):
		result = "%s %s\n" % (self.ends[0].__str__(), self.array[0].__str__())
		for b in xrange(1, len(self.ends)):
			result += "%s %s\n" % ((self.ends[b]-self.ends[b-1]).__str__(), self.array[b].__str__())
		return result

		
	# Returns an RLE array of indices of maxima
	def mode(self):
		res = stats.mode(self.array, axis=1)[0].astype(self.array.dtype).reshape(self.array.shape[0])
		return RunLengthArray(sizes=subdiff(self.ends), array=res)

	# Returns an RLE array of argmax
	def argmax(self):
		res = np.argmax(self.array, axis=1).astype(self.array.dtype)#.reshape(self.array.shape[0])
		return RunLengthArray(sizes=subdiff(self.ends), array=res)


	def __getitem__(self, key):
		return self.array[find_gt(self.ends, key)]
		
	
	# Returns a RunLengthArray object corresponding to the array in the slice  [start:end] (in uncompressed coordinates)
	def __getslice__(self, start = None, end = None, step = None):
		if step is not None and step != 1:
			assert False, "RunLengthArray does not support steps in slices!"
		if end is None:
			end = start+self.size
		assert start >=0
		assert end <= self.size, "End index %d > size of compressed array %d!" % (end, self.size)
		
		L = find_gt(self.ends, start)	# get left index
		R = find_ge(self.ends, end)	# get right index
		R+=1
		segmentSizes = copy.deepcopy(self.ends[L:R])
		# transform right boundaries to block sizes
		segmentSizes[-1]=end	# set right boundary directly
		i = len(segmentSizes)-1
		while i>0:
			segmentSizes[i] -= segmentSizes[i-1]
			i-=1
		segmentSizes[0] -= start	# adjust first blocksize
		
		assert segmentSizes.sum()==end-start, "Sum of segment sizes (%d) does not match requested range (%d)!" %(segmentSizes.sum(), end-start)
		return RunLengthArray(sizes=segmentSizes, array=self.array[L:R])
		
	
	def getSegment(self, index):
		assert index < len(self.array)
		if index==0:
			return (self.ends[index], self.array[index])
		else:
			return (self.ends[index]-self.ends[index-1], self.array[index])
		
		
	# Returns an uncompressed array
	def decompress(self, start=None, end=None):
		if start is None and end is None:
			return np.repeat(self.array, subdiff(self.ends), axis=0)
		else:
			# TODO more efficiently without creating new object
			return self[start:end].decompress()
			
	def blocksizes(self):
		return subdiff(self.ends)
	
	def nrSegments(self):
		return len(self.ends)
	
	def __len__(self):
		return self.ends[-1]
	
	
	
	
def shatter(A, B):
	"""Takes two RLE and joins their break points so that they share the same block structure in-place."""
	ends = set()
	assert A.size == B.size
	ends.update(A.ends)
	ends.update(B.ends)
	ends = np.array(sorted(ends), dtype=int)
	L = len(ends)
	valA = [None for l in xrange(L) ]
	valB = [None for l in xrange(L) ]
	a=0
	b=0
	for i in xrange(L):
		if A.ends[a]<ends[i]:
			a+=1
		valA[i] = A.array[a]
		if B.ends[b]<ends[i]:
			b+=1
		valB[i] = B.array[b]
	A.ends = ends
	B.ends = ends
	A.array = np.array(valA)
	B.array = np.array(valB)
	


