//     Copyright 2014-2015 John Wiedenhoeft, Eric Brugel
//
//     This file is part of HaMMLET.
//
//     HaMMLET is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     HaMMLET is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with HaMMLET.  If not, see <http://www.gnu.org/licenses/>.


#ifndef WAVELETTREE_HPP
#define WAVELETTREE_HPP
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include<algorithm>
#include<iterator>
#include<string>
#include<sstream>
#include "float.h"
#include "limits.h"
#include<unordered_map>
using namespace std;

typedef unsigned int uint;
typedef double dataType;


class SizeRange {
	// This keeps track of the range of treshold values that induce a certain block structure (as measured by its size), in order to avoid recreating the same block structure
	double lowest;
	double highest;
public:

	SizeRange() {
		lowest = -1;
		highest = -1;
	}

	SizeRange( double x ) {
		lowest =  x;
		highest =  x;
	}

	void update( double x ) {
		if ( lowest < 0 || highest < 0 ) {
			lowest = x;
			highest = x;
			return;
		}

		if ( x < lowest ) {
			lowest = x;
		} else {
			if ( x > highest ) {
				highest = x;
			}
		}
	}

	bool contains( double x ) {
		if ( lowest < 0 || highest < 0 ) {
			return false;	// return false if the range has not been initialized with proper values; this happens if we query a non-existing range, as it will be created
		}

		if ( lowest <= x && x <= highest ) {
			return true;
		} else {
			return false;
		}
	}
};


class WaveletNode {
public:
	unsigned char level;			// 2^level is the number of nodes in subtree
	uint size;
	dataType coeff;		// the max of all absolute wavelet coefficients
	dataType sum;		/// this holds the original data in the leaves
	dataType sumOfSquares;
	vector<uint> marginals;

	WaveletNode( dataType s ) {
		sum = s;
		sumOfSquares = pow( s, 2 );
		coeff = 0;
		level = 0;
		size=1;
	} WaveletNode() {
	}

	void addCount( uint state, uint count = 1 ) {
		if ( marginals.size() <= state ) {
			marginals.reserve( state+1 );
		}

		while ( marginals.size() <= state ) {
			marginals.push_back( 0 );
		}

		marginals[state] += count;
	}

};

class Block {
public:
	long index;
	long size;
	dataType coeff;		// the max of all absolute wavelet coefficients
	dataType sum;		/// this holds the original data in the leaves
	dataType sumOfSquares;

	Block( long index, long size, dataType sum, dataType sumOfSquares ) {
		this->index = index;
		this->size = size;
		this->sum = sum;
		this->sumOfSquares = sumOfSquares;
	}

	Block( long index, WaveletNode w ) {
		this->index = index;
		this->size = w.size;
		this->sum = w.sum;
		this->sumOfSquares = w.sumOfSquares;
	}
};


class WaveletTree {
private:

	unordered_map<int, SizeRange> ranges;

	long nrDataPoints;
	long nrLeaves; // next power of 2 above nrDataPoints NOTE there aren't actually that many leaves, since we prune the tree. This is what the number of leaves would be if we didn't.
	long nrTreeNodes; // 2*nrLeaves - 1
	double MADvariance; // the estimated variance from MAD of finest detail coefficients
	long iteration;	// a counter that is increased each time createBlocks is called
	vector < WaveletNode > data;
	uint ignoreLevels;	// how many levels to ignore (removes noise and makes blocks larger, but also decreases resolution; defaults to 0)
	vector<double> minimaxBlockMeans;	// the means of the blocks corresponding to the minimax segmentation

	double dataVariance;	// this is the maximum of the total variance of the data and the variance of the minimax block means
public:

	WaveletTree( vector < double >values, uint ignoreLevels = 0 ) ;
	void createBlocks( vector< Block >& blockSequence, double var, int max_block_len = UINT_MAX, bool debug = false );	// threshold is minimum variance in theta

	double getSigma1() {
		return data[0].sum;
	}
	double getSigma2() {
		return data[0].sumOfSquares;
	}

	double getMADvariance() {
		return MADvariance;
	}

	long getNrDataPoints() {
		return nrDataPoints;
	}

	long getNrTreeNodes() {
		return nrTreeNodes;
	}

	vector<WaveletNode> getAllNodes() {
		return data;
	}

	void setIgnoreLevels( uint ignoreLevels ) {
		this->ignoreLevels = ignoreLevels;
	}

	vector<Block> minimaxBlocks() ;
	void writeBlockStructure( vector<Block> blocks,  string filename, bool append = false );
	double getDataVariance( vector<Block> blocks, bool includeDataVariance = true );


	// debugging methods
	void allStructures( vector < Block >& blockSequence, int max_block_len );


	void addCount( Block& block,  uint state, uint count = 1 );
	
	void tallyCounts( vector<vector<long unsigned int>>& marginals );

};
#endif
