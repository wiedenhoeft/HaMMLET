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


#include "WaveletTree.hpp"
#include <sstream>
#include "float.h"



WaveletTree::WaveletTree( vector < double >values, uint ignoreLevels ) {	

	// set members to their initial values
	nrDataPoints = values.size();
	iteration = 0;
	ignoreLevels = ignoreLevels;

	//pad up to the next power of 2
	
	nrLeaves = 1;

	while ( nrLeaves < nrDataPoints ) {
		nrLeaves *= 2;
	}

	int cur = 0;

	while ( values.size() < nrLeaves ) {
		values.push_back( values[cur++] );
	}

	data.reserve( 2 * nrLeaves );

	// create the leaves at resolution level 0
	for ( vector < double >::iterator i = values.begin();
			i != values.end(); ++i ) {
		data.push_back( WaveletNode( *i ) );
	}

	//append a copy in front, which will become the actual wavelet tree
	for ( uint i = 0; i < nrLeaves; ++i ) {
		data.push_back( data[i] );
	}

	nrTreeNodes = data.size() - 1; // this is twice the number of leaves minus 1




	// Calculates the coefficients of the Haar wavelet transform on the first copy of the leaves in-place.

	dataType* temp = new dataType[nrLeaves];	

	
	for ( long i = 0; i < nrLeaves; i++ ) {
		temp[i] = 0;
	}

	// set the non-leaf coefficients to the sum in the data so we can calculate coefficients in-place
	for ( long i = 0; i < nrLeaves; ++i ) {
		data[i].coeff = data[i].sum;
	}


	// calculate the wavelet transform
	uint w = nrLeaves;

	while ( w > 1 ) {
		w /= 2;

		for ( long i = 0; i < w; i++ ) {
			temp[i] =
				( data[2 * i].coeff +
				  data[2 * i + 1].coeff ) / sqrt( 2.0 );
			temp[i + w] =
				( data[2 * i].coeff -
				  data[2 * i + 1].coeff ) / sqrt( 2.0 );
		}

		for ( long i = 0; i < ( w * 2 ); i++ ) {
			data[i].coeff = temp[i];
		}
	}

	delete[]temp;

	// remove first element (scaling coefficient), and shift all elements to the left. NOTE this shifts all elements, not just the first n.
	for ( long i = 0; i < nrTreeNodes; ++i ) {
		data[i] = data[i + 1];
	}

	data.pop_back();


	// remove data from the padded leaves
	// NOTE the level information must remain intact in order to transform the tree to DFS
	for ( long i = nrLeaves + nrDataPoints - 1; i < nrTreeNodes ; ++i ) {
		data[i].sum = 0;
		data[i].coeff = 0;
		data[i].sumOfSquares = 0;
	}


	// update parents from their children recursively (the array created from the wavelet transform is in BFS order); nodes at level 1 are added to finestDetailCoeffs in order to estimate the variance
	int p;

	long i = nrTreeNodes - 1; // index of the last leaf node/data point

	vector < double >finestDetailCoeffs;
	finestDetailCoeffs.reserve( nrLeaves / 2 );



	while ( i >= 2 ) {

		if ( i % 2 == 0 ) { // use right children
			p = ( i - 1 ) / 2;	// index of parent; NOTE since i is an integer, this will be floor((i-1)/2), which is what we want
			data[p].level = data[i].level + 1;	// step up one power of 2

			if ( data[p].level == 1  && finestDetailCoeffs.size() * 2 < nrDataPoints ) {
				finestDetailCoeffs.push_back(	data[p].coeff );
			}



			data[p].sum = data[i].sum + data[i - 1].sum;
			data[p].sumOfSquares = data[i].sumOfSquares + data[i - 1].sumOfSquares;

			if ( data[p].level > ignoreLevels + 1 ) {	
				data[p].coeff = max( abs( data[p].coeff ), max( abs( data[i].coeff ), abs( data[i - 1].coeff ) ) );
			} else {
				data[p].coeff = abs( data[p].coeff );

			}
		}

		--i;
	}


	// estimate the variance from the MAD of the finest detail coefficients
	long medianIndex = finestDetailCoeffs.size() / 2;
	nth_element( finestDetailCoeffs.begin(),
				 finestDetailCoeffs.begin() + medianIndex,
				 finestDetailCoeffs.end() );

	double median = * ( finestDetailCoeffs.begin() + medianIndex );

	for ( int i = 0; i < finestDetailCoeffs.size(); ++i ) {
		finestDetailCoeffs[i] = abs( median - finestDetailCoeffs[i] );
	}



	nth_element( finestDetailCoeffs.begin(),
				 finestDetailCoeffs.begin() + medianIndex,
				 finestDetailCoeffs.end() );
	MADvariance =
		pow( * ( finestDetailCoeffs.begin() + medianIndex ) / 0.6745, 2 );






	// change the tree from BFS to DFS order for optimal access pattern
	vector < WaveletNode > BFS = data;	// this deep copy represents the old BFS order
	uint height = data[0].level;
	vector < uint > currentIndexOnLevel;	// stores the highest index in BFS order encountered at each level
	currentIndexOnLevel.push_back( 0 );

	for ( uint i = 1; i <= height; ++i ) {
		currentIndexOnLevel.push_back( currentIndexOnLevel[i - 1] +
									   pow( 2, i - 1 ) );
	}

	reverse( currentIndexOnLevel.begin(), currentIndexOnLevel.end() );	// reverse such that the level corresponding to 2^i is at position i


	// iterate through the levels in DFS order using a stack. As in BFS, the nodes are in increasing order in each stack, we keep track of the latest visited node and increase its index after adding it to DFS.
	vector < uint > levelStack;
	levelStack.push_back( height );
	uint currentLevel;
	i = 0;

	while ( levelStack.size() > 0 ) {
		currentLevel = levelStack.back();
		levelStack.pop_back();

		if ( currentLevel > 0 ) {
			levelStack.push_back( currentLevel - 1 );
			levelStack.push_back( currentLevel - 1 );
		}

		data[i] = BFS[currentIndexOnLevel[currentLevel]];
		currentIndexOnLevel[currentLevel]++;
		i++;
	}

	

	// set block sizes, and remove all unnecessary nodes
	vector<long> currentSizeOnLevel;

	for ( uint i = 0; i <= height; ++i ) {
		currentSizeOnLevel.push_back( 0 );
	}

	for ( uint lastNode = 0; lastNode < data.size(); ++lastNode ) {

		currentLevel = data[lastNode].level;
		data[lastNode].size = pow( 2, currentLevel );	// level becomes block size N
		currentSizeOnLevel[currentLevel] += data[lastNode].size;

		if ( currentSizeOnLevel[currentLevel] >= nrDataPoints ) {
			data[lastNode].size -= ( currentSizeOnLevel[currentLevel] - nrDataPoints );

			if ( currentLevel == 0 ) {	// we are at the last node in the tree, and can remove everything after this
				data.resize( lastNode + 1 );
			}
		}
	}


	nrTreeNodes = data.size();

}



vector< Block > WaveletTree::minimaxBlocks() {
	// Creates the minimax blocks
	vector<Block> minimaxBlocks;
	createBlocks( minimaxBlocks, MADvariance, 0, false );
	--iteration;	// this method is not called as part of the regular iteration calls
	return minimaxBlocks;
}




void WaveletTree::writeBlockStructure( vector<Block> blocks,  string filename, bool append ) {	
	
	ofstream blocksFile( filename.c_str() );

	for ( vector<Block>::iterator i = blocks.begin(); i != blocks.end(); ++i ) {
		blocksFile << i->size << endl;
	}

	blocksFile.close();
}





double WaveletTree::getDataVariance( vector<Block> blocks, bool includeDataVariance ) {
	// Given a block structure, this calculates the variance of their means. If <includeDataVariance>, the maximum of the block variance and the actual data variance is returned. The typical application is to use minimax blocks to get a prior estimate of the variance for mu when facing class imbalance.

	vector<double> blockMeans;
	blockMeans.reserve( blocks.size() );

	for ( vector<Block>::iterator i = blocks.begin(); i != blocks.end(); ++i ) {
		blockMeans.push_back( i->sum / i->size );
	}

	double average = 0;

	for ( vector<double>::iterator i = blockMeans.begin(); i != blockMeans.end(); ++i ) {
		average += *i;
	}


	average /= blockMeans.size();

	dataVariance = 0;
	
	for ( vector<double>::iterator i = blockMeans.begin(); i != blockMeans.end(); ++i ) {
		dataVariance += pow( (*i) - average, 2 );
	}

	dataVariance /= blockMeans.size();

	if ( includeDataVariance ) {
		dataVariance = max( dataVariance, pow( data[0].level, 2 ) / ( data[0].level * data[0].sumOfSquares - pow( data[0].sum, 2 ) ) );
	}

	if ( dataVariance > 0 ) {
		cout << "Data variance is " << dataVariance << "." << endl;
	} else {
		cout << "[ERROR] Non-positive data variance! "<<dataVariance << endl;
	}

	return dataVariance;
}






void WaveletTree::createBlocks( vector < Block >& blockSequence,
								double var, int max_block_len, bool debug ) {


	double n = ( double ) nrLeaves ;
	double threshold = 0;
	uint nrValues = 0;		// this keeps track of how many values we have covered with the blocks, and where to stop and update N

	if ( debug ) {	// in debug mode, we use the threshold directly

		threshold = var;
		blockSequence.clear();
	} else {
		threshold = sqrt( 2 * log( n ) * min( var, MADvariance ) );	// universal threshold
	}


	// check if the current threshold would create the same block sequence as before
	if ( ranges[blockSequence.size()].contains( threshold ) ) {
		return;
	}

	blockSequence.clear();

	++iteration;
	int i = 0;
	bool createBlock;

	while ( i < nrTreeNodes ) {
		createBlock = false;

		if ( data[i].coeff <= threshold ) {
			createBlock = true;
		}

		if ( data[i].level == 0 ) {	// just in case of rounding errors
			createBlock = true;
		}

		if ( createBlock && ( data[i].size <= max_block_len || max_block_len == 0 ) ) {	// prune

			blockSequence.push_back( Block( i, data[i] ) );	

			nrValues += data[i].size;


			if ( nrValues == nrDataPoints ) {
				break;
			}

			if ( nrValues > nrDataPoints ) {
				printf( "[ERROR] The total block size is larger than the data, the wavelet tree must be malformed!" );

				return;
			}

			i += 2 * data[i].size - 1;
		} else {
			++i;
		}
	}


	if ( debug ) {	// write blocks in debugging mode

		std::stringstream ss;
		ss << var;
		ofstream blocksFile( ( "blocks_at_" + ss.str() + ".csv" ).c_str() );

		blocksFile << "# size\tsum\tsumOfSquares\tcoeff" << endl;

		for ( vector<Block>::iterator b = blockSequence.begin(); b != blockSequence.end(); ++b ) {
			blocksFile <<  b->size << "\t" << b->sum << "\t" << b->sumOfSquares << "\t" << b->coeff << endl;
		}

		blocksFile.close();
	}

	// add a new range for the block sequence if none exists, else update the existing one
	ranges[blockSequence.size()].update( threshold );

}




void WaveletTree::addCount( Block& block, uint state, uint count ) {
	data[block.index].addCount( state, count );
}



void WaveletTree::tallyCounts( vector< std::vector< long unsigned int > >& marginals ) {
	// traverses all nodes in the tree and updates the marginal counts
	// NOTE this method assumes that marginals has the correct dimensions!
	unsigned char maxLevel = data[0].level;
	vector<long> levelSize;
	levelSize.reserve( maxLevel + 1 );

	while ( levelSize.size() <= maxLevel ) {
		levelSize.push_back( 0 );
	}

	for ( vector<WaveletNode>::iterator n = data.begin(); n != data.end(); ++n ) {
		uint segmentStart = levelSize[n->level];

		for ( uint s = 0; s < n->marginals.size(); ++s ) {
			if ( n->marginals[s] > 0 ) {
				for ( uint t = segmentStart; t < segmentStart + n->size; ++t ) {
					marginals[t][s] += n->marginals[s];
				}
			}
		}

		levelSize[n->level] += n->size;
	}
}



////////// DEBUGGING METHODS //////////

void WaveletTree::allStructures( vector < Block >& blockSequence,  int max_block_len ) {
	// output all block structures possible in the tree as long as their coefficients are in blocks of size >= blocksize

	vector<double> coefficients;

	for ( vector<WaveletNode>::iterator i = data.begin(); i != data.end(); ++i ) {
		if ( i->level >= 32 ) {	
			createBlocks( blockSequence, i->coeff, max_block_len, true );
		}
	}


}


