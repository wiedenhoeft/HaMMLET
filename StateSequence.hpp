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


#ifndef STATESEQUENCE_H
#define STATESEQUENCE_H

#include <vector>
#include <iostream>
#include "Transitions.hpp"
#include "Initial.hpp"
#include "WaveletTree.hpp"
#include "Theta.hpp"
#include "Model.hpp"
#include <climits>

// class Label {
// public:
// 	int originalLabel;
// 	int newLabel;
// 	bool hasObservations;
// 	bool isAvailable;	// whether we have successfully tracked it (either direction)
// 	double mean;
// 
// 	Label( int ol, double m )	{
// 		originalLabel = ol;
// 		mean = m;
// 		hasObservations = false;
// 		isAvailable = true;
// 	}
// 
// 	friend bool operator<( const Label& a, const Label& b ) {
// 		// if one has observations and the other doesn't, sort the one that does first
// 		if ( a.hasObservations != b.hasObservations ) {
// 			if ( a.hasObservations ) {
// 				return true;
// 			} else {	// b.hasObservations
// 				return false;
// 			}
// 		}
// 
// 		if ( a.mean != b.mean ) {
// 			return a.mean < b.mean;
// 		} else {
// 			return a.originalLabel < b.originalLabel;
// 		}
// 
// 		return false;
// 
// 	}
// };
// 
// 
// class LabelDistance {
// public:
// 	double distance;
// 	int oldLabel;
// 	int newLabel;
// 
// 	LabelDistance( double dist, int oldL, int  newL ) {
// 		distance  = abs( dist );
// 		oldLabel = oldL;
// 		newLabel = newL;
// 	}
// 
// 	friend bool operator<( const LabelDistance& a, const LabelDistance& b ) {
// 		return a.distance < b.distance;
// 	}
// };




class StateSequence {
public:
	size_t dimension;
	//forward variables
	std::vector< std::vector<double> > alpha;

	//compressed variables
	std::vector<Block> blocks;
	int max_blocks_len;
	bool logs;
	bool debug;
	WaveletTree* tree;

	//counts (use this to save space and time but lose diagonostic info)
	bool counts_only;
	bool compress;
	std::vector< std::vector<size_t> > counts;	// these are the marginals




	//state sequence
	std::vector<size_t> stateSequence;	// this is the state sequence, in compressed sampling this has as many entries as the number of blocks
	size_t numStates;



	StateSequence( size_t len, size_t dimension, size_t numStates, bool counts_only, bool logs, size_t max_blocks_len, bool debug, bool compression ) {
		
		
		this->stateSequence.resize( len );

		for ( size_t i = 0; i < len; i++ ) {
			stateSequence[i] = 0;
		}

		this->debug = debug;
		this->numStates = numStates;
		this->dimension = dimension;
		this->counts_only = counts_only;
		this->max_blocks_len = max_blocks_len;
		this->logs = logs;
		this->compress = compression;
		this->tree = NULL;

		if ( counts_only ) {	
			this->counts.resize( len );

			for ( size_t i = 0 ; i < len ; i++ ) {
				counts[i].resize( numStates );
			}

			clear_counts();
		}
	}

	void update_counts();
	void clear_counts();
	void setWaveletTree( vector<double>& intensities, int ignoreLevels = 0 );

	void sample( Initial& pi, Transitions& A, Theta& theta,
				 vector<double>& intensities );
	void compressed_sample( Initial& pi, Transitions& A, Theta& theta );
	void forward( Initial& pi, Transitions& A, Theta& theta, vector<double>& intensities );
	void compressed_forward( Initial& pi, Transitions& A, Theta& theta );
	void backward( Initial& pi, Transitions& A, Theta& theta,  vector<double>& intensities );
	void compressed_backward( Initial& pi, Transitions& A, Theta& theta );
	


	size_t sampleCDF( double* dist, size_t N );
	size_t sampleLogCDF( double* dist, size_t N );

	/* Saves the current states Sequence. Not to save space an output stream
	 * is passed so we actually write instead if we are not using
	 * just count
	 */
	void addHistory( ostream& os, Theta& theta );

	void write_plain_sequence( std::ostream& os );

	void writeHistory( std::ostream& os ) {
		if ( counts_only ) {
			writeCountsHelper( os );
		} else {
			//currently writing every iteration instead of saving, so nothing to do here
			return;
		}
	}

	void writeCountsHelper( std::ostream& os );


};
#endif
