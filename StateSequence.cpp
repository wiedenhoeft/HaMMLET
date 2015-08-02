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


#include "StateSequence.hpp"
#include "pdflib.hpp"
#include <cmath>
#include "float.h"


size_t StateSequence::sampleCDF( double* dist, size_t N ) {
	double total = dist[N - 1];
	double rn = r8_uniform_sample( 0, total );

	if ( rn <= dist[0] )
		return 0;

	for ( size_t i = 1; i < N; i++ ) {
		if ( dist[i - 1] < rn && rn <= dist[i] )
			return i;
	}

	//if no prob just get random	
	cout << "[WARNING] No probabilities, using uniform sample." << endl;
	return ( size_t )r8_uniform_sample( 0, N - 1 );
}


void createCdf( double* cdf, double* p, size_t n ) {
	cdf[0] = p[0];

	for ( size_t i = 1; i < n; i++ ) {
		cdf[i] = cdf[i - 1] + p[i];
	}
}



void StateSequence::compressed_forward( Initial& pi, Transitions& A, Theta& theta ) {

	alpha.resize( blocks.size() );

	for ( size_t i = 0; i < alpha.size(); i++ ) {	
		alpha[i].resize( numStates );
	}

	// store log of self transitions
	vector<double> logSelfTrans;
	logSelfTrans.reserve( numStates );

	for ( size_t j = 0; j < numStates; j++ ) {
		logSelfTrans.push_back( A.logEntry( j, j ) );
	}

	//// first intensities ////
	double maxAlpha = -DBL_MAX;

	// compute unscaled exponent
	for ( size_t j = 0; j < numStates; j++ ) {

		alpha[0][j] = theta.likelihoodExponent( j, blocks[0] ) + ( blocks[0].size - 1 ) * logSelfTrans[j];

		maxAlpha = max( alpha[0][j], maxAlpha );
	}

	double sumOfAlphas = 0;	// sum of forward variables, for scaling

	// scale component by shifting largest exponent to 0 to avoid overflow of exp(), then calculate unscaled forward variables (these have an additional scaling factor of 1/sqrt(2pi)^N, and hence are not what you'd usually encounter as unscaled FV)
	for ( size_t j = 0; j < numStates; j++ ) {
		alpha[0][j] -= maxAlpha;
		double transitionTerm = 0;

		for ( size_t i = 0; i < numStates; i++ ) {
			transitionTerm += pi.pdf( i ) * A.entry( i, j );
		}

		alpha[0][j] = exp( alpha[0][j] ) * transitionTerm;
		sumOfAlphas += alpha[0][j];
	}

	// scale forward variables
	if ( sumOfAlphas > 0 ) {
		for ( size_t i = 0; i < numStates; i++ ) {
			alpha[0][i] /= sumOfAlphas;
		}

	} else {
		if ( sumOfAlphas == 0 ) {
			for ( size_t i = 0; i < numStates; i++ ) {
				alpha[0][i] = 1 / numStates;
			}
		} else {
			std::cout << "[ERROR] Sum of forward variables is negative (t=0)." << std::endl;
		}
	}


	//// rest of intensities ////
	for ( size_t t = 1; t < blocks.size(); t++ ) {
		maxAlpha = -DBL_MAX;

		// compute unscaled exponent
		for ( size_t j = 0; j < numStates; j++ ) {


			alpha[t][j] = theta.likelihoodExponent( j, blocks[t] ) + ( blocks[t].size - 1 ) * logSelfTrans[j];

			maxAlpha = max( alpha[t][j], maxAlpha );
		}

		sumOfAlphas = 0;	// sum of forward variables, for scaling

		// scale component by shifting largest exponent to 0 to avoid overflow of exp(), then calculate unscaled forward variables (these have an additional scaling factor of 1/sqrt(2pi)^N, and hence are not what you'd usually encounter as unscaled FV)
		for ( size_t j = 0; j < numStates; j++ ) {
			alpha[t][j] -= maxAlpha;
			double transitionTerm = 0;

			for ( size_t i = 0; i < numStates; i++ ) {
				transitionTerm += alpha[t - 1][i] * A.entry( i, j );
			}

			alpha[t][j] = exp( alpha[t][j] ) * transitionTerm;
			sumOfAlphas += alpha[t][j];
		}

		// scale forward variables
		if ( sumOfAlphas > 0 ) {
			for ( size_t i = 0; i < numStates; i++ ) {
				alpha[t][i] /= sumOfAlphas;


			}


		} else {
			if ( sumOfAlphas == 0 ) {
				for ( size_t i = 0; i < numStates; i++ ) {
					alpha[t][i] = 1 / numStates;
				}
			} else {
				std::cout << "[ERROR] Sum of forward variables is negative (t=0)." << std::endl;
			}

		}
	}
}


void StateSequence::forward( Initial& pi, Transitions& A, Theta& theta,
							 vector<double>& intensities ) {
	alpha.resize( intensities.size() / dimension );

	for ( size_t i = 0; i < alpha.size(); i++ ) {
		alpha[i].resize( numStates );
	}

	for ( size_t i = 0; i < alpha.size(); i++ ) {
		for ( size_t j = 0; j < alpha[i].size(); j++ ) {
			alpha[i][j] = 0;
		}
	}

	//first intensities
	double c = 0;

	for ( size_t i = 0; i < numStates; i++ ) {
		alpha[0][i] = pi.pdf( i ) * theta.pdf( &intensities[0], dimension, i );
		c += alpha[0][i];
	}

	if ( c > 0 ) {
		c = 1 / c;

		for ( size_t i = 0; i < numStates; i++ ) {
			alpha[0][i] *= c;
		}
	} else {
		std::cout << "bad1" << std::endl;
	}

	//rest of intensities
	for ( size_t i = 1; i < intensities.size() / dimension; i++ ) {
		c = 0;

		for ( size_t j = 0; j < numStates; j++ ) {
			alpha[i][j] = 0;
			double tmp = theta.pdf( &intensities[dimension * i], dimension,  j );

			for ( size_t k = 0; k < numStates; k++ ) {
				alpha[i][j] += alpha[i - 1][k] * A.entry( k, j ) * tmp ;
			}

			c += alpha[i][j];
		}

		if ( c > 0 ) {
			c = 1 / c;

			for ( size_t j = 0; j < numStates; j++ ) {
				alpha[i][j] *= c;
			}
		} else {
			std::cout << "bad2" << std::endl;
		}
	}
}




void StateSequence::backward( Initial& pi, Transitions& A, Theta& theta,
							  vector<double>& intensities ) {
	double cdf[numStates];
	cdf[0] = alpha[stateSequence.size() - 1][0];

	for ( size_t j = 1; j < numStates; j++ ) {
		cdf[j] = cdf[j - 1] + alpha[stateSequence.size() - 1][j];
	}

	stateSequence[stateSequence.size() - 1] = sampleCDF( cdf, numStates );

	for ( int i = stateSequence.size() - 2; i >= 0; i-- ) {
		alpha[i][0] = alpha[i][ 0 ] * A.entry( 0, stateSequence[i + 1] );

		for ( size_t j = 1; j < numStates; j++ ) {
			alpha[i][j] = alpha[i][ j ] * A.entry( j, stateSequence[i + 1] );
		}

		createCdf( cdf, &alpha[i][0], numStates );
		stateSequence[i] = sampleCDF( cdf, numStates );
	}
}




void StateSequence::compressed_backward( Initial& pi, Transitions& A, Theta& theta ) {
	stateSequence.clear();	
	stateSequence.resize( blocks.size() );

	double cdf[numStates];
	createCdf( cdf, &alpha[blocks.size() - 1][0], numStates );
	stateSequence[stateSequence.size() - 1] = sampleCDF( cdf, numStates );

	//for rest of blocks
	for ( int i = blocks.size() - 2; i >= 0; i-- ) {
		//add transition factor
		alpha[i][0] = alpha[i][ 0 ] * A.entry( 0, stateSequence[i + 1] );

		for ( size_t j = 1; j < numStates; j++ ) {
			alpha[i][j] = alpha[i][ j ] * A.entry( j,  stateSequence[i + 1] );
		}

		createCdf( cdf, &alpha[i][0], numStates );
		stateSequence[i] = sampleCDF( cdf, numStates );
	}
}




void StateSequence::setWaveletTree( vector<double>& intensities, int ignoreLevels ) {

	this->tree = new WaveletTree( intensities, ignoreLevels );
	this->blocks = this->tree->minimaxBlocks();	

}



void StateSequence::sample( Initial& pi, Transitions& A, Theta& theta,
							vector<double>& intensities ) {


	if ( !compress ) {
		forward( pi, A, theta, intensities );
		backward( pi, A, theta, intensities );
	} else {
		compressed_sample( pi, A, theta );
	}

#if DEBUG
	cout << "Printing alpha (forward) variables)\n";

	for ( size_t i = 0 ; i < alpha.size(); i++ ) {
		cout << i << ": ";

		for ( size_t j = 0; j < alpha[i].size(); j++ ) {
			cout << alpha[i][j] << " ";
		}

		cout << endl;
	}

	cout << "end alphas\n";
#endif
}




void StateSequence::compressed_sample( Initial& pi, Transitions& A, Theta& theta ) {
	tree->createBlocks( blocks, theta.getMinVar(), max_blocks_len, false );
	compressed_forward( pi, A, theta );
	compressed_backward( pi, A, theta );

}




void StateSequence::write_plain_sequence( std::ostream& os ) {
	if ( compress ) {
		for ( size_t i = 0; i < stateSequence.size(); i++ ) {
			size_t num = blocks[i].size;

			for ( size_t j = 0; j < num; j++ ) {
				os << stateSequence[i] << " ";
			}
		}

		os  << "\n";
	} else {
		for ( size_t i = 0; i < stateSequence.size(); i++ ) {
			os << stateSequence[i] << " ";
		}

		os  << "\n";
	}
}




void StateSequence::addHistory( std::ostream& os, Theta& theta ) {


	if ( counts_only ) {
		update_counts();
	} else {
		write_plain_sequence( os );
	}
}



void StateSequence::update_counts() {
	// update the marginal state counts

	if ( compress ) {

		for ( size_t i = 0; i < this->stateSequence.size(); i++ ) {
			tree->addCount( blocks[i], stateSequence[i] );
		}


	} else {
		for ( size_t i = 0; i < this->stateSequence.size(); i++ ) {

			counts[i][stateSequence[i]]++;
		}
	}

}




void StateSequence::writeCountsHelper( std::ostream& os ) {
	// this method writes the marginal counts to file

	if ( compress ) {
		clear_counts();
		tree->tallyCounts( counts );
	}

	for ( size_t j = 0; j < numStates; j++ ) {
		for ( size_t i = 0; i < counts.size(); i++ ) {
			os << counts[i][j] << " ";
		}

		os << "\n";
	}

}




void StateSequence::clear_counts() {
	for ( size_t i = 0; i < this->counts.size(); i++ ) {
		for ( size_t j = 0; j < numStates; j++ ) {
			counts[i][j] = 0;
		}
	}
}

