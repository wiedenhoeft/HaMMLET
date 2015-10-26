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


#include <fstream>
#include <string>
#include <string.h>
#include "pdflib.hpp"
#include "rnglib.hpp"
#include <sstream>
#include "StateSequence.hpp"
#include "Categorical.hpp"
#include "stdlib.h"
#include "Initial.hpp"
#include "Transitions.hpp"
#include <cmath>
#include "Normal.hpp"
#include "NormalTransformed.hpp"
#include "Model.hpp"


void gibbs( string model_filename, int iterations, bool compression,
			int max_block_len, bool logs, bool counts_only, bool write_params,
			int burn_in, int thinning, int skipLevels, char* out, bool debug, bool blocks, char* seq ) {

	//read observations
	std::vector<long long> pos;	// positions along chromosome (x-coordinate)	TODO
	std::vector<double> obs;	// measured values (y-coordinate)
	double y;
	long long x;
	std::ifstream obs_file( seq );

	if ( !obs_file.is_open() ) {
		std::cout << "[ERROR] Cannot read observations from file " << seq << endl;
		return;
	}

	while ( true ) {
		obs_file >> x;
		obs_file >> y;

		if ( obs_file.eof() ) {
			break;
		};

		pos.push_back( x );

		obs.push_back( y );
	}

	string tmp;

	while ( !obs_file.eof() ) {
		getline( obs_file, tmp );
	}

	Model model( write_params, thinning, burn_in );
	model.readModel( model_filename );


	// reserve memory for parameters if they are to be output
	if ( write_params ) {
		model.reserve_history( iterations );
	}


	//initialize state sequence
	StateSequence S( obs.size(), 1, model.A.A.size(), counts_only, logs, max_block_len, debug, compression );

	if ( debug ) {
		WaveletTree wt( obs, skipLevels );
		vector<Block> dummyblocks;
		wt.minimaxBlocks();
		wt.allStructures( dummyblocks, max_block_len );
	} else {
		S.setWaveletTree( obs, skipLevels );
	}



	vector<Block> minimaxBlocks = S.tree->minimaxBlocks()	;	// create minimax blocks (we've done this in setWaveletTree)


	//get hyperparameters from data
	
	if ( compression ) {
		model.autoHyperparams( S.tree->getSigma1(), S.tree->getSigma2(), obs.size(),
							   S.tree->getDataVariance( minimaxBlocks, true ) );

		if ( blocks ) {
			S.tree->writeBlockStructure( minimaxBlocks, string( out ) + "_blocks.csv" );
		}

	} else {
		double sigma = 0;
		double sigma2 = 0;

		for ( int i = 0; i < obs.size(); i++ ) {
			sigma += obs[i];
			sigma2 += obs[i] * obs[i];
		}

		model.autoHyperparams( sigma, sigma2, obs.size(), S.tree->getDataVariance( S.tree->minimaxBlocks(), true ) );
	}


	// free minimaxBlocks from memory
	vector<Block>().swap( minimaxBlocks );


	//run gibb sampling (NOTE out is passed so we could write the state sequence every iteration)
	model.gibbSample( iterations, out, S, obs );	// TODO spelling

#if DEBUG
	cout << "finished Gibbs sampling...\nwriting parameters now" << endl;
#endif

	//write parameters of all rvs
	model.writeHistory( out );
}

int main( int argc, char** argv ) {
	if (argc < 15){
		cout <<endl<< "Insufficient number of arguments. This executable should not"<<endl<<"be called directly, but using hammlet.py. To learn more, run"<<endl<<endl<<"                    python hammlet.py -h"<<endl<<endl<<"to see the manual."<<endl;
		return 1;
	}
	gibbs(
		argv[1],
		atoi( argv[2] ),
		atoi( argv[3] ),
		atoi( argv[4] ),
		atoi( argv[5] ),
		atoi( argv[6] ),
		atoi( argv[7] ),
		atoi( argv[8] ),
		atoi( argv[9] ),
		atoi( argv[10] ),
		argv[11],
		atoi( argv[12] ),
		atoi( argv[13] ),
		argv[14] );	
	return 0;
}

