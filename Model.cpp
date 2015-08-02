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


#include "Model.hpp"
#include <map>

void Model::autoHyperparams( double sigma1, double sigma2, int n, double sqrs ) {
	for ( int i = 0; i < continuous_indepent.size(); i++ ) {
		continuous_indepent[i]->autoHyperparams( sigma1, sigma2, n, sqrs );
	}
}

void Model::reserve_history( size_t iterations ) {
	for ( int i = 0; i < continuous_indepent.size(); i++ ) {
		continuous_indepent[i]->reserve_history( iterations );
	}
}

int Model::readModel( string filename ) {
// #if DEBUG
	cout << "Reading model from " << filename << "...";
// #endif
	std::ifstream file( filename.c_str() );

	if ( !file.is_open() ) {
		std::cout << "ERROR cant open model file\n";
		return 1;
	}

	string line;
	std::map<string, Distribution*> distributions;
	std::map<string, Categorical*> discrete_distributions;

	//not optimized parsing
	while ( getline( file, line ) ) {
		line.erase( std::remove( line.begin(), line.end(), ' ' ), line.end() );
		std::vector<string> tmp = split( line, ',' );

		if ( tmp[0].compare( "Normal" ) == 0 ) {
			Normal* dist = new Normal( atof( tmp[1].c_str() ), atof( tmp[2].c_str() ), atof( tmp[3].c_str() ), atof( tmp[4].c_str() ), tmp[6] );
			dist->knownMean = atoi( tmp[5].c_str() );

			if ( dist->knownMean ) {
				dist->mean = atof( tmp[1].c_str() );
			}

			distributions[dist->name] = dist;
			continuous_indepent.push_back( dist );
		} else if ( tmp[0].compare( "NormalTransformed" ) == 0 ) {
			Normal* p_normal = ( Normal* )distributions[tmp[1]];

			if ( p_normal == NULL ) {
				std::cout << "[ERROR] parsing normal transformed\n";
				return 1    ;
			}

			NormalTransformed* dist = new NormalTransformed( p_normal, atof( tmp[2].c_str() ), atof( tmp[3].c_str() ), tmp[4] );
			distributions[dist->name] = dist;
		} else if ( tmp[0].compare( "AutoNormal" ) == 0 ) {
			Normal* dist = new Normal( atof( tmp[1].c_str() ), atof( tmp[2].c_str() ), tmp[3] );
			distributions[dist->name] = dist;
			continuous_indepent.push_back( dist );
		} else if ( tmp[0].compare( "Mixture" ) == 0 ) {
			int num_components = atoi( tmp[1].c_str() );
			Distribution* dist;
			std::vector<Distribution*> dists;
			Categorical* comps;
			comps = discrete_distributions[tmp[2]];

			if ( comps == NULL ) {
				std::cout << "[ERROR] parsing mixutre\n";
				return 1;
			}

			for ( int i = 0; i < num_components; i++ ) {
				dist = distributions[tmp[i + 3]];

				if ( dist == NULL ) {
					std::cout << "[ERROR] parsing mixture\n";
					return 1;
				}

				dists.push_back( dist );
			}

			Mixture* mix = new Mixture( num_components, dists, comps, tmp[3 + num_components] );
			distributions[mix->name] = mix;
			continuous_indepent.push_back( mix );

		} else if ( tmp[0].compare( "Categorical" ) == 0 ) {
			int ncat = atoi( tmp[1].c_str() );
			double t[ ncat ];

			for ( int i = 0; i < ncat; i++ ) {
				t[i] = atof( tmp[i + 2].c_str() );
			}

			Categorical* dist = new Categorical( t, ncat, tmp[ncat + 2] );
			discrete_indepent.push_back( dist );
			discrete_distributions[tmp[ncat + 2]] = dist;
		} else if ( tmp[0].compare( "Uniform" ) == 0 ) {
			int a = atoi( tmp[1].c_str() );
			int b = atoi( tmp[2].c_str() );
			Uniform* u = new Uniform( a, b );
			u->name = tmp[3].c_str();
			continuous_indepent.push_back( u );
			distributions[u->name] = u;
		} else if ( tmp[0].compare( "Product" ) == 0 ) {
			int dimension = atoi( tmp[1].c_str() );
			int num_factors = atoi( tmp[2].c_str() );
			std::vector<Distribution*> dists;

			for ( int i = 0; i < num_factors; i++ ) {
				dists.push_back( distributions[tmp[i + 3]] );
			}

			Product* prod = new Product( dists, num_factors, dimension, tmp[num_factors + 3] );
			distributions[prod->name] = prod;
		} else if ( tmp[0].compare( "HMM" ) == 0 ) {
			int numStates = atoi( tmp[1].c_str() );
			Categorical* cat = discrete_distributions[tmp[2]];

			if ( cat == NULL ) {
				std::cout << "[ERROR] parsing hmm\n";
				return 1;
			}

			pi.Pi = cat;

			for ( int i = 0; i < numStates; i++ ) {
				Categorical* cat = discrete_distributions[tmp[3 + i]];

				if ( cat == NULL ) {
					std::cout << "[ERROR] parsing hmm\n";
					return 1;
				}

				A.A.push_back( cat );
			}

			for ( int i = 0; i < numStates; i++ ) {
				Distribution* dist = distributions[ tmp[3 + numStates + i] ];

				if ( dist == NULL ) {
					std::cout << "[ERROR] parsing hmm\n";
					return 1 ;
				}

				theta.emissions.push_back( dist );
			}
		} else {
			std::cout << "[ERROR] parsing " << tmp[0] << std::endl;
			return 1;
		}
	}

	file.close();
	cout << "done." << endl;
	return 0;
}

void Model::gibbSample( int iter, string filename, StateSequence& s, std::vector<double>& obs ) {
	//output

#if DEBUG
	cout << "Starting Gibbs sampling" << endl;
#endif
	int iterations_since_last_write = 0;
	int iterations_ran = 0;

	ofstream out;
	string output_name = filename + "_statesequence.csv";
#if DEBUG
	cout << output_name << endl;
#endif
	out.open( output_name.c_str() );
	out << "# state sequence\n";

	//sample initial params
	clearObs();
	sampleParameters();

#if DEBUG
	cout << "writing initial sample params\n";
	writeVars( cout );
#endif


	int nrBlocks = 0;

	//fb gibbs

	clearObs();

	for ( int i = 0; i < iter; i++ ) {
		//sample state
		s.sample( pi, A, theta, obs );

		if ( s.compress ) {
			nrBlocks += s.blocks.size();
		} else {
			nrBlocks += obs.size();
		}

#if DEBUG
		std::cout << "Sampled state sequence\n";
		std::cout << s << "\n";
#endif

		//add observations to the distributions
		A.addObs( s );
		pi.addObs( s );
		theta.addObs( obs, s );



#if DEBUG
		writeObs( cout );
#endif

		//save history
		if ( iterations_ran >= burnin && iterations_since_last_write >= thinning - 1 ) {
			s.addHistory( out , theta );	// this has to be done before we sample the other parameters
			iterations_since_last_write = 0;
		}

		//sample params
		sampleParameters();
		clearObs();

#if DEBUG
		cout << "writing " << i << " sample params\n";
		writeVars( cout );
#endif


		addHistory();
		iterations_ran++;
		iterations_since_last_write++;

	}

	s.writeHistory( out );
	out.close();

	// write compression ratio if we are using compressed state sequence
	if ( s.compress ) {
#if DEBUG
		cout << nrBlocks << " blocks" << endl;
		cout << iter* obs.size() << " points" << endl;
#endif

		string output_namec = filename + "_compression_ratio.csv";
#if DEBUG
		cout << output_namec.c_str() << endl;
#endif
		float totalCompression = ( iter * obs.size() ) / ( ( float )nrBlocks );

		out.open( output_namec.c_str() );
		out << totalCompression << endl;
		out.close();
	}

	string output_name2 = filename + "_last_statesequence.csv";
#if DEBUG
	cout << output_name2.c_str() << endl;
#endif
	out.open( output_name2.c_str() );
	out << "state sequence" << endl;
	s.write_plain_sequence( out );
	out.close();




}


void Model::sampleParameters() {
	for ( size_t i = 0; i < continuous_indepent.size(); i++ ) {
		continuous_indepent[i]->sampleParameters();
	}

	for ( size_t i = 0; i < discrete_indepent.size(); i++ ) {
		discrete_indepent[i]->sampleParameters();
	}
}

void Model::addHistory() {
	if ( !write_parameters )return;

	for ( size_t i = 0; i < continuous_indepent.size(); i++ ) {
		continuous_indepent[i]->addHistory();
	}

	for ( size_t i = 0; i < discrete_indepent.size(); i++ ) {
		discrete_indepent[i]->addHistory();
	}
}

void Model::clearObs() {
	for ( size_t i = 0; i < continuous_indepent.size(); i++ ) {
		continuous_indepent[i]->clearObs();
	}

	for ( size_t i = 0; i < discrete_indepent.size(); i++ ) {
		discrete_indepent[i]->clearObs();
	}
}

std::string Model::getHeader() {
	string header;

	for ( size_t i = 0; i < discrete_indepent.size() - 1; i++ ) {
		header += discrete_indepent[i]->name + "\t";
	}

	header += discrete_indepent[discrete_indepent.size() - 1]->name;

	for ( size_t i = 0; i < continuous_indepent.size() - 1; i++ ) {
		header += continuous_indepent[i]->name + "\t";
	}

	header += continuous_indepent[continuous_indepent.size() - 1]->name;
	return header;
}


void Model::writeHistory( std::string filename ) {
	if ( !write_parameters )return;

	for ( size_t i = 0; i < continuous_indepent.size(); i++ ) {
		string fn = std::string( filename ) + "_" + continuous_indepent[i]->name + ".txt";
		std::ofstream output( fn.c_str() );
		continuous_indepent[i]->writeHistory( output );
		output.close();
	}

	for ( size_t i = 0; i < discrete_indepent.size(); i++ ) {
		std::ofstream output( ( std::string( filename ) +  "_" + discrete_indepent[i]->name + ".txt" ).c_str() );
		discrete_indepent[i]->writeHistory( output );
		output.close();
	}
}



//---------------debug===================================
void Model::writeVars( std::ostream& out ) {
	out << "Writing all paramters" << endl;

	for ( size_t i = 0; i < continuous_indepent.size(); i++ ) {
		continuous_indepent[i]->printVars( out );
	}

	for ( size_t i = 0; i < discrete_indepent.size(); i++ ) {
		discrete_indepent[i]->printVars( out );
	}
}
void Model::writeObs( std::ostream& out ) {
	out << "Writing all paramters" << endl;

	for ( size_t i = 0; i < continuous_indepent.size(); i++ ) {
		continuous_indepent[i]->printObs( out );
	}

	for ( size_t i = 0; i < discrete_indepent.size(); i++ ) {
		discrete_indepent[i]->printObs( out );
	}
}
