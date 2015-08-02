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


#ifndef READER_HPP
#define READER_HPP

#define DEBUG 0

#include <string>
#include "Initial.hpp"
#include "Theta.hpp"
#include "Transitions.hpp"
#include "Normal.hpp"
#include "NormalTransformed.hpp"
#include "Product.hpp"
#include "Uniform.hpp"
#include <fstream>
#include <sstream>
#include "StateSequence.hpp"


class Model {
	std::vector<string> split( const std::string& s, char delim ) {
		std::vector<string> tmp;
		std::stringstream ss( s );
		std::string item;

		while ( std::getline( ss, item, delim ) ) {
			tmp.push_back( item );
		}

		return tmp;
	}

public:
	//write parameters to file
	bool write_parameters;
	//contains variables we need to sample(move to model)
	std::vector<Distribution*> continuous_indepent;//TODO SPELLING
	std::vector<Categorical*> discrete_indepent;

	//use thinning and burn in to save space (look at autocovariance to pick thinning)
	int iterations_ran;
	int iterations_since_last_write;
	int thinning;
	int burnin;

	//model's rv
	Initial pi;
	Transitions A;
	Theta theta;

	void gibbSample( int iter, string filename, StateSequence& s, std::vector<double>& obs );
	int readModel( string filename );
	void autoHyperparams( double sigma1, double sigma2, int n, double sqrs );
	void sampleParameters();
	void writeHistory( std::string filename );
	void addHistory();
	void writeVars( std::ostream& out );
	void writeObs( std::ostream& out );
	void clearObs();
	string getHeader();
	void reserve_history( size_t iterations );

	Model( bool write_parameters, int thinning, int burnin ): write_parameters( write_parameters ) {
		this->thinning = thinning;
		this->burnin = burnin;
	}
};

#endif
