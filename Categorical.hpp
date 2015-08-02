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


#ifndef MULTINOMIAL_HPP
#define MULTINOMIAL_HPP

#include "Distribution.hpp"
#include <iostream>
#include "stdlib.h"
#include <string>
#include <fstream>
#include <vector>

class Categorical {
public:
	//hyperparameters
	std::vector<double> counts; //alphas for Dirichlet, prior counts
	std::string name;

	//parameters
	std::vector<double> p;
	std::vector<double> cdf;

	//used for speed up
	std::vector<double> log_p;

	int ncat;	//number of categories
	int seed;

	//history
	std::vector<double> history;
	std::vector<double> obs;	// observed counts

	Categorical() {};

	
	Categorical( double* priorCounts, size_t ncat ) {
		this->ncat = ncat;
		this->counts.resize( ncat );
		this->obs.resize( ncat );
		this->p.resize( ncat );
		this->log_p.resize( ncat );


		for ( size_t i = 0; i < ncat; i++ ) {
			counts[i] = priorCounts[i];
		}

		this->seed = 1;
	}

	
	Categorical( double* priorCounts, size_t ncat, std::string name ) {
		this->seed = 1;
		this->ncat = ncat;
		this->counts.resize( ncat );
		this->obs.resize( ncat );
		this->p.resize( ncat );
		this->log_p.resize( ncat );
		this->cdf.resize( ncat );

		for ( size_t i = 0; i < ncat; i++ ) {
			counts[i] = priorCounts[i];
		}

		this->name = name;
	}

	void addObs( size_t i, int num );
	void clearObs() ;
	void sampleParameters();
	size_t sample();
	double pdf( size_t cat );
	double logPdf( size_t cat );
	double pdf_power( size_t cat, size_t power );
	std::string getHeader();
	void addHistory();
	inline void isample( std::vector<double>& obs );
	void writeHistory( std::ofstream& os );
	void printVars( std::ostream& out );
	void printObs( std::ostream& out );

};
#endif
