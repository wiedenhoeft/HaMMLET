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


#ifndef NORMAL_HPP
#define NORMAL_HPP

#include "string"
#include "Distribution.hpp"
#include <iostream>

class Normal: public Distribution {
	void sampleParametersKnownMean();
public:

	//params
	double mean, var;
	double log_var;
	double sd; //used t
	bool knownMean;//Bool

	double emissionTerm;	// this is the precomputed term in block emissions, i.e. log(sigma) + mu^2/(2 sigma^2)

	//autogen params
	bool autogen_params;

	//auto_s is s^2 in p(sigma^2 < s^2)
	double auto_s, auto_p;

	//hyperparams
	double p_mu0, p_nu, p_alpha, p_beta;

	Normal( double p_mu0, double p_nu, double p_alpha, double p_beta, std::string name ):
		p_mu0( p_mu0 ), p_nu( p_nu ), p_alpha( p_alpha ), p_beta( p_beta ) {
		autogen_params = false;
		knownMean = false;
		mean = 0;
		var = 0;
		sd = 0;
		this->name = name;
		sufficentStats = false;
		sampleParameters();
	}

	Normal( double auto_s, double auto_p, std::string name ): auto_s( auto_s ), auto_p( auto_p ) {	
		autogen_params = true;
		knownMean = false;
		mean = 0;
		var = 0;
		sd = 0;
		p_mu0 = 0;
		p_nu = 0;
		p_alpha = 0;
		p_beta = 0;
		sufficentStats = false;
		this->name = name;
		//cant do autgen until read data
	}

	virtual double likelihoodExponent( Block& block );

	double pdf( double* obs, size_t n );
	void sampleParameters();
	void addObs( double* obs, size_t n );
	void clearObs();
	void addHistory();
	std::string getHeader();
	void writeHistory( std::ostream& out );
	double getVar() {
		return var;
	}
	void printVars( std::ostream& out );
	void addObs( Block* obs );
	void autoHyperparams( double sigma1, double sigma2, size_t n, double sqrs );
	void reserve_history( size_t num );
	double getMean();


};
const double log_frac = 1.8378770664093453;//std::log( 2.0 * M_PI );
#endif
