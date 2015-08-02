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



#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include <string>
#include <fstream>
#include <vector>
#include "WaveletTree.hpp"
#include <float.h>

class StateSequence;

class Distribution {
public:

	//name of distribution
	std::string name;

	//data stored as sufficent statistics
	bool sufficentStats;

	//data seen by this distribution
	std::vector<double> data;

	//previously sampled parameters
	std::vector<double> history;

	//obs is the observation, n is the dimension
	virtual double pdf( double* obs, size_t n ) = 0;
	
	virtual double logPdf( double* obs, size_t n ) {
		cout << name << "logPdf not implemented\n";
		return 0;
	}
	
	virtual double pdf( Block* obs, size_t n ) {
		cout << name << "pdf not implemented\n";
		return 0;
	}
	
	virtual double logPdf( Block* obs, size_t n ) {
		cout << name << "logPdf not implemented\n";
		return 0;
	}

	virtual void reserve_history( size_t number_of_iterations ) {
		cout << name << " reserve_history not implemented\n";
	}

	virtual double likelihoodExponent( Block& block ) {
		cout << name << "likelihoodExponent not implemented\n";
		return 0;
	}

	//samples parameters from conjugate prior
	virtual void sampleParameters() {
		cout << name << "sampleParameters not implemented\n";
		return;
	}

	//stores observation for posterior
	virtual void addObs( double* obs, size_t n ) {
		cout << name << "addObs not implemented";
		return;
	}

	virtual void addObs( Block* obs ) {
		cout << name << "addobs waveletnode not implemented";
		return;
	};

	//clears data for posterior
	virtual void clearObs() {
		cout << name << "clearobs not implemented";
		return;
	}

	// creates hyperparameters automatically from data
	virtual void autoHyperparams( double sigma1, double sigma2, size_t n, double sqrs ) {
		cout << name << "automaticHyperparameters() not implemented" << endl;
		return;
	}

	//header for output
	virtual std::string getHeader() {
		cout << name << "getHeader() not implemented" << endl;
		return "";
	}

	//saved parameters
	virtual void addHistory() {
		cout << name << "addHistory() not implemented" << endl;
		return;
	}


	virtual double getMean() {
		cout << name << "getMean not implemented" << endl;
		return 0;
	}

	//returns variance (max  double if not implemented)
	virtual double getVar() {
		cout << name << "getVar() not implemented" << endl;
		return DBL_MAX;
	}

	virtual double getLogVar() {
		cout << name << "getLogVar not implemented" << endl;
		return log( getVar() );
	}

	//writes saved parameters
	virtual void writeHistory( std::ostream& out ) {
		cout << name << "writeHistory not implemented";
		return;
	}

	virtual void printVars( std::ostream& out ) {
		out << name << endl;
		return;
	}

	void printObs( std::ostream& out ) {
		out << "observations for " << name << ": ";

		for ( size_t i = 0; i < data.size(); i++ ) {
			out << data[i] << " ";
		}

		out << endl;
	}

	virtual double pdf( WaveletNode* obs, size_t n, double exponent, double scale ) {
		cout << name << "pdf not implemented";
		return 0;
	}
	
	virtual double getScale( WaveletNode* obs, size_t n ) {
		cout << name << "getScale not implemented";
		return 0;
	}

};
#endif
