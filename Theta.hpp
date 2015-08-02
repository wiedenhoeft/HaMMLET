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


#ifndef THETA_HPP
#define THETA_HPP

#include <vector>
#include "Distribution.hpp"
#include "stdlib.h"
#include "asa266.hpp"
#include "Mixture.hpp"

class StateSequence;

class Theta {
public:
	size_t dimension;
	bool write_parameters;
	//contains a distribution for each state
	std::vector<Distribution*> emissions;

	Theta() {
		dimension = 1;
	}

	void addObs( std::vector<double>& obs, StateSequence& seq );
	double pdf( double* x, size_t n, size_t s );
	double pdf( WaveletNode* x, size_t n, size_t s );
	double getLogVar(size_t s);
	double getMinVar();
	
	double likelihoodExponent( size_t s, Block& block );
};
#endif
