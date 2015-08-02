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


#ifndef PRODUCT_HPP
#define PRODUCT_HPP

#include "Distribution.hpp"
#include <vector>

class Product: public Distribution {
public:
	size_t dimension;
	size_t num_factors;
	std::vector<Distribution*> factors;


	Product( std::vector<Distribution*>& factors, size_t num_factors, size_t dimension, std::string name ):
		factors( factors ), num_factors( num_factors ), dimension( dimension ) {
		this->name = name;
	};

	double pdf( double* obs, size_t n );
	double logPdf( double* obs, size_t n ) {
		return 0;   
	}
	double pdf( WaveletNode* obs, size_t n ) {
		return 0;   
	}
	double logPdf( WaveletNode* obs, size_t n ) {
		return 0;   
	}
	void sampleParameters() {};
	void addObs( double* obs, size_t n );
	void clearObs();
	std::string getHeader() {};
	void addHistory();
	double getVar() {
		return 0;   
	}
};

#endif
