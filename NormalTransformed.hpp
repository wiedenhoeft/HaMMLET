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


#ifndef NORMALTRANSFORMED_HPP
#define NORMALTRANSFORMED_HPP

#include "Distribution.hpp"
#include "Normal.hpp"
#include <string>

class NormalTransformed: public Distribution {
public:
	double a, b;
	Normal* norm;
	double pdf( double* obs, size_t n );
	double logPdf( double* obs, size_t n ) {
		return 0;   
	}
	double logPdf( WaveletNode* obs, size_t n ) {
		return 0;   
	}
	double pdf( WaveletNode* obs, size_t n ) {
		return 0;   
	}
	void sampleParameters();
	void addObs( double* obs, size_t n );
	void clearObs();
	void addHistory();
	void write( std::ostream& out ) {
		return;
	};
	std::string getHeader();

	NormalTransformed( Normal* norm, double a, double b, std::string name ): a( a ), b( b ), norm( norm ) {
		this->name = name;
	}
	double getVar() {
		return norm->getVar();
	}
};

#endif
