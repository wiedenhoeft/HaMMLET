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


#ifndef MIXTURE_HPP
#define MIXTURE_HPP

#include "Distribution.hpp"
#include "Categorical.hpp"
#include <string>

class Mixture: public Distribution {
public:
	//params
	size_t num;
	std::vector<Distribution*> distributions;
	Categorical* components;

	//stats
	std::vector< double > counts;

	Mixture( size_t num, std::vector<Distribution*> distributions, Categorical* components, std::string name ):
		num( num ), distributions( distributions ), components( components ) {
		counts.resize( num );

		for ( size_t i = 0; i < num; i++ ) {
			counts[i] = 0;
		}

		this->name = name;
	}

	double pdf( Block* obs, size_t n );
	double pdf( double* obs, size_t n );
	double logPdf( double* obs, size_t n );
	void sampleParameters();
	void addObs( double* obs, size_t n );
	void clearObs();
	void addHistory();
	size_t sampleComponent( double obs );
	std::string getHeader();
};

#endif
