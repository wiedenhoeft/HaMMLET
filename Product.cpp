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


#include "Product.hpp"
#include "string"
#include <iostream>

double Product::pdf( double* obs, size_t n ) {
	double p = 1;

	for ( size_t i = 0; i < num_factors; i++ ) {
		p *= factors[i]->pdf( obs + i, dimension );
	}

	return p;
}

void Product::addObs( double* obs, size_t n ) {
	for ( size_t i = 0; i < num_factors; i++ ) {
		factors[i]->addObs( obs + i, dimension );
	}
}

void Product::clearObs() {
	for ( size_t i = 0; i < num_factors; i++ ) {
		factors[i]->clearObs();
	}

	return;
}
void Product::addHistory() {
	return;
}
