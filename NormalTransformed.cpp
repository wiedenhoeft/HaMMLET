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


#include "NormalTransformed.hpp"

void NormalTransformed::sampleParameters() {
	return;
}

double NormalTransformed::pdf( double* obs, size_t n ) {
	double tmp = a * obs[0] + b;
	return norm->pdf( &tmp, n );
}


void NormalTransformed::addObs( double* obs, size_t n ) {
	double tmp = a * obs[0] + b;
	norm->addObs( &tmp, n );
}

void NormalTransformed::clearObs() {
	data.clear();
}

void NormalTransformed::addHistory() {
	history.push_back( norm->mean );
	history.push_back( norm->var );
}

std::string NormalTransformed::getHeader() {
	return "";
}
