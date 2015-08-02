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


#include "Theta.hpp"
#include "StateSequence.hpp"
#include "pdflib.hpp"
#include <cmath>

double Theta::pdf( double* x, size_t n, size_t s ) {
	return emissions[s]->pdf( x, n );
}

void Theta::addObs( std::vector<double>& obs, StateSequence& seq ) {
	if ( !seq.compress ) {
		for ( size_t i = 0; i < obs.size() / dimension; i++ ) {
			emissions[seq.stateSequence[i]]->addObs( &( obs[dimension * i] ), dimension );
		}
	} else {
		for ( size_t i = 0; i < seq.stateSequence.size(); i++ ) {
			emissions[seq.stateSequence[i]]->addObs( &seq.blocks[i] );
		}

	}
}

double Theta::likelihoodExponent( size_t s, Block& block) {
	return emissions[s]->likelihoodExponent(block);
}


double Theta::getLogVar(size_t s){
	return emissions[s]->getLogVar();
}

double Theta::getMinVar() {
	double min = emissions[0]->getVar();
	double tmp;

	for ( size_t i = 0; i < emissions.size(); i++ ) {
		tmp = emissions[i]->getVar();

		if ( tmp < min ) {
			min = tmp;
		}
	}

	return min;
}
