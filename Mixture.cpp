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


#include "Mixture.hpp"
#include <vector>
#include <ostream>

double Mixture::pdf( double* obs, size_t n ) {
	double p = 0;

	for ( size_t i = 0; i < distributions.size(); i++ ) {
		p += components->pdf( i ) * distributions[i]->pdf( obs, 1 );
	}

	return p;
}

double Mixture::pdf( Block* obs, size_t n ) {
	double p = 0;

	for ( size_t i = 0; i < distributions.size(); i++ ) {
		p += components->pdf( i ) * distributions[i]->pdf( obs, 1 );
	}

	return p;
}

double Mixture::logPdf( double* obs, size_t n ) {
	return log( pdf( obs, n ) );
}

void Mixture::clearObs() {
	for ( size_t i = 0; i < distributions.size(); i++ ) {
		distributions[i]->clearObs();
		components->clearObs();
	}
}

void Mixture::sampleParameters() {
	components->sampleParameters();
}

void Mixture::addObs( double* obs, size_t n ) {
	size_t com = sampleComponent( obs[0] );
	counts[com]++;
	distributions[ com ]->addObs( obs, n );
}

size_t Mixture::sampleComponent( double obs ) {
	double c[distributions.size()];
	Categorical tmp( c, distributions.size() );

	for ( size_t i = 0; i < distributions.size(); i++ ) {
		tmp.p[i] = components->pdf( i );
		tmp.p[i] *= distributions[i]->pdf( &obs, 1 );
		c[i] = 0;
	}

	return tmp.sample();
}

void Mixture::addHistory() {
	components->addHistory();

	for ( size_t i = 0; i < distributions.size(); i++ ) {
		distributions[i]->addHistory();
	}
}

std::string Mixture::getHeader() {
	return components->getHeader();
}
