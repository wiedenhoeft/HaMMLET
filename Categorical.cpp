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


#include "Categorical.hpp"
#include "asa266.hpp"
#include "pdflib.hpp"
#include <cmath>

void Categorical::addHistory() {
	for ( size_t i = 0; i < ncat; i++ ) {
		history.push_back( p[i] );
	}
}

std::string Categorical::getHeader() {
	std::string s;

	for ( size_t i = 0; i < ncat - 1; i++ ) {
		stringstream ss;
		ss << i;
		s += ss.str();
		s += " ";
	}

	stringstream ss;
	ss << ncat - 1;
	s += ss.str();
	return s + "\n";
}

void Categorical::writeHistory( std::ofstream& os ) {
	os << getHeader();

	for ( size_t i = 0; i < history.size(); i++ ) {
		os << history[i];

		if ( i % ncat == ncat - 1 && i != 0 ) {
			os << std::endl;
		} else {
			os << " ";
		}
	}
}

inline void Categorical::isample( std::vector<double>& obs ) {
	dirichlet_sample( ncat, &( obs[0] ), seed, &p[0] );
}


void Categorical::sampleParameters() {
	std::vector<double> tmp( ncat );

	for ( size_t i = 0; i < ncat; i++ ) {
		tmp[i] = obs[i] + counts[i];
	}

	// PDF
	isample( tmp );

	//log(PDF)
	log_p.resize( ncat );

	for ( int i = 0; i < ncat; i++ ) {
		log_p[i] = log( p[i] );
	}

	//CDF
	cdf[0] = p[0];

	for ( size_t i = 1; i < ncat; i++ ) {
		cdf[i] = cdf[i - 1] + p[i];
	}
}


size_t Categorical::sample() {

	double total = cdf[ncat - 1];
	double rn = r8_uniform_sample( 0.0, total );

	for ( size_t i = 0; i < ncat; i++ ) {
		if ( rn < cdf[i] ) {
			return i;
		}
	}

	cout << "[ERROR] Discrete sampling array malformed (probably all zeros or negative numbers)." << endl;
	return ncat - 1;
}



double Categorical::pdf( size_t cat ) {
	return p[cat];
}



double Categorical::logPdf( size_t cat ) {
	return log_p[cat];
}



void Categorical::printVars( std::ostream& out ) {
	out << name << ": ";

	for ( size_t i = 0; i < p.size(); i++ ) {
		out << p[i] << " ";
	}

	out << "   Priors ";

	for ( size_t i = 0; i < counts.size(); i++ ) {
		out << counts[i] << " ";
	}

	out << endl;
}



void Categorical::printObs( std::ostream& out ) {
	out << "printing observations for " << name << ": ";

	for ( size_t i = 0; i < obs.size(); i++ ) {
		out << obs[i] << " ";
	}

	out << endl;

}



void Categorical::addObs( size_t i, int num ) {
	obs[i] = obs[i] + num;
}


void Categorical::clearObs() {
	for ( size_t i = 0; i < ncat; i++ ) {
		obs[i] = 0;
	}
}

