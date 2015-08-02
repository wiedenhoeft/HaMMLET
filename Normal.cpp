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


#include "Normal.hpp"
#include "pdflib.hpp"
#include <assert.h>
#include <math.h>

void Normal::clearObs() {
	if ( !sufficentStats ) {
		data.clear();
	} else {
		assert( data.size() == 3 );
		data[0] = data[1] = data[2] = 0;
	}
}

double Normal::likelihoodExponent( Block& block ) {
	// this calculates the unscaled exponent of of the likelihood function, together with the terms for selftransition and scaling by the N-th power of the standard deviation
	double E = ( 2 * mean * block.sum - block.sumOfSquares ) / ( 2 * var ) - block.size * emissionTerm;
	return E;
}




void Normal::sampleParameters() {


	double moment = 0;
	double moment2 = 0;
	size_t len = 0;

	if ( sufficentStats ) {
		assert( data.size() == 3 );
		//calculate using stats
		len = data[0];

		if ( len != 0 ) {
			moment = data[1] / len;
			moment2 = data[2] - len * moment * moment;
		}

#if DEBUG
		else {
			cout << "no samples for " << name << endl;
		}

#endif

	}
	//not sufficent stats
	else {
		//calculate using obs
#if DEBUG
		if ( data.size() == 0 ) {
			std::cout << " no data for " + this->name << endl;
		}

#endif

		for ( size_t i = 0; i < data.size(); i++ ) {
			moment += data[i];
		}

		if ( data.size() > 0 ) {
			moment /= data.size();
		}

		for ( size_t i = 0; i < data.size(); i++ ) {
			moment2 += std::pow( data[i] - moment, 2 );
		}

		len = data.size();
	}

	
	double new_mu0, new_nu, new_alpha, new_beta, r;
	new_nu = p_nu + len;
	new_mu0 = ( p_nu * p_mu0 + len * moment ) / new_nu;
	new_alpha = p_alpha + len / 2;
	r = len * p_nu / new_nu;
	new_beta = p_beta + .5 * moment2 + r * ( moment - p_mu0 ) * ( moment - p_mu0 ) / 2;
	var =  r8_invgam_sample( new_beta, new_alpha );
	sd = sqrt( var );
	mean =  r8_normal_sample( new_mu0, sqrt( var / new_nu ) );

	// set the term used in the forward variables
	emissionTerm = log( sd ) + mean * mean / ( 2 * var );
}

double Normal::pdf( double* obs, size_t n ) {
	return r8_normal_pdf( mean, var, *obs );
}

void Normal::addObs( double* obs, size_t n ) {
	data.push_back( *obs );
}

void Normal::reserve_history( size_t num ) {
	history.reserve( num * 2 );

}

void Normal::addHistory() {
	history.push_back( mean );
	history.push_back( var );
}

void Normal::writeHistory( std::ostream& out ) {
	out << getHeader();

	for ( size_t i = 0; i < history.size() - 1; i += 2 ) {
		out << history[i] << " " << history[i + 1] << "\n";
	}
}

void Normal::printVars( std::ostream& out ) {
	out << name << " ";
	out << "Mean " << mean << " Var " << var << " mean prior " << p_mu0 << " mean obs "
		<< p_nu << " var obs " << p_alpha << " var prior " << p_beta << endl;
}

std::string Normal::getHeader() {
	return "Mean Variance\n";
}

void Normal::addObs( Block* obs ) {
	if ( !sufficentStats ) {
		data.resize( 3 );
		sufficentStats = true;
		data[0] = data[1] = data[2] = 0;
	}

	
	data[0] += obs->size;
	data[1] += obs->sum;
	data[2] += obs->sumOfSquares;
}

//if auto_s < 0 we set it to sqrs
void Normal::autoHyperparams( double sigma1, double sigma2, size_t n, double sqrs ) {
	if ( !autogen_params ) {
		return;
	}
	
	if ( auto_s < 0 ) {
		auto_s = sqrs;
	}

	p_alpha = 2;
	double b = -log( auto_p );
	double m1 = .3361;
	double m2 = -.0042;
	double m3 = -.0201;
	double tmp_num = 2 * sqrt( b );
	double tmp_denum = sqrt( b ) * m1 + sqrt( 2 ) * ( b * m2 * exp( sqrt( b ) * m3 ) + 1 );
	p_beta = auto_s * ( tmp_num / tmp_denum + b );
	p_mu0 = sigma1 / n;
	p_nu = p_beta *min(1/sqrs,  n * n / ( n * sigma2 - sigma1 * sigma1 ));
	sampleParameters();
}

double Normal::getMean() {
	return mean;
}

