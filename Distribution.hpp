#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include "includes.hpp"
#include "Tags.hpp"


#include <random>
using std::default_random_engine;
using std::mt19937;
using std::normal_distribution;
using std::gamma_distribution;
using std::discrete_distribution;

typedef mt19937 rng_t;





template <typename DistType>
class Distribution {

		rng_t& mRNG;

	public:

		Distribution( rng_t& RNG ) : mRNG( RNG ) {};


		// sample new observation
		template< typename ParamType>
		Observation<DistType> sample(
		    const Observation< ParamType >& param ) {
			Observation<DistType> result;
			resample( result, param );
			return result;
		}


		// replace existing observation with a new sample
		template < typename ParamType>
		void resample(
		    Observation<DistType>& obs,
		    const Observation<ParamType>& param );

};





////////////////////////////////////////////////// TEMPLATE SPECIALIZATIONS //////////////////////////////////////////////////
#include "SufficientStatistics.hpp"	// required implementation
#include "Observation.hpp"	// required implementation



//////////////////// NORMAL ////////////////////

template<>	template<>
void Distribution<Normal>::resample(
    Observation<Normal>& obs,
    const Observation<NormalParam>& param ) {

	normal_distribution<real_t> dist( param.mean(),  param.stdev() );
	obs.setValue( dist( mRNG ) );
	dist.reset();
}



//////////////////// NORMAL INVERSE GAMMA ////////////////////


template<> template<>
void Distribution<NormalInverseGamma>::resample(
    Observation<NormalInverseGamma>& obs,
    const Observation<NormalInverseGammaParam>& param ) {
	gamma_distribution<real_t> gamma( param.alpha(), 1.0 / param.beta() );
	real_t var = 1.0 / gamma( mRNG );
	gamma.reset();
	normal_distribution<real_t> normal( param.mu0(),  sqrt( var / param.nu() ) );	
	real_t mean = normal( mRNG );
	normal.reset();
	obs.setValue( mean, var );
}



//////////////////// BETA ////////////////////

// Sample Beta distribution using two independent Gamma RV.
template<>	template<>
void Distribution<Beta>::resample(
    Observation<Beta>& obs,
    const Observation<BetaParam>& param ) {

	gamma_distribution<real_t> distA( param.alpha(),  1 );
	gamma_distribution<real_t> distB( param.beta(),  1 );
	real_t a = distA( mRNG );
	distA.reset();
	real_t b = distB( mRNG );
	distB.reset();
	obs.setValue( a / ( a + b ) );

}




//////////////////// DIRICHLET ////////////////////


// Sample Dirichlet distribution from normalized vector of Gamma RVs
void dirichlet_sample(
    vector<real_t>& probs,
    const vector<real_t>& alphas,
    rng_t& RNG ) {
	if ( probs.size() != alphas.size() ) {
		throw runtime_error( "Number of parameters must match the domain size of the Dirichlet RV!" );
	}
	size_t d = 0;
	real_t rand = 0;
	real_t sum = 0;
	for ( const auto & alpha : alphas ) {
		gamma_distribution<real_t> dist( alpha, 1.0 );
		rand = dist( RNG );
		dist.reset();
		probs[d] = rand;
		sum += rand;
		d++;
	}


	for ( auto & p : probs ) {
		p /= sum;
	}
}


template<> template<>
void Distribution<Dirichlet>::resample(
    Observation<Dirichlet>& obs,
    const Observation<DirichletParam>& param ) {


	if ( obs.domainSize() != param.domainSize() ) {
		throw runtime_error( "Domain sizes of Dirichlet random variable (" + to_string( obs.domainSize() ) + ") and the parameters requested for sampling (" + to_string( param.domainSize() ) + ") do not match!" );
	}

	dirichlet_sample( obs.probs(), param.alphas(), mRNG );
	// Dirichlet can be sampled using gamma(alpha_i, 1) distributions, with subsequent normalization

}



//////////////////// DIRICHLET VECTOR ////////////////////

template<> template<>
void Distribution<DirichletVector>::resample(
    Observation<DirichletVector>& obs,
    const Observation<DirichletParamVector>& param ) {

	// NOTE this class is declared friend of Observation<Dirichlet> to set obs.mValues directly

	if ( obs.nrDim() != param.nrDim() ) {
		throw runtime_error( "Dimensions of Dirichlet random variable (" + to_string( obs.nrDim() ) + ") and the parameters requested for sampling (" + to_string( param.nrDim() ) + ") do not match!" );
	}

	for ( size_t d = 0; d < obs.nrDim(); ++d ) {
		if ( obs[d].domainSize() != param[d].domainSize() ) {
			throw runtime_error( "Domain sizes of Dirichlet random variable (" + to_string( obs[d].domainSize() ) + ") and the parameters requested for sampling (" + to_string( param[d].domainSize() ) + ") do not match!" );
		}
		dirichlet_sample( obs[d].probs(), param[d].alphas(), mRNG );
	}
}



// TODO this interface differs from the others for a reason. Make this more elegant and consistent somehow.

template <>
class Distribution<Categorical> {

		rng_t& mRNG;

	public:

		Distribution( rng_t& RNG ) : mRNG( RNG ) {};


		// sample new observation
		int sample(
		    const vector<real_t>& probs,
		    const size_t begin,
		    const size_t end ) {
			discrete_distribution<int> dist( probs.begin() + begin, probs.begin() + end );
			int result = dist( mRNG );
			dist.reset();
			return result;
		}

};



#endif















