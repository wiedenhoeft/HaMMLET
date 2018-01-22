#ifndef THETA_HPP
#define THETA_HPP

#include "includes.hpp"
#include "Tags.hpp"
#include "StateSequence.hpp"
#include "ThetaHyperParam.hpp"
#include "Mapping.hpp"
#include "EFD.hpp"
#include "SufficientStatistics.hpp"
#include "KahanAggregator.hpp"

template <typename ParamType>
class Theta {

		const size_t mNrDataDim;
		vector<Observation<ParamType>> mParams;
		Mapping mMapping;
		Distribution<ParamType> mDist;

		// NOTE It would be possible to have different parameters depend on the same hyperparameters.

	public:

		// delete copy constructor
		Theta( const Theta& that ) = delete;


		template< typename HyperParamType>
		Theta(
		    ThetaHyperParam<HyperParamType>& tau_theta,	// initializers for parameters
		    const size_t nrdatadim,
		    const MappingType mappingType,
		    rng_t& RNG
		);


		template< typename HyperParamType>
		Theta(
		    ThetaHyperParam<HyperParamType>& tau_theta,	// initializers for parameters
		    const size_t nrdatadim,
		    const Mapping mapping,
		    rng_t& RNG
		);


		// TODO currently we compute log-normalizers for each parameter multiple times in the case of shared parameters. This could be memoized.
		real_t logNormalizer( size_t state ) const;


		const vector<Observation<ParamType>>& value() const;


		size_t nrDataDim() const;


		size_t nrStates() const;


		size_t nrParams() const;


		const vector<size_t>& mapping(
		    size_t state ) const ;


		template<typename ThetaParamType>
		void sample(
		    ThetaHyperParam<ThetaParamType>& tau_theta ) ;


		template<typename ThetaParamType>
		void sample(
		    ThetaParamType& tau_theta );


		// returns a distribution-specific value for threshold computation, such as the minimum variance for Gaussian emissions
		real_t thresholdValue() const;


		string str(
		    const string& sep = "\t",
		    const string& finalSep = "" ) const ;


		// TODO move outside?
		template <typename U>
		friend ostream& operator<<(
		    ostream& output,
		    const Theta<U>& D );

};



template <typename ParamType>
ostream& operator<<(
    ostream& output,
    const Theta<ParamType>& D )  {
	output << D.str( );
	return output;
};










template <typename ParamType>
template< typename HyperParamType>
Theta<ParamType>::Theta(
    ThetaHyperParam<HyperParamType>& tau_theta,	// initializers for parameters
    const size_t nrdatadim,
    const MappingType mappingType,
    rng_t& RNG
) :
	mNrDataDim( nrdatadim ),
	mParams( tau_theta.nrParams() ),
	mMapping( nrdatadim, tau_theta.nrParams(), mappingType ) ,
	mDist( RNG ) {

	// initialize by sampling from the prior
	sample( tau_theta );
}




template <typename ParamType>
template< typename HyperParamType>
Theta<ParamType>::Theta(
    ThetaHyperParam<HyperParamType>& tau_theta,	// initializers for parameters
    const size_t nrdatadim,
    const Mapping mapping,
    rng_t& RNG
) :
	mNrDataDim( nrdatadim ),
	mParams( tau_theta.nrParams() ),
	mMapping( mapping ) ,
	mDist( RNG ) {

	// initialize by sampling from the prior
	sample( tau_theta );

}

template <typename ParamType>
// TODO currently we compute log-normalizers for each parameter multiple times in the case of shared parameters. This could be memoized.
// TODO adapt for multivariate
real_t Theta<ParamType>::logNormalizer(
    size_t state
) const {
	real_t result = 0;

	for ( const auto & m : mMapping[state] ) {
		result += ::logNormalizer( mParams[m] );
	}

	return result;
};

////////// accessors //////////
//NOTE round parentheses access the data through the mapping, square brackets access the parameters directly

template <typename ParamType>
// TODO more elegant, e.g. implicit conversion?
const vector<Observation<ParamType>>& Theta<ParamType>::value() const {
	return mParams;
}


template <typename ParamType>
size_t Theta<ParamType>::nrDataDim() const {
	return mNrDataDim;
}

template <typename ParamType>
size_t Theta<ParamType>::nrStates() const {
	return mMapping.nrStates();
}


template <typename ParamType>
size_t Theta<ParamType>::nrParams() const {
	return mMapping.nrParams();
}


template <typename ParamType>
const vector<size_t>& Theta<ParamType>::mapping(
    size_t state ) const {
	return mMapping[state];
}


template <typename ParamType>
// sample each parameter from its posterior and reset the posterior afterwards
template<typename ThetaParamType>
void Theta<ParamType>::sample(
    ThetaHyperParam<ThetaParamType>& tau_theta ) {

	for ( auto d = 0; d < mParams.size(); ++d ) {
		mDist.resample( mParams[d], tau_theta.posterior( d ) );
	}
	tau_theta.reset();
	
}


template <typename ParamType>
string Theta<ParamType>::str(
    const string& sep,
    const string& finalSep ) const {
	return concat( mParams, sep, finalSep );
}






template <>
real_t Theta<NormalParam>::thresholdValue() const {

	// min
	real_t result = inf;
	for ( const auto & param : mParams ) {
		result = min( result, param.var() );
	}
	return result;

// TODO add option for averaging version for data with bad compression characteristics
// 	avg
// 	real_t result = 0;
// 	for ( const auto & param : mParams ) {
// 		result += param.var();
// 	}
// 	return result / mParams.size();

}


template <>
real_t Theta<Beta>::thresholdValue() const {
	real_t result = inf;

	for ( const auto & param : mParams ) {
		real_t p = param.value();
		result = min( result, ( 1 - p ) / ( p * p ) );
	}

	return result;
}


#endif



