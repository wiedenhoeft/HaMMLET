#ifndef EFD_HPP
#define EFD_HPP

// Functionality related to exponential family distributions

#include "SufficientStatistics.hpp"
#include "Observation.hpp"


////////// Inner Products of parameters and sufficient statistics //////////

template<class SuffStatType, class ParamType>
real_t innerProduct(
    const SufficientStatistics<SuffStatType>& suffstat,
    const Observation<ParamType>& param );





////////// Normal //////////

real_t innerProduct(
    const SufficientStatistics<Normal>& suffstat,
    const Observation<NormalParam>& param
)  {
	real_t result = ( 2.0 * param.mean() * suffstat.sum() - suffstat.sumSq() ) / ( 2.0 * param.var( ) );
	if ( !isfinite( result ) ) {
		throw runtime_error( "Result of Normal inner product is not finite!" );
	}
	return result;
}


real_t logNormalizer(
    const Observation<NormalParam>& param ) {
	return log( param.stdev() ) + param.mean() * param.mean() / ( 2 * param.var() );
}


real_t sampleMean(
    const SufficientStatistics<Normal>& suffstat,
    size_t N	) {
	if ( N <= 0 ) {
		throw runtime_error( "Cannot calculate mean from zero observations!" );
	}
	double n = N;
	return suffstat.sum() / n;
}

real_t sampleVariance(
    const SufficientStatistics<Normal>& suffstat,
    size_t N	) {
	if ( N <= 0 ) {
		throw runtime_error( "Cannot calculate variance from zero observations!" );
	}
	double n = N;
	double avg = sampleMean( suffstat, N );
	return suffstat.sumSq() / n - ( avg * avg );
}



////////// Geometric distribution //////////

real_t innerProduct(
    const SufficientStatistics<Geometric>& suffstat,
    const Observation<Beta>& param )  {
	real_t result = suffstat.sum() * param.value();
	return result;
}


real_t logNormalizer(
    const Observation<Beta>& param )  {
	return log( param.value() );
}




// calculates the inner product in the PDF of an EFD between the current sufficient statistics and a set of parameters under a current mapping
template<class EmissionObject, class ParamType>
real_t innerProduct(
    const EmissionObject& y,
    const vector<Observation<ParamType>>& param,	// TODO ParamType in template definition?
    const vector<size_t>& mapping )  {
	real_t result = 0;
	for ( auto dim = 0; dim < y.nrDim(); dim++ ) {
		result += innerProduct( y.suffStat( dim ),  param[mapping[dim]] );	
	}
	return result;
}


template<class EmissionObject, class ParamType>
real_t innerProduct(
    const EmissionObject& y,
    const vector<Observation<ParamType>>& param	// TODO ParamType in template definition?
)  {
	real_t result = 0;
	for ( auto dim = 0; dim < y.nrDim(); dim++ ) {
		result += innerProduct( y.suffStat( dim ),  param[dim] );	
	}
	return result;
}









#endif


