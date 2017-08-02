#ifndef AUTOPRIORS_HPP
#define AUTOPRIORS_HPP

#include <stdexcept>
using std::runtime_error;

#include <vector>
using std::vector;

#include "includes.hpp"

#include "Emissions.hpp"
#include "SufficientStatistics.hpp"
#include "Tags.hpp"


// Normal, breakpoint array
vector<real_t> NormalInverseGammaAutoPrior (
    real_t s2 ,	// desired variance
    real_t p ,	// desired probability of sampling a variance below s^2
    real_t dataMean,
    real_t dataVar	// in wavelet tree autopriors, this is max(sample variance of data, variance of block means)
) {

    if ( p < 0 || p > 1 ) {
        throw runtime_error ( "Parameter p for automatic priors is a probability and must be in [0,1]!" );
    }

    if ( s2 <= 0 ) {
        throw runtime_error ( "Parameter s2  for automatic priors is a variance and must be positive!" );
    }

    if ( dataVar <= 0 ) {
        throw runtime_error ( "Data variance provided to autoprior must be positive!" );
    }


    const real_t M1 = 0.3361;
    const real_t M2 = -0.0042;
    const real_t M3 = -0.0201;

    const real_t b = -log ( p );

    const real_t alpha = 2.0;
    const real_t beta = s2 * ( ( 2.0 * sqrt ( b ) ) / ( M1 * sqrt ( b ) + sqrt ( 2.0 ) * ( M2 * b * exp ( M3 * sqrt ( b ) ) + 1 ) ) + b );
    const real_t mu0 = dataMean;
    const real_t nu = beta / dataVar;

    if ( alpha <= 0 ) {
        throw runtime_error ( "Autoprior yields non-positive alpha!" );
    }

    if ( beta <= 0 ) {
        throw runtime_error ( "Autoprior yields non-positive beta!" );
    }

    if ( nu <= 0 ) {
        throw runtime_error ( "Autoprior yields non-positive nu!" );
    }

    if ( !isfinite ( alpha ) ) {
        throw runtime_error ( "Autoprior yields non-finite alpha!" );
    }

    if ( !isfinite ( beta ) ) {
        throw runtime_error ( "Autoprior yields non-finite beta!" );
    }

    if ( !isfinite ( mu0 ) ) {
        throw runtime_error ( "Autoprior yields non-finite mu0!" );
    }

    if ( !isfinite ( nu ) ) {
        throw runtime_error ( "Autoprior yields non-finite nu!" );
    }

    vector<real_t> v {alpha, beta, mu0, nu};
    return v;

}


// auto prior for Gaussian breakpoint array
vector<real_t> autoPrior (
    real_t s2 ,	// desired variance
    real_t p ,	// desired probability of sampling variance below s2
    Emissions<BreakpointArray, Normal>& y ) {

    const real_t avgWeight = y.avgWeight();
    y.initForward();
    y.createBlocks ( avgWeight );
    SufficientStatistics<Normal> muStats;	// sufficient statistics of observed block means

    while ( y.next() ) {
        for ( size_t dim = 0; dim < y.nrDim(); ++dim ) {
            muStats.addObs ( y.suffStat ( dim ).sum() / y.N() );	// NOTE sometimes, computations involve dim sum, without increasing the average weight
        }
    }

    size_t N = y.size() * y.nrDim();
    return NormalInverseGammaAutoPrior ( s2, p, sampleMean ( muStats, N ), sampleVariance ( muStats, N ) );
}







#endif
