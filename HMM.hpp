#ifndef HMM_HPP
#define HMM_HPP

#include "includes.hpp"
#include "Tags.hpp"
#include "ThetaHyperParam.hpp"
#include "Mapping.hpp"
#include "Emissions.hpp"
#include "Theta.hpp"
#include "Transitions.hpp"
#include "Initial.hpp"
#include "StateSequence.hpp"
#include "InitialHyperParam.hpp"
#include "TransitionHyperParam.hpp"
#include "StateMarginals.hpp"
#include "Records.hpp"

template < typename StateSequenceType,
         typename EmissionsType, // e.g. Normal
         typename ThetaType,	// e.g. NormalInverseGammaVector
         typename ThetaParamType,	// e.g. NormalInverseGammaParamVector
         typename TransitionType, // e.g. DirichletVector
         typename TransitionParamType,	// e.g. DirichletParamVector
         typename InitialType,	// e.g. Dirichlet
         typename InitialParamType
         >
void sampleHMM(
    EmissionsType& y,
    ThetaParamType& tau_theta,
    ThetaType& theta,
    TransitionParamType& tau_A,
    TransitionType& A,
    InitialParamType& tau_pi,
    InitialType& pi,
    StateSequenceType& q,
    const bool dynamic = true,
    const bool useSelfTransitions = true
) {

	if ( dynamic ) {
		y.createBlocks( theta );
	}
	q.sample( y, theta, A, pi, useSelfTransitions );
	theta.sample( q, y, tau_theta );
	pi.sample( q, A, y, tau_pi );
	A.sample( q, y, pi, tau_A );
}







template < typename StateSequenceType,
         typename EmissionsType, // e.g. Normal
         typename ThetaType,	// e.g. NormalInverseGammaVector
         typename ThetaParamType,	// e.g. NormalInverseGammaParamVector
         typename TransitionType, // e.g. DirichletVector
         typename TransitionParamType,	// e.g. DirichletParamVector
         typename InitialType,	// e.g. Dirichlet
         typename InitialParamType
         >
void sampleHMM(
    EmissionsType& y,
    ThetaParamType& tau_theta,
    ThetaType& theta,
    TransitionParamType& tau_A,
    TransitionType& A,
    InitialParamType& tau_pi,
    InitialType& pi,
    StateSequenceType& q,
    const Mapping& mapping,
    // insert records for parameters etc.
    const size_t iterations,
    const size_t thinning,
    Records& records,
    const bool dynamic = true,
    bool samplePrior = true,
    const bool useSelfTransitions = true
                                    // TODO RNG
) {


	if ( iterations < 0 ) {
		throw runtime_error( "Number of iterations must not be negative!" );
	}

	if ( thinning > iterations ) {
		cout << "[WARNING] Thinning parameter is larger than number of iterations. No data will be recorded!" << endl;
	}

	const size_t nrStates = mapping.nrStates();
	const size_t nrParams = mapping.nrParams();
	const size_t nrDataDim = mapping.nrDataDims();


	// TODO size checks go here

	// sample priors
	if ( samplePrior ) {
		theta.sample( tau_theta );
		pi.sample( tau_pi );
		A.sample( tau_A );
	}

	for ( auto i = 0; i <  iterations; ++i ) {
		sampleHMM( y, tau_theta, theta, tau_A, A, tau_pi, pi, q, dynamic, useSelfTransitions );

		if ( thinning > 0 ) {
			if ( ( i + 1 ) % thinning == 0 ) {
				records.record( q, y, theta );
			}
		}
	}



}




#endif





