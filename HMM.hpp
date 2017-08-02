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
// #include "Options.hpp"
#include "Records.hpp"

template <
typename StateSequenceType,
         typename EmissionDataStructure,	// e.g. WaveletTree
         typename EmissionDistType, // e.g. Normal
         typename ThetaDistType,	// e.g. NormalInverseGammaVector
         typename ThetaParamType,	// e.g. NormalInverseGammaParamVector
         typename TransitionDistType, // e.g. DirichletVector
         typename TransitionParamType,	// e.g. DirichletParamVector
         typename InitialDistType,	// e.g. Dirichlet
         typename InitialParamType 	// e.g. DirichletParam
         >
void sampleHMM(
    Emissions<EmissionDataStructure, EmissionDistType>& y,
    ThetaHyperParam<ThetaParamType>& tau_theta,
    Theta<ThetaDistType>& theta,
    TransitionHyperParam<TransitionParamType>& tau_A,
    Transitions<TransitionDistType>& A,
    InitialHyperParam<InitialParamType>& tau_pi,
    Initial<InitialDistType>& pi,
    StateSequence<StateSequenceType>& q,
    const bool useSelfTransitions = true
) {

	y.createBlocks( theta );
	q.sample( y, theta, A, pi, useSelfTransitions );
	theta.sample( q, y, tau_theta );
	pi.sample( q, A, y, tau_pi );
	A.sample( q, y, pi, tau_A );
}






template <
typename StateSequenceType,
         typename EmissionDataStructure,	// e.g. WaveletTree
         typename EmissionDistType, // e.g. Normal
         typename ThetaDistType,	// e.g. NormalInverseGammaVector
         typename ThetaParamType,	// e.g. NormalInverseGammaParamVector
         typename TransitionDistType, // e.g. DirichletVector
         typename TransitionParamType,	// e.g. DirichletParamVector
         typename InitialDistType,	// e.g. Dirichlet
         typename InitialParamType 	// e.g. DirichletParam
         >
void sampleHMM(
    Emissions<EmissionDataStructure, EmissionDistType>& y,
    ThetaHyperParam<ThetaParamType>& tau_theta,
    Theta<ThetaDistType>& theta,
    TransitionHyperParam<TransitionParamType>& tau_A,
    Transitions<TransitionDistType>& A,
    InitialHyperParam<InitialParamType>& tau_pi,
    Initial<InitialDistType>& pi,
    StateSequence<StateSequenceType>& q,
    const Mapping& mapping,
    // insert records for parameters etc.
    const size_t iterations,
    const size_t thinning,
    Records& records,
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
		sampleHMM( y, tau_theta, theta, tau_A, A, tau_pi, pi, q, useSelfTransitions );

		if ( thinning > 0 ) {
			if ( ( i + 1 ) % thinning == 0 ) {
				records.record( q, y, theta );
			}
		}
	}
}




#endif





