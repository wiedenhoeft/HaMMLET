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

// template < typename StateSequenceType,
//          typename EmissionsType, // e.g. Normal
//          typename ThetaType,	// e.g. NormalInverseGammaVector
//          typename ThetaParamType,	// e.g. NormalInverseGammaParamVector
//          typename TransitionType, // e.g. DirichletVector
//          typename TransitionParamType,	// e.g. DirichletParamVector
//          typename InitialType,	// e.g. Dirichlet
//          typename InitialParamType
//          >
// void sampleHMM(
//     EmissionsType& y,
//     StateSequenceType& q,
//     ThetaType& theta,
//     ThetaParamType& tau_theta,
//     TransitionType& A,
//     TransitionParamType& tau_A,
//     InitialType& pi,
//     InitialParamType& tau_pi,
//     const Mapping& mapping,
//     Records& records,
//     const bool doRecord,
//     const bool dynamic = true,
//     const bool useSelfTransitions = true
// ) {
//
// }
//
//
//




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
    StateSequenceType& q,
    ThetaType& theta,
    ThetaParamType& tau_theta,
    TransitionType& A,
    TransitionParamType& tau_A,
    InitialType& pi,
    InitialParamType& tau_pi,
    const Mapping& mapping,
    // insert records for parameters etc.
    const size_t iterations,
    const size_t thinning,
    Records& records,
    const bool dynamic = true,
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
	const size_t nrDataDim = mapping.nrDataDim();


	// TODO size checks go here

	// sample priors


	bool doRecord;
	for ( auto i = 0; i <  iterations; ++i ) {
		if ( dynamic ) {
			y.createBlocks( theta );
		}

		doRecord = false;
		if ( thinning > 0 ) {
			doRecord = ( ( i + 1 ) % thinning == 0 );
		}
		
		
		q.sample( y, theta, tau_theta, A, tau_A, pi, tau_pi, mapping, records, doRecord, useSelfTransitions );
		theta.sample( tau_theta );
		
		pi.sample( tau_pi ); //, records, doRecord );
		
		A.sample( tau_A ); //, records, doRecord );
		
		
		if ( doRecord ) {
			records.record( theta );
		}
	}



}




#endif





