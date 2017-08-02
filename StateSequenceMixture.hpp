#ifndef STATESEQUENCEMIXTURE_HPP
#define STATESEQUENCEMIXTURE_HPP

#include "Tags.hpp"
#include "StateSequence.hpp"
#include "Emissions.hpp"
#include "Theta.hpp"
#include "Transitions.hpp"
#include "Initial.hpp"

#include <vector>
using std::vector;



template<> template<typename EmissionsDataStructure, typename EmissionsType, typename TransitionsType, typename ThetaType, typename InitialType>
void StateSequence<Mixture>::sample(
    Emissions<EmissionsDataStructure, EmissionsType>& Y,	// TODO this cannot be const due to the use of next(), work around that somehow
    const Theta<ThetaType>& theta,
    const Transitions<TransitionsType>& A,
    const Initial<InitialType>& pi,
    const bool useSelfTransitions	// NOTE this has no effect
) {
	size_t nrStates = A.nrStates();
	mStates.clear();


	//TODO precompute carrier measure upon implementation of non-normal distributions

	vector<real_t> logNormalizers;
	logNormalizers.reserve( nrStates );



	for ( auto s = 0; s < nrStates; ++s ) {
		logNormalizers.push_back( theta.logNormalizer( s ) );
	}

	vector<real_t> weights( nrStates, 0 );




	// forward variables
	Y.initForward();
	while ( Y.next() ) {
		real_t maxE = numeric_limits<real_t>::lowest();

		real_t N = Y.N();	// typecasting to avoid integer division TODO maxBlockSize should be restricted by range of real_t (data_t)

		// TODO assertions like in StateSequenceDirectGibbs
		for ( auto s = 0; s < nrStates; ++s ) {
			auto E = innerProduct( Y, theta.value(), theta.mapping( s ) )  - N * logNormalizers[s];	// TODO carrier measure for the general EFD case
			weights[s] = E ;
			maxE = max( E, maxE );
		}


		for ( auto s = 0; s < nrStates; ++s ) {
			weights[s] = exp( weights[s] - maxE );
		}
		discrete_distribution<size_t> dist( weights.begin(), weights.end() );
		mStates.push_back( dist( mRNG ) );
	}
}






#endif
