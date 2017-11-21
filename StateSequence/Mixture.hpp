#ifndef STATESEQUENCEMIXTURE_HPP
#define STATESEQUENCEMIXTURE_HPP

#include "../Tags.hpp"
// #include "StateSequence.hpp"
#include "../Emissions.hpp"
#include "../Theta.hpp"
#include "../Transitions.hpp"
#include "../Initial.hpp"

#include <vector>
using std::vector;


// TODO if nothing is to be recorded, storing the entire sequence is wasteful, since the posterior can be updated directly for each block



template<> template<typename EmissionsType,  typename ThetaType, typename TransitionsType, typename InitialType>
void StateSequence<Mixture>::sample(
    EmissionsType& y,	// TODO this cannot be const due to the use of next(), work around that somehow
    const ThetaType& theta,
    const TransitionsType& A,
    const InitialType& pi,
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
	y.initForward();
	while ( y.next() ) {
		real_t maxE = numeric_limits<real_t>::lowest();

		real_t N = y.blockSize();	// typecasting to avoid integer division TODO maxBlockSize should be restricted by range of real_t (data_t)

		// TODO assertions like in StateSequenceDirectGibbs
		for ( auto s = 0; s < nrStates; ++s ) {
			auto E = innerProduct( y, theta.value(), theta.mapping( s ) )  - N * logNormalizers[s];	// TODO carrier measure for the general EFD case
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
