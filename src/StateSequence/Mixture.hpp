#ifndef STATESEQUENCEMIXTURE_HPP
#define STATESEQUENCEMIXTURE_HPP

#include "../Statistics.hpp"
#include "../Tags.hpp"
#include "../SufficientStatistics.hpp"
#include "../Emissions.hpp"
#include "../Theta.hpp"
#include "../Transitions.hpp"
#include "../Initial.hpp"
#include "../KahanAggregator.hpp"
#include "../Statistics.hpp"

#include <vector>
using std::vector;




// sample the state sequence and record the necessary
template<> template <
typename StatsStructure,
         typename StatsType,
         typename BlocksType,
         typename ThetaType,
         typename TauThetaType,
         typename TransitionsType,
         typename TauAType,
         typename InitialType,
         typename TauPiType >
void StateSequence<Mixture>::sample(
    Emissions<Statistics<StatsStructure, StatsType>, Blocks<BlocksType>>& y,	// TODO this cannot be const due to the use of next(), work around that somehow
    const ThetaType& theta,
    TauThetaType& tau_theta,
    const TransitionsType& A,
    TauAType& tau_A,
    const InitialType& pi,
    TauPiType& tau_pi,
    const Mapping& mapping,
    Records& records,
    const bool doRecord,
    const bool useSelfTransitions	// NOTE this has no effect for mixtures
) {
	size_t nrStates = A.nrStates();
	const size_t nrDim = y.nrDim();
	const size_t nrParams = tau_theta.nrParams();
	// TODO size checks go here


	// NOTE initial state distribution and transitions have noeffect in mixture sampling

	//TODO precompute carrier measure upon implementation of non-normal distributions

	vector<real_t> logNormalizers;
	logNormalizers.reserve( nrStates );



	for ( auto s = 0; s < nrStates; ++s ) {
		logNormalizers.push_back( theta.logNormalizer( s ) );
	}

	vector<real_t> weights( nrStates, 0 );



	vector<KahanAggregator<SufficientStatistics<StatsType>>> stats;
	stats.resize( nrParams );



	// count transitions
	SufficientStatistics<CategoricalVector> transitions( nrStates );

	// count states
	SufficientStatistics<Categorical> stateCounts( nrStates );

	// forward variables
	y.initForward();

	// TODO initial
	size_t prevState = 0;
// 	const auto stateProbs = pi.valueVector();

// 	auto logStateProbs = pi.valueVector();
// 	for ( size_t d = 0; d < nrStates; ++d ) {
// 		logStateProbs[d] = log( logStateProbs[d] );
// 	}

	while ( y.next() ) {
		real_t maxE = numeric_limits<real_t>::lowest();

		const size_t N = y.blockSize();


		// TODO assertions like in StateSequenceDirectGibbs
		for ( auto s = 0; s < nrStates; ++s ) {
			auto E = innerProduct( y, theta.value(), theta.mapping( s ) )  - N * logNormalizers[s]; // + N * logStateProbs[s];	// TODO carrier measure for the general EFD case TODO this is mixture sampling for burn-in, it does not take state probabilities into account
			weights[s] = E ;
			maxE = max( E, maxE );
		}


		for ( auto s = 0; s < nrStates; ++s ) {
			
// 			weights[s] = pow( (double)exp( weights[s] - maxE ), 1.0/(double)N);
			weights[s] = exp( weights[s] - maxE );

		}

		discrete_distribution<size_t> dist( weights.begin(), weights.end() );
		const size_t state = dist( mRNG );

		stateCounts[state] += N;
		transitions[state][state] += N - 1;
		transitions[prevState][state] += 1;	// TODO Initial

		for ( auto d = 0; d < nrDim; ++d ) {
			// TODO assert range
			stats[mapping[state][d]].add( y.suffStat( d ), N );
		}


		if ( doRecord ) {
			records.record( state, N );
		}

		prevState = state;
	}

	for ( auto p = 0; p < nrParams; ++p ) {
		const size_t N = stats[p].nrTerms();
		if ( N > 0 ) {
			tau_theta.addObservation( stats[p].sum(), N,  p );
		}
	}

	tau_A.addObservation( transitions );
	tau_pi.addObservation( stateCounts );



	// NOTE no state sequence is recorded
}






#endif

