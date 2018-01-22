#ifndef STATESEQUENCEFORWARDBACKWARD_HPP
#define STATESEQUENCEFORWARDBACKWARD_HPP



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
void StateSequence<ForwardBackward>::sample(
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
	const size_t nrParams = tau_theta.nrParams();
	const size_t nrDim = y.nrDim();
	const size_t nrStates = A.nrStates();
	mTrellis.clear();
	mTrellis.setNrStates( nrStates );



	//TODO precompute carrier measure upon implementation of non-normal distributions

	// precompute log of self-transitions etc.
	vector<real_t> logA;
	logA.reserve( nrStates );
	vector<real_t> logNormalizers;
	logNormalizers.reserve( nrStates );



	for ( auto s = 0; s < nrStates; ++s ) {
		if ( useSelfTransitions ) {
			logA.push_back( log( A( s, s ) ) );
		}
		logNormalizers.push_back( theta.logNormalizer( s ) );
	}


	// set the first element of the trellis to the initial state distribution
	// NOTE t-th index of the state sequence is t+1 in the trellis
	mTrellis.push_back( pi.valueVector() );

	size_t t = 0;	// TODO rename to tt?


	real_t prevN = 1;

	// FORWARD FILTERING
	y.initForward();
	vector<real_t> forward( nrStates, 0 );
	while ( y.next() ) {
		++t;
		real_t maxE = numeric_limits<real_t>::lowest();


		real_t N = y.blockSize();	// typecasting to avoid integer division TODO maxBlockSize should be restricted by range of real_t (data_t)

		for ( auto s = 0; s < nrStates; ++s ) {
			auto E = innerProduct( y, theta.value(), theta.mapping( s ) )  - N * logNormalizers[s];	// TODO carrier measure for the general EFD case
			if ( useSelfTransitions ) {
				E += ( N - 1 ) * logA[s];	// include self-transitions
			}
			forward[s] = E;
			maxE = max( E, maxE );
		}
		for ( auto s = 0; s < nrStates; ++s ) {
			forward[s]  = exp( forward[ s ] - maxE );
		}


		// calculate transition term and include in forward variables
		real_t forwardSum = 0;

		for ( auto j = 0; j < nrStates; ++j ) {
			real_t transitionTerm = 0;

			for ( auto i = 0; i < nrStates; ++i ) {
				transitionTerm += mTrellis( t - 1, i ) * A( i, j );
			}

			forward[ j ] *= transitionTerm;
			forwardSum += forward[j];
		}

		//normalize forward variables
		if ( forwardSum != 0 ) {
			for ( auto j = 0; j < nrStates; ++j ) {
				forward[ j ] /= forwardSum;

			}
		} else {
			cout << "[WARNING] Uniform sampling of forward variables!" << endl;
			for ( auto j = 0; j < nrStates; ++j ) {
				forward[ j ] = 1.0 / ( ( real_t )nrStates );
			}
		}


		if ( useSelfTransitions ) {	// we are done calculating the next forward variables. In the backward step, we have to scale them by their block size, which we can do now already.
			for ( auto s = 0; s < nrStates; ++s ) {
				mTrellis.back( s ) *= exp( ( prevN - 1 ) * logA[s] );
			}
		}

		mTrellis.push_back( forward );

		prevN = N;

	}

	size_t T = mTrellis.size() - 1;		// -1 because pi is in the trellis, but not part of the state sequence


	// BACKWARD SAMPLING
	// each forward variable is multiplied by the probability to transition into the sampled state

	mStates.resize( T );

	auto j =  mTrellis.sample( mTrellis.size() - 1 ) ;	// the sampled state
	mStates[T - 1] = j;

	real_t N = 0;

	for ( auto tt = T - 1; tt > 0; --tt ) {	// index in the trellis

		t = tt - 1;

		// update forward variable based on sampled state
		for ( auto i = 0; i < nrStates; ++i ) {
			mTrellis( tt, i ) = mTrellis( tt, i ) * A( i, j ) ;
			if ( mTrellis( tt, i ) < 0 ) {
				throw runtime_error( "Negative backward variable!" );
			}
		}

		// sample
		j = mTrellis.sample( tt );

		// set sampled state or initial value accordingly

		mStates[t] = j;

		// NOTE the value for the initial state distribution is NOT sampled by FB, but within pi itself, depending on the type of distribution.


	}

	mTrellis.clear();


	// POSTERIOR RECORDING
	// TODO is there any way to do the second iteration without recomputing stats?

	vector<KahanAggregator<SufficientStatistics<StatsType>>> stats;
	stats.resize( nrParams );
	SufficientStatistics<Categorical> stateCounts( nrStates );

	// count transitions
	SufficientStatistics<CategoricalVector> transitions( nrStates );

	size_t prevState = 0;
	t = 0;
	marginal_t state;
	y.initForward();
	
	while ( y.next() ) {
		N = y.blockSize();
		state = mStates[t];
		transitions[state][state] += N - 1;
		transitions[prevState][state] += 1;	// TODO Initial
		stateCounts[state] += N;

		for ( auto d = 0; d < nrDim; ++d ) {
			// TODO assert range
			stats[mapping[state][d]].add( y.suffStat( d ), N );
		}

		if ( doRecord ) {
			records.record( mStates[t], N );
		}
		
		prevState = state;
		t++;
	}
	

	for ( auto p = 0; p < nrParams; ++p ) {
		const size_t N = stats[p].nrTerms();
		if ( N > 0 ) {
			tau_theta.addObservation( stats[p].sum(), stats[p].nrTerms(),  p );
		}
	}

	tau_A.addObservation( transitions );
	tau_pi.addObservation( stateCounts );
	// TODO initial
}




#endif
