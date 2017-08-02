#ifndef INITIAL_HPP
#define INITIAL_HPP


#include "includes.hpp"
#include "Tags.hpp"
#include "Distribution.hpp"	
#include "InitialHyperParam.hpp"
#include "StateSequence.hpp"
#include "Transitions.hpp"

template <typename  DistType>	// e.g. Dirichlet
class Initial {

		Observation<DistType> mValue;
		SufficientStatistics< Categorical > mCounts;
		Distribution<DistType> mDist;

	public:

		// delete copy constructor
		Initial( const Initial& that ) = delete;

		Initial( size_t nrStates, rng_t& RNG ) :
			mValue( nrStates ),
			mCounts( nrStates ),
			mDist( RNG ) {}

		Initial( vector< real_t>& vec, rng_t& RNG ) :
			mValue( vec ),
			mDist( RNG ) {};


		template<typename InitialHyperParamType>
		void sample(
		    const InitialHyperParam<InitialHyperParamType>& tau_pi ) {
			mDist.resample( mValue,  tau_pi.posterior() );
		}

		template<typename StateSequenceType, typename TransitionsType, typename EmissionDataStructure, typename EmissionType, typename InitialHyperParamType>
		void sample(
		    const StateSequence<StateSequenceType>& q,
		    const Transitions<TransitionsType>& A,
		    Emissions<EmissionDataStructure, EmissionType>& y,	// NOTE pi does not statistically depend on y; y is only passed because it contains the block sizes which we need for self-transitions
		    InitialHyperParam<InitialHyperParamType>& tau_pi, // NOTE tau_pi cannot be constant since we update the parameters
		    bool considerFullSequence = true	) {	//TODO there are two ways to integrate the observed state sequence: count only the first state, or the complete sequence to get the stable distribution (considerFullSequence=true for the latter)

			mCounts.clear();

// 		TODO replace this functionality somehow
			if ( q.size() != y.size() ) {
				throw runtime_error( "Emissions and state sequence have different sizes!" );
			}


			size_t t = 0;

			if ( considerFullSequence ) {
				y.initForward();
				while ( y.next() ) {

					if ( t >= q.size() ) {
						throw runtime_error( "Emissions and state sequence have different sizes!" );
					}

					mCounts[q[t]] += y.N( );
					++t;
				}

				tau_pi.addObservation( mCounts );	// TODO more efficiently?
			} /*else {

			tau_pi.addObservation( q[0] );
			//TODO use A to sample backward from q[0]
		}*/

			mDist.resample( mValue, tau_pi.posterior() );

			tau_pi.reset();
		}



		vector<real_t> valueVector() const {	// NOTE this is intermediate level is necessary, since dist might be a more complicated structure than a simple probability vector itself, e.g. when using Dirichlet process priors
			return mValue.probs();
		};

		size_t nrStates() const {
			return mValue.domainSize();
		}


		string str() const {
			return mValue.str();
		}
};

#endif
