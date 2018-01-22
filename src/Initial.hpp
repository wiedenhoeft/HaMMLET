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
		    InitialHyperParamType& tau_pi // NOTE tau_pi cannot be constant since we
		) {
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
