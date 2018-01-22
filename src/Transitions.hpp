#ifndef TRANSITIONS_HPP
#define TRANSITIONS_HPP

#include "includes.hpp"
#include "Tags.hpp"
#include "TransitionHyperParam.hpp"
#include "Distribution.hpp"
#include "Observation.hpp"
#include "StateSequence.hpp"



template <typename ParamType>
ostream& operator<<(
    ostream& output,
    const Transitions<ParamType>& D )  {
	output << D.str();
	return output;
}


template <typename DistType> // e.g. DirichletVector
class Transitions {

		size_t mNrStates;
		Distribution<DistType> mDist;
		Observation<DistType> mValue;
		SufficientStatistics<CategoricalVector> mCounts;	// the count matrix TODO are there any cases where this is not CategoricalVector?


	public:

		// delete copy constructor
		Transitions( const Transitions& that ) = delete;

		////////// constructors //////////

		Transitions( size_t nrStates,
		             rng_t& RNG
		           ) :
			mNrStates( nrStates ),
			mDist( RNG ),
			mValue( nrStates ),
			mCounts( nrStates ) {};

		inline const real_t& operator()(
		    const size_t from,
		    const size_t to ) const {
			return mValue( from, to );
		};

		inline real_t& operator()(
		    const size_t from,
		    const size_t to ) {
			return mValue( from, to );
		};




		//////////  const methods //////////

		size_t nrStates() const {
			return mNrStates;
		}

		string str() const {
			return mValue.str();
		}


		//////////  non-const methods //////////

		template<typename TransitionParamType>
		void sample(
		    TransitionHyperParam<TransitionParamType>& tau_A ) {	// NOTE tau_A cannot be const since we update the parameters
			mDist.resample( mValue, tau_A.posterior() );
			tau_A.reset();
		}

// 		template<typename TransitionParamType>
// 		void sample(
// 		    TransitionParamType& tau_A
// 		) {
// 			mDist.resample( mValue,  tau_A.posterior() );
// 			tau_A.reset();
// 		}


};



// dummy specializations

template<> template<typename TransitionParamType>
void Transitions<Dummy>::sample(
    TransitionHyperParam<TransitionParamType>& tau_A ) {	// NOTE tau_A cannot be const since we update the parameters
}
/*
template<> template< typename StateSequenceType,  typename EmissionType, typename InitialType,  typename TransitionParamType>
void Transitions<Dummy>::sample(
    const StateSequenceType& q,
    EmissionType& y,	// NOTE A does not stochastically depend on y; y is only passed because it contains the block sizes which we need for self-transitions
    const InitialType& pi,
    TransitionParamType& tau_A ) {
}*/


#endif

