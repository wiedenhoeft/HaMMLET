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

		const real_t& operator()(
		    const size_t from,
		    const size_t to ) const {
			return mValue( from, to );
		};

		real_t& operator()(
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
		}

		template<typename StateSequenceType, typename EmissionDataStructure, typename EmissionType, typename InitialType,  typename TransitionParamType>
		void sample(
		    const StateSequence<StateSequenceType>& q,
		    Emissions<EmissionDataStructure, EmissionType>& y,	// NOTE A does not stochastically depend on y; y is only passed because it contains the block sizes which we need for self-transitions
		    const Initial<InitialType>& pi,
		    TransitionHyperParam<TransitionParamType>& tau_A
// 		    , size_t ignoreBlockSize = 0	// ignore blocks <= ignoreBlockSize, to handle salt-and-pepper noise TODO implement this, efficiently, minimizing queries to next()
		) {	// NOTE tau_A cannot be const since we update the parameters

			// tau_A.addObservation( pi.value(), q[0] );


			if ( q.size() != y.size() ) {
				throw runtime_error( "Sizes of state sequence (" + to_string( q.size() ) + ") and emissions (" + to_string( y.size() ) + ") do not match!" );
			}

			auto T = q.size();

			mCounts.clear();

			// add between-block transitions
			for ( auto t = 0; t < T - 1; ++t ) {
				mCounts[ q[t] ][ q[t + 1] ] += 1;
			}

			// add within-block transitions
			y.initForward();	// move to front
			size_t t = 0;
			while ( y.next() ) {
				mCounts[ q[t] ][ q[t] ] += ( y.N( ) - 1 );
				t++;
			}

			if ( t != T ) {
				throw runtime_error( "State sequence has unexpected length!" );
			}

			tau_A.addObservation( mCounts );
			mDist.resample( mValue,  tau_A.posterior() );
			tau_A.reset();
		}


};



// dummy specializations

template<> template<typename TransitionParamType>
void Transitions<Dummy>::sample(
    TransitionHyperParam<TransitionParamType>& tau_A ) {	// NOTE tau_A cannot be const since we update the parameters
}

template<> template< typename StateSequenceType, typename EmissionDataStructure, typename EmissionType, typename InitialType,  typename TransitionParamType>
void Transitions<Dummy>::sample(
    const StateSequence<StateSequenceType>& q,
    Emissions<EmissionDataStructure, EmissionType>& y,	// NOTE A does not stochastically depend on y; y is only passed because it contains the block sizes which we need for self-transitions
    const Initial<InitialType>& pi,
    TransitionHyperParam<TransitionParamType>& tau_A ) {
}


#endif

