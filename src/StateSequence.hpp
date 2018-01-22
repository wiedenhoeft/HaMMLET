#ifndef STATESEQUENCE_HPP
#define STATESEQUENCE_HPP

#include "includes.hpp"
#include "Emissions.hpp"
#include "Theta.hpp"
#include "Transitions.hpp"
#include "Initial.hpp"
#include"Blocks.hpp"
#include "Statistics.hpp"
#include "KahanAggregator.hpp"
#include "Trellis.hpp"
#include "Records.hpp"




template <typename Type>
class StateSequence {

		vector<marginal_t> mStates;
		rng_t& mRNG;
		vector<size_t> mPrevStateSequence;	// for direct Gibbs
		Trellis mTrellis;	// implementation as member avoids frequent allocations

	public:

		// delete copy constructor
		StateSequence( const StateSequence& that ) = delete;

		StateSequence( rng_t& RNG ) : mTrellis( RNG ), mRNG( RNG ) {};

// 		template<typename EmissionsType, typename ThetaType, typename TransitionsType, typename InitialType>
// 		void sample(
// 		    EmissionsType& Y,	// TODO cannot be const due to next()
// 		    const ThetaType& theta,
// 		    const TransitionsType& A,
// 		    const InitialType& pi,
// 		    const bool useSelfTransitions = true
// 		);



		template <
		typename StatsStructure,
		         typename StatsType,
		         typename BlocksType,
		         typename ThetaType,
		         typename TauThetaType,
		         typename TransitionsType,
		         typename TauAType,
		         typename InitialType,
		         typename TauPiType >
		void sample(
		    Emissions<Statistics<StatsStructure, StatsType>, Blocks<BlocksType>>& y,
		    const ThetaType& theta,
		    TauThetaType& tau_theta,
		    const TransitionsType& A,
		    TauAType& tau_A,
		    const InitialType& pi,
		    TauPiType& tau_pi,
		    const Mapping& mapping,
		    Records& records,
		    const bool doRecord,
		    const bool useSelfTransitions );


		size_t size() const {
			return mStates.size();
		}

		const vector<marginal_t>& states() const  {
			return mStates;
		}

		const marginal_t& operator[](
		    const size_t s	) const  {
			if ( s >= mStates.size() ) {
				throw runtime_error( "State sequence index " + to_string( s ) + " out of bounds!" );
			}
			return mStates[s];
		}

		marginal_t operator[](
		    const size_t s	) {
			if ( s >= mStates.size() ) {
				throw runtime_error( "State sequence index " + to_string( s ) + " out of bounds!" );
			}
			return mStates[s];
		}

		string str() const {
			stringstream ss;
			copy( mStates.begin(), mStates.end(), ostream_iterator<marginal_t>( ss, " " ) );
			string s = ss.str();
			s = s.substr( 0, s.length() - 1 );	// TODO inefficient, copies all but the last character
			return s;
		};

		void clear() {
			deleteVector( mStates );
			deleteVector( mPrevStateSequence );
			mTrellis.clear();
		}

};





////////////////////////////////////////////////// TEMPLATE SPECIALIZATIONS //////////////////////////////////////////////////


#include "StateSequence/ForwardBackward.hpp"
#include "StateSequence/Mixture.hpp"
// #include "StateSequence/DirectGibbs.hpp"


















#endif



