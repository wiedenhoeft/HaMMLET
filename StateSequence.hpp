#ifndef STATESEQUENCE_HPP
#define STATESEQUENCE_HPP

#include "includes.hpp"
#include "Emissions.hpp"
#include "Theta.hpp"
#include "Transitions.hpp"	
#include "Initial.hpp"	
// #include "StateMarginals.hpp"
#include "Trellis.hpp"





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

		template<typename EmissionsDataStructure, typename EmissionsType, typename TransitionsType, typename ThetaType, typename InitialType>
		void sample(
		    Emissions< EmissionsDataStructure, EmissionsType >& Y,	// TODO cannot be const due to next()
		    const Theta< ThetaType >& theta,
		    const Transitions< TransitionsType >& A,
		    const Initial< InitialType >& pi,
		    const bool useSelfTransitions = true
		);

		size_t size() const {
			return mStates.size();
		}

		const vector<marginal_t>& states() const  {
			return mStates;
		}

		const marginal_t& operator[](
		    const size_t s	) const  {
			if ( s >= mStates.size() ) {
				throw runtime_error( "State sequence index "+to_string(s)+" out of bounds!" );
			}
			return mStates[s];
		}

		marginal_t operator[](
		    const size_t s	) {
			if ( s >= mStates.size() ) {
				throw runtime_error( "State sequence index "+to_string(s)+" out of bounds!" );
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

};





////////////////////////////////////////////////// TEMPLATE SPECIALIZATIONS //////////////////////////////////////////////////


#include "StateSequenceForwardBackward.hpp"
#include "StateSequenceMixture.hpp"
// #include "StateSequenceDirectGibbs.hpp"	


















#endif



