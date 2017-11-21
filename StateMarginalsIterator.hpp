#ifndef STATEMARGINALITERATOR_HPP
#define STATEMARGINALITERATOR_HPP

#include "includes.hpp"


template<typename T>
class StateMarginalIterator {

		// internal states of a finite state machine parsing the marginals
		enum fsa_state {start, afterPos, afterState, afterCount};
		fsa_state mInternal;

		size_t mIndex;
		size_t mState;
		size_t mCount;
		size_t mPos;
		const T& mArray;





	public:

		StateMarginalIterator( const T& array ):
			mArray( array ),
			mIndex( 0 ),
			mState( 0 ),
			mCount( 0 ),
			mPos( 0 ),
			mInternal( start ) {
// 			if ( array[0] != 0 ) {
// 				throw runtime_error( "Container holding compressed marginal counts must contain 0 at the first position!" );
// 			}
			initForward();
		}

		void initForward() {
			mInternal = start;
		}


		// TODO we do not check the sum of counts for each position. Technically speaking, it would be possible for the array to contain runs of zeros.
		bool next() {
			while ( mIndex + 1 < mArray.size() ) {
				mIndex++;
				marginal_t value = mArray[mIndex];

				switch ( mInternal ) {

					case start:

						// first position, must be 0, TODO implementation of MarginalRecords does this slightly differently, with 0 added to the end
						mIndex = 0;
						value = mArray[0];
						mPos = 0;
						if ( value >= 0 ) {
							mState = value;
							mInternal = afterState;
							break;
						} else {
							mState = 0;
							mCount = -value;
							mInternal = afterCount;
							return true;
						}
						//throw runtime_error( "Malformed marginals, first entry must be zero!" );

					case afterPos:

						// counts for state 0
						if ( value < 0 ) {
							mCount = -value;
							mInternal = afterCount;
							return true;
						}

						// new state
						if ( value > 0 ) {
							mState = value;
							mInternal = afterState;
							break;
						}

						throw runtime_error( "Malformed marginals, expected count or state label!" );

					case afterCount:

						// another count, for next state
						if ( value < 0 ) {
							mState++;
							mCount = -value;
							return true;
						}

						// a new state
						if ( value > 0 ) {
							if ( value <= mState ) {
								throw runtime_error( "Malformed marginals, new state label must be larger than previous one!" );
							}
							mState = value;
							mInternal = afterState;
							break;
						}

						// a new pos
						if ( value == 0 ) {
							mPos++;
							mState = 0;
							mInternal = afterPos;
							break;
						}

					case afterState:

						// a new state count
						if ( value < 0 ) {
							mCount = - value;
							mInternal = afterCount;
							return true;
						}
						throw runtime_error( "Malformed marginals, expected state count." );
				}
			}

			// TODO change MarginalRecords
// 			if ( mInternal != afterCount ) {
// 				throw runtime_error( "Malformed marginals, last element must be a state count!" );
// 			}
			return false;
		}




		size_t count() const {
			return mCount;
		}

		size_t state() const {
			return mState;
		}

		size_t pos() const {
			return mPos;
		}


		void print() const {
			cout << pos() << " \t" << state() << " \t" << count() << endl;
		}
};



#endif



