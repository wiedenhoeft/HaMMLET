#ifndef STATEMARGINALS_HPP
#define STATEMARGINALS_HPP

#include "includes.hpp"


#include <deque>
using std::deque;



// NOTE Strictly speaking, there is the possibility that two adjacent segments will have the same marginals even though they were different in previous iterations. They could be merged, but this is not worth the effort in practice.




// template<typename T>
class StateMarginals {
		deque<marginal_t> mCountQ;	// store marginal counts as negative values, and use positive values to denote the state of the next element of states are skipped due to having zero count. The size of this is guaranteed to be t most #segments*(1+max(#states, #iterations)), plus 1/32 bit overhead per element due to using std::deque.
		deque<size_t> mSizeQ;
		size_t mCurrentindex;
		size_t mSize;
		size_t mNrIterations;
		size_t mNrStates;	// the total number of states, i.e. the highest added state label + 1
		size_t mNrSegments;

		void addRecord(
		    const marginal_t& state, // the state to be recorded for the current mCountQ.front()
		    size_t blockSize ,	// the block size associated with it
		    const marginal_t count = 1	// how many counts to add for this state
		) {

			if ( state >= mNrStates ) {
				mNrStates = state + 1;
			}

			bool done = false;
			marginal_t s = 0;	// holds the state associated with the current mCountQ.front()
			while ( blockSize > 0 ) {
				s = 0;
				done = false;
				if ( mCountQ.size() == 0 ) {
					throw runtime_error( "Empty count queue, this is a bug!" );
				}
				for ( size_t t = 0; t < mCountQ.size(); ++t ) { // push back the front segment and include the curent count
					marginal_t entry = mCountQ[t];
					// end of the record
					if ( entry == 0 ) {
						if ( !done) { // if the highest state seen before is less than the one we want to add, we have yet to do so
							if ( s < state ) { // do we need index information?
								mCountQ.push_back( state );
							}
							mCountQ.push_back( -count );
						}
						mCountQ.push_back( 0 );
						break;
					}
					if ( done ) {
						mCountQ.push_back( entry );
						continue;
					}
					// indexing information, current entry will be the next state
					if ( entry > 0 ) {
						if ( state < entry ) {	// we need to insert the state before we proceed
							if ( s < state ) {	// gap to previous state?
								mCountQ.push_back( state );
							}
							mCountQ.push_back( -count );
							if ( state + 1 < entry ) {	// gap to next state?
								mCountQ.push_back( entry );
							}
							done = true;
						} else {
							mCountQ.push_back( entry );
						}
						s = entry; //next state
					}
					// this is the negative count for the current state
					else {
						if ( s == state ) {
							mCountQ.push_back( entry - count );
							done = true;
						} else {
							mCountQ.push_back( entry );
						}
						s++;	// next state
					}
				}

				if ( blockSize < mSizeQ.front() ) {		// residual front segment remains, decrease its size
					mSizeQ.push_back( blockSize );
					mSizeQ.front() -= blockSize;
					mNrSegments++;
					break;
				} else {	// front segment was completely absorbed, pop it
					mSizeQ.push_back( mSizeQ.front() );
					blockSize -= mSizeQ.front();
					mSizeQ.pop_front();	// pop_front(  )
					while ( mCountQ.front() != 0 ) {
						mCountQ.pop_front();
					}
					mCountQ.pop_front();	// pop zero
				}
			}
		}

	public:

		// delete copy constructor
		StateMarginals( const StateMarginals& that ) = delete;

		StateMarginals( size_t size ): mCurrentindex( 0 ), mSize( size ), mNrIterations( 0 ), mNrStates( 0 ), mNrSegments( 1 ) {
			mSizeQ.push_back( size );
			mCountQ.push_back( 0 );
		}



		template<typename StateSequenceType, typename EmissionsDataStructure, typename EmissionsDistType >
		void record(
		    const StateSequence<StateSequenceType>& q ,
		    Emissions< EmissionsDataStructure, EmissionsDistType>& y	// TODO use a non-const wrapper around y?
		) {
			++mNrIterations;

			size_t totalBlockSize = 0;
			size_t currentState = 0;
			size_t blockSize = 0;


			size_t t = 0;
			y.initForward();	// this uses the previous threshold etc.
			while ( y.next() ) {
				if ( t >= q.size() ) {
					throw runtime_error( "The state sequence is shorter (" + to_string( q.size() ) + ") than the number of blocks (" + to_string( t ) + ")!" );
				}
				if ( t >= q.size() ) {
					throw runtime_error( "State sequence index out of bounds!" );
				}
				if ( q[t] == currentState ) {
					blockSize += y.N();
				} else {

					addRecord( currentState, blockSize, 1 );
					totalBlockSize += blockSize;
					currentState = q[t];
					blockSize = y.N();
				}
				++t;
			}

			if ( t < q.size() ) {
				throw runtime_error( "The state sequence is longer (" + to_string( q.size() ) + ") than the number of blocks (" + to_string( t ) + ")!" );
			}

			if ( t > q.size() ) {
				throw runtime_error( "The state sequence is shorter (" + to_string( q.size() ) + ") than the number of blocks (" + to_string( t ) + ")!" );
			}


			addRecord( currentState, blockSize, 1 );

			totalBlockSize += blockSize;

			if ( totalBlockSize != y.T() ) {
				throw runtime_error( "Total number of chunks recorded (" + to_string( totalBlockSize ) + ") in marginals does not match the number of positions in the data(" + to_string( y.T() ) + ")!" );
			}

			if ( mNrSegments != mSizeQ.size() ) {
				throw runtime_error( "The number of segments (" + to_string( mNrSegments ) + ") does not match the count of segment sizes (" + to_string( mSizeQ.size() ) + ")" );
			}
		}


		// return the number of distinct segments
		size_t nrSegments() const {
			return mNrSegments;
		}

		// return the number of values used to store the compressed marginals
		size_t internalSize() const {
			return mCountQ.size();
		}

		void save(
		    string filename,
		    size_t chunkSize = 1, // multiply the segment size by this number in the output file
		    bool verbose = false
		) const {
			ofstream ofs( filename.c_str() );
			save( ofs, chunkSize, verbose );
		}


		void save(
		    ofstream& ofs,
		    size_t chunkSize = 1, // multiply the segment size by this number in the output file
		    bool verbose = false
		) const  {	// NOTE since marginals are stores only up to the largest state that has non-zero counts for space efficiency reasons, nrStates can be used to fill the remaining counts with zero. Otherwise, these states are not written.
			if ( chunkSize <= 0 ) {
				throw runtime_error( "Chunk size must be at least one!" );
			}
			if ( mCurrentindex != 0 ) {
				throw runtime_error( "Cannot output incomplete marginals!" );
			}

			size_t i = 0;
			marginal_t s = 0;
			for ( size_t segSize : mSizeQ ) {
				size_t iterations = 0;

				s = 0;
				ofs << segSize;
				while ( mCountQ[i] != 0 ) {
					if ( mCountQ[i] > 0 ) {
						while ( s < mCountQ[i] ) {
							ofs << "\t" << 0;
							s++;
						}
					} else {
						s++;
						ofs << "\t" << -mCountQ[i];
						iterations -= mCountQ[i];
					}
					i++;
				}
				while ( s < mNrStates ) {
					ofs << "\t" << 0;
					s++;
				}
				i++;	// skip over zero
				ofs << endl;
				if ( iterations != mNrIterations ) {
					throw runtime_error( "Sum of marginals (" + to_string( iterations ) + ") does not match the number of iterations (" + to_string( mNrIterations ) + ")!" );
				}
			}
		}
};

#endif










