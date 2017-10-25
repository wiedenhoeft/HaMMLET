#ifndef BREAKPOINTARRAY_HPP
#define BREAKPOINTARRAY_HPP

#include "includes.hpp"
#include "SufficientStatistics.hpp"
#include "Tags.hpp"
#include "Theta.hpp"
#include "uintmath.hpp"


#include <deque>
using std::deque;


typedef uint16_t PointerType;

template<typename SuffStatType>
class Emissions<BreakpointArray, SuffStatType > {

		// number of input data points that the breakpoints are derived from, mNrDim* mSize is the size of the stats arrays
		const size_t mSize;

		// the maximum size that an iterator can jump forward (pruning limit)
		const size_t mMaxJumpSize;

		Direction mDirection;

		// mWeights[i] represents the weight of breakpoint [i-1,i]. This also means that mWeights[0] is essentially ignored.
		vector<real_t> mWeights;


		// mStats[i] represents the sufficient statistics of [i,..,i+mPointers[i]-1] (inclusive)
		vector<SufficientStatistics<SuffStatType>> mStats;

		// mPointers[i] means that for all j in [i+1, i+mPointers[i]-1] (inclusive), mWeights[j] < mWeights[i]
		vector<PointerType> mPointers;

		// number of input dimensions (for stats arrays) TODO current implementation only works for 1D
		const size_t mNrDim;

		// current state during iteration
		vector<SufficientStatistics<SuffStatType>> mCurrentSuffStat;
		real_t mThreshold;
		size_t mBlockCounter;

		// the boundaries of the current block
		size_t mBlockStart;
		size_t mBlockEnd;
		size_t mCurrentBlockSize;


		// this creates a pointer target to <right> based on the top of the stacks
		inline void reduceStacks(
		    deque<size_t>& indexStack,
		    const size_t& right ) {

			// set pointer for stretch
			size_t left = indexStack.back();
			mPointers[left] =  right - left;
			indexStack.pop_back();
		}

	public:

		// delete copy constructor
		Emissions( const Emissions& that ) = delete;


		// NOTE this constructor swaps its input vectors, i.e. they are empty outsize of this class
		Emissions(
		    vector<real_t>& weights,
		    vector<SufficientStatistics< SuffStatType>>& stats
		) :
			mSize( weights.size() ),
			mMaxJumpSize( min( mSize, ( size_t )numeric_limits<PointerType>::max() ) ),
			mPointers( mSize, mMaxJumpSize ),	// initialize all pointers to their maximum allowed value
			mNrDim( stats.size() / weights.size() ),
			mDirection( unset ),
			mBlockCounter( 0 ),
			mCurrentSuffStat( mNrDim ) {

			mWeights.swap( weights );
			mStats.swap( stats );


			// check that weights contain data
			if ( mSize <= 0 ) {
				throw runtime_error( "Input vector for breakpoint weights is empty!" );
			}

			//check that stats contain data
			if ( mStats.size() <= 0 ) {
				throw runtime_error( "Input vector for breakpoint weights is empty!" );
			}

			if ( !divides( mStats.size(), mSize ) ) {
				throw runtime_error( "Error constructing breakpoint array: number of sufficient statistics must be an integer multiple of the number of weights!" );
			}


			// check that size of stats array is integer mutliple of size of weights array
			if ( !( divides( mStats.size(), mSize ) && mStats.size() >= mSize ) ) {
				throw runtime_error( "Cannot infer data dimension, size of statistics vector (" + to_string( mStats.size() ) + ") must be multiple of size of weight vector(" + to_string( mSize ) + ")!" );
			}


			// calculate pointers
			// NOTE stacks are implemented without container adapters, since we require random access to indexStack[0]
			deque<size_t> indexStack;
			indexStack.push_back( 0 );
			size_t left = 0;
			for ( size_t right = 1; right < mSize; ++right ) {

				// check if the furthest element in the deque  has reached its maximum jump size, and set its pointer if necessary
				if ( !indexStack.empty() ) {
					size_t furthestIndex = indexStack.front();
					if ( right - furthestIndex == mMaxJumpSize ) {
						mPointers[furthestIndex] = mMaxJumpSize;
						indexStack.pop_front();
					}
				}

				while ( !indexStack.empty() ) {
					left = indexStack.back();

					if ( mWeights[left] <= mWeights[right] ) {
						reduceStacks( indexStack, right );
					} else {
						break;	// weights only get larger further down the stack
					}
				}
				indexStack.push_back( right );
			}
			// elements still on the stack all point past the end
			while ( indexStack.size() > 0 ) {
				reduceStacks( indexStack, mSize );
			}




			// move sufficient statistics to right and compute cumulative sums; Having the first entry be zero means we don't have to check for t=0 start positions and also don't worry about underflow of t
			mStats.reserve( mStats.size() + mNrDim );
			for ( size_t d = 0; d < mNrDim; ++d ) {
				mStats.push_back( SufficientStatistics<SuffStatType>( 0 ) );
			}
			SufficientStatistics<SuffStatType> currStats;
			vector<SufficientStatistics<SuffStatType>> sumStats;
			sumStats.resize( mNrDim, SufficientStatistics<SuffStatType>( 0 ) );
			size_t index = 0;
			for ( size_t t = 0; t <= mSize; ++t ) {	// NOTE t <= T is correct due to resizing
				for ( size_t d = 0; d < mNrDim; ++d ) {
					currStats = mStats[index];
					mStats[index] = sumStats[d];
					sumStats[d] += currStats;
					index++;
				}
			}
		};




		void createBlocks( real_t threshold ) {
			mThreshold = threshold;
		}

		template<typename ParamType>
		void createBlocks( Theta<ParamType>& param );

		void initForward() {
			mDirection = forward;
			mBlockStart = 0;
			mBlockEnd = 0;
			mCurrentBlockSize = 0;
			mBlockCounter = 0;
			// initialize block at 0
		}


		// get the next block under the current threshold
		bool next() {
			if ( mBlockEnd >= mSize ) {
				mDirection = unset;
				return false;
			}
			mBlockCounter++;

			//// determine the end of the block and set its size ////
			mBlockStart = mBlockEnd;
			mBlockEnd = mBlockStart + 1;
			while ( mBlockEnd < mSize ) {
				// NOTE to future self: you might think that aggregating very short blocks here will help with increasing compression for high-variance components; it does, but leads to dramatic overcompression and is a very bad idea...
				if ( mWeights[mBlockEnd] < mThreshold ) {
					mBlockEnd += mPointers[mBlockEnd];	// NOTE this involves typecasting by necessity
				} else {
					break;
				}
			}


			// set current statistics to the ones observed at the block start these might be too large or too small, depending on whether and how the target differs from the block end
			size_t suffStatStartIndex = mNrDim * mBlockStart;
			size_t suffStatEndIndex = mNrDim * mBlockEnd;


			for ( size_t dim = 0; dim < mNrDim; ++dim ) {
				mCurrentSuffStat[dim] = mStats[suffStatEndIndex] - mStats[suffStatStartIndex];
				suffStatStartIndex++;
				suffStatEndIndex++;
			}
			mCurrentBlockSize =  mBlockEnd - mBlockStart ;

			return true;
		}


		const SufficientStatistics<SuffStatType>& suffStat( size_t dim ) const {
			return mCurrentSuffStat[dim];

		}

		size_t nrDim() const {
			return mNrDim;
		}

		size_t nrBlocks() const {
			if ( mDirection == unset ) {
				throw runtime_error( "Cannot obtain number of blocks before block structure is finished!" );
			}
			return mBlockCounter;
		}


		// average weight of breakpoints, can be used to derive a block structure for automatic priors for instance
		real_t avgWeight() const {
			real_t mean = 0;
			for ( const auto & w : mWeights ) {
				if ( isfinite( w ) ) {
					mean += w;
				}
			}
			if ( isfinite( mWeights[0] ) ) {
				mean -= mWeights[0];	// the first element isn't really a true weight, as there is always a breakpoint before the first element
			}
			return mean / ( ( real_t )mWeights.size() - 1 );
		}


		size_t N() const {
			return mCurrentBlockSize ;
		}

		size_t T() const {
// 			return mNrInputPositions;
			return mSize;
		}

		void printBlock() const {
			cout << "[" << mBlockStart << ":" << mBlockEnd << ") " << mCurrentBlockSize << " ";
		}


		// Return the number of blocks induced by the current threshold. This can only be called once the iteration is complete.
		size_t size() const {
			if ( mDirection != unset ) {
				throw runtime_error( "Cannot determine size of block structure before all blocks have been seen!" );
			}
			return mBlockCounter;
		}

};


// TODO find more elegant and flexible way to compute thresholds
template<> template<>
void Emissions<BreakpointArray, Normal>::createBlocks(
    Theta<NormalParam>& param ) {
	createBlocks( sqrt( 2 * log( ( real_t )mSize ) *param.thresholdValue() ) );
}

#endif



