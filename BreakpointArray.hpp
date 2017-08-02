#ifndef BREAKPOINTARRAY_HPP
#define BREAKPOINTARRAY_HPP

#include "includes.hpp"
#include "SufficientStatistics.hpp"
#include "Tags.hpp"
#include "Theta.hpp"
#include "uintmath.hpp"

#include <type_traits>
using std::is_integral;
using std::is_unsigned;



// TODO support optional backward iteration? (Not necessary at the moment)

typedef uint16_t PointerType;	// TODO make template argument?

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
		vector<SufficientStatistics< SuffStatType>> mTempStats;	// used to aggregate sum of stats to subtract (numerically, adding is more stable than subtracting!)

		// mPointers[i] means that for all j in [i+1, i+mPointers[i]-1] (inclusive), mWeights[j] < mWeights[i]
		vector<PointerType> mPointers;

		// number of input dimensions (for stats arrays) TODO current implementation only works for 1D
		const size_t mNrDim;

// 		const size_t mNrInputPositions; // number of data points that created this data


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
		    vector<size_t>& indexStack,
		    vector<SufficientStatistics<SuffStatType>>& statsStack,
		    const size_t& right,
		    SufficientStatistics<SuffStatType>& tempStats ) {

			size_t left = indexStack.back();
			mPointers[left] =  right - left;
			indexStack.pop_back();

			tempStats = statsStack.back();
			statsStack.pop_back();
			mStats[left] = tempStats;
			if ( statsStack.size() > 0 ) {
				statsStack.back() += tempStats;
			}
		}

	public:

		// delete copy constructor
		Emissions( const Emissions& that ) = delete;

		Emissions(
		    vector<real_t>& weights,
		    vector<SufficientStatistics< SuffStatType>>& stats,
		    size_t nrDim = 1
		) :
			mSize( weights.size() ),
			mMaxJumpSize( min( mSize, ( size_t )numeric_limits<PointerType>::max() ) ),
			mPointers( mSize, mMaxJumpSize ),	// initialize all pointers to their maximum allowed value
			mNrDim( nrDim ),
// 			mNrInputPositions( weights.size() ),
			mDirection( unset ),
			mBlockCounter( 0 ),
			mCurrentSuffStat( mNrDim ) {

			mWeights.swap( weights );
			mStats.swap( stats );

			if ( mNrDim != 1 ) {
				throw runtime_error( "Multivariate emissions not implemented yet!" );	// TODO
			}

			mTempStats.resize( mNrDim );



			if ( mWeights.size() <= 0 ) {
				throw runtime_error( "Input vector for breakpoint weights is empty!" );
			}

			if ( !( divides( mStats.size(), mWeights.size() ) && mStats.size() >= mWeights.size() ) ) {
				throw runtime_error( "Cannot infer data dimension, size of statistics vector (" + to_string( mStats.size() ) + ") must be multiple of size of weight vector(" + to_string( mWeights.size() ) + ")!" );
			}

			if ( mStats.size() / mWeights.size() != mNrDim ) {
				throw runtime_error( "Number of statistics (" + to_string( mStats.size() ) + ") and number of breakpoint weights (" + to_string( mWeights.size() ) + ") do not result in the requested number of data dimenions (" + to_string( mNrDim ) + ")!" );
			}


			if ( mSize * mNrDim != mStats.size() ) {
				throw runtime_error( "Sufficient statistics array has the wrong size!" );
			}



			if ( !( ( is_integral<PointerType>::value ) && ( is_unsigned<PointerType>::value ) ) ) {
				throw runtime_error( "Pointer type must be unsigned and a type of integer (char, int, size_t etc.)" );
			}



			// calculate pointers and block sufficient statistics
			// NOTE stacks are implemented without container adapters, since we require random access to indexStack[0]
			vector<size_t> indexStack;
			vector<SufficientStatistics<SuffStatType>> statsStack;
			SufficientStatistics<SuffStatType> tempStats;
			indexStack.push_back( 0 );
			statsStack.push_back( mStats[0] );	// TODO multidimensional
			size_t left = 0;

			for ( size_t right = 1; right < mSize; ++right ) {
				size_t furthestIndex = indexStack[0];	// the index of the furthest unprocessed element still in the stack, used to determine whether we would jump too far
				while ( !indexStack.empty() ) {
					left = indexStack.back();

					// TODO we repeatedly check max jump size,  there are more elegant ways to implement this, but it is correct
					// if the furthest active index is as far away as we can jump

					if ( right - furthestIndex == mMaxJumpSize || mWeights[left] <= mWeights[right] ) {
						reduceStacks( indexStack, statsStack, right, tempStats );
					} else {
						break;	// weights only get larger further down the stack
					}
				}
				indexStack.push_back( right );
				statsStack.push_back( mStats[right] );
			}
			size_t right = mSize;
			// elements still on the stack all point past the end
			while ( statsStack.size() > 0 ) {
				reduceStacks( indexStack, statsStack, right, tempStats );
			}
		};




		void createBlocks( real_t threshold ) {
			// TODO should we ignore states with no observation for minimum variance calculation?
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


		// get the end of a block for a given start <pos>ition under the current threshold
		inline size_t blockEnd( size_t pos ) {
			pos++;
			while ( pos < mSize ) {
				// NOTE to future self: you might think that aggregating very short blocks here will help with increasing compression for high-variance components; it does, but leads to dramatic overcompression and is a very bad idea...
				if ( mWeights[pos] < mThreshold ) {
					pos += mPointers[pos];
				} else {
					break;
				}
			}
			return pos;	// TODO enforce length?
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
			mBlockEnd = blockEnd( mBlockStart );


			// set current statistics to the ones observed at the block start these might be too large or too small, depending on whether and how the target differs from the block end
			size_t suffStatIndex = mNrDim * mBlockStart;
			for ( size_t dim = 0; dim < mNrDim; ++dim ) {
				mCurrentSuffStat[dim] = mStats[suffStatIndex + dim];
			}
			mCurrentBlockSize =  mBlockEnd - mBlockStart ;


			// target represents the position the block start points to. This may or may not coincide with the block end. We first set the current stats to the one between mBlockStart and target, and then use blockiter to remove whatever was too much
			const size_t target = mBlockStart + mPointers[mBlockStart];

			if ( target == mBlockEnd ) {	// if both match, the statistics are correct
				return true;
			} else {	// the statistics might be off, so we accumulate the difference in mTempStats

				// clear mTempStats from previous use
				for ( size_t dim = 0; dim < mNrDim; ++dim ) {
					mTempStats[dim].clear();
				}

				// define the range for mTempStats
				size_t left = min( mBlockEnd, target );
				const size_t right = max( mBlockEnd, target );

				// accumulate difference
				while ( left < right ) {
					suffStatIndex = mNrDim * left;
					for ( size_t dim = 0; dim < mNrDim; ++dim ) {
						mTempStats[dim] += mStats[suffStatIndex + dim];
					}
					left += mPointers[left];
				}
				if ( left != right ) {
					throw runtime_error( "Left (" + to_string( left ) + ") does not match right (" + to_string( right ) + ") index in block adjustment!" );
				}

				// add or subtract mTempStats, depending on whether the jump was too short or too far
				if ( right == mBlockEnd ) {
					// add additional stats
					for ( size_t dim = 0; dim < mNrDim; ++dim ) {
						mCurrentSuffStat[dim] += mTempStats[dim];
					}
				} else {	// right == target
					// subtract additional stats
					for ( size_t dim = 0; dim < mNrDim; ++dim ) {
						mCurrentSuffStat[dim] -= mTempStats[dim];
					}
				}
				return true;
			}
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

		void printData( real_t tolerance = 0.0001 )const {
			for ( size_t start = 0; start < mSize; ++start ) {
				size_t jump = mPointers[start];
				size_t end = start + jump;

				real_t sum = 0;
				real_t sumSq = 0;
				for ( size_t i = start; i < end; ++i ) {
					sum += mWeights[i];
					sumSq += mWeights[i] * mWeights[i];
				}
				if ( abs( sum - mStats[start].sum() ) > tolerance || abs( sumSq - mStats[start].sumSq() ) > tolerance ) {
					cout << start << "-" << jump << "->" << end << " " << endl;
					cout << sum << " " << sumSq << endl;
					cout  << mStats[start].sum() << " " << mStats[start].sumSq() << " " << endl;
					cout << endl;
				}
			}
			cout << "max jump size: " << mMaxJumpSize << endl;
		}

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

