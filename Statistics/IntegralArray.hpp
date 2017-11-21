#ifndef INTEGRALARRAY_HPP
#define INTEGRALARRAY_HPP

#include "../includes.hpp"
#include "../SufficientStatistics.hpp"
#include "../Tags.hpp"
#include "../Theta.hpp"
#include "../uintmath.hpp"

#include <algorithm>
using std::rotate;

#include <deque>
using std::deque;

#include "../utils.hpp"

#include "../MultiVector.hpp"

typedef uint16_t PointerType;
const size_t CELLSIZE = 65535;	// for numeric reasons, the cumulative sum array is divided into cells of this size, and the reverse cumulative sum is calculated within that cell

// TODO rename
template<typename SuffStatType>
class Statistics<IntegralArray, SuffStatType > {

		// number of input data points that the breakpoints are derived from, mNrDim* mSize is the size of the stats arrays
		const size_t mSize;

		// number of input dimensions (for stats arrays) TODO current implementation only works for 1D
		const size_t mNrDim;

		// mStats[i] represents the sufficient statistics of [i,..,i+mPointers[i]-1] (inclusive)
		MultiVector<SufficientStatistics<SuffStatType>> mStats;


		// current state during iteration
		vector<SufficientStatistics<SuffStatType>> mCurrentSuffStat;


		// to compute the current sufficient statistics, we aggregate the positive and negative parts (as positive value) in a Kahan aggregator for numerical stability
		vector<KahanAggregator<SufficientStatistics<SuffStatType>>> mStatsSums;

// 		const size_t mBlockLimit;

		// Given a start and end position as well as two Kahan aggregators, add the positive and negative statistics for the segment.
		void addBlockStats(
		    size_t start,
		    size_t end,	// TODO check bounds?
		    size_t d,
		    KahanAggregator<SufficientStatistics<SuffStatType>>& stats
		) {
			// we use the Kahan aggregator to record the number of positions included in the sum, which we need to set manually
			size_t N = stats.nrTerms() + end - start;


			stats.add( mStats( start , d ) );
			for ( start = floorMult( start + CELLSIZE, CELLSIZE ); start < end; start += CELLSIZE ) {
				stats.add( mStats( start, d ) );
			}
			if ( end % CELLSIZE != 0 ) {
				stats.subtract( mStats( end, d ) );
			}

			// manually set number of terms
			stats.setNrTerms( N );
		}



	public:

		// delete copy constructor
		Statistics( const Statistics& that ) = delete;


		// NOTE this constructor swaps its input vectors, i.e. they are empty outsize of this class
		// TODO Multivector?
		Statistics(
		    vector<SufficientStatistics< SuffStatType>>& stats,
		    const size_t nrDim
		) :
			mSize( stats.size() / nrDim ),
			mNrDim( nrDim ),
			mCurrentSuffStat( nrDim, 0 ),
			mStats( nrDim ) ,
			mStatsSums( nrDim, KahanAggregator<SufficientStatistics<SuffStatType>>( ) ) {
			// TODO make parameter?

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

			// move sufficient statistics to right and compute cumulative sums; Having the first entry be zero means we don't have to check for t=0 start positions and also don't worry about underflow of t

			SufficientStatistics<SuffStatType> zero( 0 );
			mStats.reserve( mStats.size() + 1 );
			mStats.push_back( zero );


			// compute the partial cumulative sums in subarrays, for all dimensions

			const size_t skip = mNrDim * CELLSIZE;
			for ( size_t start = 0; start < mStats.size(); start += skip ) {
				for ( size_t d = 0; d < mNrDim; ++d ) {
					KahanCumulativeSum( mStats, start + d, start + skip + d, mNrDim, true );
				}
			}
		};


		// get the next block under the current threshold
		// NOTE: for a cumulative sum array A shifted one position to the right, with A[0]=0, the sum of [start, end) = [start, end-1] is A[end]-A[start]
		void setStats( const Blocks<BreakpointArray>& blocks ) {

			for ( size_t dim = 0; dim < mNrDim; ++dim ) {
				mStatsSums[dim].reset();
// 				cout << blocks.start() << " " << blocks.end() << endl;
				addBlockStats( blocks.start(), blocks.end(), dim, mStatsSums[dim] );
				mCurrentSuffStat[dim] = mStatsSums[dim].sum();
			}

		}

		const SufficientStatistics<SuffStatType>& suffStat( size_t dim ) const {
			return mCurrentSuffStat[dim];

		}

		size_t nrDim() const {
			return mNrDim;
		}


		size_t size() const {
			return mSize;
		}

		// given a state sequence and mappings, aggregate the sufficient statistics for each pemission posterior. NOTE This is handled here instgead of Theta, since it can be done in a much more numerically stable way
		template<typename StateSequenceType, typename HyperParamType, typename BlocksType>
		void aggregateStatistics(
		    const StateSequenceType& q,
		    BlocksType& blocks,
		    HyperParamType& tau_theta,
		    const Mapping& mapping,
		    const size_t ignoreBlockSize = 1
		) {
			// NOTE we collect start and end separately. Instead of having the sum of up to T subtractions, we have two Kahan summations and a single subtraction of numbers which can be expected to differ considerably
			vector<KahanAggregator<SufficientStatistics<SuffStatType>>> aggregates;


			blocks.initForward();
			while ( blocks.next() ) {
// 				setStats( blocks );	// NOTE we don't need to set statistics here
				if ( blocks.blockSize() > ignoreBlockSize ) {

					for ( auto d = 0; d < mNrDim; ++d ) {
						const size_t index = mapping[q[blocks.pos()]][d];
						if ( index >= aggregates.size() ) {
							aggregates.resize( index + 1 );
						}
						addBlockStats( blocks.start(), blocks.end(), d, aggregates[index] );
					}
				}
			}

			// add aggregated observations to posterior
			for ( size_t i = 0; i < aggregates.size(); ++i ) {
				size_t N = aggregates[i].nrTerms();
				if ( N > 0 ) {
					tau_theta.addObservation( aggregates[i].sum(), N, i );
				}
			}
		}

};


#endif

// TODO ::solidify(): given maximum margins, transform mStats into fixed stats array, and mPointers into blockSizes

