#ifndef INTEGRALARRAY_HPP
#define INTEGRALARRAY_HPP

#include "../Statistics.hpp"
#include "../Blocks/BreakpointArray.hpp"
#include "../includes.hpp"
#include "../SufficientStatistics.hpp"
#include "../Tags.hpp"
#include "../uintmath.hpp"
#include "../KahanAggregator.hpp"


#include <algorithm>
using std::rotate;

#include <deque>
using std::deque;

#include "../utils.hpp"

#include "../MultiVector.hpp"

typedef uint16_t PointerType;
const size_t CELLSIZE = 65535;	// for numeric reasons, the cumulative sum array is divided into cells of this size, and the reverse cumulative sum is calculated within that cell


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
		vector<KahanAggregator<SufficientStatistics<SuffStatType>>> mCurrentStats;


		// Given a start and end position as well as two Kahan aggregators, add the positive and negative statistics for the segment.
		void addBlockStats(
		    size_t start,
		    size_t end,	// TODO check bounds?
		    size_t d,
		    KahanAggregator<SufficientStatistics<SuffStatType>>& stats );


	public:

		Statistics( const Statistics& that ) = delete;

		Statistics(
		    vector<SufficientStatistics< SuffStatType>>& stats,
		    const size_t nrDim	);

		template<typename T>
		void setStats(
		    const Blocks<T>& blocks );

		const SufficientStatistics<SuffStatType>& suffStat(
		    size_t dim ) const;

		size_t nrDim() const;

		size_t size() const;

};


























// Given a start and end position as well as two Kahan aggregators, add the positive and negative statistics for the segment.
template<typename SuffStatType>
void Statistics<IntegralArray, SuffStatType >::addBlockStats(
    size_t start,
    size_t end,	// TODO check bounds?
    size_t d,
    KahanAggregator<SufficientStatistics<SuffStatType>>& stats
) {
	// we use the Kahan aggregator to record the number of positions included in the sum, which we need to set manually
	size_t N = stats.nrTerms() + end - start;


	stats.add( mStats( start , d ) );
	for ( start = higher_mult( start, CELLSIZE ); start < end; start += CELLSIZE ) {
		stats.add( mStats( start, d ) );
	}
	if ( end % CELLSIZE != 0 ) {
		stats.subtract( mStats( end, d ) );
	}

	// manually set number of terms
	stats.setNrTerms( N );
}



// delete copy constructor
// template<typename SuffStatType>
// Statistics<IntegralArray, SuffStatType >::Statistics( const Statistics& that ) = delete;


// NOTE this constructor swaps its input vectors, i.e. they are empty outsize of this class
// TODO Multivector?
template<typename SuffStatType>
Statistics<IntegralArray, SuffStatType >::Statistics(
    vector<SufficientStatistics< SuffStatType>>& stats,
    const size_t nrDim
) :
	mSize( stats.size() / nrDim ),
	mNrDim( nrDim ),
	mCurrentSuffStat( nrDim, 0 ),
	mStats( nrDim ) ,
	mCurrentStats( nrDim, KahanAggregator<SufficientStatistics<SuffStatType>>( ) ) {
	// TODO make parameter?


	
	// check that weights contain data
	if ( mSize <= 0 ) {
		throw runtime_error( "Input vector for breakpoint weights is empty!" );
	}

	//check that stats contain data
	if ( stats.size() <= 0 ) {
		throw runtime_error( "Input vector for breakpoint weights is empty!" );
	}

	if ( !divides( stats.size(), mSize ) ) {
		throw runtime_error( "Error constructing breakpoint array: number of sufficient statistics must be an integer multiple of the number of weights!" );
	}


	// check that size of stats array is integer mutliple of size of weights array
	if ( !( divides( stats.size(), mSize ) && stats.size() >= mSize ) ) {
		throw runtime_error( "Cannot infer data dimension, size of statistics vector (" + to_string( stats.size() ) + ") must be multiple of size of weight vector(" + to_string( mSize ) + ")!" );
	}

	// move sufficient statistics to right and compute cumulative sums; Having the first entry be zero means we don't have to check for t=0 start positions and also don't worry about underflow of t

	// TODO might cause reallocation
	SufficientStatistics<SuffStatType> zero( 0 );
	stats.reserve( stats.size() + mNrDim );
	for ( size_t d = 0; d < mNrDim; ++d ) {
		stats.push_back( zero );
	}


	// compute the partial cumulative sums in subarrays, for all dimensions

	const size_t skip = mNrDim * CELLSIZE;
	for ( size_t start = 0; start < stats.size(); start += skip ) {
		for ( size_t d = 0; d < mNrDim; ++d ) {
			KahanCumulativeSum( stats, start + d, start + skip + d, mNrDim, true );
		}
	}


	mStats.swap( stats );

};




template<typename SuffStatType>
template<typename T>
void Statistics<IntegralArray, SuffStatType >::setStats(
    const Blocks<T>& blocks ) {

	/*
	 Get the next block under the current threshold.
	NOTE: for a cumulative sum array A shifted one position to the right, with A[0]=0, the sum of [start, end) = [start, end-1] is A[end]-A[start]
	*/

	for ( size_t dim = 0; dim < mNrDim; ++dim ) {
		mCurrentStats[dim].reset();
		addBlockStats( blocks.start(), blocks.end(), dim, mCurrentStats[dim] );
		mCurrentSuffStat[dim] = mCurrentStats[dim].sum();
	}

}




template<typename SuffStatType>
const SufficientStatistics<SuffStatType>& Statistics<IntegralArray, SuffStatType >::suffStat( size_t dim ) const {
	return mCurrentSuffStat[dim];

}


template<typename SuffStatType>
size_t Statistics<IntegralArray, SuffStatType >::nrDim() const {
	return mNrDim;
}


template<typename SuffStatType>
size_t Statistics<IntegralArray, SuffStatType >::size() const {
	return mSize;
}







#endif


