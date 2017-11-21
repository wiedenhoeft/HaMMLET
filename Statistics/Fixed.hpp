#ifndef FIXEDEMISSIONS_HPP
#define FIXEDEMISSIONS_HPP

// The fixed version of data compression, different thresholds have no effect.


#include "../includes.hpp"
#include "../SufficientStatistics.hpp"
#include "../Tags.hpp"
#include "../Theta.hpp"
#include "../uintmath.hpp"

#include "../utils.hpp"

#include "../MultiVector.hpp"



template<typename SuffStatType>
class Statistics<Fixed, SuffStatType > {

		MultiVector<SuffStatType> mStats;
		vector<size_t> mBlockSizes;
		const size_t mSize;	// number of data points
		const size_t mNrDim;
		size_t mSizeIndex;	// current position in data
		size_t mStatsIndex;	// current position of dimension 0 in stats array
		size_t mBlockStart;
		size_t mBlockEnd;


	public:


		// delete copy constructor
		Statistics( const Statistics& that ) = delete;


		// NOTE this constructor swaps its input vectors, i.e. they are empty outsize of this class
		// TODO subdivide larger blocks to fit uint16_t
		Statistics(
		    vector<SufficientStatistics< SuffStatType>>& stats,
		    vector<size_t>& sizes
		) :
			mSize( accumulate( sizes.begin(), sizes.end(), 0 ) ),
			mNrDim( stats.size() / sizes.size() ) {

			mBlockSizes.swap( sizes );
			mStats.swap( stats );

			if ( mBlockSizes.size() == 0 ) {
				throw runtime_error( "Cannot create empty array for sufficient statistics!" );
			}
		}




		// this is a dummy, it doesn't do anything
		void createBlocks( real_t threshold ) {	}

		template<typename ParamType>
		void createBlocks( Theta<ParamType>& param );

		void initForward() {
			mSizeIndex = 0;
			mStatsIndex = 0;
			mBlockStart = 0;
			mBlockEnd = 0;
		}


		// get the next block under the current threshold
		// NOTE: for a cumulative sum array A shifted one position to the right, with A[0]=0, the sum of [start, end) = [start, end-1] is A[end]-A[start]
		bool next() {
			if ( mBlockEnd == 0 ) {
				mBlockEnd = mBlockSizes[0];
				return true;
			}

			if ( mSizeIndex + 1 < mSize ) {
				mSizeIndex++;
				mStatsIndex += mNrDim;
				mBlockStart = mBlockEnd;
				mBlockEnd += mBlockSizes[mSizeIndex];
				return true;
			} else {
				return false;
			}
		}

		const SufficientStatistics<SuffStatType>& suffStat( size_t dim ) const {
			return mStats[mStatsIndex + dim];

		}

		size_t nrDim() const {
			return mNrDim;
		}

		size_t nrBlocks() const {
			return mBlockSizes.size();
		}

		// Return the number of blocks induced by the current threshold. This can only be called once the iteration is complete.
		// TODO why do we have this twice?
		size_t size() const {
			return mBlockSizes.size();
		}

		size_t N() const {
			return mBlockSizes[mSizeIndex] ;
		}

		size_t T() const {
			return mSize;
		}

		void printBlock() const {
			cout << "[" << mBlockStart << ":" << mBlockEnd << ") " << mBlockEnd - mBlockStart << " ";
		}






		// given a state sequence and mappings, aggregate the sufficient statistics for each pemission posterior. NOTE This is handled here instgead of Theta, since it can be done in a much more numerically stable way
		template<typename StateSequenceType, typename HyperParamType>
		void aggregateStatistics(
		    const StateSequence<StateSequenceType>& q,
		    ThetaHyperParam<HyperParamType>& tau_theta,
		    const Mapping& mapping,
		    const size_t ignoreBlockSize = 1
		) {

			vector<KahanAggregator<SufficientStatistics<SuffStatType>>> aggregates;

			size_t i = 0;
			for ( size_t t = 0; t < nrBlocks(); ++t ) {
				for ( auto d = 0; d < mNrDim; ++d ) {
					const size_t index = mapping[q[t]][d];
					if ( index >= aggregates.size() ) {
						aggregates.resize( index + 1 );
					}
					aggregates[index].add( mStats[i] );
					++i;
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
