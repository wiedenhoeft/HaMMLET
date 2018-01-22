#ifndef EMISSIONS_HPP
#define EMISSIONS_HPP

#include "includes.hpp"
#include "Tags.hpp"
#include "Blocks.hpp"
#include "Statistics.hpp"


// A wrapper around a combination of a data structure holding data points/sufficient statistics and an associated block structure
template<typename S, typename T, typename B>
class Emissions<Statistics<S, T>, Blocks<B>> {

		Statistics<S, T>& mStats;
		Blocks<B>& mBlocks;

	public:

		Emissions(
		    Statistics<S, T>& stats,
		    Blocks<B>& blocks
		):
			mStats( stats ),
			mBlocks( blocks ) {
			if ( mStats.size() != mBlocks.size() ) {
				throw runtime_error( "Block structure and statistics have different number of data points!" );
			}
		}


		Statistics<S, T>& stats()  {
			return mStats;
		}

		Blocks<B>& blocks()  {
			return mBlocks;
		}


		const Statistics<S, T>& stats() const {
			return mStats;
		}

		const Blocks<B>& blocks() const {
			return mBlocks;
		}



		void createBlocks( real_t thresh ) {
			mBlocks.createBlocks( thresh );
		}

		template <typename ParamType>
		void createBlocks( const Theta<ParamType>& theta ) {
			mBlocks.createBlocks( theta );
		}


		size_t nrBlocks() const {
			return mBlocks.nrBlocks();
		}

		size_t nrDim() const {
			return mStats.nrDim();
		}

		size_t start() const {
			return mBlocks.start();
		}

		size_t end() const {
			return mBlocks.end();
		}

		size_t blockSize() const {
			return mBlocks.blockSize();
		}

		size_t size() const {
			return mBlocks.size();
		}

		void initForward() {
			mBlocks.initForward();
		}

		bool next() {
			if ( mBlocks.next() ) {
				mStats.setStats( mBlocks );
				return true;
			} else {
				return false;
			}
		}


		const SufficientStatistics<T>& suffStat( size_t dim ) const {
			return mStats.suffStat( dim );
		}


		template<typename StateSequenceType, typename HyperParamType>
		void aggregateStatistics(
		    const StateSequenceType& q,
		    HyperParamType& tau_theta,
		    const Mapping& mapping,
		    const size_t ignoreBlockSize = 1
		) {
			mStats.aggregateStatistics( q, mBlocks, tau_theta, mapping, ignoreBlockSize );
		}
};






#endif



