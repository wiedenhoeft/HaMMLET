#ifndef BREAKPOINTARRAY_HPP
#define BREAKPOINTARRAY_HPP

#include "../Blocks.hpp"

#include "../includes.hpp"
#include "../Tags.hpp"
#include "../Theta.hpp"
#include "../uintmath.hpp"
#include "../utils.hpp"

#include <algorithm>
using std::rotate;

#include <deque>
using std::deque;



typedef uint16_t PointerType;

// generates a block structure for any go
template<>
class Blocks<BreakpointArray> {

		// number of input positions
		const size_t mSize;

		// the maximum size that an iterator can jump forward (pruning limit)
		Direction mDirection;

		// mWeights[i] represents the weight of breakpoint [i-1,i]. This also means that mWeights[0] is essentially ignored.
		vector<real_t> mWeights;

		// mPointers[i] means that for all j in [i+1, i+mPointers[i]-1] (inclusive), mWeights[j] < mWeights[i]
		vector<PointerType> mPointers;

		real_t mThreshold;
		size_t mBlockCounter;

		// the boundaries of the current block
		size_t mBlockStart;
		size_t mBlockEnd;
		size_t mBlockSize;


		// this creates a pointer target to <right> based on the top of the stacks
		inline void reduceStacks(
		    deque<size_t>& indexStack,
		    const size_t& right );

	public:

		// delete copy constructor
// 		Blocks( const Blocks& that ) = delete;


		// NOTE this constructor swaps its input vectors, i.e. they are empty outsize of this class
		Blocks(
		    vector<real_t>& weights );


		void createBlocks( real_t threshold );

		template<typename ParamType>
		void createBlocks( const Theta<ParamType>& param );

		void initForward();



		// get the end of a block starting at <start> for a given threshold
		// return false if the block end is the last possible value
		inline bool next();



		// average weight of breakpoints, can be used to derive a block structure for automatic priors for instance
		real_t avgWeight() const;


		size_t start() const;

		size_t end() const;

		size_t pos() const;

		// Return the size of the current block.
		size_t blockSize() const;

		// Return the total size, i.e. the sum of all block sizes.
		size_t size() const;

		size_t nrBlocks() const;


		void printBlock() const;


};












// this creates a pointer target to <right> based on the top of the stacks

inline void Blocks<BreakpointArray>::reduceStacks(
    deque<size_t>& indexStack,
    const size_t& right ) {

	// set pointer for stretch
	size_t left = indexStack.back();
	mPointers[left] =  right - left;
	indexStack.pop_back();
}




// NOTE this constructor swaps its input vectors, i.e. they are empty outsize of this class

Blocks<BreakpointArray>::Blocks(
    vector<real_t>& weights
) :
	mSize( weights.size() ),
	mDirection( unset ),
	mBlockCounter( 0 ) {
	// TODO make parameter?

	mWeights.swap( weights );



	// check that weights contain data
	if ( mSize <= 0 ) {
		throw runtime_error( "Input vector for breakpoint weights is empty!" );
	}

	// calculate pointers
	// NOTE stacks are implemented without container adapters, since we require random access to indexStack[0]
	// the maximum value the pointers can take
	const PointerType maxJumpSize = min( mSize, ( size_t )numeric_limits<PointerType>::max() );

	// initialize all pointers to their maximum allowed value
	mPointers.assign( mSize, maxJumpSize );
	deque<size_t> indexStack;
	indexStack.push_back( 0 );
	size_t left = 0;
	for ( size_t right = 1; right < mSize; ++right ) {

		// check if the furthest element in the deque  has reached its maximum jump size, and set its pointer if necessary
		if ( !indexStack.empty() ) {
			size_t furthestIndex = indexStack.front();
			if ( right - furthestIndex == maxJumpSize ) {
				mPointers[furthestIndex] = maxJumpSize;
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

};




void Blocks<BreakpointArray>::createBlocks( real_t threshold ) {
	mThreshold = threshold;
}


// TODO find more elegant and flexible way to compute thresholds
template<>
void Blocks<BreakpointArray>::createBlocks(
    const  Theta<NormalParam>& param ) {
	createBlocks( sqrt( 2 * log( ( real_t )mSize ) *param.thresholdValue() ) );
}



void Blocks<BreakpointArray>::initForward() {
	mDirection = forward;
	mBlockStart = 0;
	mBlockEnd = 0;
	mBlockSize = 0;
	mBlockCounter = 0;
}



// get the end of a block starting at <start> for a given threshold
// return false if the block end is the last possible value

inline bool Blocks<BreakpointArray>::next() {
	if ( mBlockEnd >= mSize ) {
		mDirection = unset;
		return false;
	} else {
		mBlockCounter++;
		mBlockStart = mBlockEnd;
		mBlockEnd = mBlockStart + 1;
		while ( mBlockEnd < mSize ) {
			// TODO how to handle overflow. Maximum block size?
			if ( mWeights[mBlockEnd] < mThreshold ) {
				mBlockEnd += mPointers[mBlockEnd];	// NOTE this involves typecasting by necessity
			} else {
				break;
			}
		}
		mBlockSize = mBlockEnd - mBlockStart;
		return true;
	}
}



// average weight of breakpoints, can be used to derive a block structure for automatic priors for instance

real_t Blocks<BreakpointArray>::avgWeight() const {
	real_t sum = 0;
	for ( const auto & w : mWeights ) {
		if ( isfinite( w ) ) {
			sum += w;
		}
	}
	if ( isfinite( mWeights[0] ) ) {
		sum -= mWeights[0];	// the first element isn't really a true weight, as there is always a breakpoint before the first element
	}
	return sum / ( ( double )mWeights.size() - 1 );
}



size_t Blocks<BreakpointArray>::start() const {
	return mBlockStart;
}



size_t Blocks<BreakpointArray>::end() const {
	return mBlockEnd;
}



size_t Blocks<BreakpointArray>::pos() const {
	if ( mBlockCounter > 0 ) {
		return mBlockCounter - 1;
	} else {
		throw runtime_error( "No blocks created yet, position is undefined!" );
	}
}



// Return the size of the current block.

size_t Blocks<BreakpointArray>::blockSize() const {
	return mBlockSize ;
}



// Return the total size, i.e. the sum of all block sizes.

size_t Blocks<BreakpointArray>::size() const {
	return mSize;
}



// Return the number of blocks induced by the current threshold. This can only be called once the iteration is complete.

size_t Blocks<BreakpointArray>::nrBlocks() const {
	if ( mDirection != unset ) {
		throw runtime_error( "Cannot determine size of block structure before all blocks have been seen!" );
	}
	return mBlockCounter;
}



void Blocks<BreakpointArray>::printBlock() const {
	cout << "[" << mBlockStart << ":" << mBlockEnd << ") " << mBlockSize << " ";
}

#endif
