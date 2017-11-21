#ifndef FIXEDBLOCKS_HPP
#define FIXEDBLOCKS_HPP

template<>
class Blocks<Fixed> {
		// number of input positions
		const size_t mSize;
		const vector<size_t>& mSizes;

		size_t mBlockCounter;
		// the maximum size that an iterator can jump forward (pruning limit)
		Direction mDirection;

		// the boundaries of the current block
		size_t mBlockStart;
		size_t mBlockEnd;
		size_t mBlockSize;

	public:

		// delete copy constructor
		Blocks( const Blocks& that ) = delete;


		// NOTE this constructor swaps its input vectors, i.e. they are empty outsize of this class
		Blocks(
		    vector<size_t>& sizes
		) :
			mSizes( sizes ),
			mBlockCounter( 0 ),
			mSize( accumulate( sizes.begin(), sizes.end(), 0 ) ),
			mDirection( unset ),
			mBlockStart( 0 ),
			mBlockEnd( 0 ) {

			// check that weights contain data
			if ( mSize <= 0 ) {
				throw runtime_error( "Input vector for breakpoint weights is empty!" );
			}
		};


		template<typename ParamType>
		void createBlocks( Theta<ParamType>& param );

		void initForward() {
			mDirection = forward;
			mBlockStart = 0;
			mBlockEnd = 0;
			mBlockSize = 0;
			mBlockCounter = 0;
			// initialize block at 0
		}



		// get the end of a block starting at <start> for a given threshold
		// return false if the block end is the last possible value
		inline bool next() {
			if ( mBlockEnd >= mSize ) {
				mDirection = unset;
				return false;
			} else {
				mBlockStart = mBlockEnd;
				mBlockEnd = mBlockStart + mSizes[mBlockCounter];
				mBlockCounter++;
				mBlockSize = mBlockEnd - mBlockStart;
				return true;
			}
		}


		size_t nrBlocks() const {
			return mSizes.size();
		}

		size_t start() const {
			return mBlockStart;
		}

		size_t end() const {
			return mBlockEnd;
		}

		size_t pos() const {
			if ( mBlockCounter > 0 ) {
				return mBlockCounter - 1;
			} else {
				throw runtime_error( "No blocks created yet, position is undefined!" );
			}
		}

		size_t N() const {
			return mBlockSize ;
		}

		size_t T() const {
			return mSize;
		}

		void printBlock() const {
			cout << "[" << mBlockStart << ":" << mBlockEnd << ") " << mBlockSize << " ";
		}


};

#endif
