#ifndef SPLITTABLEBLOCKS_HPP
#define SPLITTABLEBLOCKS_HPP

template<>
class Blocks<Splittable> {
		const size_t mSize;
		deque<size_t> mSizes;
		size_t mStart;
		size_t mEnd;
		size_t mPos;	// the current position as seen from the outside.
		size_t mIndex; // the index in mSizes corresponding to mPos
		bool mFirstIter;

		void rotate() {
			mSizes.push_back( mSizes.front() );
			mSizes.pop_front();
			if ( mIndex > 0 ) {
				mIndex--;
			} else {
				mIndex = mSizes.size() - 1;
			}
		}

	public:


		Blocks( size_t size ):
			mSize( size ),
			mPos( 0 ),
			mStart( 0 ),
			mEnd( 0 ),
			mIndex( 0 ),
			mFirstIter( true ) {
			if ( size == 0 ) {
				throw runtime_error( "Cannot create block structure for 0 positions!" );
			}
			mSizes.push_back( size );
		}

		// TODO is there any way to swap this?
		Blocks( const vector<size_t>& sizes ):
			mSize( accumulate( sizes.begin(), sizes.end(), 0 ) ),
			mPos( 0 ),
			mStart( 0 ),
			mEnd( 0 ),
			mIndex( 0 ),
			mFirstIter( true ) {
			for ( auto & x : sizes ) {
				mSizes.push_back( x );
			}
		}

		// Split the current block such that the first new block is of size s, and move to the first block; this preserves pos(). If s is >= the size of the block, an exception is thrown
		void split( size_t s ) {
			if ( mSizes[mIndex] <= s ) {
				throw runtime_error( "Cannot split block into this size!" );
			}

			// rotate until current position is at the front
			while ( mIndex != 0 ) {
				rotate();
			}

			mSizes.push_back( s );
			mSizes[0] = mSizes[0] - s;
			mIndex = mSizes.size() - 1;
		};


		// Move to the first position, so pos()==0.
		void initForward() {
			mFirstIter = true;
			if ( mIndex >= mPos ) {
				mIndex -= mPos;
			} else {
				mIndex += ( mSizes.size() - mPos );
			}
			mPos = 0;
			mStart = 0;
			mEnd = 0;
		}

		// move to the next block, in modular fashion (wrap around)
		// returns false if wrapped around, otherwise true.
		bool next() {
			if ( mFirstIter ) {
				mFirstIter = false;
				mEnd = mStart + mSizes[mIndex];
			} else {
				mIndex = ( mIndex + 1 ) % mSizes.size();
				mPos++;
				mStart = mEnd;
				mEnd += mSizes[mIndex];
			}
			if ( mPos == mSizes.size() ) {
				mPos = 0;
				return false;
			} else {
				return true;
			}
		};

		size_t start() const {
			return mStart;
		}

		size_t end() const {
			return mEnd;
		}


		// return the number of blocks
		size_t nrBlocks()const {
			return mSizes.size();
		};


		// return the current block size
		size_t size() const {
			return mSizes[mIndex];
		};


		// return the index of the current block
		size_t pos() const {
			return mPos % mSizes.size();
		}
};


#endif
