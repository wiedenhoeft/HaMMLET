#ifndef MULTIVECTOR_HPP
#define MULTIVECTOR_HPP

#include <vector>
using std::vector;

#include <stdexcept>
using std::runtime_error;

template <typename T>
class MultiVector {
		vector<T> mVec;
		const size_t mNrDim;

	public:


		MultiVector(
		    size_t nrDim
		):
			mNrDim( nrDim ) {
			if ( mNrDim <= 0 ) {
				throw runtime_error( "Number of dimensions in multivector must be positive!" );
			}
		}

		MultiVector(
		    T entry,
		    size_t size,
		    size_t nrDim
		):
			mVec( entry, nrDim* size ),
			mNrDim( nrDim ) {
			if ( mNrDim <= 0 ) {
				throw runtime_error( "Number of dimensions in multivector must be positive!" );
			}
		}


		// direct access to underlying vector
		inline T& operator[]( size_t i ) {
			if ( i >= mVec.size() ) {
				throw runtime_error( "Direct index out of bounds for multivector!" );
			}
			return mVec[i];
		}

		inline const T& operator[]( size_t i ) const {
			if ( i >= mVec.size() ) {
				throw runtime_error( "Direct index out of bounds for multivector!" );
			}
			return mVec[i];
		}


		inline T& operator()( size_t pos, size_t dim ) {
			if ( dim >= mNrDim ) {
				throw runtime_error( "Multivector dimension index out of bounds!" );
			}
			const size_t i = pos * mNrDim + dim;
			if ( i < mVec.size() ) {
				return mVec[i];
			} else {
				throw runtime_error( "Multivector index out of bounds!" );
			}
		}

		size_t size() const {
			return mVec.size() / mNrDim;
		}

		void reserve( size_t size ) {
			mVec.reserve( mNrDim * size );
		}


		void resize( size_t N ) {
			mVec.resize( N * mNrDim );
		}


		void push_back( T& entry ) {
			mVec.reserve( mVec.size() + mNrDim );
			for ( size_t d = 0; d < mNrDim; ++d ) {
				mVec.push_back( entry );
			}
		}


		void swap( vector<T>& vec ) {
			const size_t s = mVec.size();
			if ( s != ( s / mNrDim )*mNrDim ) {
				throw runtime_error( "Cannot swap into multivector, size is not a multiple of dimensions!" );
			}
			mVec.swap( vec );
		}

		size_t nrDim()const {
			return mNrDim;
		}
};

#endif
