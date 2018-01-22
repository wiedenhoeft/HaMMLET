#ifndef KAHANAGGREGATOR_HPP
#define KAHANAGGREGATOR_HPP

// This class aggregates values into a sum, while keeping track of the error using Kahan's algorithm. Subtractions are stored separately so that only a single subtraction is used when calling sum()
template <typename T>
class KahanAggregator {
		T mPosSum;
		T mNegSum;
		T mPosError;
		T mNegError;
		size_t mN;	// number of terms


	public:

		// TODO operator+, operator()

		KahanAggregator():
			mPosSum( 0 ),
			mNegSum( 0 ),
			mPosError( 0 ),
			mNegError( 0 ),
			mN( 0 ) {}


		void add( T x, size_t N = 1 ) {
			T y = x - mPosError;
			T temp = mPosSum + y;
			mPosError = ( temp - mPosSum ) - y;
			mPosSum = temp;
			mN += N;
		}


		void subtract( T x, size_t N = 1 ) {
			T y = x - mNegError;
			T temp = mNegSum + y;
			mNegError = ( temp - mNegSum ) - y;
			mNegSum = temp;
			mN += N;
		}

		T sum() const {
			return mPosSum - mNegSum;
		}

		T error() const {
			return mPosError + mNegError;	// TODO minus?
		}

		size_t nrTerms() const {
			return mN;
		}

		void setNrTerms( size_t N ) {
			mN = N;
		}

		void reset() {
			mPosSum = T( 0 );
			mPosError = T( 0 );
			mNegSum = T( 0 );
			mNegError = T( 0 );
			mN = 0;
		}
};

#endif

