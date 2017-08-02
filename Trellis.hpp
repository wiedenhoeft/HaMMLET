#ifndef TRELLIS_HPP
#define TRELLIS_HPP

#include <stdexcept>

#include "uintmath.hpp"

class Trellis {

		vector<real_t> mVec;
		size_t mNrStates;
		rng_t& mRNG;

		void assertRange( size_t d ) const {
			if ( d >= mNrStates ) {
				throw runtime_error( "Trellis dimension index out of bounds!" );
			}
		}
	public:

		// delete copy constructor
		Trellis( const Trellis& that ) = delete;

		Trellis( rng_t& RNG ):
			mNrStates( 2 ) ,
			mRNG( RNG )
		{};

		Trellis( size_t nrStates, rng_t& RNG ):
			mNrStates( nrStates ) ,
			mRNG( RNG )
		{};


		real_t& operator()( size_t t, size_t d ) {
			assertRange( d );
			return mVec[t * mNrStates + d];
		}


		// return reference to last element at dimension d
		real_t& back( size_t d ) {
			return mVec[mVec.size() - mNrStates + d];
		}


		// this interface is for cases when the number of states is not known beforehand, e.g. for Dirichlet Process Priors
		void setNrStates( size_t K ) {
			mNrStates = K;
		}

		size_t size() const {
			return divide( mVec.size(), mNrStates );
		}


		void push_back( const vector<real_t>& vec ) {
			mVec.insert( mVec.end(), vec.begin(), vec.end() );
		}

		size_t sample( size_t t ) const {
			discrete_distribution<size_t> dist( mVec.begin() + ( t * mNrStates ), mVec.begin() + ( ( t + 1 )*mNrStates ) );
			size_t result = dist( mRNG );
			dist.reset();
			return result;
		}

		void reserve( size_t N ) {
			mVec.reserve( N * mNrStates );
		}


		void clear() {
			mVec.clear();
		}
};


#endif
