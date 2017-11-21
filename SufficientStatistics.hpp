#ifndef SUFFICIENTSTATISTICS_HPP
#define SUFFICIENTSTATISTICS_HPP

#include "Tags.hpp"
#include "includes.hpp"
#include "Observation.hpp"

//NOTE empty constructors should create objects that are meaningful for not having any observatgions, e.g. all 0 for <Normal>


template <typename DistType>
class SufficientStatistics;


template<typename Type>
inline SufficientStatistics<Type> operator+(
    SufficientStatistics<Type> lhs,
    const SufficientStatistics<Type>&  rhs ) {
	if ( lhs.nrDim() != rhs.nrDim() ) {
		throw runtime_error( "Cannot add sufficient statistics for categorical of different dimensions!" );
	}

	lhs += rhs;
	return lhs;
}


template<typename Type>
inline SufficientStatistics<Type> operator-(
    SufficientStatistics<Type> lhs,
    const SufficientStatistics<Type>&  rhs ) {
	if ( lhs.nrDim() != rhs.nrDim() ) {
		throw runtime_error( "Cannot add sufficient statistics for categorical of different dimensions!" );
	}

	lhs -= rhs;
	return lhs;
}



// TODO define common methods outside?



//////////////////// UNIVARIATE NORMAL ////////////////////

template <>
class SufficientStatistics<Normal> {
		real_t mSum;
		real_t mSumSq;



	public:

		SufficientStatistics() {
			clear();
		};

		SufficientStatistics(
		    real_t singleValue ) {
			mSum = singleValue;
			mSumSq = singleValue * singleValue;
		}

		SufficientStatistics(
		    vector<real_t> vec,
		    size_t begin = 0,
		    size_t end = 0 ) {
			if ( end == 0 || end > vec.size() ) {
				end = vec.size();
			}
			clear();
			for ( auto i = begin; i < end; ++i ) {
				addObs( vec[i] );
			}
		}

		SufficientStatistics(
		    const real_t sum,
		    const real_t sumSq
		) {
			mSum = sum;
			mSumSq = sumSq;
		}

		inline void addObs( const real_t x ) {
			mSum += x;
			mSumSq += x * x;
		}


		size_t nrDim() const {
			return 1;
		}

		real_t sum() const {
			return mSum;
		}

		real_t sumSq() const {
			return mSumSq;
		}


		SufficientStatistics<Normal>&  operator+=(
		    const SufficientStatistics<Normal>&  rhs )  {
			mSum += rhs.sum();
			mSumSq += rhs.sumSq();
// 			if (mSumSq <0){
// 				throw runtime_error( "Sum of squares is negative after addition  (" + to_string( mSumSq - rhs.sumSq() ) + "+" + to_string( rhs.sumSq() ) + "=" + to_string( mSumSq ) + "), possibly cause by overflow!");
// 			}
			return *this;
		}

		SufficientStatistics<Normal>&  operator-=(
		    const SufficientStatistics<Normal>&  rhs )  {
			mSum -= rhs.sum();
			mSumSq -= rhs.sumSq();
// 			if ( mSumSq < 0 ) {
// 				throw runtime_error( "Sum of squares is negative after subtraction  (" + to_string( rhs.sumSq() + mSumSq ) + "-" + to_string( rhs.sumSq() ) + "=" + to_string( mSumSq ) + ")!" );	// TODO numerics?
// 			}
			return *this;
		}


		void clear() {
			mSum = 0;
			mSumSq = 0;
		}


// TODO move outside?
		friend ostream& operator<<(
		    ostream& output,
		    const SufficientStatistics<Normal>& D )  {
			output << D.mSum << "\t" << D.mSumSq;
			return output;
		};

};









//////////////////// Categorical ////////////////////


template <>
class SufficientStatistics<Categorical> {
		vector<size_t> mCounts;

	public:

		SufficientStatistics( size_t domainsize ) {
			mCounts.assign( domainsize, 0 );
		}

		SufficientStatistics( vector<size_t>& counts ) {
			mCounts = counts;
		}

		SufficientStatistics<Categorical>&  operator+=(
		    const SufficientStatistics<Categorical>&  rhs )  {
			if ( domainSize() != rhs.domainSize() ) {
				throw runtime_error( "Cannot add sufficient statistics for categorical of different domain sizes!" );
			}

			for ( size_t i = 0; i < mCounts.size(); ++i ) {
				mCounts[i] += rhs[i];
			}

			return *this;
		}

		SufficientStatistics<Categorical>&  operator-=(
		    const SufficientStatistics<Categorical>&  rhs )  {
			if ( domainSize() != rhs.domainSize() ) {
				throw runtime_error( "Cannot add sufficient statistics for categorical of different domain sizes!" );
			}

			for ( size_t i = 0; i < mCounts.size(); ++i ) {
				mCounts[i] -= rhs[i];
			}

			return *this;
		}

		size_t nrDim() const {
			return 1;	// NOTE This is correct, these are the sufficient statistics for a 1-D categorical variable.
		}


		// how many different values can the underlying categorical variable take?
		size_t domainSize() const {
			return mCounts.size();
		}

		const size_t& operator[]( size_t i ) const {
			return mCounts[i];
		}

		size_t& operator[]( size_t i ) {
			return mCounts[i];
		}

		void clear() {
			mCounts.assign( mCounts.size(), 0 );
		}

		string str() const {
			return concat( mCounts );
		}
};



//////////////////// Categorical Vector ////////////////////
// NOTE if nrDim() == domainSize(), then this is a symmetric matrix that can be used for counting transitions


template <>
class SufficientStatistics<CategoricalVector> {
		vector < SufficientStatistics<Categorical>> mCounts;

	public:

		SufficientStatistics(
		    size_t nrdim,
		    size_t domainsize = 0 ) {
			if ( domainsize == 0 ) {
				domainsize = nrdim;
			}

			for ( size_t d = 0; d < nrdim; ++d ) {
				mCounts.push_back( SufficientStatistics< Categorical >( domainsize ) );
			}
		}

		size_t nrDim() const {
			return mCounts.size();
		}



		const SufficientStatistics<Categorical>& operator[]( int i ) const {
			return mCounts[i];
		}


		SufficientStatistics<Categorical>& operator[]( int i )  {
			return mCounts[i];
		}


		SufficientStatistics<CategoricalVector>&  operator+=(
		    const SufficientStatistics<CategoricalVector>&  rhs )  {

			for ( size_t i = 0; i < nrDim(); ++i ) {
				mCounts[i] += rhs[i];
			}

			return *this;
		}

		SufficientStatistics<CategoricalVector>&  operator-=(
		    const SufficientStatistics<CategoricalVector>&  rhs )  {

			for ( size_t i = 0; i < nrDim(); ++i ) {
				mCounts[i] -= rhs[i];
			}

			return *this;
		}



		string str() const {
			string result = "";
			for ( auto & c : mCounts ) {
				result += c.str() + "\n";
			}
			return result;
		}


		void clear() {
			for ( auto & c : mCounts ) {
				c.clear();
			}
		}

};





//////////////////// Geometric ////////////////////


template <>
class SufficientStatistics<Geometric> {
		size_t mSum;

	public:

		SufficientStatistics() {
			mSum = 0;
		}

		SufficientStatistics( real_t& sum ) {	// TODO parsing real_t for count data is not the best option
			mSum = sum;
		}

		SufficientStatistics( size_t& sum ) {
			mSum = sum;
		}


		SufficientStatistics(
		    vector<real_t> vec,
		    size_t begin = 0,
		    size_t end = 0 ) {
			if ( end == 0 || end > vec.size() ) {
				end = vec.size();
			}
			clear();
			for ( auto i = begin; i < end; ++i ) {
				mSum += vec[i] ;
			}
		}


		SufficientStatistics<Geometric>&  operator+=(
		    const SufficientStatistics<Geometric>&  rhs )  {
			if ( domainSize() != rhs.domainSize() ) {
				throw runtime_error( "Cannot add sufficient statistics for categorical of different domain sizes!" );
			}

			mSum += rhs.sum();

			return *this;
		}


		SufficientStatistics<Geometric>&  operator-=(
		    const SufficientStatistics<Geometric>&  rhs )  {
			if ( domainSize() != rhs.domainSize() ) {
				throw runtime_error( "Cannot add sufficient statistics for categorical of different domain sizes!" );
			}

			mSum -= rhs.sum();

			return *this;
		}


		size_t nrDim() const {
			return 1;	// NOTE This is correct, these are the sufficient statistics for a 1-D categorical variable.
		}


		// how many different values can the underlying categorical variable take?
		size_t domainSize() const {
			return 1;
		}

		size_t sum() const {
			return mSum;
		}

		void clear() {
			mSum = 0;
		}

		string str() const {
			return to_string( mSum );
		}

};

#endif













