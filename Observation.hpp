#ifndef OBSERVATION_HPP
#define OBSERVATION_HPP

#include "Tags.hpp"
#include "includes.hpp"
#include "utils.hpp"
#include <boost/concept_check.hpp>



template <typename DistType>
class Observation {
	template<typename Dummy> friend class Distribution;
	template<typename Dummy> friend class Conjugate;

public:
	Observation() {};

	Observation( const vector<real_t>& vec );


	friend ostream& operator<<(
		ostream& output,
		const Observation<DistType>& D )  {
		output << D.str();
		return output;
	};




	void checkFinite() const;

	bool operator==(
		const Observation< DistType >& other ) const;

	bool operator!=(
		const Observation< DistType >& other ) const {
		return !( ( *this ) == other );
	};

};




////////////////////////////////////////////////// TEMPLATE SPECIALIZATIONS //////////////////////////////////////////////////


//////////////////// NORMAL ////////////////////

template<>
class Observation<Normal> {

	real_t mValue;

public:

	// implicit conversion
// 	operator real_t() {
// 		return mValue;
// 	}

	Observation<Normal>() {
		setValue( NAN );
	}

	Observation<Normal>(
		real_t val ) {
		setValue( val );
	}

	void setValue(
		real_t val ) {
		mValue = val;
	};

	real_t value() const {
		return mValue;
	}

	size_t nrDim() const {
		return 1;
	}

	// TODO str()?

	friend ostream& operator<<(
		ostream& output,
		const Observation<Normal>& D )  {
		output << D.mValue;
		return output;
	}

	void checkFinite() const {
		if ( !isfinite( mValue ) ) {
			throw runtime_error( "Value of Normal observation is not finite!" );
		}
	}
};





//////////////////// NORMAL PARAMETER (mu, sigma^2) ////////////////////



template<>
class Observation<NormalParam> {
	template<typename> friend class Distribution;

// NIG is typically used as a prior for Normal parameters, hence the names if these values:
	real_t mMean;
	real_t mVar;
	real_t mStdev;
	real_t mPrec;

public:

	Observation() {
		mMean = NAN;
		mVar = NAN;
		mStdev = NAN;
		mPrec = NAN;
//         setValue( NAN, NAN );
	}

	Observation( const vector<real_t>& vec ) {
		// vec = {mean, variance}
		if ( vec.size() != 2 ) {
			throw runtime_error( "Need exactly 2 parameters for Normal distribution!" );
		}


		setValue( vec[0], vec[1] );
	};



	Observation(
		real_t mean,
		real_t var ) {
		setValue( mean, var );
	}


	size_t nrDim() const {
		return 1;
	}


	void setValue(
		real_t mean,
		real_t var ) {

		if ( !isfinite( mean ) ) {
			throw runtime_error( "Mean (" + to_string( mean ) + ") must be set to a finite value!" );
		}

		setMean( mean );

		if ( !isfinite( var ) ) {
			throw runtime_error( "Variance(" + to_string( var ) + ") must be set to a finite value!" );
		}

		setVar( var );
	}

	void setMean(
		real_t mean ) {
		mMean = mean;
	}

	void setVar(
		real_t var ) {
		if ( var <= 0 ) {
			throw runtime_error( "Variance ("+to_string(var)+") must be positive!" );
		}


		mVar = var;
		mPrec = 1 / var;
		mStdev = sqrt( var );
	}

// 	real_t logNormalizer() const {
// 		return log( mStdev ) + mMean * mMean / ( 2 * mVar );
// 	}

	real_t mean() const {
		return mMean;
	}
	real_t var() const {
		return mVar;
	}
	real_t stdev() const {
		return mStdev;
	}
	real_t prec() const {
		return mPrec;
	}


	string str(
		const string sep = "\t",
		const string finalSep = ""
	)  const {
		return to_string( mMean ) + sep + to_string( mVar ) + finalSep;
	}


	friend ostream& operator<<(
		ostream& output,
		const Observation<NormalParam>& D )  {
		output << D.str();
		return output;
	};

	void checkFinite() const {
		if ( !isfinite( mMean ) ) {
			throw runtime_error( "Mean of Normal parameters is not finite!" );
		}
		if ( !isfinite( mVar ) ) {
			throw runtime_error( "Variance of Normal parameters is not finite!" );
		}
		if ( !isfinite( mStdev ) ) {
			throw runtime_error( "Standard deviation of Normal parameters is not finite!" );
		}
		if ( !isfinite( mPrec ) ) {
			throw runtime_error( "Precision of Normal parameters is not finite!" );
		}
	}
};




//////////////////// INVERSE GAMMA ////////////////////

template <>
class Observation<InverseGamma> {
	real_t mVariance;

public:

	Observation() {
		setValue( NAN );
	}

	Observation( real_t x ) {
		setValue( x );
	}

	void setValue( real_t x ) {
		mVariance = x;
	}
};







//////////////////// INVERSE GAMMA PARAMETERS ////////////////////
template <>
class Observation<InverseGammaParam> {
	real_t mAlpha;
	real_t mBeta;

public:



	Observation(
		real_t alpha,
		real_t beta ) {
		setValue( alpha, beta );
	}

	void setValue(
		real_t alpha,
		real_t beta ) {
		mAlpha = alpha;
		mBeta = beta;
		checkFinite();
	}

	real_t alpha() const {
		return mAlpha;
	}

	real_t beta() const {
		return mBeta;
	}

// 	real_t logNormalizer() const;

	void checkFinite()const {
		if ( !isfinite( mAlpha ) ) {
			throw runtime_error( "Alpha of Inverse Gamma parameters is not finite!" );
		}
		if ( !isfinite( mBeta ) ) {
			throw runtime_error( "Beta of Inverse Gamma parameters is not finite!" );
		}
	}
};








//////////////////// NORMAL INVERSE GAMMA ////////////////////

template <>
class Observation<NormalInverseGammaParam> {

	real_t mAlpha;
	real_t mBeta;
	real_t mMu0;
	real_t mNu;

public:

	Observation() {
		setValue( NAN, NAN, NAN, NAN );
	}

	Observation(
		const vector<real_t>& v ) {
		if ( v.size() != 4 ) {
			throw runtime_error( "Parameter vector for Normal-Inverse Gamma must have 4 elements!" );
		}
		setValue( v[0], v[1], v[2], v[3] );
	}

	Observation(
		const vector<vector<real_t>>& v ) {
		if ( v.size() != 1 ) {
			throw runtime_error( "Trying to initialize univariate observation with multivariate values!" );
		}

		if ( v[0].size() != 4 ) {
			throw runtime_error( "Initialization for Normal-Inverse gamma must have 4 parameters!" );
		}

		setValue( v[0][0], v[0][1], v[0][2], v[0][3] );
	}

	Observation(
		real_t alpha,
		real_t beta,
		real_t mu0,
		real_t nu ) {

		setValue( alpha, beta, mu0, nu );
	}

	size_t nrDim() const {
		return 1;
	}

	void setValue(
		real_t alpha,
		real_t beta,
		real_t mu0,
		real_t nu
	) {
// 		cout<<alpha<<" "<<beta<<" "<<mu0<<" "<<nu<<endl;
		if ( alpha <= 0 ) {
			throw runtime_error( "Alpha (" + to_string(alpha) + ") must be positive!" );
		}
		if ( beta <= 0 ) {
			throw runtime_error( "Beta (" + to_string(beta) + ") must be positive!" );
		}
		if ( nu <= 0 ) {
			throw runtime_error( "Nu (" + to_string(nu) + ")must be positive!" );
		}
		if ( !isfinite( mu0 ) ) {
			throw runtime_error( "Mu0 (" + to_string(mu0) + ")  must be finite!" );
		}

		mAlpha = alpha;
		mBeta = beta;
		mMu0 = mu0;
		mNu = nu;

	}

	real_t alpha( size_t d = 0 ) const {
		return mAlpha;
	}

	real_t beta( size_t d = 0 ) const {
		return mBeta;
	}

	real_t mu0( size_t d = 0 ) const {
		return mMu0;
	}

	real_t nu( size_t d = 0 ) const {
		return mNu;
	}

// 	real_t logNormalizer() const;

	string str(
		const string sep = "\t",
		const string finalSep = ""
	) const {
		string s;
		s += to_string( mAlpha ) + sep + to_string( mBeta ) + sep + to_string( mMu0 ) + sep + to_string( mNu ) + finalSep;
		return s;
	}


	bool operator==( const Observation< NormalInverseGammaParam >& other ) const {
		if ( alpha() != other.alpha() ) {
			return false;
		}

		if ( beta() != other.beta() ) {
			return false;
		}

		if ( mu0() != other.mu0() ) {
			return false;
		}

		if ( nu() != other.nu() ) {
			return false;
		}

		return true;
	}

	bool operator!=( const Observation< NormalInverseGammaParam >& other ) const {
		return !( ( *this ) == other );
	};

	friend ostream& operator<<(
		ostream& output,
		const Observation< NormalInverseGammaParam >& D )  {
		output << D.str();
		return output;
	}

	size_t domainSize() const {
		return 4;
	}
};











template<>
class Observation<Beta> {

	real_t mValue;

public:

	// implicit conversion
// 	operator real_t() {
// 		return mValue;
// 	}

	Observation<Beta>() {
		setValue( NAN );
	}

	Observation<Beta>(
		real_t val ) {
		setValue( val );
	}

	void setValue(
		real_t val ) {
		mValue = val;
		checkFinite();
	};

	real_t value() const {
		return mValue;
	}

	size_t nrDim() const {
		return 1;
	}

	// TODO str()?

	friend ostream& operator<<(
		ostream& output,
		const Observation<Beta>& D )  {
		output << D.mValue;
		return output;
	}

	void checkFinite() const {
		if ( !isfinite( mValue ) ) {
			throw runtime_error( "Value of Normal observation is not finite!" );
		}
	}
};





template<>
class Observation<BetaParam> {

	real_t mAlpha;
	real_t mBeta;

public:

	// implicit conversion
// 	operator real_t() {
// 		return mValue;
// 	}

	Observation<BetaParam>() {
		setValue( NAN, NAN );
	}

	Observation<BetaParam>(
		real_t alpha,
		real_t beta ) {
		setValue( alpha, beta );
	}

	Observation<BetaParam>(
		vector<real_t> args ) {
		if ( args.size() != 2 ) {
			throw runtime_error( "Input vector for Beta parameters must be of size 2!" );
		}
		setValue( args[0], args[1] );
	}

	void setValue(
		real_t alpha, real_t beta ) {
		mAlpha = alpha;
		mBeta = beta;
		checkFinite();
	};


	size_t nrDim() const {
		return 1;
	}

	// TODO str()?

	real_t alpha() const {
		return mAlpha;
	}

	real_t beta() const {
		return mBeta;
	}

	friend ostream& operator<<(
		ostream& output,
		const Observation<BetaParam>& D )  {
		output << D.alpha() << " " << D.beta();
		return output;
	}

	void checkFinite() const {
		if ( !isfinite( mAlpha ) ) {
			throw runtime_error( "Alpha of Beta distribution parameter is not finite!" );
		}
		if ( !isfinite( mBeta ) ) {
			throw runtime_error( "Beta of Beta distribution parameter is not finite!" );
		}
	}
};

//////////////////// CATEGORICAL ////////////////////

template <>
class Observation<Categorical> {
	int mValue;	// TODO size_t?

public:

	// implicit conversion
	operator int() {
		return mValue;
	}

	Observation() {};

	Observation( int val ) {
		mValue = val;
	}

	int value() const {
		return mValue;
	}

	size_t nrDim() const {
		return 1;
	}

};


//////////////////// CATEGORICAL VECTOR ////////////////////
//  implement if necessary






//////////////////// DIRICHLET ////////////////////


template<>
class Observation<Dirichlet> {
	template<typename> friend class Distribution;
	template<typename> friend class Conjugate;
	vector<real_t> mProbs;


public:

	Observation( ) {	};

	// zero-instantiation by size
	Observation( size_t size ) {
		mProbs.assign( size, NAN );
	};

	Observation( size_t size, real_t x ) {
		mProbs.assign( size, x );
	};

	Observation( const vector<real_t>& vec ) {
		mProbs = vec;
	}

	size_t domainSize() const {
		return mProbs.size();
	}

	size_t nrDim() const {
		return 1;
	}

	real_t& operator[]( size_t i ) {
		if ( i >= mProbs.size() ) {
			throw runtime_error( "Probability index out of bounds!" );
		}
		return mProbs[i];
	}

	const real_t& operator[]( size_t i ) const {
		if ( i >= mProbs.size() ) {
			throw runtime_error( "Probability index out of bounds!" );
		}
		return mProbs[i];
	}


	const vector<real_t>& probs() const {
		return mProbs;
	}

	vector<real_t>& probs()  {
		return mProbs;
	}

	void setValue( const vector<real_t>& v ) {
		mProbs = v;
	}


	
	string str(
		const string sep = "\t",
		const string finalSep = ""
	) const {
		return concat( mProbs, sep, finalSep );
	}


	friend ostream& operator<<(
		ostream& output,
		const Observation<Dirichlet>& D )  {
		output << D.str();
		return output;
	};

};




//TODO this should be removed
////////// DIRICHLET VECTOR  //////////


template <>
class Observation<DirichletVector> {
	template<typename> friend class Distribution;
	template<typename> friend class Conjugate;
	vector<Observation<Dirichlet>> mRows;


public:

	Observation() {};

	// quadratic NAN-instantiation by size
	Observation( size_t size ) {
		mRows.assign( size, Observation<Dirichlet>( size ) );
	};

	Observation( size_t size, size_t domainSize ) {
		mRows.assign( size, Observation<Dirichlet>( domainSize) );
	};
	
	Observation( const vector<vector<real_t>>& vec ) {
		mRows.reserve( vec.size() );

		for ( const auto & v : vec ) {
			mRows.push_back( Observation<Dirichlet>( v ) );
		}
	}


	size_t nrDim() const {
		return mRows.size();
	}


	const Observation<Dirichlet>& operator[]( size_t i ) const {
		if ( i >= mRows.size() ) {
			throw runtime_error( "Dirichlet index out of bounds!" );
		}
		return mRows[i];
	}

	Observation<Dirichlet>& operator[]( size_t i ) {
		if ( i >= mRows.size() ) {
			throw runtime_error( "Dirichlet index out of bounds!" );
		}
		return mRows[i];
	}


	const real_t& operator()( size_t from, size_t to ) const {
		if ( from >= mRows.size() ) {
			throw runtime_error( "Dirichlet vector row index out of bounds!" );
		}
		if ( to >= mRows[from].domainSize() ) {
			throw runtime_error( "Dirichlet vector column index out of bounds!" );
		}
		return mRows[from][to];
	}

	real_t& operator()( size_t from, size_t to ) {
		if ( from >= mRows.size() ) {
			throw runtime_error( "Dirichlet vector row index out of bounds!" );
		}
		if ( to >= mRows[from].domainSize() ) {
			throw runtime_error( "Dirichlet vector column index out of bounds!" );
		}
		return mRows[from][to];
	}

	string str(
		const string sep = "\t",
		const string finalSep = ""
	) const  {
		return concat( mRows, sep, finalSep );
	}
};






////////// DIRICHLET PARAMETERS //////////


template <>
class Observation<DirichletParam> {
	template<typename> friend class Distribution;
	template<typename> friend class Conjugate;
	vector<real_t> mAlphas;


public:

	Observation( ) {	};


	// zero-instantiation by size
	Observation( size_t size ) {
		mAlphas.assign( size, NAN );
	};


	Observation(
		size_t size,
		real_t value ) {
		if ( value <= 0 ) {
			throw runtime_error( "All values in Dirichlet parameters must be greater than zero!" );
		}
		mAlphas.assign( size, value );
	};



	Observation(
		size_t size,
		real_t value,
		size_t specialIndex,
		real_t specialValue ) {

		if ( specialIndex >= size ) {
			throw runtime_error( "Special index out of bounds for Dirichlet parameter." );
		}
		if ( value <= 0 ) {
			throw runtime_error( "All values in Dirichlet parameters must be greater than zero!" );
		}
		if ( specialValue <= 0 ) {
			throw runtime_error( "All values in Dirichlet parameters must be greater than zero!" );
		}

		mAlphas.assign( size, value );
		if ( specialIndex >= mAlphas.size() ) {
			throw runtime_error( "Special index out of bounds!" );
		}
		mAlphas[specialIndex] = specialValue;
	};


	Observation( const vector<real_t>& vec ) {
		mAlphas = vec;
	}

	size_t domainSize() const {
		return mAlphas.size();
	}

	size_t nrDim() const {
		return 1;
	}

	const vector<real_t>& alphas() const {
		return mAlphas;
	}


	const real_t& operator[]( size_t i ) const {
		if ( i >= mAlphas.size() ) {
			throw runtime_error( "Alpha index out of bounds!" );
		}
		return mAlphas[i];
	}

	real_t& operator[]( size_t i )  {
		if ( i >= mAlphas.size() ) {
			throw runtime_error( "Alpha index out of bounds!" );
		}
		return mAlphas[i];
	}


	string str(
		const string sep = "\t",
		const string finalSep = ""	) const {
		return concat( mAlphas, sep, finalSep );
	}


	friend ostream& operator<<(
		ostream& output,
		const Observation<DirichletParam>& D )  {
		output << D.str();
		return output;
	};


	bool operator==( const Observation< DirichletParam>& other ) const {
		return ( mAlphas == other.alphas() );
	};

};








////////// DIRICHLET PARAMETER VECTOR //////////


template <>
class Observation<DirichletParamVector> {
	template<typename> friend class Distribution;
	template<typename> friend class Conjugate;
	vector<Observation<DirichletParam>> mRows;

public:

	Observation() {};

	Observation(
		size_t size,
		real_t value ) {
		for ( size_t r = 0; r < size; ++r ) {
			Observation<DirichletParam> p( size, value );
			mRows.push_back( p );
		}
	};

	Observation(
		size_t size,
		real_t value,
		real_t diagonalValue ) {
		size_t d = 0;

		for ( size_t r = 0; r < size; ++r ) {
			Observation<DirichletParam> p( size, value, d, diagonalValue );
			mRows.push_back( p );
			++d;
		}
	}

	// quadratic zero-instantiation by size
	Observation( size_t size ) {
		mRows.assign( size, Observation<DirichletParam>( size ) );
	};

	Observation( const vector<vector<real_t>>& vec ) {
		mRows.reserve( vec.size() );

		for ( const auto &  v : vec ) {
			mRows.push_back( Observation<DirichletParam>( v ) );
		}
	}


	size_t nrDim() const {
		return mRows.size();
	}


	const Observation<DirichletParam>& operator[]( size_t i ) const {
		if ( i >= mRows.size() ) {
			throw runtime_error( "Dirichlet parameter index out of bounds!" );
		}
		return mRows[i];
	}

	Observation<DirichletParam>& operator[]( size_t i ) {
		if ( i >= mRows.size() ) {
			throw runtime_error( "Dirichlet parameter index out of bounds!" );
		}
		return mRows[i];
	}


	string str(
		const string sep = "\t",
		const string finalSep = "" ) const {
		return concat( mRows, sep, finalSep );
	}


};








// TODO checkFinite in constructors



#endif











