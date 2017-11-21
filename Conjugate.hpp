// Class file for conjugate parameter pairs (prior and posterior)


#ifndef CONJUGATE_HPP
#define CONJUGATE_HPP

#include "includes.hpp"
#include "Tags.hpp"
#include "Observation.hpp"
#include "SufficientStatistics.hpp"




// parameters of ParamType will be updated using sufficient statistics of ObsType
template <typename ParamType>
class Conjugate {

		// parameters
		Observation<ParamType> mPrior;
		Observation<ParamType> mPosterior;

	public:

		////////// constructors //////////


		// A conjugate is initialized using the same parameters as an observation of the same type
		template<typename ... Types>
		Conjugate( Types ... args ) : mPrior( args ... ), mPosterior( args ... ) {};


		Conjugate( Observation<ParamType> prior ): mPrior( prior ), mPosterior( prior ) {}



		const Observation<ParamType>& prior() const {
			return mPrior;
		}

		Observation<ParamType>& prior() {
			return mPrior;
		}

		const Observation<ParamType>& posterior() const {
			return mPosterior;
		}

		Observation<ParamType>& posterior() {
			return mPosterior;
		}


		size_t domainSize() const {
			return mPrior.domainSize();
		}


		size_t nrDim() const {
			return mPrior.nrDim();
		}




		template<typename ObsType>
		inline void addObservation(
		    const SufficientStatistics<ObsType>& obs );

		template<typename ObsType>
		inline void addObservation(
		    const SufficientStatistics<ObsType>& obs,
		    const size_t N	);


		template<typename ObsType>
		inline void addObservation(
		    const SufficientStatistics<ObsType>& suffstat,
		    const size_t N,
		    const vector<int>& mapping );


		template<typename ObsType>
		inline void addObservation(
		    const Observation<ObsType>& obs );


		void reset() {
			mPosterior = mPrior;
		}

		string str( string sep = " " ) const {	// TODO sep
			return " " + mPosterior.str() + " (prior: " + mPrior.str() + ")";
		}

		friend ostream& operator<<(
		    ostream& output,
		    const Conjugate<ParamType>& D )  {
			output << D.str();
			return output;
		};

};







////////////////////////////////////////////////// TEMPLATE SPECIALIZATIONS //////////////////////////////////////////////////
#include "SufficientStatistics.hpp"	// required implementation
#include"Observation.hpp"	// required implementation



////////// Normal-Inverse Gamma parameters //////////


template<> template<>
void Conjugate<NormalInverseGammaParam>::addObservation(
    const SufficientStatistics< Normal>& obs,
    const size_t counts ) {

	const real_t sum = obs.sum();
	const real_t sumSq = obs.sumSq();

	if ( sumSq < 0 ) {
		throw runtime_error( "Sum of squares is negative (" + to_string( sumSq ) + ")!" );
	}

	if ( counts == 0 ) {
		if ( sumSq > 0 ) {
			throw runtime_error( "Sufficient statistics contain values, but no observation count!" );
		}
		cout << "[WARNING] No observation count for sufficient statistics, there might be an index error!" << endl;
		// TODO this shouldn't occur anyway...
		return;
	}

	const double N = ( double )counts;
	const real_t xbar = sum / N;	// sample mean


	const real_t alpha = mPosterior.alpha( );
	const real_t beta = mPosterior.beta( );
	const real_t mu0 = mPosterior.mu0( );
	const real_t nu = mPosterior.nu( );



	// NOTE sometimes (sum*sum)/N > sumSq, i.e. sample variance is negative. This is akin to the numerically unstable way to calculate the sample variance naively.
	real_t ssN = ( sum * sum ) / N;
	if ( ssN > sumSq ) {
		ssN = sumSq;
		cerr << "[WARNING] Encountered numerical instability in sample variance!" << endl;
	}

	// TODO better numerics, especially for beta?
	mPosterior.setValue(
	    alpha + N / 2.0,	// alpha	NOTE using 2 instead of 2.0 previously caused the strange "label-switching" bug in maxBlockLen=1
	    // NOTE This involves a term which makes the naive calculation of sample variance numerically unstable (sumSq-ssN); we subtract ssN last in hopes that the term added to sumSq is large enough to alleviate catastrophic cancellation.
	    beta + ( ( sumSq   + ( N * nu / ( N + nu ) ) * ( ( xbar - mu0 ) * ( xbar - mu0 ) ) )  -  ssN ) / 2.0, 	// beta
	    ( nu * mu0 + sum ) / ( nu + N ), 	// mu0
	    nu + N	// nu	// TODO size_t?
	);
}






////////// Dirichlet parameters //////////

template<> template<>
void Conjugate<DirichletParamVector>::addObservation(
    const SufficientStatistics<CategoricalVector>& countMatrix ) {

	if ( countMatrix.nrDim() != mPosterior.nrDim() ) {
		throw runtime_error( "Dimensions of count matrix (" + to_string( countMatrix.nrDim() ) + ") and posterior observations (" + to_string( mPosterior.nrDim() ) + ") do not match!" );
	}

	for ( size_t d = 0; d < countMatrix.nrDim() ; ++d ) {
		for ( size_t c = 0; c < countMatrix[d].domainSize(); ++c ) {
			mPosterior[d][c] += countMatrix[d][c];
		}
	}
}



template<> template<>
void Conjugate<DirichletParam>::addObservation(
    const SufficientStatistics<Categorical>& obs ) {

	if ( obs.domainSize() != mPosterior.domainSize() ) {
		throw runtime_error( "Domain size of observations (" + to_string( obs.domainSize() ) + ") does not match that of posterior (" + to_string( mPosterior.domainSize() ) + ")!" );
	}

	for ( size_t i = 0; i < mPosterior.domainSize(); ++i ) {
		mPosterior.mAlphas[i] += obs[i];
	}
}



template<> template<>
void Conjugate<BetaParam>::addObservation(
    const SufficientStatistics<Geometric>& obs,
    const size_t N
) {
	mPosterior.setValue( mPosterior.alpha() + N, mPosterior.beta() + obs.sum() );
}



#endif

