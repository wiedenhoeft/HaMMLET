#ifndef THETAHYPERPARAM_HPP
#define THETAHYPERPARAM_HPP

#include "Tags.hpp"
#include "Mapping.hpp"
#include "Conjugate.hpp"
#include "Observation.hpp"
#include "SufficientStatistics.hpp"




template <typename ParamType>	// e.g. NormalInverseGammaParam
class ThetaHyperParam {

		const size_t mNrParams;

		vector<Conjugate<ParamType>> mParams;

	public:


		// automatic priors
		ThetaHyperParam(
		    const vector< vector< real_t > >& hyperparams
		) : mNrParams( hyperparams.size() ) {

			if ( mNrParams <= 0 ) {
				throw runtime_error( "Number of hyperparameters must be positive" );
			}

			if ( mNrParams <= 0 ) {
				throw runtime_error( "Number of emission hyperparameters must be positive! Did you forget to provide them, or to use -a?" );
			}


			for ( const auto & hp : hyperparams ) {
				mParams.push_back( Conjugate<ParamType>( hp ) );
			}
		}

		ThetaHyperParam(
		    const vector < Observation<ParamType>>& hyperparams
		) : mNrParams( hyperparams.size() ) {

			if ( mNrParams <= 0 ) {
				throw runtime_error( "Number of hyperparameters must be positive" );
			}

			for ( const auto & hp : hyperparams ) {
				mParams.push_back( Conjugate<ParamType>( hp ) );
			}
		}


		ThetaHyperParam(
		    const  Observation<ParamType>& hyperparams,
		    const size_t nrDim
		) : mNrParams( nrParams ),
			mParams( hyperparams, nrDim ) {

			if ( mNrParams <= 0 ) {
				throw runtime_error( "Number of hyperparameters must be positive" );
			}

		}

		size_t nrParams() const {
			return mNrParams;
		}



		////////// accessors //////////
		//NOTE round parentheses access the data through the mapping, square brackets access the parameters directly

		template<typename EmissionsType>
		inline void addObservation(
		    const SufficientStatistics<EmissionsType>& suffStat,
		    const size_t N,
		    const size_t dim ) {

			mParams[dim].addObservation( suffStat, N );
		}


		// TODO some objects use posterior(), make consistent
		// TODO implicit conversion?
		const Observation<ParamType>& posterior(
		    const size_t d	) const {
			return mParams[d].posterior();
		}


		const Observation<ParamType>& prior(
		    const size_t d	) const {
			return mParams[d].prior();
		}

		void reset() {
			for ( auto & p : mParams ) {
				p.reset();
			}
		}


		string str() const {
			return concat( mParams, "\t", "\n" );
		}

		const Conjugate<ParamType> operator[]( size_t d ) const {
			return mParams[d];
		}
};




////////////////////////////////////////////////// TEMPLATE SPECIALIZATIONS //////////////////////////////////////////////////
// #include "WaveletTree.hpp"	// required implementation



#endif



