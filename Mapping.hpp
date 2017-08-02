#ifndef MAPPING_HPP
#define MAPPING_HPP

#include "Tags.hpp"
#include "includes.hpp"
#include "Parser.hpp"

// Mapping[s][d] is the parameter dimension for state s and data dimension d
// TODO this must only be a many-to-one mapping; there should be one from data dimensions to parameter dimensions, one from parameter to hyperparameter, and also a method to create a compound mapping to update hyperparameters from data directly. States are just indices for a vector of parameter mappings. The mapping to hyperparameters makes it easy to have one prior per parameter or a joint prior such as a Dirichlet Process Prior.

// for Parser
template<>
MappingType convertType(
    const string& s ) {

	if ( s == "combinations" || s == "C" ) {
		return combinations;
	} else {
		throw runtime_error( "Unknown mapping type " + s + "!" );
	}
}



// helper function to get the number of states, to be used in the initializer list of the constructor so the member can be const
size_t nrOfStates(
    size_t nrDataDim,
    size_t nrParam,
    MappingType mappingType ) {
	size_t result;

	switch ( mappingType ) {
		case combinations:
			result = pow( nrParam, nrDataDim );
			break;

		default:
			throw runtime_error( "Mapping type not implemented!" );
	}

	if ( result <= 1 ) {
		throw runtime_error( "Requested parameters would yield an HMM with less than 2 states!" );
	}

	return result;
}

// maps a state s to the indices in theta and tau_theta




class Mapping {
		vector<vector<size_t>> mValue;

		const size_t mNrDataDim;
		const size_t mNrParams;
		const size_t mNrStates;

	public:

		// delete copy constructor
		Mapping( const Mapping& that ) = delete;

		Mapping(
		    const size_t nrdatadim,
		    const size_t nrparams,	// number of parameters
		    const MappingType mappingType )
			:
			mNrDataDim( nrdatadim ),
			mNrParams( nrparams ),
			mNrStates( nrOfStates( nrdatadim, nrparams, mappingType ) ) {


			if ( mNrDataDim <= 0 ) {
				throw runtime_error( "Number of data dimensions must be positive!" );
			}

			if ( mNrParams <= 0 ) {
				throw runtime_error( "Number of parameters must be positive!" );
			}

			if ( mNrStates <= 0 ) {
				throw runtime_error( "Number of states must be positive!" );
			}


			// setup the mapping using vectors of reference wrappers
			switch ( mappingType ) {

				case combinations:	// all combinations of shared parameters

					// states are assigned by generating reversed nrParam-ary numbers of nrDataDim digits
					for ( size_t x = 0; x < mNrStates; ++x ) {
						size_t n = x;
						vector<size_t> mapping;
						mapping.reserve( nrdatadim );

						for ( size_t d = 0; d < nrdatadim; ++d ) {
							auto residue = n % mNrParams;
							n /= mNrParams;
							mapping.push_back( residue );
						}

						mValue.push_back( mapping );
					}

					break;

				case independent:
					throw runtime_error( "Mapping type \"independent\" not implemented yet!" );
					break;

				default:
					throw runtime_error( "Unknown mapping type!" );
					break;
			}
		};


		const vector<size_t>& operator[](
		    size_t state ) const {
			return mValue[state];
		}

		size_t nrStates() const {
			return mNrStates;
		}

		size_t nrParams() const {
			return mNrParams;
		}

		size_t nrDataDims() const {
			return mNrDataDim;
		}
};

#endif


