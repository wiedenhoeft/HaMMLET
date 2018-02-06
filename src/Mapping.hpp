#ifndef MAPPING_HPP
#define MAPPING_HPP

#include "Tags.hpp"
#include "includes.hpp"
#include "Parser.hpp"
#include "uintmath.hpp"

// Mapping[s][d] is the parameter dimension for state s and data dimension d
// TODO this must only be a many-to-one mapping; there should be one from data dimensions to parameter dimensions, one from parameter to hyperparameter, and also a method to create a compound mapping to update hyperparameters from data directly. States are just indices for a vector of parameter mappings. The mapping to hyperparameters makes it easy to have one prior per parameter or a joint prior such as a Dirichlet Process Prior.

// for Parser
template<>
MappingType convertType(
    const string& s ) {

	if ( s == "shared" || s == "S" ) {
		return shared;
	} else if ( s == "combinations" || s == "C" ) {
		return combinations;
	} else if ( s == "manual" || s == "M" ) {
		return manual;
	} else {
		throw runtime_error( "Unknown mapping type " + s + "!" );
	}
}



// maps a state s to the indices in theta and tau_theta
class Mapping {
		vector<vector<size_t>> mValue;

		size_t mNrDataDim;
		size_t mNrParams;
		size_t mNrStates;

	public:

		Mapping() {};

		Mapping(
		    const size_t nrparams,	// number of parameters
		    const size_t nrdatadim,
		    const vector<size_t>& pointers )
			:
			mNrDataDim( nrdatadim ),
			mNrParams( nrparams ),
			mNrStates( pointers.size() / nrdatadim ) {

			if ( mNrDataDim <= 0 ) {
				throw runtime_error( "Number of data dimensions must be positive!" );
			}

			if ( mNrParams <= 0 ) {
				throw runtime_error( "Number of parameters must be positive!" );
			}

			if ( mNrStates <= 0 ) {
				throw runtime_error( "Number of states must be positive!" );
			}

			if ( !divides( pointers.size(), mNrDataDim ) ) {
				throw runtime_error( "Number of mapping pointers must be a multiple of the number of data dimensions!" );
			}

			for ( const auto & p : pointers ) {
				if ( p >= mNrParams ) {
					throw runtime_error( "Invalid mapping pointer " + to_string( p ) + ", there are only " + to_string( mNrParams ) + " emission distributions!" );
				}
			}

			mValue.reserve( mNrStates );
			for ( size_t offset = 0; offset < pointers.size(); offset += mNrDataDim ) {
				mValue.push_back( vector<size_t>( pointers.begin() + offset, pointers.begin() + offset + mNrDataDim ) );
			}
		}


		Mapping(
		    const size_t nrparams,	// number of parameters
		    const size_t nrdatadim,
		    const MappingType mappingType )
			:
			mNrDataDim( nrdatadim ),
			mNrParams( nrparams ) {


			if ( mNrDataDim <= 0 ) {
				throw runtime_error( "Number of data dimensions must be positive!" );
			}

			if ( mNrParams <= 0 ) {
				throw runtime_error( "Number of parameters must be positive!" );
			}




			// setup the mapping using vectors of reference wrappers
			switch ( mappingType ) {

				case combinations:	// all combinations of shared parameters
					mNrStates = pow( mNrParams, mNrDataDim );

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

				case shared:
					mNrStates = mNrParams;
					for ( size_t s = 0; s < mNrStates; ++s ) {
						mValue.push_back( vector<size_t>( mNrDataDim, s ) );
					}
					break;

				default:
					throw runtime_error( "Unknown mapping type!" );
					break;
			}

			if ( mNrStates <= 0 ) {
				throw runtime_error( "Number of states must be positive!" );
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

		size_t nrDataDim() const {
			return mNrDataDim;
		}

		void print() const {
			for ( size_t s = 0; s < mNrStates; ++s ) {
				cout << "state " << s << ":";
				for ( const auto & p : mValue[s] ) {
					cout << "\t" << p;
				}
				cout << endl;
			}
			cout << flush;
		}
};

#endif


