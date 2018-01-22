#ifndef WAVELET_HPP
#define WAVELET_HPP

#include "includes.hpp"
#include "uintmath.hpp"



// TODO for implementation of chunks etc.: be carefull, we don't maintain scale coefficients!
// In-place Haar transform, for arbitrary data sizes. Data is sorted by (position, dimension), and dimensions are processed in pseudo-parallel fashion. TODO acually make this parallel?
// Data is treated like a greedy concatenation of vectors of sizes that are powers of two. In place where there would be the scale coefficient, the result is infinity.
void HaarDetailCoeffs(
    vector<real_t>& y,
    size_t dim = 1 ) {


	const size_t Tdim = y.size();	// T*dim
	if ( !divides( Tdim, dim ) ) {
		throw runtime_error( "Cannot compute Haar detail coefficients, array size is not a multiple of the number of dimensions!" );
	}
	const size_t T = Tdim / dim;	// number of input positions

	real_t yL;
	real_t yR;
	size_t L;
	size_t R;
	size_t Ldim;
	size_t Rdim;
	const size_t Nup = ceilPow2( T );
	for ( size_t N = 2; N <= Nup; N *= 2 ) {	// N is size of non-zero support interval
		const real_t s = 1 / sqrt( N );
		L = 0;
		R = N / 2;
		while ( L < T ) {
			if ( R < T ) {
				Ldim = L * dim;
				Rdim = R * dim;
				for ( size_t d = 0; d < dim; ++d ) {
					yL = y[Ldim];
					yR = y[Rdim];
					y[Ldim] = yL + yR;
					y[Rdim] = s * ( yL - yR );
					Ldim ++;
					Rdim++;
				}
			} else {
				Ldim = L * dim;
				for ( size_t d = 0; d < dim; ++d ) {
					y[Ldim] = inf;
					Ldim++;
				}
			}
			L += N;
			R += N;
		}
	}

	// set scale coefficient at first position force breakpoint at first element
	for ( size_t d = 0; d < dim; ++d ) {
		y[d] = inf;
	}
}




// Takes a maxlet transform, and computes the breakpoint weights, i.e. for each position t it computes the maximum absolute coefficient of all wavelets which have a discontinuity at t. Complexity is in-place in linear time.
void HaarBreakpointWeights(
    vector< real_t >& weights	// absolute Haar wavelet coefficients
) {
	const size_t size = weights.size();
	if ( size <= 0 ) {
		throw runtime_error( "Cannot compute Haar breakpoint weights, vector is empty!" );
	}
	size_t index;
	size_t L;
	size_t R;
	for ( size_t interval = ceilPow2( size ) / 2; interval >= 1; interval = interval / 2 ) {	// size of the positive (or negative) support interval (half the support)
		const size_t shift = 2 * interval;	// how far to advance to the next index
		for ( size_t index = interval; index < size; index += shift ) {
			L = index - interval;
			R = index + interval;
			if ( R < size ) {
				weights[R] = max( weights[R], weights[index] );
			} else {
				// NOTE just a precaution we expect this to be the case in the input
				weights[L] = inf;
				weights[index] = inf;
			}
			weights[L] = max( weights[L], weights[index] );
		}
	}
}


// Computes the maxlet transform (absolute Haar wavelet transform  for each dimension, then maximum of corresponding values across dimensions) from streaming input (dimensions first, then position), using only space T for coefficients and nrDim*T for statistics, plus nrDim*log2(T) for a stack. Output: coeffs.size()=T, suffstats.size() = nrDim*T
template< typename T>
void MaxletTransform(
    istream& input,
    vector<real_t>& coeffs,
    vector< SufficientStatistics<T> >& suffstats,
    const size_t nrDim = 1,
    const size_t reserveT = 0	// an estimate of the number of data points to avoid reallocation
) {

	if ( nrDim <= 0 ) {
		throw runtime_error( "Number of dimensions must be positive!" );
	}


	if ( coeffs.size() > 0 ) {
		throw runtime_error( "Coefficient array must be empty!" );
	}

	if ( suffstats.size() > 0 ) {
		throw runtime_error( "Statistics array must be empty!" );
	}

	if ( input ) {


		coeffs.reserve( ( reserveT + nrDim ) / nrDim + nrDim );
		suffstats.reserve( reserveT + nrDim );

// 	stack<real_t, vector<real_t> > S;	// stack never gets larger than nrDim*log2(T), so we don't expect a lot of reallocation, and save a lot of push and pop operations due to random access
		vector<real_t> S;
		size_t i = 0;
		real_t v = 0;
		size_t dim = 0;

		while ( input >> v ) {
			S.push_back( v );
			suffstats.push_back( SufficientStatistics<T>( v ) );
			dim++;	// set dimension of next value
			if ( dim == nrDim ) {	// filled all dimensions at index i
				dim = 0;	// next value will be first dimension again


				coeffs.push_back( inf );


				size_t j = i;	// points to node indices on an upward-left path (i.e. DFS post-order)
				size_t m = 1;	// mask to determine whether j is an index of a left child
				real_t normalizer = sqrt2half;

				while ( ( j & m ) > 0 ) {	// while j is on a left-upward path (DFS post-order)

					real_t maxCoeff = 0;	// the maximum detail coefficient across dimensions at j; NOTE we cannot take the maximum with coeffs because it contains infinity

					size_t L = S.size() - 2 * nrDim;	// index of left element in stack, get incremented to iterate over dimensions
					size_t R = L + nrDim;	// likewise, index of right element in stack


					// compute maximum of detail coefficients across dimensions
					for ( size_t d = 0; d < nrDim; ++d ) {
						maxCoeff = max( maxCoeff, normalizer * abs( S[L] - S[R] ) );
						S[L] += S[R];	// add right values to left values, so only the right values need to be popped
						L++;		// go to next dimension
						R++;
					}
					coeffs[j] = maxCoeff;


					// pop the right values
					for ( size_t d = 0; d < nrDim; ++d ) {
						S.pop_back();
					}


					j = j - m;	// move to left parent (if current position is not a right child, the loop will exit)
					m *= 2;	// move bit-mask to the left, i.e. check if i is still on a left-up path)
					normalizer *= sqrt2half;	// moving up one level changes normalization factor
				}
				i++;
			}
		}


		if ( dim != 0 ) {
			throw runtime_error( "Input stream did not contain enough values to fill all dimensions at last position!" );
		}

		coeffs[0] = inf;

	} else {
		throw runtime_error( "Cannot read input file or stream!" );
	}
}



#endif


