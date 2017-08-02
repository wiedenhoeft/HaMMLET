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


// Takes a vector, computes its absolute values and merges all dimension by taking the maximum of their absolute values. This changes the vector size by a factor of 1/dim.
void mergeDimensions(
    vector<real_t>& y,
    const size_t dim = 1
) {
	const size_t T = y.size();
	if ( !divides( T, dim ) ) {
		throw runtime_error( "Cannot merge dimensions, array size is not a multiple of the number of dimensions!" );
	}
	if ( T > 0 ) {
		y[0] = abs( y[0] );
	}
	if ( dim > 1 ) {
		if ( T > 1 ) {
			size_t L = 0;
			size_t R = 1;
			size_t d = 1;
			while ( R < T ) {
				if ( d == dim ) {
					d = 0;
					L++;
					y[L] = abs( y[R] );
				} else {
					y[L] = max( y[L], abs( y[R] ) );
					d++;
				}
				R++;
			}
		}
		y.resize( T / dim );
	} else {
		for ( size_t t = 0; t < T; ++t ) {
			y[t] = abs( y[t] );
		}
	}
}



// takes a vector of absolute wavelet coefficients in DFS in-order layout (possibly the maxima over several dimensions), and combines all the breakpoint weights for all levels into the appropriate position.
void HaarBreakpointWeights(
    vector< real_t >& weights,
    const size_t nrDim = 1 ) {	// TODO for sizes larger than 1

	// IDEA we can use positive/negative values to mark the boundary to which the largest weight in a chunk is closer
	const size_t size = weights.size();
	if ( size <= 0 ) {
		throw runtime_error( "Cannot compute Haar breakpoint weights, vector is empty!" );
	}


	HaarDetailCoeffs( weights );
	mergeDimensions( weights, nrDim );	// merge the data dimensions

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


#endif
