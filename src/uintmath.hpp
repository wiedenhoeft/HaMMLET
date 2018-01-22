////////// Mathematical operations and bit manipulations for unsigned integer types //////////


#ifndef UINTMATH_HPP
#define UINTMATH_HPP

#include "includes.hpp"

// keep lowest bit:
// set all except the rightmost 1-bit (LSB) to 0.
// this is 2^ctz(x, nrDigits)
inline size_t klb( size_t x ) {
	if ( x > 0 ) {
		return ( x & ( ( ~x ) + 1 ) );
	} else {
		throw runtime_error( "klb(0) is undefined!" );
	}
}


// Calculate the number of trailing zeros.
// Given an integral type, this is log2(klb(x)).
inline size_t ctz( size_t x ) {
	if ( x > 0 ) {
		size_t c = 0;
		x = ( x ^ ( x - 1 ) ) >> 1;
		for ( c = 0; x; c++ ) {
			x >>= 1;
		}
		return c;
	} else {
		throw runtime_error( "ctz(0) is undefined" );
	}
}


// Set the lowest 1-bit to 0.
inline size_t ulb( size_t x ) {
	return x & ( x - 1 );
}


// set trailing zero bits to 1
inline size_t stb( size_t x ) {
	return x | ( x - 1 );
}





// check if a is a multiple of b
template <typename T>
inline bool divides( T a, T b ) {
	return a == b * ( a / b );
}

// divide a/b if it is divisible, throw error otherwise
template <typename T>
inline T divide( T a, T b ) {
	T d = a / b;
	if ( b * d == a ) {
		return d;
	} else {
		throw runtime_error( "Truncated integer division: " + to_string( b ) + " does not divide " + to_string( a ) + "!" );
	}
}


// for unsigned integer types, round division a/b to the closest integer (.5 round up) without casting to float. round(a/b) = floor((a+ floor(b/2))/b)
template<typename T>
inline T rounddiv( T a, T b ) {
	return ( a + ( b / 2 ) ) / b;
}





// integer 2**x for non-negative x
inline size_t iexp2( size_t x ) {
	if ( x > numeric_limits<size_t>::digits ) {
		throw runtime_error( "Exponent too large for iexp2()!" );
	}
	return ( 1 << x );
}





// round up to the next multiple of m
// Returns the lowest multiple of m greater-equal n
inline size_t ceil_mult( size_t n, size_t m ) {
	return ( ( ( n + m ) - 1 ) / m ) * m;
}

// Return the lowest multiple of m strictly greater than n.
inline size_t higher_mult( size_t n, size_t m ) {
	return ( ( n + m ) / m ) * m;
}


// Returns the highest multiple of m less-equal than n
inline size_t floor_mult( size_t n, size_t m ) {
	return ( n / m ) * m;
}


// Returns the highest multiple of m strictly smaller than n
// NOTE not defined if n==0
inline size_t lower_mult( size_t n, size_t m ) {
	if ( n > 0 ) {
		return ( ( n - 1 ) / m ) * m;
	} else {
		throw runtime_error( "lower_mult(0) is undefined!" );
	}
}


// test wether a number is a power of 2
inline bool isPow2( size_t x ) {
	return ( x > 0 ) && ( ulb( x ) == 0 );
}


// the next higher power of two
inline size_t ceilPow2( size_t n ) {
	// TODO fails if first bit and some other bit are set (overflow)
	if ( isPow2( n ) ) {
		return n;
	}
	size_t p = 1;
	while ( p < n ) {
		p <<= 1;
	}
	return p;
}



// round to the smaller power of two
// essentially: keep highest bit
inline size_t floorPow2( size_t n ) {
	if ( isPow2( n ) ) {
		return n;
	}
	size_t p = 1;
	while ( p <= n ) {
		p <<= 1;
	}
	return p / 2;
}


// Check if a number is even.
inline bool isEven( const size_t& x ) {
	return ( 1 & x ) == 0;
}



#endif
