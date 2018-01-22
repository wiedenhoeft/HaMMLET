#ifndef INCLUDES_HPP
#define INCLUDES_HPP

#include <cstdint>
// using int_16_t;

// This file contains header inclusions, because we are lazy (and consistent), as well as some basic constants and typedefs.

// Type to use for real numbers. Note that we don't introduce separate types for data, wavelet coefficients and probabilities, as one might be tempted, because that would lead to many implicit type conversions in likelihood computations, wavelet thresholding etc.
typedef float real_t;

// Type to use for the marginal counts. Has to be a signed integer type. Its maximum is the maximum number of counts per state and position as well as the maximum number of recorded states.
typedef int16_t marginal_t;	// the type used to record marginal counts


#include <type_traits>
using std::is_integral;
using std::is_unsigned;

#include <cstddef>
using std::size_t;


#include <vector>
using std::vector;

#include <queue>
using std::queue;

#include <deque>
using std::deque;

#include <array>
using std::array;


#include <string>
using std::string;
using std::to_string;


#include <sstream>
using std::stringstream;
using std::istringstream;


#include <iostream>
using std::istream;
using std::ostream;
using std::endl;
using std::cin;
using std::cout;
using std::cerr;
using std::clog;
using std::wcout;
using std::flush;
using std::boolalpha;
using std::ios;


#include <fstream>
using std::ifstream;
using std::ofstream;


#include <stdexcept>
using std::runtime_error;	// TODO throw the appropriate errors, like logic_error etc.
using std::exception;

#include <cmath>
using std::pow;
using std::exp;	// e^x
using std::exp2;
using std::log;	// natural log
using std::log2;
using std::log10;
using std::sqrt;
using std::ceil;
using std::floor;
using std::abs;

using std::isfinite;



#include <algorithm>
using std::min;
using std::max;
using std::nth_element;
using std::reverse;
using std::fill;


#include <numeric>
using std::partial_sum;
using std::plus;
using std::accumulate; //e.g. sum of vector


#include <iterator>
using std::istream_iterator;
using std::ostream_iterator;
using std::back_inserter;


#include <unordered_map>
using std::unordered_map;


#include <stack>
using std::stack;


#include <climits>
using std::numeric_limits;


#include <ctime>
using std::time;


#include <iomanip>
using std::setprecision;


const real_t inf = numeric_limits<real_t>::infinity();
const real_t sqrt2 = sqrt( 2.0 );
const real_t sqrt2half = sqrt2 / 2.0;	// sqrt(2)/2 = 1/sqrt(2)


#endif

