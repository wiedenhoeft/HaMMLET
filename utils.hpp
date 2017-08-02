////////// This file contains all those small utility functions we need. //////////

#ifndef UTILS_HPP
#define UTILS_HPP

#include "includes.hpp"



enum Direction {forward, backward, unset};


// Return a reference to the lower median of a vector.
template<typename T>
T& lowerMedian( vector<T>& vec ) {
	nth_element( vec.begin(), vec.begin() + ( vec.size() / 2 ), vec.end() );
	return vec[vec.size() / 2];
}


// Delete a vector and release its memory.
template <typename T>
void deleteVector( vector<T>& vec ) {
	vector<T>().swap( vec );
}





// concatenate elements of a vector and output to stream
template <typename T>
string concat(
    const vector<T>& vec, 	// vector of elements
    const string sep = "\t", 	// inner separator
    const string finalSep = "",	// outer separator, e.g. "\n"
    const size_t minSize = 0,	// minimum size, if vector is smaller, defaultValue will be appended
    const T defaultValue = T()
) {
	stringstream ss;
	size_t size = max( vec.size(), minSize );
	for ( const auto & v : vec ) {
		ss << v;
		size--;
		if ( size > 0 ) {
			ss << sep;
		}
	}
	while ( size > 0 ) {
		ss << defaultValue << sep;
		size--;
	}
	ss << finalSep;
	return ss.str();
}


// check whether a file exists
bool fileExists( const string& path ) {
	ifstream f( path.c_str() );
	bool status = f.good();
	f.close();
	return status;
}

// count the number of lines in a file based on the occurrence of newline
size_t nrLinesInFile( ifstream& infile ) {
	infile.unsetf( std::ios_base::skipws );
	size_t line_count = count( istream_iterator<char> ( infile ), istream_iterator<char>(), '\n' );
	infile.clear();
	infile.seekg( 0, ios::beg );
	infile.setf( std::ios_base::skipws );
	return line_count;
}




// append the values obtained from an input strea to a vector. If <ignoreInvalid> is true, this will throw an exception if some input cannot be converted to the proper type.
template<class T>
void istreamToVector(
    istream& stream,
    vector<T>& vec,
    const bool ignoreInvalid = false ) {

	T val;
	if ( ignoreInvalid ) {
		for ( ;; ) {
			stream >> val;
			if ( stream.eof() || stream.bad() ) {
				break;
			} else if ( stream.fail() ) {
				stream.clear(); // unset failbit
				stream.ignore( 1 ); // skip next char
			} else {
				vec.push_back( val );
			}
		}
	} else {
		while ( stream >> val ) {
			vec.push_back( val );
		}
		if ( !( stream.eof() || stream.bad() ) ) {
			throw runtime_error( "Invalid input encountered!" );
		}
	}
}


#endif
