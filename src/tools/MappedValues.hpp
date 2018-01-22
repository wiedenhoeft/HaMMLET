// Values mapped to a genome, such as read-depth, start counts etc.

#include <vector>
using std::vector;


#include <string>
using std::string;

#include <map>
using std::map;

#include <stdexcept>
using std::runtime_error;

#include <algorithm>
using std::sort;


template <typename T>
class MappedValueEntry {
	public:
		size_t pos;
		T entry;


		MappedValueEntry(): pos( 0 ), entry() {}

		MappedValueEntry( size_t p, T v ): pos( p ), entry( v ) {}


		bool operator<( const MappedValueEntry<T>& b ) const {
			return pos < b.pos;
		}


		bool operator>( const MappedValueEntry<T>& b ) const {
			return pos > b.pos;
		}

		MappedValueEntry<T>& operator+=( const MappedValueEntry<T>& rhs ) {
			if ( pos != rhs.pos ) {
				throw runtime_error( "Cannot add values, positions don't match!" );
			}
			entry += rhs.entry;
			return *this;
		}

		const MappedValueEntry<T> operator+( const MappedValueEntry<T>& other ) const {
			return MappedValueEntry<T>( *this ) += other;
		}

		MappedValueEntry<T>& operator-=( const MappedValueEntry<T>& rhs ) {
			if ( pos != rhs.pos ) {
				throw runtime_error( "Cannot add values, positions don't match!" );
			}
			entry -= rhs.entry;
			return *this;
		}

		const MappedValueEntry<T> operator-( const MappedValueEntry<T>& other ) const {
			return MappedValueEntry<T>( *this ) -= other;
		}
};


template<typename T>
void sortAddAndCompress(
    vector<MappedValueEntry<T>>& vec
) {
	sort( vec.begin(), vec.end() );
	size_t L = 0;
	size_t R = 1;
	while ( R < vec.size() ) {
		if ( vec[L].pos == vec[R].pos )  {
			vec[L].entry += vec[R].entry;
		} else {
			++L;
			vec[L] = vec[R];
		}
		++R;
	}
	vec.resize( L + 1 );
}

template<typename T>
void sortMultiplyAndCompress(
    vector<MappedValueEntry<T>>& vec
) {
	sort( vec.begin(), vec.end() );
	size_t L = 0;
	size_t R = 1;
	while ( R < vec.size() ) {
		if ( vec[L].pos == vec[R].pos )  {
			vec[L].entry *= vec[R].entry;
		} else {
			++L;
			vec[L] = vec[R];
		}
		++R;
	}
	vec.resize( L + 1 );
}



template <typename T>
class MappedValues {
		map<string, vector<MappedValueEntry<T>>> mEntries;	// map refseq ID to vector of (pos, T) tuples


		void update( const MappedValues<T>& other ) {};

		// TODO add(), subtract(bool removeZero=false)
};







