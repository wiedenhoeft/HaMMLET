// reads data stream and puts average of non-overlapping windows to stdout
#include <iostream>
using std::ostream;
using std::endl;
using std::cin;
using std::cout;
using std::flush;


#include <stdexcept>
using std::runtime_error;
using std::exception;

#include <sstream>
using std::istringstream;

int main( int argc, const char* argv[] ) {


	if ( argc <= 1 ) {
		throw runtime_error( "Not enough arguments!" );
	}

	istringstream ss( argv[1] );
	size_t windowSize;
	ss >> windowSize ;
	float v;
	float sum = 0;
	size_t pos = 0;
	while ( cin >> v ) {
		sum += v;
		pos++;
		if ( pos == windowSize ) {
			cout << sum / pos << endl;
			pos = 0;
			sum = 0;
		}
	}
	if ( pos != 0 ) {
		cout << sum / pos << endl;
	}
}
