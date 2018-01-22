// Given a file of state marginals, compute the maximum posterior margins.

#include "../Parser.hpp"


#include <iostream>
using std::cin;
using std::cout;
using std::endl;
using std::istream;
using std::ostream;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <string>
using std::stoi;
using std::getline;

#include <stdexcept>
using std::runtime_error;

#include <sstream>
using std::stringstream;

int main( int argc, const char* argv[] ) {

	Parser args( argc, argv );
	args.registerFlags( {"-i", "-infile"}, "" );
	args.registerFlags( {"-h", "--help", "-help"}, "" );
	args.parseArgs();

	if ( args.isSet( "-h" ) ) {
		cout << "Given a marginals file (-i) or input from STDIN, computes the maximum posterior margins segmentation, combining adjacent segments whenever possible." << endl;
		return 0;
	}

	const bool readFromFile = args.isSet( "-i" );
	ifstream realInFile;
	if ( readFromFile ) {
		realInFile.open( args.parse<string>( "-i" ), std::ios::in );
	}
	istream& inFile = ( readFromFile ? realInFile : cin );

	string line;
	size_t count;
	size_t RLE = 0;
	size_t totalRLE = 0;
	size_t col = 0;
	size_t maxCol = 0;
	size_t index = 0;
	size_t maxIndex = 0;
	size_t prevIndex = 0;

	while ( getline( inFile, line ) ) {
		stringstream ss( line );

		index = 0;
		maxIndex = 0;
		col=0;
		maxCol = 0;
		ss >> RLE;
		while ( ss >> col ) {
			if ( col > maxCol ) {
				maxIndex = index;
				maxCol = col;
			}
			index++;
		}

		if ( maxIndex == prevIndex ) {
			totalRLE += RLE;
		} else {
			cout << totalRLE << "\t" << prevIndex << endl;
			totalRLE = RLE;
			prevIndex = maxIndex;
		}

	}
	cout << totalRLE << "\t" << maxIndex << endl;

}
