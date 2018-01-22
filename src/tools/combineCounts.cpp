// Combine counts from different files, by adding, subtracting.

#include "MappedValues.hpp"

#include <iostream>
using std::cout;
using std::endl;
using std::flush;

#include <algorithm>
using std::sort;

#include "../Parser.hpp"

#include <string>
using std::string;

#include <fstream>
using std::ifstream;
using std::ofstream;
using std::getline;

#include <sstream>
using std::stringstream;

#include <unordered_map>
using std::unordered_map;

#include "gzstream.h"

int main( int argc, const char* argv[] ) {

	Parser args( argc, argv );

	// OPTIONS
	// file containing
	args.registerFlags( {"-i", "-input-prefices"}, "" );	// prefix list for files that should be counted positive
	args.registerFlags( {"-p", "-pos-suffix"}, "-pos.csv.gz" );
	args.registerFlags( {"-c", "-count-suffix"}, "-count.csv.gz" );
	args.registerFlags( {"-s", "-size-suffix"}, "-size.csv" );

	args.registerFlags( {"-n", "-normalization-prefix"}, "mappability" );
	
	args.registerFlags( {"-o", "-out-prefix"} );


	args.registerFlags( {"-h", "--help", "-help"}, "" );
	args.parseArgs();

	const string posSuffix = args.parse<string>( "-pos-suffix", 0 );
	const string countSuffix = args.parse<string>( "-count-suffix", 0 );
	const string sizeSuffix = args.parse<string>( "-size-suffix", 0 );
	const string outPrefix = args.parse<string>( "-out-prefix" );
	vector<string> prefices = args.parseVector<string>( "-input-prefices" );

	if ( args.isSet( "-h" ) ) {
		cout << "Takes lists of file prefices (-i) and adds their counts (use + and - before lists of files), and adds them. Shared suffices for files can be set using -p, -c, and -s, for position, count and size. The output prefix is set using -o. If -n is provided, its argument is used as a prefix for normalization, i.e. counts are multiplied if a position exists (e.g. for mappability correction)." << endl;
		return 0;
	}


	unordered_map<string, vector<MappedValueEntry<long int>> > refseqToData;
	vector<string> observedRefSeqs;	// refseqs should be stored in the order observed in files, not by string comparison in an ordered map. Also, unordered_map is faster.
	vector<MappedValueEntry<long int>> currentData;


	// TODO implement with GenomeGetter
	
	long int sign = 1;
	if ( prefices[0] != "+" && prefices[0] != "-" ) {
		throw runtime_error( "First token of -i must be + or -!" );
	}
	for ( const auto & prefix : prefices )	{
		if ( prefix == "+" ) {
			sign = 1;
			continue;
		}
		if ( prefix == "-" ) {
			sign = -1;
			continue;
		}

		if ( sign > 0 ) {
			cout << "Adding";
		} else {
			cout << "Subtracting";
		}
		cout  << " counts for " << prefix << "*" << endl << flush;

		const string sizefilename = prefix + sizeSuffix;
		const string posfilename = prefix + posSuffix;
		const string countfilename = prefix + countSuffix;

		ifstream sizeFile( sizefilename );
		igzstream posFile( posfilename.c_str() );
		igzstream countFile( countfilename.c_str() );



		// string to hold line read from file
		string line;

		// variable to stream unised parts of string to
		string dump;

		// variable to hold current refseq
		string refseq;


		// variable for number of entries in a refseq
		size_t nrEntries;

		// variable to hold current genome position
		size_t pos;

		// variable to hold current count
		long int count;

		if ( !sizeFile.good() ) {
			throw runtime_error( "Cannot open " + sizefilename + "!" );
		}
		while ( getline( sizeFile, line ) ) {

			// string stream for line
			stringstream ss( line );

			// parse line in sizeFile
			ss >> refseq >> nrEntries >> dump;

// 			cout << "\t" << refseq << endl;

			// check whether refseq has been observed before
			if ( refseqToData.find( refseq ) == refseqToData.end() ) {

				// insert new refseq into records
				observedRefSeqs.push_back( refseq );

				// add empty vector to new refseq
				refseqToData.insert( {refseq, {}} );
			}


			// swap out data for refseq to currentData to be processed
			refseqToData[refseq].swap( currentData );

			currentData.reserve( currentData.size() + nrEntries );

			for ( size_t i = 0; i < nrEntries; ++i ) {

				// read position
				getline( posFile, line );
				pos = atoi( line.c_str() );

				// read count
				// TODO negative numbers
				getline( countFile, line );
				count = sign * atoi( line.c_str() );

				currentData.push_back( MappedValueEntry<long int>( pos, count ) );
			}

			// compress currentData
			sortAddAndCompress( currentData );

			// put the current data into the refseq that ends here, currentData is now an empty vector
			refseqToData[refseq].swap( currentData );

		}
		sizeFile.close();
		posFile.close();
		countFile.close();

	}


	// normalize counts if necessary
	if (args.isSet("-n")){
		throw runtime_error("-n not implemented yet!");
	}


	// output
	cout << "Writing output to " << outPrefix << "*" << endl << flush;
	ofstream sizeFile( outPrefix + sizeSuffix );
	ogzstream posFile( (outPrefix + posSuffix).c_str() );
	ogzstream countFile( (outPrefix + countSuffix).c_str() );
	size_t totalSize = 0;
	size_t size = 0;
	for ( const string & refseq : observedRefSeqs ) {
		size = refseqToData[refseq].size();
		totalSize += size;
		sizeFile << refseq << "\t" << size << "\t" << totalSize << endl;
		for ( const auto & entry : refseqToData[refseq] ) {
			posFile << entry.pos << endl;
			countFile << entry.entry << endl;
		}
	}
	sizeFile.close();
	posFile.close();
	countFile.close();
}




