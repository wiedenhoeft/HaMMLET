#include "GenomeGetter.hpp"

#include "../Parser.hpp"

#include <limits>
using std::numeric_limits;

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

#include "gzstream.h"


int main( int argc, const char* argv[] ) {

	Parser args( argc, argv );

	// OPTIONS
	args.registerFlags( {"-g", "-genome-prefix"}, "" );
	args.registerFlags( {"-c", "-coordinates"}, "" );
	args.registerFlags( {"-w", "-window-size"}, "1" );
	args.registerFlags( {"-r", "-range"} );	// print start and end column for each refseq, instead of individual lines for each position
	args.registerFlags( {"-i", "-infile"}, "" );	// read from cin, unless -i is defined
	args.registerFlags( {"-o", "-outfile"}, "" );	// write to cout, unless -o is defined
	args.registerFlags( {"-b", "-blocks"}, "" );	// TODO implement: there is no size column, treat each position as block size 1
	args.registerFlags( {"-h", "--help", "-help"}, "" );
	args.parseArgs();

	if ( args.isSet( "-h" ) ) {
		cout << "Prepend genomic coordinates to lines, separated by tabs. If -c/-coordinates is set, coordinates are refseq:start:inclusiveend, otherwise they are separated by tabs. The input file name is set using -i, otherwise lines are read from STDIN. The PREFIX for the genomes size and position files is set using -g/-genome-prefix, and genomic coordinates are read from PREFIX-size.csv and PREFIX-pos.csv. The window size -w/-window size specifies the number of genome coordinates corresponding to one line in the data, for example if mapping counts have been averaged over adjacent, non-overlapping windows. The last data line may map to less genome positions than the window size, and no warning is issued. If -b/-blocks is set, the first entry in each input line (up to the first tab) is considered the segment size for a run-length encoding, and thus specifies that the lines should be repeated this many times; in that case, the window size is multiplied by that number. If -r/-range is specified, only the first and last genome position (in columns) is printed for each input segment per refseq, instead of mapping to all genome positions. If an INT is provided as argument to -r, this specifies the maximum distance between adjacent positions within a range, otherwise a new range is started. If -o/-outfile is specified, output is written to that path, otherwise it is written to STDOUT." << endl;
		return 0;
	}


	string prefix = args.parse<string>( "-genome-prefix" );

	// specify input stream (cin or -i)
	const bool readFromFile = args.isSet( "-i" );
	ifstream realInFile;
	if ( readFromFile ) {
		realInFile.open( args.parse<string>( "-i" ), std::ios::in );
	}
	istream& inFile = ( readFromFile ? realInFile : cin );

	// specify output stream (cout or -o)
	const bool writeToFile = args.isSet( "-o" );
	ofstream realOutFile;
	if ( writeToFile ) {
		realOutFile.open( args.parse<string>( "-o" ), std::ios::out );
	}
	ostream& outFile = ( writeToFile ? realOutFile : cout );


	const bool isRLE = args.isSet( "-b" );
	const size_t windowSize = args.parse<size_t>( "-window-size" );
	const bool outputRanges = args.isSet( "-range" );

	GenomeGetter gg( prefix );

	string sep1 = "\t";
	string sep2 = "\t";
	if ( args.isSet( "-coordinates" ) ) {
		sep1 = ":";
		sep2 = "-";
	}


	size_t start = 0;	// start position of a segment if isRLE, data position otherwise
	size_t end = 0;

	string refseq;
	string line;
	size_t segmentSize = 1;
	const size_t windowsize = args.parse<size_t>( "-w" );
	size_t maxMergeDist = numeric_limits<size_t>::max();	// maximum distance of entries to still be considered for merging into range
	if ( outputRanges ) {
		if ( args.nrTokens( "-range" ) > 0 ) {
			maxMergeDist = args.parse<size_t>( "-range", 0 );
		}
	}

	// the number of lines in the pos file corresponding to a line in the data file. This is the window size times the segment size
	size_t nrGenomeLines = 0;
	while ( getline( inFile, line ) ) {	// TODO check status etc.

		// if the first input column represents the segment size, update it
		if ( isRLE ) {
			// get the size of the segment
			segmentSize = stoi( line.substr( 0, line.find_first_of( "\t" ) ) ) ;	// get segment size from string up to the first tab
			line = line.substr(line.find_first_of( "\t" )+1 ) ;
			if ( segmentSize == 0 ) {
				throw runtime_error( "Segment size must be positive!" );
			}
		}

		// specify the number of genome lines (positions) that this data line represents
		nrGenomeLines = windowsize * segmentSize;

		if ( outputRanges ) {

			// get values for new segment
			if ( gg.next() ) {
				if ( gg.refseqChanged() ) {
					refseq = gg.refseq();
				}
				start = gg.pos();
				end = start;

			} else {
				throw runtime_error( "Genome ended before all data was processed!" );
			}
			nrGenomeLines--;

			// read as many lines as the run-length information specifies
			while ( nrGenomeLines > 0 ) {
				if ( gg.next() ) {

					// if the current line starts a new refseq
					if ( gg.refseqChanged() || gg.pos() - end > maxMergeDist ) {	// TODO what if this happens in first line

						// print previous refseq, start and end
						cout << refseq << sep1 << start << sep2 << end << "\t" << line << endl;

						refseq = gg.refseq();
						start = gg.pos();
					}
					end = gg.pos();
				} else {
					break;	// silently ignore window size going past the end of the data
				}
				nrGenomeLines--;
			}

			cout  << refseq << sep1 << start << sep2 << end << "\t" << line << endl;

		} else { // output each line
			while ( nrGenomeLines > 0 ) {
				if ( gg.next() ) {
					cout << gg.refseq() << sep1 << gg.pos() << "\t" << line << endl;
				} else {
					break;
				}
				nrGenomeLines--;
			}
		}
		// we only allow for the last window to be incomplete, since the genome size might not be a multiple of the window size, but otherwise we enforce correct size
		if ( nrGenomeLines >= windowsize ) {
			throw runtime_error( "Data too long for genome!" );
		}
	}

	// see if there are unprocessed parts of the genome
	if ( gg.next() ) {
		throw runtime_error( "Data ended before genome!" );
	}

	return 0;
}




