// Class that reads compressed genome representations and serves as a kind of iterator.

#include <string>
using std:: string;

#include <vector>
using std::vector;

#include <fstream>
using std::ifstream;
using std::getline;

#include <sstream>
using std::stringstream;

#include <stdexcept>
using std::runtime_error;

#include "gzstream.h"

class GenomeGetter {
		string mRefSeq;
		string mPrevRefSeq;
		size_t mPos;
		size_t mPrevPos;

		bool mIsNewRefSeq;

		// size of the current refseq, and the index of the line within that refseq
		size_t mRefSeqSize;
		size_t mRefSeqIndex;

		// same as above, but cumulative
		size_t mTotalSize;
		size_t mTotalIndex;

		igzstream posfile;
		igzstream sizefile;

		string line;
	public:

		// if only a prefix is provided, the input is assumed to consist of 3 files: size, pos, and count
		GenomeGetter( const string& prefix ):
			mRefSeq( "" ),
			mPrevRefSeq( "" ),
			mPos( 0 ),
			mPrevPos( 0 ),
			mRefSeqSize( 0 ),
			mRefSeqIndex( 0 ),
			mTotalSize( 0 ),
			mTotalIndex( 0 ),
			mIsNewRefSeq( false ) {
			sizefile.open( ( prefix + "-size.csv" ).c_str() );
			if ( !sizefile ) {
				throw runtime_error( "Cannot read " + prefix + "-size.csv!" );
			}
			posfile.open( ( prefix + "-pos.csv.gz" ).c_str() );
			if ( !posfile ) {
				throw runtime_error( "Cannot read " + prefix + "-pos.csv.gz!" );
			}
		}


		bool next() {
			stringstream ss( line );
			if ( mRefSeqIndex == mRefSeqSize ) { // start a new refseq
				mIsNewRefSeq = true;
				mPrevRefSeq = mRefSeq;
				if ( getline( sizefile,  line ) ) {
					ss.str( line );
					ss >> mRefSeq;
					ss >> mRefSeqSize;
					ss >> mTotalSize;	// TODO check that this increased, adds up etc.
					mRefSeqIndex = 0;
					ss.clear();

				} else {	// no more lines
					// TODO assert there is nothing left in the other files
					mRefSeq = "";
					mPos = 0;
					return false;
				}
			} else {
				mIsNewRefSeq = false;
			}
			if ( getline( posfile, line ) ) {
				ss.str( line );
				mPrevPos = mPos;
				ss >> mPos;
				ss.clear();
			} else {
				throw runtime_error( "Not enough entries in position file!" );
			}

			mRefSeqIndex++;
			mTotalIndex++;
			return true;
		}

		const string& refseq() const {
			return mRefSeq;
		}

		const string& prevRefseq() const {
			return mPrevRefSeq;
		}

		const size_t& pos() const {
			return mPos;
		}

		const size_t& prevPos() const {
			return mPrevPos;
		}

		bool refseqChanged() const {
			return mIsNewRefSeq;
		}


};
