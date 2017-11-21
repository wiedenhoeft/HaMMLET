// Records data to file for each iteration.

#ifndef RECORDS_HPP
#define RECORDS_HPP

#include "StateMarginals.hpp"
#include "Theta.hpp"
#include "StateMarginalsIterator.hpp"


#include <string>
using std::string;

#include <fstream>
using std::ofstream;




class Records {

		string mPrefix;
		string mSuffix;

		StateMarginals mMarginals;

		bool mRecordMarginals;
		bool mRecordBlocks;
		bool mRecordCompression;
		bool mRecordSequences;
		bool mRecordTheta;
		bool mRecordSegments;

		// TODO allow individual names for each file?

		ofstream mMarginalsFile;
		ofstream mSequenceFile;
		ofstream mBlocksFile;
		ofstream mThetaFile;
		ofstream mCompressionsFile;
		ofstream mSegmentFile;

		// helper method to avoid copying of code
		void setRecordX(
		    ofstream& mFile, 	// reference to member ofstream
		    const string& type, 	// type of record, also used as suffix
		    bool& mFlag,  	// reference to bool member (whether to record or not)
		    const bool flag,	// whether to record or not
		    const bool overwrite	// whether to allow overwriting existing files
		) {
			mFlag = flag;
			if ( mFlag && !mFile.is_open() ) {
				string filename = mPrefix +  type + mSuffix;
				if ( fileExists( filename ) && !overwrite ) {
					throw runtime_error( "File " + filename + " already exists! Use -w to allow overwrite!" );
				}
				mFile.open( ( filename ).c_str() );
				if ( !mFile.is_open() ) {
					throw runtime_error( "Cannot write to file " + filename + "!" );
				}
			}
		}

	public:
		// delete copy constructor
		Records( const Records& that ) = delete;


		Records(
		    size_t T,
		    string prefix,
		    string suffix,
		    const size_t nrStates
		) :
			mPrefix( prefix ),
			mSuffix( suffix ),
			mMarginals( T ),
			mRecordMarginals( true ),
			mRecordSegments( false ),
			mRecordBlocks( false ),
			mRecordCompression( false ),
			mRecordSequences( false ),
			mRecordTheta( false ) {}

		~Records() {
			close();
		}

		void close() {
			if ( mRecordMarginals ) {
				mMarginals.save( mMarginalsFile );
				mMarginalsFile.close();
			}
			if ( mRecordSequences ) {
				mSequenceFile.close();
			}
			if ( mRecordBlocks ) {
				mBlocksFile.close();
			}
			if ( mRecordTheta ) {
				mThetaFile.close();
			}
			if ( mRecordSegments ) {
				mSegmentFile.close();
			}
		}

		void setRecordMarginals( bool b, bool overwrite = false ) {
			setRecordX( mMarginalsFile, "marginals", mRecordMarginals, b, overwrite );
		}

		void setRecordBlocks( bool b, bool overwrite = false ) {
			setRecordX( mBlocksFile, "blocks", mRecordBlocks, b, overwrite );
		}

		void setRecordCompression( bool b, bool overwrite = false ) {
			setRecordX( mCompressionsFile, "compression", mRecordCompression, b, overwrite );
		}

		void setRecordStateSequence( bool b, bool overwrite = false ) {
			setRecordX( mSequenceFile, "sequences", mRecordSequences, b, overwrite );
		}

		void setRecordTheta( bool b, bool overwrite = false ) {
			setRecordX( mThetaFile, "parameters", mRecordTheta, b, overwrite );	// TODO rename?
		}

		void setRecordSegments( bool b, bool overwrite = false ) {
			setRecordX( mSegmentFile, "segments", mRecordSegments, b, overwrite );	// TODO rename?
		}


		template<typename StateSequenceType, typename EmissionsType,  typename ThetaType >
		void record(
		    StateSequenceType& q,
		    EmissionsType& y,
		    ThetaType& theta ) {

			mMarginals.record( q, y.blocks() );
			y.initForward();
			size_t w = 0;
			bool firstSegmentDone = false;	// whether the curent segment to be output is the first one, to decide whether to put a tab
			size_t segSize = 0;
			int prevState;
			int state;
			if ( mRecordSegments ) {
				mSegmentFile << mMarginals.nrSegments() << "\t" << mMarginals.internalSize() << endl;
			}
			while ( y.next() ) {
				if ( mRecordSequences ) {
					if ( w == 0 ) {
						segSize = y.blockSize();
						prevState = q[0];
					} else {
						state = q[w];
						if ( prevState != state ) {
							if ( firstSegmentDone ) {
								mSequenceFile << " ";
							} else {
								firstSegmentDone = true;
							}
							mSequenceFile << segSize << ":" << prevState;
							segSize = y.blockSize();
							prevState = state;
						} else {
							segSize += y.blockSize();
						}
					}

				}
				if ( mRecordBlocks ) {
					mBlocksFile << y.blockSize() << " ";
				}
				w++;
			}
			if ( mRecordBlocks ) {
				mBlocksFile << endl;
			}
			if ( mRecordSequences ) {
				if ( firstSegmentDone ) {
					mSequenceFile << " ";
				}
				mSequenceFile << segSize << ":" << state;
				mSequenceFile << endl;
			}

			if ( mRecordCompression ) {
				mCompressionsFile << ( ( real_t ) y.size() ) / ( ( real_t ) y.nrBlocks() ) << endl;
			}
			// write parameters
			if ( mRecordTheta ) {
				mThetaFile << theta << endl;
			}
		}



		vector<size_t> maxMarginSegmentation() {
			return mMarginals.maxMarginSegmentation();
		}

};

#endif


