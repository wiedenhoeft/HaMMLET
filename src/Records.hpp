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

		size_t mNrObservedPos;
		size_t mNrBlocks;
		size_t mNrSegments;
		size_t mSegmentState;	// the state being aggregated
		size_t mSegmentSize;	// the current aggregated segment size
		const size_t mSize;
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
			mNrBlocks( 0 ),
			mSize( T ),
			mNrObservedPos( 0 ),
			mPrefix( prefix ),
			mSuffix( suffix ),
			mMarginals( T ),
			mRecordMarginals( true ),
			mRecordSegments( false ),
			mRecordBlocks( false ),
			mRecordCompression( false ),
			mRecordSequences( false ),
			mRecordTheta( false ),
			mSegmentSize( 0 ),
			mSegmentState( 0 ),
			mNrSegments( 0 ) {}

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

		template<typename ThetaType>
		void record(
		    const Theta<ThetaType>& theta ) {
			// write parameters
			if ( mRecordTheta ) {
				mThetaFile << theta << endl;
			}
		}

		void record(
		    const size_t state,
		    const size_t N ) {

			string blockPreChar = "\t";
			string segmentPreChar = "";
			if (mNrSegments>0){
				segmentPreChar = "\t";
			}
			string postChar = "";

			// is this the first block?
			if ( mNrBlocks == 0 ) {
				mSegmentState = state;
				mSegmentSize = N;
				blockPreChar="";
			} else {	// this is not the first block
				// record the previous segment if we encountered a new state
				if ( state != mSegmentState ) {
					if ( mRecordMarginals ) {
						mMarginals.addRecord( mSegmentState, mSegmentSize );
					}
					if ( mRecordSequences ) {
						mSequenceFile << segmentPreChar << mSegmentSize << ":" << mSegmentState << postChar;
					}
					mSegmentState = state;
					mSegmentSize = N;
					mNrSegments++;
					segmentPreChar = "\t";
				} else {
					// do nothing except extend segment
					mSegmentSize += N;
				}
			}


			mNrBlocks++;

			mNrObservedPos += N;
			

			// end of line? In that case we need to record everything that's left
			if ( mNrObservedPos >= mSize ) {
				if ( mNrObservedPos > mSize ) {
					throw runtime_error( "Cannot record block, exceeding data size!" );
				}

				postChar = "\n";

				if ( mRecordCompression ) {
					mCompressionsFile << ( ( double ) mSize ) / ( ( double ) mNrBlocks ) << endl;
				}

				if ( mRecordSegments ) {
					mSegmentFile << mMarginals.nrSegments() << "\t" << mMarginals.internalSize() << endl;
				}



				// record the last block
				if ( mRecordMarginals ) {
					mMarginals.addRecord( mSegmentState, mSegmentSize );
				}
				if ( mRecordSequences ) {
					mSequenceFile << segmentPreChar << mSegmentSize << ":" << mSegmentState << postChar;
				}


				mNrObservedPos = 0;
				mNrBlocks = 0;
				mSegmentSize = 0;
				mNrSegments = 0;
			}


			if ( mRecordBlocks ) {
				mBlocksFile << blockPreChar << N << postChar;
			}


		}



		vector<size_t> maxMarginSegmentation() {
			return mMarginals.maxMarginSegmentation();
		}

};

#endif








