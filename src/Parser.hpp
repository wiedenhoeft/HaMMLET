#ifndef PARSER_HPP
#define PARSER_HPP

// Implements a rudimentary but flexible parser. The command line is split at whitespaces (multiple whitespaces are merged) into a vector of string tokens. Parser::registerFlags() is used to register groups of tokens as representing the same command line option, and an optional string of default parameter tokens, e.g. Parser::registerFlags({"-f", "--files"}, "file1.txt file2.txt"). It does not enforce the -/-- convention; instead, arbitrary strings can be defined as flags. This also means that stacking like -abc will NOT be treated as -a -b -c. However, it is possible to define subparsers on the argument lists associated with each token, for instance, define Parser subargs, register flags that have not been registered with args, then use args.tokens(flag) to obtain a token list for a given flag and parse it into subargs.

// These basic functions are available:
// isSet(flag) to check flag has been set by the user
// parse<type>(flag, i) to parse the i-th parameter token into a variable of a given type
// parseVector<type>(flag, begin, end) to do the same, but yield a vector as a result
// nrTokens(flag) to get the number of tokens passed to a given flag


// NOTE Though default tokens can be defined, there is a method that checks whether a flag has been provided by the user. Depending on the design choices of a software, providing the flag on the command line could change the behavior of the software even if the default tokens are provided (presence and absence of a flag can toggle variables).

#include <vector>
using std::vector;

#include <string>
using std::string;
using std::to_string;

#include <iostream>
using std::cout;
using std::endl;
using std::flush;

#include <sstream>
using std::istringstream;


#include<stdexcept>
using std::runtime_error;


#include<map>
using std::map;

#include <initializer_list>
using std::initializer_list;

#include <set>
using std::set;


template<typename T>
T convertType(
    const string& s ) {
	T x;
	if ( istringstream( s ) >> x ) {
		return x;
	} else {
		throw runtime_error( "Conversion failed for string \"" + s + "\"!" );
	}
}



// tokenize a string
vector<string> inline tokenize(
    const string& source,
    const char* delimiter = " ",
    bool keepEmpty = false ) {
	vector<string> results;

	size_t prev = 0;
	size_t next = 0;

	while ( ( next = source.find_first_of( delimiter, prev ) ) != string::npos ) {
		if ( keepEmpty || ( next - prev != 0 ) ) {
			results.push_back( source.substr( prev, next - prev ) );
		}
		prev = next + 1;
	}

	if ( prev < source.size() ) {
		results.push_back( source.substr( prev ) );
	}

	return results;
}


class Parser {

		vector< string > args;

		vector<string> mTokenVector;	// the vector of input tokens


		map<string, size_t> mFlagToIndex;

		vector<vector<string>> mTokenGroups;
		vector<bool> mFlagsProvided;	// information about whether or not a flag has been parsed already
		vector<string> mHelpTexts;	// holds the help text for the flags
		set<string> mBlockedFlags;	// list of flags that are not allowed, useful when creating subparsers
		bool mTokensParsed;


		// assert that the command line has been parsed
		void assertParsed() const {
			if ( !mTokensParsed ) {
				throw runtime_error( "Command line has not been parsed yet!" );
			}
		}

		// assert that a flag is registered, otherwise throw exception
		void assertFlag( string flag ) {
			if ( !isFlag( flag ) ) {
				throw runtime_error( flag + " is not registered as a flag!" );
			}
		}

		// check whether a string is a registered flag 
		bool isFlag( string flag ) {
			return ( mFlagToIndex.count( flag ) > 0 );
		}


		// check whether a string is prohibited as a flag
		bool isBlocked( string flag ) {
			return mBlockedFlags.count( flag ) == 1;
		}





	public:

		Parser( int argc, const char* argv[] ) :
			mTokenVector( argv + 1, argv + argc ),
			mTokensParsed( false ) {
		}

		Parser( vector<string> tokens ) :
			mTokenVector( tokens ),
			mTokensParsed( false ) {
		}


		// register a group of strings (such as -f, --file) as a flag for the same group of arguments
		void registerFlags( initializer_list<string> flags, const string& defaults = "" ) {


			if ( mTokensParsed ) {
				throw runtime_error( "Cannot register flags, tokens have already been parsed!" );
			}

			for ( string flag : flags ) {
				if ( isFlag( flag ) ) {
					throw runtime_error( "Flag " + flag + " has already been registered!" );
				}
				if ( isBlocked( flag ) ) {
					throw runtime_error( "Flag " + flag + " is blocked!" );
				}
				mFlagToIndex[flag] = mTokenGroups.size();
			}
			mTokenGroups.push_back( tokenize( defaults ) );
			mFlagsProvided.push_back( false );
		}


		void parseArgs() {
			mTokensParsed = true;
			string currentFlag = "";
			vector<string> currentTokens;
			if ( mTokenVector.size() == 0 ) {
				return;
			}
			if ( !isFlag( mTokenVector[0] ) ) {
				throw runtime_error( "First input token (" + mTokenVector[0] + ") is not a registered flag; parser does not support positional arguments!" );
			}
			size_t currentIndex;
			for ( const string & token : mTokenVector ) {
				if ( isFlag( token ) ) {	// token is a registered flag

					// make sure the flag was not already provided
					if ( isSet( token ) ) {
						throw runtime_error( "Duplicate flag " + token + "!" );
					}

					// get index of flag
					currentIndex = mFlagToIndex[token];
					mFlagsProvided[currentIndex] = true;	// set token as being set by the user
					mTokenGroups[currentIndex].clear();	// remove default tokens

				} else {	// token is parameter to the previous flag
					mTokenGroups[currentIndex].push_back( token );
				}
			}
			mTokenVector.clear();
			mTokensParsed = true;
		}



		template< class T>
		T parse( const string& flag, size_t index=0 )  {
			assertParsed();
			assertFlag( flag );
			const vector<string>& tokens = mTokenGroups[mFlagToIndex[flag]];
			if ( index >= tokens.size() ) {
				throw runtime_error( "Not enough arguments for flag " + flag + "!" );
			} else {
				T result = convertType<T> ( tokens[index] );
				return result;
			}
		}


		template< class T>
		vector<T> parseVector( const string& flag, size_t begin = 0, size_t end = 0 )  {
			assertParsed();
			assertFlag( flag );
			const vector<string>& tokens = mTokenGroups[mFlagToIndex[flag]];
			if ( end == 0 ) {
				end = tokens.size();
			}
			if ( end <= begin ) {
				throw runtime_error( "Invalid range for flag " + flag + "!" );
			}
			if ( end > tokens.size() ) {
				throw runtime_error( "Not enough arguments for flag " + flag + "!" );
			}
			vector<T> result;
			result.reserve( end - begin );
			for ( size_t i = begin; i < end; ++i ) {
				result.push_back( convertType<T> ( tokens[i] ) );
			}
			return result;
		}


		// check whether a flag from a given group has been provided in the command line
		bool isSet( string flag ) {
			assertFlag( flag );
			return mFlagsProvided[mFlagToIndex[flag]];
		}



		void print() {
			assertParsed();
			vector<vector<string>> flags;
			size_t nrFlags = mTokenGroups.size();
			flags.resize( nrFlags );
			vector<bool> setbits( nrFlags, false );
			for ( auto args : mFlagToIndex ) {
				string flag = args.first;
				size_t index = args.second;
				flags[index].push_back( flag );
				setbits[index] = isSet( flag );
			}
			for ( auto i = 0; i < nrFlags; ++i ) {
				if ( setbits[i] ) {
					cout << "[*]";
				} else {
					cout << "[ ]";
				}
				for ( auto flag : flags[i] ) {
					cout << " " << flag ;
				}
				cout << " :";
				for ( auto token : mTokenGroups[i] ) {
					cout << " " << token;
				}
				cout << endl;
			}
		}

		void blockFlag( string flag ) {
			mBlockedFlags.insert( flag );
		}

		// returns the number of arguments for a given flag
		size_t nrTokens( string flag ) {
			assertParsed();
			assertFlag( flag );
			return mTokenGroups[mFlagToIndex[flag]].size();
		}

		void blockFlags( vector<string> flags ) {
			for ( string flag : flags ) {
				blockFlag( flag );
			}
		}

		// returns the tokens for a given flag
		vector<string> tokens( string flag ) {
			assertParsed();
			assertFlag( flag );
			return mTokenGroups[mFlagToIndex[flag]];
		}


		Parser subparser( string flag ) {
			assertParsed();
			Parser result( tokens( flag ) );
			for ( auto & item : mFlagToIndex ) {
				result.blockFlag( item.first );
			}
			return result;
		}

};


#endif
