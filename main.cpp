#include "hammlet-manpage.hpp"	// a hexdump of hammlet-manpage.txt wrapped into a variable declaration, so that -h/--help can be used without external file access.
#include "Tags.hpp"
#include "HMM.hpp"
// #include "Options.hpp"
#include "Parser.hpp"
#include "Emissions.hpp"
#include "Blocks.hpp"
#include "AutoPriors.hpp"
#include "Records.hpp"
#include "wavelet.hpp"
#include "StateSequence.hpp"
#include "Statistics.hpp"
#include "includes.hpp"

#include "utils.hpp"


// TODO generalize static sampling by setting minimum and maximum threshold to aavoid huge trellis as well as bias through overcompression?


int main( int argc, const char* argv[] ) {


	try {

		Parser args( argc, argv );



		// OPTIONS
		args.registerFlags( {"-v", "-verbose"} );
		args.registerFlags( {"-g", "-arguments"} );
		args.registerFlags( {"-h", "-help" , "--help"} );	// we leave --help an undocumented convenience


		// I/O
		args.registerFlags( {"-f", "-input-file"} );
		args.registerFlags( {"-o", "-output-pattern"}, "hammlet .csv" );
		args.registerFlags( {"-O", "-output-data"}, "marginals" );
		args.registerFlags( {"-w", "-overwrite"} );


		// MODEL
		args.registerFlags( {"-s", "-states"}, "3" );
		args.registerFlags( {"-e", "-emissions"}, "normal 0.2 0.9" );
		args.registerFlags( {"-a", "-auto-priors"} );
		args.registerFlags( {"-t", "-transitions"}, "0.5 0.5" );
		args.registerFlags( {"-S", "-no-self-transitions"} );
		args.registerFlags( {"-I", "-initial-dist"}, "0.5" );


		// SAMPLING SCHEME
		args.registerFlags( {"-R", "-random-seed"}, to_string( time( 0 ) ) );
		args.registerFlags( {"-i", "-iterations"}, "M 500 0 P S 150 0 S 250 1" );


		// COMPRESSION
// 		args.registerFlags( {"-y", "-data-structure"}, "B" );	// TODO provide alternative data structures
		args.registerFlags( {"-b", "-block-limits"}, "0 0" );
		args.registerFlags( {"-m", "-weight-multiplier"}, "1" );	// multiply weights by this factor, to avoid overcompression

		args.parseArgs();


		// Print the internal state of the parser if requested
		if ( args.isSet( "-g" ) ) {
			args.print();
		}

		// Set verbosity
		const bool verbose = args.isSet( "-v" );

		// Allow to overwrite existing files?
		const bool overwrite = args.isSet( "-w" );


		// Display manpage (as .txt) if -h or --help are provided, and exit
		if ( args.isSet( "-h" ) ) {
			cout << endl << doc_hammlet_manpage_txt << endl;
			return 0;
		}

		//// Input/output ////

		// Get the prefix for the output files
		const string outputPrefix = args.parse<string> ( "-o", 0 );
		const string outputSuffix = args.parse<string> ( "-o", 1 );





		//// Random number generator
		size_t rng_seed = args.parse<size_t> ( "-R", 0 );
		rng_t RNG( rng_seed );



		//// State descriptions, mappings, and numbers of states ////

		size_t p;	// number of emission parameters
		size_t d;	// data dimensions
		MappingType m;
		if ( args.nrTokens( "-s" ) == 1 ) {
			p = args.parse<size_t> ( "-s", 0 );
			d = 1;
			m = combinations;
		} else {
			m = args.parse<MappingType> ( "-s", 0 );
			p = args.parse<size_t> ( "-s", 1 );
			d = 1;
			if ( args.nrTokens( "-s" ) >= 3 ) {
				d = args.parse<size_t> ( "-s", 2 );
			}
		}

		const MappingType mappingType = m;
		const size_t nrDataDim = d;
		const size_t nrParams = p;
		Mapping mapping(
		    nrDataDim,
		    nrParams,
		    mappingType );
		const size_t nrStates = mapping.nrStates();




		//// Transition parameters ////

		// TODO accept full vector as well
		// default values of 0.5 correspond to Jeffreys prior
		const real_t selfTrans = args.parse<real_t> ( "-t", 1 );
		real_t t = selfTrans;
		if ( args.nrTokens( "-t" ) > 1 ) {
			t = args.parse<real_t> ( "-t", 0 );
		}
		const real_t trans = t;


		Transitions<DirichletVector> A( nrStates, RNG );
		TransitionHyperParam<DirichletParamVector> tau_A( nrStates, trans, selfTrans );

		const bool useSelfTrans = !args.isSet( "-S" );	// use self-transitions within block? otherwise 1 is used




		//// initial state distribution parameters ////
		// TODO accept full vector
		const real_t initialAlpha = args.parse<real_t> ( "-I", 0 );		// 0.5 is Jeffreys prior for Dirichlet
		Initial<Dirichlet> pi( nrStates, RNG );
		InitialHyperParam<DirichletParam> tau_pi( nrStates, initialAlpha );



		//// Compression ////

		const size_t chunkSize = max( ( size_t ) 1, args.parse<size_t> ( "-b", 0 ) );
		const size_t maxBlockSize = args.parse<size_t> ( "-b", 1 );

// 		const string dataStructure =  args.parse<string> ( "-y", 0 );	// TODO make command line option
		const string dataStructure = "B";
		if ( dataStructure == "B" ||  dataStructure == "breakpointarray" ) {
			// TODO any parameters to breakpoint array would go here




// 		} else
// 				if ( dataStructure == "wavelettree" ) {

			// TODO parameters to wavelet tree
// 			const size_t nrSkipLevels = args.parse<size_t>( "-w", 0, 0 );
// 		const size_t nrAdditionalLevels = args.parse<size_t>( "-w", 1, 0 );
// 		const bool useFullCompression = args.parse<bool>( "-w", 2, true );

		} else {
			throw runtime_error( "Unknown data structure \"" + dataStructure + "\", or not implemented yet!" );
		}

		const real_t weightMultiplier = args.parse<real_t>( "-m" );


		// Emissions
		vector<vector<real_t>> thetaParams;
		thetaParams.reserve( nrParams );
		const string emissionType = "normal";	// TODO parse this once more distribution types are implemented
		const bool autoPriors = args.isSet( "-a" );

		// parse parameters for creating automatic priors if -a is set
		if ( autoPriors ) {
			// parse var and p
			vector<real_t> thp = args.parseVector<real_t> ( "-e", 1, 3 );

			for ( size_t i = 0; i < nrParams; ++i ) {
				thetaParams.push_back( thp );
			}

			// TODO implement auto priors: general method for block means etc., leave specifics to template specialization
		} else {
			throw runtime_error( "Manual theta priors not implemented!" );
		}




		//// Parse the input data from a list of files or stdin. For multivariate data, dimensions are filled before filling the next position.




		//// Sampling scheme ////

		if ( verbose ) {
			cout << "Data dimensions: " << nrDataDim << endl;
			cout << "Emission distributions: " << nrParams << endl;
			cout << "States: " << nrStates << endl;
			cout << "Sampling scheme: " << concat( args.tokens( "-i" ), " " ) << endl;
			cout << "Random seed: " << rng_seed << endl;
		}



		// Check if additional output files should be created
		Parser outputArgs = args.subparser( "-output-data" );
		outputArgs.registerFlags( {"M", "marginals"} );	// marginals: segmentsize, counts for each state
		outputArgs.registerFlags( {"S", "sequences"} );
		outputArgs.registerFlags( {"P", "parameters"} );
		outputArgs.registerFlags( {"B", "blocks"} );
		outputArgs.registerFlags( {"C", "compression"} );
		outputArgs.registerFlags( {"D", "mapping"} );	// output the emission mappings for each state
		outputArgs.registerFlags( {"G", "segments"} );	// in each iteration: number of marginal segments, number of values used to store marginals (for diagnostics)
		outputArgs.parseArgs();


		// output state mappings to file if required
		if ( outputArgs.isSet( "D" ) ) {
			//TODO
		}


		// inputValues holds things like breakpoint weights, depending on the data structure being used
		vector<real_t> inputValues;


		// TODO allow to ignore invalid input?
		if ( emissionType == "normal" ) {	// univariate normal with automatic priors

			// create sufficient statistics for input data
			vector<SufficientStatistics<Normal>> stats;

			// TODO right now, individual files are concatenated. We should also allow multiple files to contain multiple dimensions.
			if ( args.isSet( "-f" ) ) { // read from input files
				for ( string fname : args.parseVector<string>( "-f" ) ) {	// iterate over input file names
					if ( verbose ) {
						cout << "Reading " + fname + "..." << endl << flush;
					}
					ifstream fin( fname );

					if ( fin ) {
						// TODO this can still lead to reallocation, fix later
						// TODO MaxletTransform does not work for multiple files in its current state
						MaxletTransform( fin, inputValues, stats, nrDataDim, inputValues.size() + nrLinesInFile( fin ) + 1 );
						// NOTE Reserving +1 is really important here! In the integral array, an element is appended to stats, and not reserving space for that element can lead to reallocations in the gigabyte range!

					} else {
						throw runtime_error( "Cannot read from input file " + fname + "!" );
					}
				}
			} else {	// read from STDIN
				if ( verbose ) {
					cout << "Reading from standard input..." << endl << flush;
				}
				MaxletTransform( cin, inputValues, stats, nrDataDim );
			}

			const size_t T = inputValues.size();

			if ( verbose ) {
				cout << "Number of data points: " + to_string( T ) << endl << flush;
			}

			Records records( T, outputPrefix, outputSuffix, nrStates );
			records.setRecordStateSequence( outputArgs.isSet( "sequences" ), overwrite );
			records.setRecordTheta( outputArgs.isSet( "parameters" ), overwrite );
			records.setRecordBlocks( outputArgs.isSet( "blocks" ), overwrite );
			records.setRecordCompression( outputArgs.isSet( "compression" ), overwrite );
			records.setRecordMarginals( outputArgs.isSet( "marginals" ), overwrite );
			records.setRecordSegments( outputArgs.isSet( "segments" ), overwrite );


			// inputValues holds the maxlet transform, now transform it to breakpoint weights
			if ( verbose ) {
				cout << "Calculating Haar breakpoint weights..." << endl << flush;
			}
			HaarBreakpointWeights( inputValues );

			// multiply weights
			for ( auto & w : inputValues ) {
				w *= weightMultiplier;
			}

			if ( dataStructure == "B" || dataStructure == "breakpointarray" ) {

				typedef Statistics<IntegralArray, Normal> S;
				typedef Blocks<BreakpointArray> B;
				S ia( stats, nrDataDim );
				B waveletBlocks( inputValues );

				Emissions<S, B> y( ia, waveletBlocks );



				// TODO this version calculates the same autopriors for all dimensions, adapt for flexible mapping
				thetaParams[0] = autoPrior( thetaParams[0][0], thetaParams[0][1], y );
				for ( auto & param : thetaParams ) {
					param = thetaParams[0];
				}

				ThetaHyperParam<NormalInverseGammaParam> tau_theta(
				    thetaParams );


				Theta<NormalInverseGamma> theta(
				    tau_theta,
				    nrDataDim,
				    mappingType,
				    RNG	// TODO pass to sampler instead of making it a member?
				);




				// TODO run a general check on the tokens to avoid running the sampler if there are parsing errors
				size_t nrTokens = 0;
				for ( auto c : args.tokens( "-i" ) ) {
					if ( c != "P" ) {
						nrTokens++;
					}
				}
				// check that iterations are grouped in triples
				if ( nrTokens % 3 != 0 ) {
					throw runtime_error( "Parameters for -i, excluding \"P\", must be multiples of 3!" );
				}
				nrTokens = args.nrTokens( "-i" );


				// get iteration types
				bool samplePrior = true;
				bool reuseBlocks = false;	// in a sequence of static smapling iterations, we want to reuse the previous block tructure, not resample it
				string method;
				size_t iterations, thinning;
				for ( size_t i = 0; i < nrTokens; ) {
					string method = args.parse<string> ( "-i", i );	// D=direct gibbs, M=mixture, F = forward-backward gibbs
					if ( method == "P" ) {
						samplePrior = true;
						i ++;
						continue;
					} else {
						if ( i + 2 >= nrTokens ) {
							throw runtime_error( "Incomplete command line for -i!" );
						}
						iterations = args.parse<size_t> ( "-i", i + 1 );
						thinning = args.parse<size_t> ( "-i", i + 2 );
						i += 3;
					}

					if ( verbose && samplePrior ) {
						cout << "Sampling prior..." << endl << flush;
					}

					if ( method == "S" ) {
						// sample static FBG. If thinning is 0, the priors are sampled and a burnin is performed so as to not use the previous parameters from mixture, for example.

						// create blocks from the last sampling
						if ( !reuseBlocks ) {
							cout << "Setting static block structure..." << endl << flush;
							y.createBlocks( theta );
							reuseBlocks = true;
						}
						
						if ( verbose ) {
							cout << "Sampling static Forward-Backward..." << endl << flush;
						}
						StateSequence< ForwardBackward > q( RNG );
						sampleHMM( y, tau_theta, theta, tau_A, A, tau_pi, pi, q, mapping, iterations, thinning, records, false, samplePrior, useSelfTrans );

// 						StateSequence< ForwardBackward > q( RNG );

					} else {
						reuseBlocks = false;
						if ( method == "F" ) {	// Forward-Backward sampling
							if ( verbose ) {
								// TODO more detailed output
								cout << "Sampling Forward-Backward...";
							}
							StateSequence< ForwardBackward > q( RNG );
							sampleHMM( y, tau_theta, theta, tau_A, A, tau_pi, pi, q, mapping, iterations, thinning, records, true, samplePrior, useSelfTrans );

						} else if ( method == "M" ) {	// Mixture sampling
							StateSequence< Mixture > q( RNG );
							if ( verbose ) {
								cout << "Sampling mixture..." << endl;
							}
							sampleHMM( y, tau_theta, theta, tau_A, A, tau_pi, pi, q, mapping, iterations, thinning, records, true, samplePrior,  useSelfTrans );

							/*} else if ( method == "D" ) {	// Direct Gibbs sampling
								StateSequence< DirectGibbs > q( RNG );
								sampleHMM( y, tau_theta, theta, tau_A, A, tau_pi, pi, q, mapping, iterations, thinning, records, i == 0, useSelfTrans );*/

						} else {
							throw runtime_error( "Unknown sampling type " + method + "!" );
						}
					}
					samplePrior = false;
				}
				// NOTE if marginals are to be saved, the output routine is automatically triggered by the destructor of records
			}
		}

		else {
			throw runtime_error( "Emission type " + emissionType + " unknown or not implemented yet!" );
		}


		return 0;

	} catch
		( exception& e ) {
		cerr << endl << flush << "\e[97;41;1m [ERROR] " << e.what()  << " \e[0m" << endl << endl << flush;

		return 1;
	}


}




