#!/usr/bin/env python
# -*- coding: utf8 -*-

###     Copyright 2014-2015 John Wiedenhoeft, Eric Brugel
###
###     This file is part of HaMMLET.
### 
###     HaMMLET is free software: you can redistribute it and/or modify
###     it under the terms of the GNU General Public License as published by
###     the Free Software Foundation, either version 3 of the License, or
###     (at your option) any later version.
### 
###     HaMMLET is distributed in the hope that it will be useful,
###     but WITHOUT ANY WARRANTY; without even the implied warranty of
###     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###     GNU General Public License for more details.
### 
###     You should have received a copy of the GNU General Public License
###     along with HaMMLET.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division, print_function
import sys
if sys.version_info < (2,6,0):
	print("[ERROR] HaMMLET requires at least Python 2.6.")
	sys.exit()
import argparse
import os


parser = argparse.ArgumentParser(
	description='Dynamically compressed Forward-Backward Gibbs sampling using Haar wavelets (http://bioinformatics.rutgers.edu/Software/HaMMLET/, URL is case-sensitive). NOTE: This version of HaMMLET was pulled from the biorxiv branch, and is provided for reproducibility purposes only; it is a little rough around the edges, might be missing features, and contain unfixed bugs. We do not recommend using this version in a production environment. For a maintained, up-to-date version, please checkout the master branch.',
	epilog="HaMMLET is ©2014-2015 John Wiedenhoeft, Eric Brugel, and distributed under GNU GPL. It uses rnglib by Pierre L'Ecuyer and Serge Cote, implemented in C++ by John Burkardt (http://people.sc.fsu.edu/~jburkardt/f_src/rnglib/rnglib.html), which is distributed under GNU LGPL. \n\nWhen using HaMMLET, please cite the following paper: John Wiedenhoeft, Eric Brugel and Alexander Schliep: \"Fast Bayesian Inference of Copy Number Variants using Hidden Markov Models with Wavelet Compression\". bioRxiv preprint (2015). http://dx.doi.org/10.1101/023705. ")
parser.add_argument("modelFile",
		   metavar="model",
		   type=str,
		   help="A file containing the model description. If automatic priors are used, this is the file HaMMLET writes the model to.")
parser.add_argument("files",
		    metavar="input",
		    nargs="+",
		    type=str,
		   help="A list of input files. Each file is assumed to represent one dimension in the data. Columns have to be sorted by position. If files contain two columns, they are assumed to represent xy-coordinates (e.g. positions on the chromosome and hybridization levels). One-column files are assumed to represent y-values. NOTE: at the moment, only univariate data (i.e. a single file) is supported.")	# TODO rephrase for multivariate support
parser.add_argument('-a', '--auto-priors',
		    dest="autopriors",
		    metavar=("NRSTATES", "VAR", "P", "SELFTRANS", "TRANS", "PI"),
		    nargs=6,
		    default=None,
		    help="Parameters for autopriors (desired variance, p, self-transitions, transitions, pi). A model file will be created and written to <model>, possibly overwriting an existing model file of the same name."
		    )
parser.add_argument('-A', '--auto-priors-mad',
		    dest="autopriorsmad",
		    metavar=("NRSTATES", "P", "SELFTRANS", "TRANS", "PI"),
		    nargs=5,
		    default=None,
		    help="Same as -a, but the variance will be estimated from the data itself, by taking the MAD of the finest detail coefficients (experimental)."
		    )
parser.add_argument('-s', '--states',
		    dest="nrStates",
		    metavar="INT",
		    type=int,
		    help="Number of states, must be supplied when using -P (only in developer version, deprecated)."
		    )


outputGroup = parser.add_argument_group("Output options", "The general format for any output file is <PREFIX>_<SUFFIX>, where PREFIX is set manually and SUFFIX is specific to each option.")
outputGroup.add_argument('-o', '--output',
		    dest="output",
		    metavar="PREFIX",
		    default=None,
		    type=str,
		    help="Prefix for output path (default is the first input file without its file suffix). Note that this is a filename, not a directory, so use outdir/prefix format. If no output prefix is provided, the path of the first input file is used."
		    )
outputGroup.add_argument('-S', '--sequences',
		    dest="sequences",
		    action="store_true",
		    help="By default, only the marginal counts for each state per position will be output as <PREFIX>_statesequence.csv, and the last sampled state sequence as <PREFIX>_last_statesequence.csv. Use this option if you want to write every single state sequence sample to <PREFIX>_statesequence.csv (Warning: this involves a lot of write operations and will be very slow)."
		    )
outputGroup.add_argument('-P', '--parameters',
		    dest="parameters",
		    action="store_true",
		    help="Output and plot the posterior samples for all model parameters as well (slow). Each sample will be written to a file <PREFIX>_<ID>.csv."
		    )
outputGroup.add_argument('-k', '--blocks',
		   dest="blocks",
		   action="store_true",
                   help="Output the coarsest possible block structure as induced by the estimated noise variance to <PREFIX>_blocks.csv."
                   )
outputGroup.add_argument('-p', '--plot',
		   dest="plot",
		   action="store_true",
                   help="Plot the results, depending on what other output options have been selected. If -P is selected, convergence plots for the parameter samples are also included."
                   )
outputGroup.add_argument('-m', '--multiple-plots',
		   dest="multiplots",
		   metavar="FILE",
		   type=str,
		   default=None,
                   help="If the plot is to be split up into multiple plots, e.g. one for each chromosome, this file contains the lengths in the first and the names in the second column."
                   )
outputGroup.add_argument('-x', '--xlabel',
		   dest="xlabel",
		   metavar="XLABEL",
		   type=str,
		   default="Position",
                   help="Label for the x-coordinate."
                   )
outputGroup.add_argument('-y', '--ylabel',
		   dest="ylabel",
		   metavar="YLABEL",
		   type=str,
		   default="Position",
                   help="Label for the y-coordinate."
                   )
samplingGroup = parser.add_argument_group("Sampling parameters", "Various parameters for the sampling process.")
samplingGroup.add_argument("-i", "--iterations",
		   dest="iterations",
		   metavar="INT",
		   default=1000,
		   type=int,
		   help="The number of sampling iterations (default 1000).")
samplingGroup.add_argument('-b', '--burnin',
		    dest="burnin",
		    metavar="INT",
		    type=int,
		    default=0,
		    help="Number of burn-in iterations (default=0). The total number of iterations for the sampler is the sum of -b and -i." 
		    )
samplingGroup.add_argument('-t', '--thinning',
		    dest="thinning",
		    metavar="INT",
		    type=int,
		    default=1,
		    help="Thinning parameter (t>0): only the iterations that divide this number will be considered (default=1, i.e. no thinning)."
		    )



treeGroup = parser.add_argument_group("Wavelet tree parameters", "These options affect the wavelet tree data structure itself.")
treeGroup.add_argument('-c', '--compression-off',
		   dest="compression",
		   action="store_false",
                   help="By default, we use dynamic Haar wavelet compression for faster sampling. Use -c to sample points individually (this might take a very long time).")
treeGroup.add_argument('-B', '--max-block-size',
		    dest="maxblocksize",
		    metavar="INT",
		    type=int,
		    default=0,
		    help="Maximum size of compression block, in case of numerical instabilities. If this parameter is zero, this is ignored and blocks can be arbitrarily large (default=0). Typically, this option will not be necessary, as the implementation uses various techniques to be numerically stable, but it might be a good idea if -l is used. -B 1 is equivalent to -c, i.e. instead of using blocks of size 1, uncompressed sampling is used to avoid the overhead from the data structure."
		    )
treeGroup.add_argument('-L', '--skip-levels',
		    dest="skiplevels",
		    metavar="INT",
		    type=int,
		    default=0,
		    help="The number of levels to skip when recursively calculating the maximum subtree coefficient. This is used in cases where vastly different variances lead to oversegmentation of high-variance components, thus reducing compression. Setting this to 1 is usually safe, as the first level contains almost only noise, and setting it to larger values is usually not required. Default is 0."
		    )







debugGroup = parser.add_argument_group("Developer options", "The following options are for debugging and development purposes and are usually not required.")
debugGroup.add_argument('-l', '--logs',
		    dest="logs",
		    action="store_true",
		    help="If this option is used, probabilities will be calculated in log-space, as one would usually do (for testing and debugging purposes only)."
		    )
debugGroup.add_argument('-d', '--debug',
		    dest="debug",
		    action="store_true",
		    help="Debuggin mode, saves various things to files for regression tests."
		    )
debugGroup.add_argument('-T', '--time',
		    dest="time",
		    action="store_true",
		    help="On Linux systems, use /usr/bin/time to benchmark time and memory consumption of the software (without plotting etc.) and write this to file. Use the time command outside this Python script to include the total runtime with plotting etc."
		    )
debugGroup.add_argument('-V', '--valgrind',
		    dest="valgrind",
		    action="store_true",
		    help="Use valgrind for debugging. This is VERY slow, use only a handful of iterations!"
		    )
debugGroup.add_argument('-G', '--gdb',
		    dest="gdb",
		    action="store_true",
		    help="Use gdb for debugging."
		    )
debugGroup.add_argument('-g', '--git-hash',
		    dest="git",
		    action="store_true",
		    help="Append the current git hash to the filename prefix."
		    )



def autoPriors(nrStates, variance, p, selfTransitions, transitions, pi, filename):
	modelString = ""
	for i in xrange(nrStates):	
		modelString += "AutoNormal,%f,%f,s%d\n" % (variance, p, i+1,)
	
	for i in xrange(nrStates):
		row = [transitions for x in xrange(nrStates)]
		row[i] = selfTransitions
		modelString += "Categorical,%d,%s,a%d\n" % (nrStates, ",".join([str(x) for x in row]), i+1)
	modelString += "Categorical,%d,%s,pi\n" % (nrStates, ",".join([str(pi) for i in xrange(nrStates)]))
	modelString += "HMM,%d,pi,%s,%s\n" % (nrStates, ",".join(["a%d"%i for i in xrange(1, nrStates+1)]), ",".join(["s%d"%i for i in xrange(1, nrStates+1)]))
	f = file(filename, "w")
	print("Writing automatic priors to {0}".format(filename))
	f.write(modelString)
	f.close()	





args = parser.parse_args()

assert args.thinning>0, "[ERROR] Thinning parameter -t must be a positive integer."



if args.autopriors != None:
	nrStates, variance, p, selfTransitions, transitions, pi = args.autopriors
	autoPriors(int(nrStates), float(variance), float(p), float(selfTransitions), float(transitions), float(pi), args.modelFile)
	args.nrStates = int(nrStates)

if args.autopriorsmad != None:
	nrStates, p, selfTransitions, transitions, pi = args.autopriorsmad
	autoPriors(int(nrStates), -1, float(p), float(selfTransitions), float(transitions), float(pi), args.modelFile)
	args.nrStates = int(nrStates)


if args.output == None:
	args.output = os.path.splitext(args.files[0])[0]
	
# append git hash to file prefix if -g option is set
if args.git:
	from subprocess import Popen, PIPE
	dataDir = os.path.abspath(os.getcwd())
	hammletDir = os.path.split(sys.argv[0])[0]
	
	os.chdir(hammletDir)
	gitproc = Popen(['git', 'rev-parse', 'HEAD'], stdout = PIPE)
	(stdout, stderr) = gitproc.communicate()
	for row in stdout.split('\n'):
		hash = row.strip(" \t\n")
		print("git hash is {0}".format(hash))
		args.output += "_" + hash
		os.chdir(dataDir)
		break
	


args.output = os.path.abspath(args.output)
args.files = [os.path.abspath(x) for x in args.files]
outdir = os.path.dirname(args.output)
if not os.path.exists(outdir):
	os.makedirs(outdir)


if args.maxblocksize == 1:
	print("[NOTE] Maximum block size is set to 1, using uncompressed option (-c).")
	args.compression = False
	args.maxblocksize = 0


if args.debug:	# debug mode, only one iteration
	print("===== DEBUG MODE =====")
	callstring = "'%s' '%s' %d %d %d %d %d %d %d %d %d '%s' %s %d '%s'" % (
		os.path.join(os.path.dirname(os.path.abspath(__file__)), "./hammlet"),
		args.modelFile,
		1,
		args.compression,
		args.maxblocksize,
		int(args.logs),
		1-int(args.sequences),
		int(args.parameters),
		args.burnin,
		args.thinning,
                int(args.skiplevels),
		args.output,
		int(args.debug),
		int(args.blocks),
		"' '".join(args.files))
else:
	callstring = "'%s' '%s' %d %d %d %d %d %d %d %d %d '%s' %s %d '%s'" % (
		os.path.join(os.path.dirname(os.path.abspath(__file__)), "./hammlet"),
		args.modelFile,
		args.iterations+args.burnin,
		args.compression,
		args.maxblocksize,
		int(args.logs),
		1-int(args.sequences),
		int(args.parameters),
		args.burnin,
		args.thinning,
                int(args.skiplevels),
		args.output,
		int(args.debug),
		int(args.blocks),
		"' '".join(args.files))



if args.valgrind:
	callstring = ("valgrind --tool=callgrind --callgrind-out-file='%s_callgrind.out' " % args.output )+callstring

elif args.gdb:
	callstring = "gdb --args "+callstring
if args.time:
	callstring = ("/usr/bin/time -f \"command\t%s\ncall\t%%C\nexit status\t%%x\nmax resident size [kB]\t%%M\navg resident size [kB]\t%%t\nuser CPU time [s]\t%%U\nprocess CPU time [s]\t%%S\nwall clock time [s]\t%%e\nmajor page faults\t%%F\nminor page faults\t%%R\n\" -o '%s' " % (" ".join(sys.argv), args.output + ".time") ) + callstring
os.system(callstring)



args.outdir = os.path.dirname(args.output)

# plot parameters for convergence test if option -P was provided
if args.parameters:
	print("Plotting parameters...")
	import matplotlib
	from hammlet_plotting import *
	# Force matplotlib to not use any Xwindows backend.
	#matplotlib.use('Agg')
	import numpy as np
	
	files = [os.path.join(args.outdir, ("%s_s%d.txt" % (args.output, s))) for s in xrange(1, args.nrStates+1)]
	trueFile = args.output+"_truetheta.csv"
	trueTheta=None
	trueMeans=None
	trueVariances=None
	if os.path.exists(trueFile):
		trueTheta = np.loadtxt(trueFile)
		trueMeans=trueTheta[:,0]
		trueVariances=trueTheta[:,1]
	means = np.zeros((args.nrStates, args.iterations+args.burnin))
	variances = np.zeros((args.nrStates, args.iterations+args.burnin))
	for i in xrange(len(files)):
		f = files[i]
		theta = np.loadtxt(f, skiprows=1)
		means[i,:]=theta[:,0]
		variances[i,:]=theta[:,1]
	saveParameters(args.output+"_parameters.png", means, variances, trueMeans, trueVariances, args.burnin)

if args.plot or args.debug:
	
	print("Plotting results...")
	from hammlet_plotting import *

	blocksFile = None
	if args.blocks:
		blocksFile = args.output+"_blocks.csv"
	saveResultsFromFile(
		outFileName=args.output+".png", 
		dataFiles=args.files, 
		sequenceFile=args.output+"_statesequence.csv",  
		nrStates=args.nrStates, 
		blocksFile=blocksFile, 
		maxBlockSize = args.maxblocksize, 
		burnin=args.burnin, 
		allSequences=args.sequences,
		multipleOutputs=args.multiplots,
		xlabel=args.xlabel,
		ylabel=args.ylabel)
