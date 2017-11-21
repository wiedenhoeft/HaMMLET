#import argparse
import sys
import os
from pyhammlet.plotting import *
from pyhammlet.io import *
from matplotlib.ticker import MaxNLocator  
from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm, Colormap
import bisect

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument(
	'-i', 
	'--input-pattern', 
	dest='inPattern', 
	metavar="PATTERN",
	type=str,
	default=None,
	nargs=2,
	help="A pattern for the output files created by HaMMLET. If '-o PREFIX .EXT' is provided, the data files are assumed to be PREFIXmarginals.EXT, PREFIXparameters.EXT etc., and the output is PREFIX<start>-<end>.ext, where <start> and <end> are the first and last data position in the data file.")
parser.add_argument(
	'-o', 
	'--output-pattern', 
	dest='outPattern', 
	metavar="PATTERN",
	type=str,
	default=None,
	nargs=2,
	help="A pattern for the output files to be created. If '-o PREFIX SUFFIX' is provided, the output will be PREFIXstart-endSUFFIX. If -o is not provided, the pattern from -i will be used.")
parser.add_argument(
	'-R', 
	'--range', 
	dest='range', 
	metavar="PATTERN",
	type=int,
	default=[0, None],
	nargs=2,
	help="The range of positions to pe plotted: [START END).")
parser.add_argument(
	'-f', 
	'--data-file', 
	dest='datafile', 
	metavar="FILENAME",
	type=str,
	default="sample.csv",
	help='Name of the input data file (default: sample.csv). If -o is not provided, the name of the input file is split into basename and extension, and used as the values for -o.')
parser.add_argument(
	'-d', 
	'--dimensions', 
	dest='dimensions', 
	metavar=("WIDTH", "HEIGHT"),
	type=float,
	nargs=2,
	default=(10, 10),
	help='The width and height of the figure, in inches.')
parser.add_argument(
	'-r', 
	'--resolution', 
	metavar="DPI",
	dest='resolution', 
	type=int,
	default=300,
	help='The resolution, in dpi.')
parser.add_argument(
	'-s', 
	'--subfigures', 
	nargs='+',
	dest='subfigures', 
	default=["b", "dm", "msp"],
	help='List of strings describing subfigures.')	# TODO details
parser.add_argument(
	'-S', 
	'--split', 
	type=int,
	dest='splitsize', 
	default=None,
	help='Create one file for this many data points each (default: everything in one file).')	# TODO details
parser.add_argument(
	'-x', 
	'--xlabel', 
	dest='xlabel', 
	metavar="xlabel",
	type=str,
	default="Position",
	help='Label -f x-axis.')
parser.add_argument(
	'-y', 
	'--ylabels', 
	metavar="XLABELS",
	nargs='+',
	dest='ylabels', 
	default=["Iterations", "Data", "Marginal probabilities"],
	help='List of labels for the y-axes.')	# TODO details
parser.add_argument(
	'-p', 
	'--palette', 
	metavar="FILENAME",
	dest='palette', 
	default=os.path.join(sys.path[0], "palette.txt"),
	help='Name of a file for the color palette, with one #RRGGBB entry per line. If there are fewer colors in the palette than there are states, those missing states will be plotted in black.')	# TODO details
parser.add_argument(
	'-c', 
	'--chunksize', 
	metavar="CHUNKSIZE",
	dest='chunksize', 
	type=int,
	default=1,
	help='Size of input chunks.')

# parse and process command line arguments
args = parser.parse_args()
args.width=args.dimensions[0]
args.height=args.dimensions[1]
args.nrFigures = len(args.subfigures)
#assert args.nrFigures <= len(args.ylabels), "Number of y-labels is smaller than the number of figures!"
if args.nrFigures > len(args.ylabels):
	print("[WARNING] Number of y-labels is smaller than the number of figures!")
dataFilename = args.datafile
if args.inPattern is None:
	args.inPattern = os.path.splitext(args.datafile)
if args.outPattern is None:
	args.outPattern = [args.inPattern[0],""]
	args.outPattern[1] = args.inPattern[1]+".png"

sequencesFilename = "%ssequences%s" % tuple(args.inPattern)
marginalsFilename = "%smarginals%s" % tuple(args.inPattern)
blocksFilename = "%sblocks%s" % tuple(args.inPattern)



data=None
marginals = None
blocks=None
sequences=None

# TODO check which ones we need, using args.subfigures
T=None
NR_STATES=None		
NR_ITERATIONS=None

# set the above variables and check that they are not conflicting
# TODO make assignment within function work somehow
def setT(val):
	global T
	if T is not None:
		assert T == val, "Conflicting data sizes detected: %s, %s!" %(T, val)
	T=val
	
def setNrStates(val):
	global NR_STATES
	if NR_STATES is not None:
		assert NR_STATES == val, "Conflicting number of states detected: %s, %s!" %(NR_STATES, val)
	NR_STATES=val
	
def setNrIterations(val):
	global NR_ITERATIONS
	if NR_ITERATIONS is not None:
		assert NR_ITERATIONS == val, "Conflicting number of iterations detected: %s, %s!" %(NR_ITERATIONS, val)
	NR_ITERATIONS=val
		
computeMaxMargins=False # whether we need to calculate max margins
# load data and extract values for T, NR_STATES, and NR_ITERATIONS whenever possible
for string in args.subfigures:
	if string[0]=='d' and data is None:
		print "Loading data from %s" % dataFilename
		data = np.loadtxt(dataFilename)
		setT(len(data))
		if string[1]=='m':
			computeMaxMargins=True
	elif string[0]=='b' and blocks is None:
		print "Loading block sizes from %s" % blocksFilename
		blocks = readBlockSizes(blocksFilename)
		setT(sum([int(x) for x in file(blocksFilename).readline().split()]))
		setNrIterations(blocks.shape()[1])
	elif string[0]=='m' and marginals is None:
		print "Loading compressed state marginals from %s" % marginalsFilename
		marginals = readMarginals(marginalsFilename)
		setT( marginals.shape()[0])
		setNrStates(marginals.shape()[1])
		setNrIterations(marginals[0:1].decompress().sum())
	elif string[0]=='s' and sequences is None:
		print "Loading compressed state sequences from %s" % sequencesFilename
		sequences = readCompressedStateSequences(sequencesFilename)
		setNrIterations(len(sequences))
		# set number of states manually, since it could potentially increase later when reading marginals
		if NR_STATES is None:
			for seq in sequences:
				NR_STATES = max(NR_STATES, seq.array.max()+1)
	else:
		assert False, "Unknown sampling type: %s" % string


assert T is not None, "Could not determine data size, input is incomplete!"

# if we use colors for different states anywhere, compute the colormap TODO this shouldn't be necessary, the map can stay fixed


# Get palette from file.
# The default palette contains 56 colors. The first twelve colors are ColorBrewer's Paired12, and the others have been created with HalPal so that no two colors are closer to each other than the closest colors in Paired12 (CIEDE2000 >= 13.790368).
from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm, Colormap
# create a ListedColormap with norm from a palette and a number of classes
palette=[x.strip() for x in file(args.palette).readlines()]
norm = BoundaryNorm(range(len(palette)+1), len(palette))
cmap = ListedColormap(palette, name="HaMMLET")
cmap.set_bad("k")	
cmap.set_over("k")
cmap.set_under("k")
if NR_STATES is not None:
	if len(palette) < NR_STATES:
		print "[WARNING] Not enough colors in palette (%d) to represent %d states! Missing colors will be shown in black!"




# compute the maximum state margins if necessary
if computeMaxMargins:
	if marginals is None and sequences is None:
		assert False, "Need at least one of marginals and sequences to determine most common states."
	if marginals is not None:
		maxMargins = marginals.argmax().decompress()	#TODO get this from sequences if necessary
	else:	
		maxMargins = combineRLE(sequences).argmax().decompress()	# TODO implement getting marginals from sequences as an alternative (see above)
		
# determine first start and end positions
if args.splitsize is None:
	args.splitsize=T
if args.range[1] is None:	
	args.range[1]=T
start=max(0, args.range[0])
end = min(args.range[1], start+args.splitsize)

# create individual subplots
while T==-1 or start < end:	
	fig, axes = plt.subplots(args.nrFigures, figsize=(args.width, args.height), dpi=args.resolution, sharex=True, sharey=False, squeeze=False)
	figfile="%s%d-%d%s" % (args.outPattern[0], start, end-1, args.outPattern[1])
	print "Plotting to", figfile
	for i in xrange(args.nrFigures):
		plt.sca(axes[i,0])
		string = args.subfigures[i]
		if string == "b":	# blocks
			cm = matplotlib.cm.Greys_r
			cm.set_bad("k")
			cm.set_under("k")
			cm.set_over("k")
			nm = LogNorm(vmin=1)
			plotBlockSizes(
				blocks, 
				start=start, 
				end=end, 
				chunkSize=args.chunksize, 
				ylabel=args.ylabels[i],
				cmap=cm,
				norm=nm
				)	
		elif string == "dm":	# data with maximum colors
			plotData(
				data, 
				states=maxMargins, 
				start=start, 
				end=end, 
				cmap=cmap, 
				norm=norm, 
				ylabel=args.ylabels[i]) 
		elif string == "ds":	# data with a single color
			plotData(
				data, 
				states=None, 
				start=start, 
				end=end, 
				xlabel=args.xlabel, 
				ylabel=args.ylabels[i]
				)
		elif string[0] == "m":	# marginals
			if string[1] == "s":	# sort by state
				pass	# TODO
			elif string[1] == "f":	# sort by frequency of occurence
				pass	# TODO
			else:
				assert False, "Unknown sorting type \"%s\"!" % string[1]
			normalize = False
			if string[2] == "p":	# y-axis as probabilities
				normalize = True
			elif string[2] == "c":	# y-axis as counts
				normalize = False
			else:
				assert False, "Unknown y-scaling \"%s\"!" % string[2]
			plotMarginals(
				marginals, 
				start=start, 
				end=end,  
				cmap=cmap, 
				norm=norm, 
				normalize=normalize, 
				ylabel=args.ylabels[i]
				)
		elif string == "s":	# sequences
			plotSequences(
				sequences, 
				start=start, 
				end=end, 
				nrStates=NR_STATES, 
				cmap=cmap, 
				norm=norm, 
				ylabel=args.ylabels[i]
				)
	plt.xlabel(args.xlabel)
	
	# remove all but the bottom xlabels
	#for ax in axes[1:].reshape(-1):
		#ax.set_xticklabels([])
	#xticklabels = sum([ax.get_xticklabels() for ax in axes])
	#plt.setp(xticklabels, visible=False)
	
	for ax in axes.reshape(-1):
		nbins = len(ax.get_xticklabels())
		ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both'))
		
	fig.subplots_adjust(hspace=0)
	plt.savefig( figfile, bbox_inches='tight', dpi=args.resolution)
	plt.close()
	start = end
	end = min(start+args.splitsize, T, args.range[1])
	
