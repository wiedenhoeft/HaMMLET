#!/usr/bin/env python

#creates boxplots of F-measures
#sys.argv[1] is the name of the yaml file in the subdirs, e.g. test.yaml or batch.yaml
from __future__ import division
import sys
sys.path.append("..")
import matplotlib
matplotlib.use('Agg')
from Fmeasures import *
import os
import numpy as np
import matplotlib.pyplot as plt
import yaml


try:
	hash = "_"+file("githash.txt").readlines()[1]
except:
	print "Not git hash found..."
	hash = ""




simYaml = os.path.abspath(sys.argv[1])
simFolder = os.path.dirname(simYaml)
folders = [os.path.join(simFolder, x) for x in filter(lambda x: x.isdigit(), os.listdir(simFolder))]


print "Found %d folders..." % len(folders)
yaml_content = yaml.load(file(os.path.join(simFolder, "0", os.path.basename(simYaml)), "r"))
iterations = yaml_content["iterations"] - yaml_content["burnin"]


def timesFromFiles(folders, compression, resultFilename):
	result = np.zeros(len(folders), dtype = float)
	for i in range(len(folders)):
		folder = folders[i]
		filename = os.path.join(folder, "%s.time" % (compression,))
		if os.path.exists(filename):
			for line in file(filename):
				if line.startswith("user CPU time"):
					result[i]=float(line.split("\t")[1])
					break
		else:
			print "[ERROR] No such file", filename
	np.savetxt(resultFilename, result)
	return result



timeCompressedFilename = os.path.join(simFolder, "times_compressed.csv")
timeUncompressedFilename = os.path.join(simFolder, "times_uncompressed.csv")
if os.path.exists(timeCompressedFilename):
	print timeCompressedFilename, "already exists."
	timeCompressed = np.loadtxt(timeCompressedFilename)
else:
	print "Creating", timeCompressedFilename
	timeCompressed = timesFromFiles(folders, "compressed_time", timeCompressedFilename)
if os.path.exists(timeUncompressedFilename):
	print timeUncompressedFilename, "already exists."
	timeUncompressed = np.loadtxt(timeUncompressedFilename)
else:
	print "Creating", timeUncompressedFilename
	timeUncompressed = timesFromFiles(folders, "uncompressed_time", timeUncompressedFilename)
	

compressionsFilename = os.path.join(simFolder, "compressions.csv")
if os.path.exists(compressionsFilename):
	print compressionsFilename, "already exists!"
	compressions = np.loadtxt(compressionsFilename)
else:
	print "Creating", compressionsFilename
	compressions = np.zeros(len(folders), dtype=float)
	for i in xrange(len(folders)):
		folder = folders[i]
		filename = os.path.join(folder, "compressed_fmeasure_compression_ratio.csv") 
		if os.path.exists(filename):
			compressions[i] = float(np.loadtxt(filename))
		else:
			print "[ERROR] No such file", filename
	np.savetxt(compressionsFilename, compressions)	
	


speedups = timeUncompressed/timeCompressed
print "total time compressed:", timeCompressed.sum()
print "total time uncompressed:", timeUncompressed.sum()
print "total speedup:", timeUncompressed.sum()/timeCompressed.sum()

plt.scatter(compressions, speedups, color="white", edgecolor="black")
plt.xlabel("Average compression ratio", fontsize=18)
plt.ylabel("Speedup", fontsize=18)
plt.xlim(0, plt.xlim()[1])
plt.ylim(0, plt.ylim()[1])
plt.tight_layout()
plt.savefig(os.path.join(simFolder, "compression_speedup%s.png" % hash), bbox_inches='tight')
plt.savefig(os.path.join(simFolder, "compression_speedup%s.pdf" % hash), bbox_inches='tight')
plt.close()

plt.scatter(timeCompressed, timeUncompressed, color="white", edgecolor="black")
plt.xlabel("Running time (compressed) [s]", fontsize=18)
plt.ylabel("Running time (uncompressed) [s]", fontsize=18)
maxval=max(plt.xlim()[1], plt.ylim()[1])
plt.xlim(0, maxval)
plt.ylim(0, maxval)
ax=plt.gca()
ax.set_aspect(1)
plt.tight_layout()
plt.savefig(os.path.join(simFolder, "times%s.png" % hash), bbox_inches='tight')
plt.savefig(os.path.join(simFolder, "times%s.pdf" % hash), bbox_inches='tight')
plt.close()



os.system("bash aggregateFmeasures.sh '%s'" % os.path.dirname(os.path.abspath(sys.argv[1])))

fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)

(umi, cmi, uma, cma) = ax.flatten()

plt.sca(umi)
plt.title("Uncompressed")
FmeasureQuantilePlots(os.path.join(simFolder, "uncompressed_fmeasure_microF.tsv"), iterations, hash=hash, ylabel="Micro-averaged F-measure", xlabel="")

plt.sca(uma)
FmeasureQuantilePlots(os.path.join(simFolder, "uncompressed_fmeasure_macroF.tsv"), iterations, hash=hash, ylabel="Macro-averaged F-measure")

plt.sca(cmi)
plt.title("Compressed")
FmeasureQuantilePlots(
	os.path.join(simFolder, "compressed_fmeasure_microF.tsv"), 
	iterations, 
	hash=hash, 
	ylabel="", 
	xlabel="",
	insetXlim=[0, 50], 
	insetYlim=[0.8, 1.1], 
	insetXticks=[0, 25], 
	insetYticks=[0.8, 1.0], 
	insetWidth="40%", 
	insetHeight="40%", 
	insetLoc=4
	)


plt.sca(cma)
FmeasureQuantilePlots(
	os.path.join(simFolder, "compressed_fmeasure_macroF.tsv"), 
	iterations, 
	hash=hash, 
	ylabel="",
	insetXlim=[0, 50], 
	insetYlim=[0.8, 1.1], 
	insetXticks=[0, 25], 
	insetYticks=[0.8, 1.0], 
	insetWidth="40%", 
	insetHeight="40%", 
	insetLoc=4
	)

fig.subplots_adjust(hspace=0.1)
fig.subplots_adjust(wspace=0.1)
plt.savefig(os.path.join(simFolder, "Fmeasures.png"), bbox_inches="tight")
plt.savefig(os.path.join(simFolder, "Fmeasures.pdf"), bbox_inches="tight")
plt.close()
