from __future__ import division
import sys
import os
sys.path.append("..")
hammlet = "python ../../../hammlet.py "
import conf
import yaml
import numpy as np
from Fmeasures import *
import itertools
import copy
if __name__=="__main__":
			
	
	import matplotlib
	# Force matplotlib to not use any Xwindows backend.
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	
	
		
	### parse the config file
	conf = conf.Config()
	conf.simSampleFile = os.path.abspath(conf.simSampleFile)
	conf.modelFile = os.path.abspath(conf.modelFile)


	
	## set means uniformly spaced between min and max of data
	modelString = ""
	for i in xrange(conf.nrCNstates):	
		modelString += "AutoNormal,0.2,0.9,s%d\n" % (i+1,)
	
	for i in xrange(conf.nrCNstates):
		row = [conf.transitions for x in xrange(conf.nrCNstates)]
		row[i] = conf.selfTransitions
		modelString += "Categorical,%d,%s,a%d\n" % (conf.nrCNstates, ",".join([str(x) for x in row]), i+1)
	modelString += "Categorical,%d,%s,pi\n" % (conf.nrCNstates, ",".join([str(conf.alphas) for i in xrange(conf.nrCNstates)]))
	modelString += "HMM,%d,pi,%s,%s\n" % (conf.nrCNstates, ",".join(["a%d"%i for i in xrange(1, conf.nrCNstates+1)]), ",".join(["s%d"%i for i in xrange(1, conf.nrCNstates+1)]))
	f = file(conf.modelFile, "w")
	f.write(modelString)
	f.close()	
	
	
	### get simulated state sequence
	trueStateSeq = np.loadtxt(conf.simSampleFile, dtype="int", usecols=(0,))
	trueStateSeqByClass = splitLabels(trueStateSeq, conf.nrCNstates)
	
	
	
	for trial in conf.trials:
		for benchmark in conf.benchmarks:
			
			prefix = trial+"_"+benchmark
			
			# write the true means and variances to a file format NumPy understands
			np.savetxt(prefix+"_truetheta.csv", np.column_stack((conf.means, conf.variances)))
			
			if trial not in ["uncompressed", "compressed"]:
				print "[ERROR] Unknown trial", trial
			
			if benchmark not in ["time", "fmeasure"]:
				print "[ERROR] Unknown benchmark", benchmark
				
			
			if benchmark == "time":
				if trial=="compressed":
					callstring = "%s    --max-block-size %d    --blocks    --burnin %d    --states %d    --iterations %d    --time    --output '%s'  '%s' '%s'" % (hammlet, conf.maxBlockSize,   conf.burnin, conf.nrCNstates, conf.iterations, prefix, conf.modelFile, conf.simSampleFile)
				elif trial=="uncompressed":
					callstring = "%s    --compression-off    --burnin %d    --states %d    --iterations %d    --time    --output '%s' '%s' '%s'" % (hammlet, conf.burnin, conf.nrCNstates, conf.iterations, prefix, conf.modelFile, conf.simSampleFile)
			
			elif benchmark == "fmeasure":
				if trial=="compressed":
					callstring = "%s    --max-block-size %d    --blocks        --sequences    --burnin %d      --states %d    --iterations %d    --output '%s' '%s' '%s'" % (hammlet, conf.maxBlockSize,   conf.burnin, conf.nrCNstates, conf.iterations, prefix, conf.modelFile, conf.simSampleFile)
				elif trial=="uncompressed":
					callstring = "%s    --compression-off        --sequences    --burnin %d    --states %d    --iterations %d    --output '%s' '%s' '%s'" % (hammlet, conf.burnin, conf.nrCNstates, conf.iterations, prefix, conf.modelFile, conf.simSampleFile)
			
			
			try:
				os.system(callstring)				
			except:
				print "[ERROR] Failed to run", callstring
				pass
			
			
			if benchmark == "fmeasure":
				FmeasuresFromSamples(conf.simSampleFile, prefix, conf.nrCNstates, conf.iterations, burnin=conf.burnin, plot=False)
				
		
				# cleanup
				if benchmark == "fmeasure":
					try:
						os.remove(prefix+"_statesequence.csv")
					except:
						pass
					try:
						os.remove(prefix+"_last_statesequence.csv")
					except:
						pass
					try:
						os.remove(prefix+"_statesequence.tsv")
					except:
						pass
					try:
						os.remove(prefix+"_last_statesequence.tsv")
					except:
						pass
					