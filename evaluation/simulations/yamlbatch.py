#!/usr/bin/env python

import yaml
import sys
import os
import copy
import itertools
import conf
import multiprocessing as mp


"""Unrolls a YAML file of variable: valuelist entries for batch processing. If valuelist is a key, the following entry is again considered a variable:valuelist. This is extremely useful for generating config files for batch evaluation in which different parameter combinations are used. For instance, if different classes are used, the parameters for these classes depend on which class is used.""" 

def expandYAML(d):
	# make sure the entries are lists (even constants must be represented as 1-element lists)
	for variable, values in d.iteritems():
		if not isinstance(values, list) and not isinstance(values, dict) and not isinstance(values, tuple):
			d[variable]=[values]
	#variables, values =  zip(*sorted(d.items()))
	variables = d.keys()
	values = d.values()
	for valuecombi in itertools.product(*values):
		d = dict(zip(variables, valuecombi))
		branching = False	# check if there are any variables below
		for i in xrange(len(valuecombi)):
			value = valuecombi[i]
			variable = variables[i]
			if isinstance(value, dict):	# if value is a dict {k:...}, then k is the actual value and ... is a dict of dependent variables
				branching=True
				assert len(value.keys())==1, "Values are the only key in their dict, something's wrong!"
				k = value.keys()[0]	# k is the value
				#d[variable] = [k]
				d[variable] = k
				for x in value[k].keys():
					assert not d.has_key(x), "Key %s already in use!" % str(x)
				d.update(value[k])
			else:
				d[variable]=[value]
		if branching:
			for e in expandYAML(d):
				yield e
		else:
			# remove [ ] from values
			for k, v in d.iteritems():
				assert len(d[k])==1
				d[k] = v[0]
			yield d

def runCalls(f):
	python = "python"
	print "%s: thread started" % f
	r = os.system(" && ".join(["%s '%s' '%s'" % (python, os.path.expanduser(call), f) for call in conf.__calls__]))
	print "%s: thread finished" % f
	return 0
	
	
if __name__=="__main__":
	conf = conf.Config()
	YAML = conf.dictCopy(excludeKeys=["__yamlname__", "__lists__", "__groups__", "__calls__", "__maximum__", "__nrThreads__"])
	
	d='.'
	filenames = sorted(set([os.path.join(d,o, os.path.basename(conf.__internals__["yamlFile"])) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o)) and not o.startswith(".")]))
	
	nrFiles = len(filenames)
	assert nrFiles <= conf.__maximum__, "This batch file would run on %d existing output files, which is larger than the limit of %d." % (nrFiles, conf.__maximum__)
	
	if nrFiles>0:
		print "Directory already contains folders with YAML files. Assuming they have been created during previous runs; delete them if you'd like to create the YAML files again!"
	else:
		results = []
		for x in expandYAML(YAML):
			results.append(x)
			nrFiles += 1
			assert nrFiles <= conf.__maximum__, "This batch file would produce more output files than the limit of %d." % (conf.__maximum__)
		

		print "Producing %d YAML files in enumerated subdirectories..." % nrFiles
		for c in xrange(len(results)):
			YAML = results[c]
			if not os.path.exists(str(c)):
				os.makedirs(str(c))
			filename = os.path.join(str(c), os.path.basename(conf.__internals__["yamlFile"]))
			filenames.append(filename)
			yaml.dump(copy.deepcopy(YAML), file(filename, "w"), default_flow_style=False)
		print "...done."
	
	if conf.__dict__.has_key("__calls__"):
		pool = mp.Pool(conf.__nrThreads__)
		pool.map(runCalls, filenames)	
		
		
		
		
# variables: even level, no dashes
# values: odd level, always dashes
#
# 01v:
#         - 4
#         - 5
#         - 6
# 02v:
#         - 7
#         - 8:
#                 11v:
#                         - 0
#                 12v:
#                         - 15:
#                                 17v:
#                                         - 19
#                                 18v:
#                                         - 20
#                                         - 21
#                         - 16
#                 13v:
#                         - 100
#         - 9:
#                 14v:
#                         - 200
# 03v:
#         - 10
# 
# yields:
#[('01v', 4), ('02v', 7), ('03v', 10)]
#[('01v', 4), ('02v', 8), ('03v', 10), ('11v', 0), ('12v', 15), ('13v', 100), ('17v', 19), ('18v', 20)]
#[('01v', 4), ('02v', 8), ('03v', 10), ('11v', 0), ('12v', 15), ('13v', 100), ('17v', 19), ('18v', 21)]
#[('01v', 4), ('02v', 8), ('03v', 10), ('11v', 0), ('12v', 16), ('13v', 100)]
#[('01v', 4), ('02v', 9), ('03v', 10), ('14v', 200)]
#[('01v', 5), ('02v', 7), ('03v', 10)]
#[('01v', 5), ('02v', 8), ('03v', 10), ('11v', 0), ('12v', 15), ('13v', 100), ('17v', 19), ('18v', 20)]
#[('01v', 5), ('02v', 8), ('03v', 10), ('11v', 0), ('12v', 15), ('13v', 100), ('17v', 19), ('18v', 21)]
#[('01v', 5), ('02v', 8), ('03v', 10), ('11v', 0), ('12v', 16), ('13v', 100)]
#[('01v', 5), ('02v', 9), ('03v', 10), ('14v', 200)]
#[('01v', 6), ('02v', 7), ('03v', 10)]
#[('01v', 6), ('02v', 8), ('03v', 10), ('11v', 0), ('12v', 15), ('13v', 100), ('17v', 19), ('18v', 20)]
#[('01v', 6), ('02v', 8), ('03v', 10), ('11v', 0), ('12v', 15), ('13v', 100), ('17v', 19), ('18v', 21)]
#[('01v', 6), ('02v', 8), ('03v', 10), ('11v', 0), ('12v', 16), ('13v', 100)]
#[('01v', 6), ('02v', 9), ('03v', 10), ('14v', 200)]
