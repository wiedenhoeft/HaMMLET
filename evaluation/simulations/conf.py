#!/usr/bin/env python

import sys
import os
import os.path
import time
import shutil
from datetime import datetime
from dateutil.tz import tzlocal
import yaml
import copy

class Config:

	"""
	This is a class for serialization. All input is stored in a YAML file, e.g. config.yaml. In your main.py, add

		import conf
		conf = conf.Config()

	and call

		python main.py config.yaml
		
	Parameters can then be accessed as members, e.g. conf.inputfile. Also, for each call a hidden log directory .config.yaml_files is created, in which the calls are registered and main.py as well as its imports are saved for reproducibility.
	
	It contains the following internals:
	self.__internals__["pythonFile"]: absolute path of main.py
	self.__internals__["classFile"]: absolute path of this file (conf.py)
	self.__internals__["yamlFile"]: absolute path of the YAML file being processed
	There is also a ...Dir version of the above, i.e. without the base name
	"""

	def __init__(self):
		self.__internals__ = dict()
		
		python_file, config_file = sys.argv[0:2]
		if len(sys.argv) > 2:
			comment = " ".join(sys.argv[2:])
		else:
			comment = "--no comment--"


		if not python_file.endswith(".py"):
			python_file += ".py"
		python_file = os.path.abspath(python_file)
		python_dir = os.path.dirname(python_file)	# directory of the python file
		self.__internals__["pythonFile"] = python_file
		self.__internals__["pythonDir"] = python_dir
		self.__internals__["pythonBase"] = os.path.basename(python_file)
		#self.__python_dir__ = python_dir

		callenv_path = os.path.abspath(__file__)	# path of this file
		callenv_dir = os.path.dirname(callenv_path)	# directory of this file
		self.__internals__["classFile"] = callenv_path
		self.__internals__["classDir"] = callenv_dir
		self.__internals__["classBase"] = os.path.basename(callenv_path)
		callenv_symlink = os.path.join(python_dir, os.path.basename(callenv_path))

		config_file = os.path.abspath(config_file)
		config_dir = os.path.dirname(config_file)
		self.__internals__["yamlFile"] = config_file
		self.__internals__["yamlDir"] = config_dir
		self.__internals__["yamlBase"] = os.path.basename(config_file)
		os.chdir(config_dir)	# every output and logs from callenv.py will be written there




		# load YAML file
		if not config_file.lower().endswith(".yaml"):
				config_file += ".yaml"
		if not os.path.exists(config_file):
			print "ERROR: config file %s does not exist!" % config_file
			sys.exit()
		yaml_content = yaml.load(file(config_file, "r"))
		for key in ["__internals__", "__dict__", "dictCopy", "dumpDictCopy"]:
			assert not yaml_content.has_key(key), "YAML file contained illegal entry \"%s\", which is used by Config class!" % key

		# create output directories, use timestring if not specified in YAML file
		timestamp = datetime.now(tzlocal())
		timestring = datetime.strftime(timestamp, "%s")	# decimal part missing
		log_dir = list(os.path.split(os.path.join(config_dir, config_file+"_files")))
		log_dir = os.path.join(log_dir[0], "."+log_dir[1])	# append dot so it's hidden
		log_dir_ts = os.path.join(log_dir, timestring)
		if not os.path.exists(log_dir_ts):
			os.makedirs(log_dir_ts)
		else:
			print "[WARNING] log directory %s already exists!" % log_dir_ts
		imports_dir = os.path.join(log_dir_ts, "imports")
		if not os.path.exists(imports_dir):
			os.makedirs(imports_dir)
		else:
			print "[WARNING] imports directory %s already exists!" % imports_dir


		# copy all basic files there
		shutil.copy(config_file, log_dir_ts)

		# find all imports in all relevant files and copy them to the imports dir
		modules = [python_file]
		seen_modules = set()
		modules_to_zip = set([])
		while len(modules) > 0:
			module = modules.pop()
			
			if not module.endswith(".py"):
				module += ".py"
			if module in seen_modules:
				continue
			seen_modules.add(module)
			if os.path.exists(module):
				#print "\tadding module", module[0:-3]
				modules_to_zip.add(module)
				for line in file(module, "r"):
					line = line.split()
					if len(line) == 0:
						continue
					if line[0] == "import":
						modules.append(line[1])
						#print "\t\tfound module", line[1]
					elif line[0] == "from" and line[2] == "import":
						modules.append(line[1])
						#print "\t\tfound module", line[1]
		shutil.copy(python_file, log_dir_ts)
		for module in set(modules_to_zip):
			if module != python_file:
				shutil.copy(module, imports_dir)	

		timestring_file = file(os.path.join(log_dir, "log.txt"), "a")
		timestring_file.write("\t".join([timestring, datetime.strftime(timestamp, "%A, %Y-%m-%d %H:%M:%S (%Z)"), python_file, config_file, comment]))
		timestring_file.write("\n")
		timestring_file.close()

		self.__dict__.update(yaml.load(file(config_file, "r")))

	def setItem(self, key, value, overwrite=False):
		#assert overwrite or not self.__dict__.has_key(key)
		self.__dict__[key] = value
	
	def removeItem(self, key):
		self.__dict__.pop(key, None)
		
	def dictCopy(self, excludeKeys=[]):
		d = copy.deepcopy(self.__dict__)
		d.pop("__internals__", None)
		for k in excludeKeys:
			d.pop(k, None)
		return d
		
	def saveConf(self, filename=None, excludeKeys=[]):
		"""Store conf in YAML file (default: input YAML file)"""
		if filename==None:
			filename = self.__internals__["yamlFile"]
		yaml.dump(self.dictCopy(excludeKeys), file(filename, "w"), default_flow_style=False)
		
		
		
	def has_key(self, key):
		return self.__dict__.has_key(key)
	
