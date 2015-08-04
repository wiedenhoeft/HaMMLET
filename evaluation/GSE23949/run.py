#!/usr/bin/env python

from __future__ import division
import sys
import os

sys.path.append("..")
sys.path.append("../..")

hammlet = "python ../../hammlet.py"

from hammlet_plotting import *

from plotCBSresults import *




for prefix in ["chr20", "all_bt474"]:
	callstring = "/usr/bin/time -f \"exit status\t%%x\nmax resident size [kB]\t%%M\navg resident size [kB]\t%%t\nuser CPU time [s]\t%%U\nprocess CPU time [s]\t%%S\nwall clock time [s]\t%%e\nmajor page faults\t%%F\nminor page faults\t%%R\n\" -o '%s_CBS.time' Rscript ../DNAcopy_MergeLevels.R '%s'" % (prefix, prefix)
	os.system(callstring)

plotCBSresults("chr20", paletteFile="../palette.txt")
plotCBSresults("all_bt474", multiFile="all_bt474_multi.txt", paletteFile="../palette.txt")


os.system("%s -a 20 0.1 0.8 10 1 1 model.txt '%s.csv' -m all_bt474_multi.txt -Tp -L1 -b 900 -i 100 -o '%s'_HaMMLET -x 'Probe number' -y 'Log2 ratio'" % (hammlet, "all_bt474", "all_bt474")) 
os.system("%s -a 19 0.04 0.9 10 1 1 model.txt '%s.csv' -p -L1 -T -i 600 -o '%s'_HaMMLET_19states -x 'Probe number' -y 'Log2 ratio'" % (hammlet, "chr20", "chr20"))


