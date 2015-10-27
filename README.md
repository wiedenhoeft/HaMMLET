HaMMLET
========

Bayesian Hidden Markov Model with Wavelet Compression.

**NOTE** This version of HaMMLET was pulled from the biorxiv branch, and is provided for reproducibility purposes only; it is a little rough around the edges, might be missing features, and contain unfixed bugs. We do not recommend using this version in a production environment. For a maintained, up-to-date version, please checkout the master branch.
 
When using HaMMLET, please cite the following reference:

John Wiedenhoeft, Eric Brugel, Alexander Schliep. Fast Bayesian Inference of Copy Number Variants using Hidden Markov Models with Wavelet Compression. bioRxiv (2015). DOI: 10.1101/023705.

**Requirements** Python 2.7 (tested with 2.7.3 and 2.7.6), NumPy (tested with 1.6.2 and 1.8.2), Matplotlib (tested with 1.1.1 and 1.3.1).

How to use HaMMLET
=======

HaMMLET is called via the Python script "hammlet.py"; never call the executable "hammlet" directly. We have provided a small sample file to demonstrate the main features of HaMMLET. Calling
	
        python hammlet.py -p -P -k -s 6 -i 100 sample_model.txt sample.csv

does the following: It takes the data in "sample.csv", runs the sampler for 100 iterations (-i option), and outputs the marginal counts of each state to "sample_statesequence.csv". It plots the result to a file called "sample.png" (-p option), as well as the sampled values for emission means and variances to "sample_parameters.png" (-P option). Additionally, the -P option writes the samples for each model parameter to "sample_*.txt", where * represents the name of the variable as specified in "sample_model.txt". Notice that the data only contains 5 different states, but the model contains 6 (-s option; this MUST be provided if -P is used, but can be omitted otherwise), thus you will see that the variates for the emission parameters of one state vary widely, as they are sampled from the prior and do not contribute to any CNV call, whereas the others converge almost instantly and yield flat curves in "sample_parameters.png". The -k option writes the minimum block boundaries to "sample_blocks.csv", which 
appear as grey vertical lines in "sample.png". Note that these are NOT the segment boundaries of the calls, but rather represent the largest blocks used during compression as derived from the MAD variance (see paper for details); this is mostly for evaluation purposes, and is of limited value in actual applications.

HaMMLET comes with a variety of options. For more information, see the manual page by running

        python hammlet.py -h
	



License
=======

Copyright 2014-2015 John Wiedenhoeft, Eric Brugel

HaMMLET is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


HaMMLET uses rnglib by Pierre L'Ecuyer and Serge Cote, implemented in C++ by John Burkardt (http://people.sc.fsu.edu/~jburkardt/f_src/rnglib/rnglib.html), which is distributed under GNU LGPL.
