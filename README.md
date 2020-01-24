![HaMMLET](https://github.com/wiedenhoeft/HaMMLET/blob/dev/logo/logo-inv-noborder.png)

HaMMLET – Fast Bayesian HMM segmentation for big data
=====================================================

This software implements Forward-Backward Gibbs sampling for Bayesian segmentation in Hidden Markov Models (HMM). It uses dynamic wavelet compression to drastically improve convergence and memory consumption, making inference possible on large-scale data. 

For instance, HaMMLET can be used on a regular laptop for segmentation of genomic data, such as array-CGH or depth-of coverage from whole-genome sequencing (WGS), to find candidates for copy-number variants (CNV). For details, please refer to the doc/ directory.

For implementation details and the theory behind this approach, please refer to my [thesis](https://rucore.libraries.rutgers.edu/rutgers-lib/59275/) (DOI: [10.7282/t3-4e1k-ph18](https://doi.org/doi:10.7282/t3-4e1k-ph18)).
