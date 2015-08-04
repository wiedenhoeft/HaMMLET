import numpy as np
import pandas as pd
import os

# NOTE This script has been tested with pandas version 0.13.1.

os.system("wget -nc ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE16/matrix/GSE16_series_matrix.txt.gz && gunzip GSE16_series_matrix.txt.gz")
os.system("wget -nc ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE16/soft/GSE16_family.soft.gz && gunzip GSE16_family.soft.gz")
 

mdata = pd.read_csv("GSE16_series_matrix.txt", comment="!", sep="\t", skiprows=88)
platform = pd.read_csv("GSE16_family.soft", comment="!", sep="\t", usecols=[0, 4, 5], skiprows=147)
#data = pd.merge(platform, mdata, left_on="ID", right_on="ID_REF", how="inner")	# BUG this doesn't seem to work for some versions of pandas
data = pd.merge(platform, mdata, left_index=True, right_index=True, how="inner")
keys = list(data.columns.values)
data['CHR'] = data['CHR'].astype(float)


for name in keys:
	if name not in ["ID",  "CHR", "POS_GENOME_KB", "ID_REF"]:
		d = data[(data["CHR"]>=1) & (data["CHR"]<=23)][["CHR", name]]
		chrom = list(d["CHR"])
		d.rename(columns={name: 'log2ratio', 'CHR':'# chromosome'}, inplace=True)
		d = d[np.isfinite(d['log2ratio'])]
		d['# chromosome'] = d['# chromosome'].astype(int)
		d.to_csv(name+".csv", columns=["# chromosome", "log2ratio"], sep="\t", index=False, header=False)
		
		
		blocksFile = file(name+"_blocks.tsv" , "w")
		blocksFile.write("0\n")
		for i in xrange(1,len(chrom)):
			if chrom[i] != chrom[i-1]:
				blocksFile.write("%d\n" % i)
