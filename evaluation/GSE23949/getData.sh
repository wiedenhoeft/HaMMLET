#!/usr/bin/env bash

wget -nc ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE23nnn/GSE23949/suppl/GSE23949_RAW.tar
tar -xvf GSE23949_RAW.tar
gunzip GSM590105_bt474.130209_252152911064_S01_CGH_105_Dec08.txt.gz




echo "Preparing files..."
cut --output-delimiter=" " -f 8,11 'GSM590105_bt474.130209_252152911064_S01_CGH_105_Dec08.txt'  | grep chr | grep -v X  | grep -v Y | grep -v random | grep -v M | sed 's/ /\t/g' | sed 's/\([0-9]\)-/\1\t/g'  | sed 's/:\([0-9]\)/\t\1/g' | sed 's/chr//g' | cut -f 1,2,4 | sort -n -k1,1 -n -k2,2 | cut -f 1,3 > all_bt474.csv 
cat all_bt474.csv | grep $'20\t' > chr20.csv
