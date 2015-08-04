#!/usr/bin/env bash


for compression in 'uncompressed' 'compressed'
do
	for fmeasure in 'micro' 'macro'
	do
		filename=${compression}_fmeasure_${fmeasure}F.tsv
		if [ -e ${filename} ] 
		then
			echo "${filename} already exists."
		else
			{
			for folder in */ 
			do 
				(cat ${folder}/${filename} | tr '\n' '\t') >> $filename
				echo "" >> $filename
			done
			} &
		fi
	done
done

