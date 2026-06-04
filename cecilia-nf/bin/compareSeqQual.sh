#!/bin/bash

# Copyright © 2022 Gian M.N. Benucci, Ph.D.
# email: benucci@msu.edu URL: https://github.com/Gian77

# Adapted form https://www.biostars.org/p/56088/
# This script compare the length of the sequence and the quality to debug
# errors in sequence assemblies or rename steps. It will output the number
# of lines where the length do not match the sequence length.

# Usage: sh compareSeqQual.sh <dir>

if test -f sequence_quality_nomatch.counts ; then
	rm sequence_quality_nomatch.counts 
	echo -e "\noverwriting output file!\n"
fi


cd $1

for file in *.fastq 
do 
	awk '{if(NR%4==2) print NR"\t"$0"\t"length($0)}' $file > ${file}_seqLength.res
	awk '{if(NR%4==0) print NR"\t"$0"\t"length($0)}' $file > ${file}_qualityLength.res
	echo "$file: `echo $(awk 'NR==FNR{a[$3]++;next}!a[$3]' ${file}_seqLength.res ${file}_qualityLength.res | wc -l)`" >> sequence_quality_nomatch.counts	
	rm ${file}_seqLength.res; rm ${file}_qualityLength.res
done 

sed -i "1i#SampleID No_Match_Seq-Qual" sequence_quality_nomatch.counts 
