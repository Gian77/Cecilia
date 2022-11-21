#!/bin/bash 

# usage: combinePairs.sh file_R1.fastq file_R2.fastq
R1=$1
R2=$2

# Check if any unpaired paired reads are present.
c1=$(echo $(cat ${R1} | wc -l)/4 | bc)
c2=$(echo $(cat ${R2} | wc -l)/4 | bc)

if (($c1 == $c2)); 
then 
	echo -e "\nFiles contain the same number of reads. You re good to go.\n"; 
else 
	echo -e "\nFiles contain different number of reads. I will pair them for you!\n"; 
fi

# Extract the unique sequence names (but not the actual sequences, nor the quality scores) and 
# clean the headers to the unique original portion. Note: headers must be sorted.

grep "^@sample" $R1 | cut -f 2 -d" " | sort > ${R1%.*}.headers
grep "^@sample" $R2 | cut -f 2 -d" " | sort > ${R2%.*}.headers

# List all intersecting elements.
comm -12 ${R1%.*}.headers ${R2%.*}.headers > intersect.headers

# Extract those reads from the input files. The -v '^--$' because sometimes qaulity line can start 
# with 2 -- and that would screw the results.
grep -f intersect.headers -F -A3 $R1 | grep -v '^--$' > ${R1%.*}_pe.fastq
grep -f intersect.headers -F -A3 $R2 | grep -v '^--$' > ${R2%.*}_pe.fastq

