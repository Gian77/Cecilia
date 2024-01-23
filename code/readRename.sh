#!/bin/bash

# Copyright © 2023 Gian M.N. Benucci, Ph.D.
# email: benucci@msu.edu URL: https://github.com/Gian77

# This script crename fastq files according to a $prefix name and also
# all the reads in each file progressively from $count to N read in the file.

fastq_directory=$1  # Specify the path to the directory containing the FASTQ files
prefix="sample"  # Specify the desired prefix for the renamed files
count=1  # Specify the starting number for the renaming sequence

# Output file for storing old and new file names
output_file="$fastq_directory"/file_mapping.txt

# Loop through the files in the directory
for file in "$fastq_directory"/*.fastq; do
    # Check if the file is a regular file
    if [[ -f "$file" ]]; then
        # Pad the count number with leading zeros
        count_padded=$(printf "%03d" "$count")
	 echo $count_padded
        sample_name="${prefix}${count_padded}"
        new_name="${sample_name}.fastq"
        
	# Write old and new file names into the output file
        echo -e "$(basename "$file")\t$new_name" >> "$output_file"

        # Rename the file
        mv "$file" "$fastq_directory/$new_name"
        
	# Update read headers in the renamed file
        awk -v new_prefix="$sample_name" \
            -v count_start="$count" '
            {
                if (NR % 4 == 1) {
                    header = sprintf("%s.%d", new_prefix, (NR + 3) / 4 + count_start - 1)
                    printf "@" header "\n"
                } else {
                    print
                }
            }' "$fastq_directory/$new_name" > temp.fastq
        mv temp.fastq "$fastq_directory/$new_name"
        
        count=$((count + 1))
    fi
done

#Explanation:

#The -v new_prefix="$sample_name" option passes the value of the sample_name variable from the Bash script to the new_prefix variable in awk.
#The -v count_start="$start_number" option passes the value of the start_number variable from the Bash script to the count_start variable in awk.
#The main code block inside single quotes (') is executed for each line of the input file.
#Explanation of the main code block:

#{ if (NR % 4 == 1) { ... } else { ... } } is a conditional statement that checks if the line number (NR) modulo 4 is equal to 1. This condition identifies the read header lines since they appear as the first line in every 4-line record.
#header = sprintf("%s.%d", new_prefix, (NR + 3) / 4 + count_start - 1) constructs the updated read header by combining the new_prefix (sample name) and the progressive number. (NR + 3) / 4 + count_start - 1 calculates the progressive number by dividing the line number by 4 (since each record has 4 lines) and adding count_start - 1 to adjust for the starting number.
#printf "@" header "_\n" prints the updated read header lines, with the old header completely removed. It uses printf instead of print to ensure that the line is printed without a trailing newline character.
#print is executed for all other lines (sequence, quality header, and quality values), and it simply prints those lines without any modification.
