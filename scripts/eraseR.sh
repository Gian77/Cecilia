#!/bin/bash

# Copyright © 2023 Gian M.N. Benucci, Ph.D.
# email: benucci@msu.edu URL: https://github.com/Gian77

# This is a simple script to remove an R that appear at the beginning
# of  the file names so that Cecilia does not get confused with the
# R1 and R2 that mark the forward or revrese read

# use: sh eraseR.sh <path-to-files>

directory_path=$1 

for file in "$1"/R*.fastq.gz; do
    new_name="$1/M${file##*/R}"
    mv "$file" "$new_name"
done
