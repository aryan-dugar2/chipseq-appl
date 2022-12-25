#!/bin/sh

# Load python
ml python

# Read relevant files
ls *.bam>bams.txt

# Define relevant file name
file="bams.txt"

# Call peaks!
while read -r line; do filename="${line:0:2}"; echo "macs2 callpeak -t $line --outdir $filename -g hs --bdg -q 0.05 -f BAM"; done < $file



