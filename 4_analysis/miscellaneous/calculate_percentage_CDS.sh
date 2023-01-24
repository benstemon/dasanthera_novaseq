#!/bin/bash

# Set default values for parameters
window_size=10000
sliding_increment=2500
genome_size="genomesize_scaffolds_davidsonii_1mb.txt"
CDS="annot_Pdavidsonii_1mb.gffread.genes.bed"
outfile="test3.bed"

# Prompt user to enter new values for parameters
read -p "Enter window size (default = $window_size): " new_window_size
read -p "Enter sliding increment (default = $sliding_increment): " new_sliding_increment
read -p "Enter genome size file name (default = $genome_size): " new_genome_size
read -p "Enter CDS.bed file name (default = $CDS): " new_CDS
read -p "Enter output file name (default = $outfile): " new_outfile

# Update parameters with new values entered by the user
if [ -n "$new_window_size" ]; then
  window_size=$new_window_size
fi
if [ -n "$new_sliding_increment" ]; then
  sliding_increment=$new_sliding_increment
fi
if [ -n "$new_genome_size" ]; then
  genome_size=$new_genome_size
fi
if [ -n "$new_CDS" ]; then
  CDS=$new_CDS
fi
if [ -n "$new_outfile" ]; then
  outfile=$new_outfile
fi

# Run commands using the updated parameters
bedtools makewindows -g $genome_size -w $window_size -s $sliding_increment | bedtools coverage -a - -b $CDS | awk '{print $1, $2, $3, $7}' > $outfile

echo "Output written to $outfile"
