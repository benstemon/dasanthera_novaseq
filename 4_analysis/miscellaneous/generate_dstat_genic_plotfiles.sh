#!/bin/bash

# Set default values for variables
dstat_infile=""
CDS=""
outfile=""

# Loop through command-line arguments
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -d|--dstat_infile)
    dstat_infile="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--CDS)
    CDS="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--outfile)
    outfile="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    shift # past argument
    ;;
esac
done

# Check if all required variables have been entered
if [[ -z "$dstat_infile" || -z "$CDS" || -z "$outfile" ]]; then
    echo "Error: dstat_infile, CDS, and outfile must be entered"
    exit 1
fi

# Use user-entered values in commands
awk 'BEGIN {OFS = "\t"} NR>1 {print $1, $2-1, $3}' $dstat_infile | 
 bedtools coverage -a - -b $CDS | 
 cut -f7 | (echo -e "genic_fraction" && cat -) | 
 paste -d'\t' - $dstat_infile > $outfile
