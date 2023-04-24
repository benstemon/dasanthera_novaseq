# Set input and output file names
input_file = "features_in_shared_regions.bed"
output_file = "features_in_shared_regions-GENESONLY.bed"

# Open input and output files
with open(input_file, "r") as input, open(output_file, "w") as output:
    # Loop through lines in input file
    for line in input:
        # Split line into columns
        columns = line.strip().split("\t")
        # Check if second column contains "maker" and third column contains "gene"
        if len(columns) >= 3 and "maker" in columns[1] and "gene" in columns[2]:
            # Write matching line to output file
            output.write(line)
