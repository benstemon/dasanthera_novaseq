import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Script to match gene coordinates with functional annotations.")
parser.add_argument("--bedfile", help="Path to the .bed file")
parser.add_argument("--genesfile", help="Path to the genes of interest .txt file")
parser.add_argument("--gfffile", help="Path to the .gff file")
parser.add_argument("--outputfile", help="Path to the output file")
args = parser.parse_args()

# Read the .bed file and store the gene model names, scaffold, and coordinates in a dictionary
gene_coordinates = {}
with open(args.bedfile, "r") as bed_file:
    for line in bed_file:
        fields = line.strip().split("\t")
        scaffold = fields[0]
        gene_name = fields[3]
        start = int(fields[1])
        end = int(fields[2])
        gene_coordinates[gene_name] = (scaffold, start, end)

# Read the genes of interest from the separate .txt file
genes_of_interest = []
with open(args.genesfile, "r") as genes_file:
    for line in genes_file:
        gene = line.strip()
        genes_of_interest.append(gene)

# Read the functional annotations from the .gff file and store in a dictionary
functional_annotations = {}
with open(args.gfffile, "r") as gff_file:
    for line in gff_file:
        if not line.startswith("#"):
            fields = line.strip().split("\t")
            if len(fields) >= 9 and fields[2] == "mRNA":
                gene_id = fields[8].split(";")[0].split("=")[1]
                note = fields[8].split("Note=")[1]
                functional_annotations[gene_id] = note

# Create a new output file to store the results
output_file = open(args.outputfile, "w")

# Look up the coordinates for each gene of interest and write to the output file
for gene in genes_of_interest:
    if gene in gene_coordinates:
        scaffold, start, end = gene_coordinates[gene]
        annotation = functional_annotations.get(gene, "N/A")
        output_file.write(f"{scaffold}\t{start}\t{end}\t{gene}\t{annotation}\n")
    else:
        output_file.write(f"Gene: {gene} not found in the .bed file.\n")

# Close the output file
output_file.close()
