import re
import glob
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Parse .fa.iqtree files')
parser.add_argument('--directory', type=str, help='Directory containing .fa.iqtree files', required=True)
parser.add_argument('--outfile', type=str, help='Output file name', required=True)
args = parser.parse_args()

# Open output file for writing
with open(args.outfile, 'w') as outfile:
    # Write header to output file
    outfile.write('File name\tTotal nucleotide sites\tParsimony informative sites\tInvariant sites\n')
    # Loop over all files with .fa.iqtree suffix in specified directory
    for filename in glob.glob(args.directory + '/*.iqtree'):
        with open(filename, 'r') as infile:
            # Extract file name from input file
            file_match = re.search(r'Input file name:\s+(\S+)', infile.read())
            if file_match:
                file_name = file_match.group(1)
            else:
                print(f'Error: File name not found in file {filename}')
                continue
            # Initialize counts
            total_sites = None
            pi_sites = None
            invariant_sites = None
            # Reset file pointer to beginning of file
            infile.seek(0)
            # Loop over lines in file
            for line in infile:
                # Check for total nucleotide sites
                if not total_sites:
                    match = re.search(r'(\d+)\s+nucleotide sites', line)
                    if match:
                        total_sites = int(match.group(1))
                        continue
                # Check for parsimony informative sites
                if not pi_sites and 'parsimony informative sites' in line:
                    match = re.search(r'(\d+)$', line)
                    if match:
                        pi_sites = int(match.group(1))
                        continue
                # Check for invariant sites
                if not invariant_sites:
                    match = re.search(r'Number of invariant (?:\(constant or ambiguous constant\) )?sites:\s+(\d+)', line)
                    if match:
                        invariant_sites = int(match.group(1))
                        continue
                # If all counts have been found, exit loop
                if total_sites is not None and pi_sites is not None and invariant_sites is not None:
                    break
            else:
                # If loop completed without finding all counts, report error
                print(f'Error: Count(s) not found in file {filename}')
                continue
        # Write results to output file
        outfile.write('{}\t{}\t{}\t{}\n'.format(file_name, total_sites, pi_sites, invariant_sites))

