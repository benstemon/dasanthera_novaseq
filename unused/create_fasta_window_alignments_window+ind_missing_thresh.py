import argparse
import math
import os

# Parse command-line arguments
#If using "both", will filter individuals first, and then windows.
parser = argparse.ArgumentParser(description='Split an input file into windows and apply filters.')
parser.add_argument('-i', '--input', help='Input file name', required=True)
parser.add_argument('--windowsize', type=int, help='Size of the window', required=True)
parser.add_argument('--mode', choices=['window', 'individual', 'both'], help='Filter mode', required=True)
parser.add_argument('--individual', type=float, help='Individual missing data threshold (proportion). INDIVIDUALS are removed if they have missing data proportion > this proportion. This filter is applied BEFORE the window filter, if both filters are applied.')
parser.add_argument('--window', type=float, help='Window missing data threshold (proportion). WINDOWS are removed if any individuals in the alignment have missing data proportion > this proportion. This filter is applied AFTER the individual filter, if both filters are applied.')
parser.add_argument('--prefix', help='Prefix for output files', required=True)
args = parser.parse_args()

input_file = args.input
windowsize = args.windowsize
filter_mode = args.mode
window_missing_thresh = args.window
individual_missing_thresh = args.individual
prefix = args.prefix

# Read the input file to get alignment length and calculate the number of cycles
with open(input_file, 'r') as fasta_input:
    infile = fasta_input.read().splitlines()
    alignment_length = len(infile[1])
    cycles = int(math.ceil(alignment_length / windowsize))

# Loop through each window, apply filters, and generate output files
for i in range(1, cycles + 1):
    with open(input_file, 'r') as fasta_input:
        window_start = -windowsize + windowsize * i
        window_end = windowsize * i
        has_lines = True
        name_line = ''
        seq_line = ''

        if window_end < alignment_length:
            filename = f'{prefix}_bp_{window_start + 1}-{window_end}.fa'
        else:
            filename = f'{prefix}_bp_{window_start + 1}-{alignment_length}.fa'

        with open(filename, 'w') as fasta_output:
            delete_file = False
            while has_lines:
                fasta_line = fasta_input.readline()
                if len(fasta_line) == 0:
                    has_lines = False
                else:
                    fasta_line = fasta_line.strip()
                    if fasta_line.startswith('>'):
                        name_line = fasta_line
                    elif len(fasta_line) > 0:
                        seq_line = fasta_line
                        if name_line.startswith('>REF '):  # remove reference sequence
                            missing_count = 0
                        else:
                            individual_seq = seq_line[window_start:window_end]
                            missing_count = individual_seq.count('N')
                            missing_prop = missing_count / windowsize

                            if filter_mode == 'window' and window_missing_thresh is not None:
                                if missing_prop > window_missing_thresh:
                                    delete_file = True
                            elif filter_mode == 'individual' and individual_missing_thresh is not None:
                                if missing_prop > individual_missing_thresh:
                                    continue
                            elif filter_mode == 'both' and window_missing_thresh is not None and individual_missing_thresh is not None:
                                if missing_prop > individual_missing_thresh:
                                    continue
                                elif missing_prop > window_missing_thresh:
                                    delete_file = True

                            fasta_output.write(name_line + '\n')
                            fasta_output.write(individual_seq + '\n')

            if delete_file:
                os.remove(filename)
