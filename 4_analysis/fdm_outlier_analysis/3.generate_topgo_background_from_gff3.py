#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Extract GO terms for mRNA features from a GFF file')
parser.add_argument('-i', '--input', dest='input_file', help='Input GFF file', required=True)
parser.add_argument('-o', '--output', dest='output_file', help='Output file', required=True)

args = parser.parse_args()

mRNA_dict = {}

with open(args.input_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            line = line.rstrip().split('\t')
            if line[1] == 'maker' and line[2] == 'mRNA':
                mRNA_id = line[-1].split('ID=')[1].split(';')[0]
                go_terms = line[-1].split('Ontology_term=')[1].split(';')[0].split(',') if 'Ontology_term=' in line[-1] else []
                mRNA_dict[mRNA_id] = go_terms

with open(args.output_file, 'w') as f:
    for mRNA_id, go_terms in mRNA_dict.items():
        f.write(mRNA_id + '\t' + ','.join(go_terms) + '\n')
