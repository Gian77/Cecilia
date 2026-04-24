#!/usr/bin/env python3

# reformat_silva_for_sintax.py

import re
import sys

# SILVA rank prefixes in order
rank_prefixes = ['d', 'p', 'c', 'o', 'f', 'g', 's']

input_fasta  = sys.argv[1]
output_fasta = sys.argv[2]

with open(input_fasta) as fin, open(output_fasta, 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            # SILVA header: >AB001718.1.1521 Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus
            parts = line[1:].strip().split(None, 1)  # split on first whitespace
            seq_id = parts[0]
            tax_string = parts[1] if len(parts) > 1 else ''

            # Split taxonomy on semicolons, strip whitespace/trailing empties
            taxa = [t.strip() for t in tax_string.split(';') if t.strip()]

            # Pair with rank prefixes (only as many as available)
            tagged = ','.join(
                f"{rank_prefixes[i]}:{taxa[i]}"
                for i in range(min(len(taxa), len(rank_prefixes)))
            )

            fout.write(f">{seq_id};tax={tagged}\n")
        else:
            fout.write(line)

print("Done.")

