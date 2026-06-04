#!/usr/bin/env python3
import sys
import re

rank_prefixes = ['d', 'p', 'c', 'o', 'f', 'g', 's']

input_fasta  = sys.argv[1]
output_fasta = sys.argv[2]

with open(input_fasta) as fin, open(output_fasta, 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            parts = line[1:].strip().split(None, 1)
            seq_id = parts[0]
            tax_string = parts[1] if len(parts) > 1 else ''

            taxa = [t.strip() for t in tax_string.split(';') if t.strip()]

            cleaned = []
            for t in taxa:
                t = t.replace(',', '')      # remove commas (SINTAX delimiter)
                t = t.replace(' ', '_')     # replace spaces with underscores
                t = re.sub(r'[^\w\-\.]', '_', t)  # sanitize any other special chars
                cleaned.append(t)

            tagged = ','.join(
                f"{rank_prefixes[i]}:{cleaned[i]}"
                for i in range(min(len(cleaned), len(rank_prefixes)))
            )

            fout.write(f">{seq_id};tax={tagged}\n")
        else:
            fout.write(line)

print("Done.")

