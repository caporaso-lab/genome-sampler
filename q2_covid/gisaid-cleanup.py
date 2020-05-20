#!/usr/bin/env python

# This script has been tested with scikit-bio 0.5.6.

from sys import argv
import skbio
import skbio.io
import skbio.sequence

if len(argv) != 3:
    print("Usage: gisaid-cleanup.py input-path output-path")
    exit()

count_w_space = 0
count_too_much_N = 0
count_written = 0
count_read = 0


class GisaidDNA(skbio.sequence.GrammaredSequence):
    gap_chars = set('-. ')
    complement_map = skbio.DNA.complement_map
    definite_chars = skbio.DNA.definite_chars
    degenerate_map = skbio.DNA.degenerate_map
    default_gap_char = skbio.DNA.default_gap_char


seq_gen = skbio.io.read(argv[1],
                        format='fasta',
                        lowercase=True,
                        constructor=GisaidDNA)

with open(argv[2], 'w') as output_f:
    for s in seq_gen:
        count_read += 1
        if ' ' in s:
            count_w_space += 1
            continue

        s = s.degap()

        frac_n = s.frequencies(chars={'N'}, relative=True)['N']
        if frac_n >= 0.1:
            count_too_much_N += 1
            continue

        s.metadata['id'] = ' '.join([s.metadata['id'], s.metadata['description']]).strip().replace(' ', '_')
        del s.metadata['description']
        s = skbio.DNA(str(s), metadata=s.metadata)
        s.write(output_f, format='fasta')
        count_written += 1

print("Count read: %d" % count_read)
print("Count discared (contains spaces): %d" % count_w_space)
print("Count discarded (>10%% N): %d" % count_too_much_N)
print("Count written: %d" % count_written)
