#!/usr/bin/env python

import click
import pandas as pd
import numpy as np
import skbio

@click.command()
@click.option('--tsv', help='filepath to nextstrain formatted tsv metadata',
              type=click.File(mode='r'), required=True)
@click.argument('input', type=click.Path())
@click.argument('output', type=click.Path())
def main(input, output, tsv):
    df = pd.read_csv(tsv, sep='\t')
    df.index = df['gisaid_epi_isl']

    total = 0
    kept = 0

    def generator():
        for seq in skbio.io.read(input, format='fasta', constructor=skbio.DNA):
            nonlocal total
            nonlocal kept
            total += 1
            try:
                to_find = seq.metadata['id'].split('|')[1]
            except IndexError:
                continue

            if to_find in df.index:
                kept += 1
                yield seq

    skbio.io.write(generator(), format='fasta', into=output)

    print("For %d ids to keep, %d were found out of %d total samples"
          % (len(df), kept, total))


if __name__ == '__main__':
    main()

