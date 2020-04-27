#!/usr/bin/env python

import click
import pandas as pd
import numpy as np
import skbio

@click.command()
@click.option('--tsv', help='filepath to sequence metadata',
              type=click.File(mode='r'), required=True)
@click.option('--max-conservation', 
              help=('only show positions with conservation less'
                    ' than (exclusive) this value (default=1.0'),
              type=float, default=1.0)
@click.argument('alignment', type=click.Path())
@click.argument('output', type=click.Path())
@click.argument('output-gap-degen', type=click.Path())
def main(alignment, output, output_gap_degen, tsv, max_conservation):
    group_column = 'clade'
    df = pd.read_csv(tsv, sep='\t', header=0)
    try:
        df = df.set_index('id')
    except KeyError:
        raise KeyError('First column in tsv must have header "id", but header is: "%s"' % df.columns[0])
    
    try:
        clades = df.groupby(by=[group_column]).groups
    except KeyError:
        raise KeyError('Column named "%s" does not exist in metadata. Existing columns are: %s' 
                       % (group_column, ', '.join(df.columns)))

    msa = skbio.TabularMSA.read(alignment, format='fasta', constructor=skbio.DNA)
    msa.reassign_index(minter='id')

    # indices must be identical (we can consider relaxing this if necessary)
    metadata_index = set(df.index)
    msa_index = set(msa.index)
    if metadata_index != msa_index:
        if len(metadata_index - msa_index) != 0:
            raise KeyError("Some ids in metadata are not in the alignment: %s" % ", ".join(metadata_index - msa_index))

        if len(msa_index - metadata_index) != 0:
            raise KeyError("Some ids in alignment are not in the metadata: %s" % ", ".join(msa_index - metadata_index))

    conservation = msa.conservation(degenerate_mode='nan', gap_mode='nan')

    clade_msas = {
        clade_id: msa.loc[list(seq_ids)] for clade_id, seq_ids in clades.items()}

    results = []
    results_gap_or_degen = []
    for position, cons in enumerate(conservation):
        if np.isnan(cons):
            clade_columns = _get_clade_columns(position, clade_msas)
            results_gap_or_degen.append([position] + clade_columns)
        elif cons < max_conservation:
            clade_columns = _get_clade_columns(position, clade_msas)
            results.append([position] + clade_columns)
        else:
            continue
    
    _write_output(output, results, clade_msas)
    _write_output(output_gap_degen, results_gap_or_degen, clade_msas)

def _get_clade_columns(position, clade_msas):
    clade_columns = [
        str(clade_msa.iloc[:, position]) for clade_msa in clade_msas.values()
    ]
    return clade_columns

def _write_output(output_fp, results, clade_msas):
    results = pd.DataFrame(results, 
                           columns=['position (0-indexed)'] + list(clade_msas.keys()))
    results = results.set_index('position (0-indexed)')
    results.to_csv(output_fp, sep='\t')

if __name__ == '__main__':
    main()

