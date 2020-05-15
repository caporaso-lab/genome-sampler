#!/usr/bin/env python

import click
import pandas as pd
import numpy as np

@click.command()
@click.option('--rate', help='set sample rate, format as: <n>,<days>',
              required=True)
@click.option('--hist', help='filepath for histogram',
              type=click.Path())
@click.argument('input', type=click.File(mode='r'))
@click.argument('output', type=click.File(mode='w'))
def main(input, output, rate, hist):

    filter_before = pd.Timestamp("2019-12-1")
    n, w = rate.split(',')
    window_size = '%dD' % int(w)
    n_per_window = int(n)

    df = pd.read_csv(input, sep='\t')
    orig = len(df)
    time = pd.to_datetime(df['date'], errors='coerce')
    df.index = time
    df = df.iloc[np.where(time > pd.Timestamp("2019-12-1"))]
    gb = df.groupby(pd.Grouper(freq=window_size),  group_keys=False)
    df = gb.apply(
        lambda x: x.sample(n_per_window) if len(x) > n_per_window else x)

    print('Subsampled %d down to %d samples.' % (orig, len(df)))

    if hist is not None:
        df.index.to_series().hist(bins=200, figsize=(10, 4)).figure.savefig(str(hist))

    with output.open() as fh:
        fh.write(df.to_csv(sep='\t', index=False))


if __name__ == '__main__':
    main()
