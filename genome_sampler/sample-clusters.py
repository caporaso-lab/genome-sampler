#!/usr/bin/env python
# flake8: noqa

import click
import pandas as pd
import numpy as np
import skbio


def cigar_parse(cigar):
    if cigar == '=':
        yield (0, '=')
        return
    try:
        count = ""
        for char in cigar:
            try:
                int(char)
            except:
                if not count:
                    count = 0
                yield (int(count), char)
                count = ""
            else:
                count += char
    except ValueError:
        print(cigar)
        raise


def slice_query_target(query, target, cigar_tuple):
    """IMPORTANT: this slice out indels, leaving only match/mismatch
       It appears that the middle is almost always ~29k M, so not penalizing
       a middle indel as part of the hamming probably won't matter too much
    """
    new_query = ""
    new_target = ""
    for count, op in cigar_tuple:
        if count > min(len(query), len(target)):
            count = min(len(query), len(target))
        if op == 'M':
            new_query = new_query + query[:count]
            new_target = new_target + target[:count]
            target = target[count:]
            query = query[count:]
        elif op == 'I':
            target = target[count:]
        elif op == 'D':
            query = query[count:]
        elif op == '=':
            new_query = query
            new_target = target

    return skbio.DNA(new_query), skbio.DNA(new_target)


def snp_count(record, seqs):
    cigar_tuple = cigar_parse(record.cigar)
    q, t = slice_query_target(seqs[record.query], seqs[record.target],
                              cigar_tuple)
    try:
        mask = q.definites() & t.definites()
    except ValueError:
        print(record.cigar)
        print(repr(q))
        print(repr(t))
    return (q[mask].values != t[mask].values).sum()


def location(gisaid, df):
    return ':'.join(df.loc[gisaid][['country', 'division']].dropna())


def get_id(header):
    if type(header) is str and 'EPI_ISL' in header:
        return header.split('|')[1]
    return header


def clusters_from_uc(uc):
    hits = uc[uc['type'] == 'H']
    clusters = {}
    for r in hits.itertuples():
        if r.target not in clusters:
            clusters[r.target] = []
        clusters[r.target].append(r)
    return clusters


def make_snp_tables(clusters, seqs, df):
    snps = {}
    for target, hits in clusters.items():
        idx = [get_id(h.query) for h in hits]
        table = pd.DataFrame(
            {'snp': [snp_count(h, seqs) for h in hits],
             'location': [location(get_id(h.query), df) for h in hits],
             'date': pd.to_datetime(df['date'][df.index.isin(idx)],
                                    errors='coerce')},
            index=idx)
        table = table.sort_values(by=['snp', 'date'])
        snps[get_id(target)] = table
    return snps


def sample_table(n, table):
    if n > len(table):
        n = len(table)
    loc = table['location']
    # convert locations into number of obs, take the inverse, then scale by
    # number of unique locations so it all sums to 1. Infrequent location obs
    # will have a high weight, frequent location obs will have an individually
    # lower weight, but collectively all locations have the same weight when
    # all obs are summed within a location class.
    weights = ((1 / loc.map(loc.value_counts())) / len(loc.unique()))
    sample = table.sample(n, weights=weights)
    if any(sample['snp'] == 1):
        return sample  # already have something like a MRCA
    elif not any(table['snp'] == 1):
        return sample  # no good candidates
    else:
        # patch in a better one
        extra = sample_table(1, table[table['snp'] == 1])
        sample = sample[:(n-1)]
        return pd.concat([sample, extra])


@click.command()
@click.option('--n', help='number of additional "context" samples to include'
                          ' per cluster', required=True, type=int)
@click.option('--uc', help='UC cluster map',
              type=click.Path(), required=True)
@click.option('--target', help='target sequences as fasta',
              type=click.Path(), required=True)
@click.option('--query', help='query sequences as fasta',
              type=click.Path(), required=True)
@click.option('--tsv', help='nextstrain metadata tsv',
              type=click.Path(), required=True)
@click.argument('output', type=click.File(mode='w'))
def main(n, uc, target, query, tsv, output):
    df = pd.read_csv(tsv, sep='\t')
    df.index = df['gisaid_epi_isl']

    targets = pd.Series(skbio.io.read(target, format='fasta',
                                      constructor=skbio.DNA))
    targets = targets.rename(lambda x: targets[x].metadata['id'] + targets[x].metadata['description'])
    queries = pd.Series(skbio.io.read(query, format='fasta',
                                      constructor=skbio.DNA))
    queries = queries.rename(lambda x: queries[x].metadata['id'] + queries[x].metadata['description'])
    qidx = queries.index.to_series().apply(lambda x: get_id(x) in df.index)
    queries = queries[qidx]
    print("Metadata is missing %d samples. Any missing samples will be filtered."
          % (~qidx).sum())

    seqs = queries.combine_first(targets).apply(str)

    uc = pd.read_csv(
        uc, sep='\t', na_values='*',
        names=['type', 'cluster_id', 'length', 'perc_id', 'strand', 'BLANK1',
               'BLANK2', 'cigar', 'query', 'target'])
    avail = (uc['target'].apply(lambda x: x in seqs.index)
             & uc['query'].apply(lambda x: x in seqs.index))
    uc = uc[avail]

    clusters = clusters_from_uc(uc)
    snps = make_snp_tables(clusters, seqs, df)

    context_df = pd.DataFrame(
        columns=['gisaid_epi_isl', 'snp', 'location', 'date', 'centroid-id'])
    for centroid, table in snps.items():
        table = sample_table(n, table)

        table['gisaid_epi_isl'] = table.index
        table['centroid-id'] = centroid

        context_df = context_df.append(table, sort=False)

    with output.open() as fh:
        fh.write(context_df.to_csv(sep='\t', index=False))


if __name__ == '__main__':
    main()
