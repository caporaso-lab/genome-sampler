#!/usr/bin/env python
import datetime

import click
import pandas as pd
import numpy as np
import skbio

def tgen_id_to_datetime(s):
     # this is very specific to one data set, and probably
     # shouldn't be used widely
     fields = s.split('|')
     date_field = fields[-1]
     if '/' in date_field:
         year, month, day = date_field.split('/')
     elif '.' in date_field:
         month, day, year = date_field.split('.')
     else:
         raise ValueError("Can't parse date from id: %s" % s)
     return '-'.join([year, month, day])


def fix_ua(s):
    if s == 'UA|SARS-COV-2|171|PIMA':
        return '2020-03-19'

    raise Exception("Unknown date for UA sequence")

@click.command()
@click.option('--remove', help='fasta IDs to remove',
              type=click.File(mode='r'), required=False)
@click.option('--fix-date', help='fix tgen/ua dates to match gisaid',
              type=bool, required=False)
@click.argument('input', type=click.Path())
@click.argument('output', type=click.Path())
def main(input, output, remove, fix_date):

    if remove:
        remove = set(l.strip() for l in remove)
    else:
        remove = set()


    total = 0
    duplicates = 0
    removed = 0
    renamed = 0

    seen = set()

    def generator():
        for seq in skbio.io.read(input, format='fasta', constructor=skbio.DNA):
            nonlocal total
            nonlocal duplicates
            nonlocal removed
            nonlocal renamed
            total += 1

            if seq.metadata['id'] in seen:
                duplicates += 1
                continue
            seen.add(seq.metadata['id'])

            if seq.metadata['id'] in remove:
                removed += 1
                continue

            if fix_date:
                change = False
                try:
                    start, *mid, date = seq.metadata['id'].split('|')
                except ValueError:
                    start = ""

                if 'TGEN' in start:
                    date = tgen_id_to_datetime(seq.metadata['id'])
                    change = True
                elif 'UA' in start:
                    date = fix_ua(seq.metadata['id'])
                    change = True

                if change:
                    seq.metadata['id'] = '|'.join([start, *mid, date])
                    renamed += 1

            yield seq

    skbio.io.write(generator(), format='fasta', into=output)

    print("%d sequences were duplicated, %d were removed, and %d were renamed"
          " out of %d total." % (duplicates, removed, renamed, total))


if __name__ == '__main__':
    main()

