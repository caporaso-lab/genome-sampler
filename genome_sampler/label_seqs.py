# ----------------------------------------------------------------------------
# Copyright (c) 2020-2025, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

import qiime2


def label_seqs(seqs: pd.Series, delimiter: str,
               metadata: qiime2.Metadata = None, columns: str = None,
               missing_value: str = 'missing') \
                   -> pd.Series:
    if columns is not None and metadata is None \
            or metadata is not None and columns is None:
        raise ValueError('Columns and metadata must be passed or not passed '
                         'together.')

    if delimiter in missing_value:
        raise ValueError(f'The provided delimiter ({repr(delimiter)}) cannot '
                         'be contained in the missing value placeholder '
                         f'({repr(missing_value)}).')

    # This is necessary because QIIME 2 will not accept an empty list as an
    # argument of type List[str]
    if columns is None:
        columns = []

    # Make sure we have strings at this point not skbio DNA objects because we
    # experienced a bizarre segmentation fault while using DNA objects
    seqs = seqs.apply(str)
    seqs.index = seqs.index.map(lambda x: x.split(delimiter)[0])

    if metadata is not None:
        md_df = metadata.to_dataframe()

        for column in columns:
            if column not in md_df.columns:
                raise ValueError(f'The column {repr(column)} is not present '
                                 'in the metadata')

        missing_ids = seqs.index.difference(md_df.index)
        if len(missing_ids):
            difference = \
                ' '.join(repr(value) for value in missing_ids.values[0:10])
            additional_missing = len(missing_ids.values[10:])

            error_message = ('The following ids are present in the sequences '
                             f'but not the metadata {difference}')

            if additional_missing > 0:
                error_message += (f' ({additional_missing} additional ids are'
                                  ' missing from metadata but omitted from'
                                  ' this list)')

            raise ValueError(error_message)
    else:
        md_df = pd.DataFrame({}, index=seqs.index)

    selected = md_df[columns]
    selected = selected.fillna(missing_value)
    rename = pd.Series([delimiter.join(row) for row in selected.itertuples()],
                       index=selected.index)
    seqs.index = seqs.index.map(rename)

    return seqs
