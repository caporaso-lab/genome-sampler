import pandas as pd

import qiime2


def label_seqs(seqs: pd.Series, delimiter: str,
               metadata: qiime2.Metadata = None, columns: str = None) \
                   -> pd.Series:
    if columns is not None and metadata is None \
            or metadata is not None and columns is None:
        raise ValueError('Columns and metadata must be passed or not passed '
                         'together.')

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
        if missing_ids.values.size != 0:
            difference = ' '.join(repr(value) for value in missing_ids.values)
            raise ValueError('The following ids are present in the sequences '
                             f'but not the metadata {repr(difference)}')
    else:
        md_df = pd.DataFrame({}, index=seqs.index)

    selected = md_df[columns]
    rename = pd.Series([delimiter.join(row) for row in selected.itertuples()],
                       index=selected.index)
    seqs.index = seqs.index.map(rename)

    return seqs
