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
    else:
        md_df = pd.DataFrame({}, index=seqs.index)

    selected = md_df[columns]
    rename = pd.Series([delimiter.join(row) for row in selected.itertuples()],
                       index=selected.index)

    seqs.index = seqs.index.map(rename)

    return seqs
