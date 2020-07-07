import pandas as pd


import qiime2


def label_seqs(seqs: pd.Series, metadata: qiime2.Metadata, columns: str,
               delimiter: str) -> pd.Series:
    seqs = seqs.apply(str)
    md_df = metadata.to_dataframe()

    selected = md_df[columns]
    rename = pd.Series([delimiter.join(row) for row in selected.itertuples()],
                       index=selected.index)

    seqs.index = seqs.index.map(rename)

    return seqs
