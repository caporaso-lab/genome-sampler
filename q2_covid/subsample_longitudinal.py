import pandas as pd
import numpy as np

import qiime2

from q2_covid.common import IDSelection


def subsample_longitudinal(time: qiime2.CategoricalMetadataColumn,
                           start: str,
                           rate_n: int = 7,
                           rate_duration: int = 7) -> IDSelection:

    window_size = '%dD' % rate_duration
    filter_before = pd.Timestamp(start)

    dt_series = pd.to_datetime(time.to_series(), errors='coerce')
    df = pd.DataFrame({'ids': time.to_series().index}, index=dt_series)

    df = df.iloc[np.where(dt_series > filter_before)]
    grouped = df.groupby(pd.Grouper(freq=window_size), group_keys=False)
    filtered_df = grouped.apply(
        lambda g: g.sample(rate_n) if len(g) > rate_n else g)

    selection = pd.Series(False, index=df['ids'])
    selection[filtered_df['ids']] = True

    md = qiime2.Metadata(time.to_dataframe())
    return IDSelection(selection, md, 'subsample_longitudinal')
