import pandas as pd
import numpy as np

import qiime2

from q2_covid.common import IDSelection


def _sample_group(samples_per_interval, seed):
    state = np.random.RandomState(seed=seed)
    def _sampler(df):
        df = df.dropna(axis=0)
        if len(df) > samples_per_interval:
            return df.sample(samples_per_interval, random_state=state)
        else:
            return df
    
    return _sampler


def subsample_longitudinal(dates: qiime2.CategoricalMetadataColumn,
                           start_date: str = None,
                           samples_per_interval: int = 7,
                           days_per_interval: int = 7,
                           seed: int = None) -> IDSelection:

    window_size = '%dD' % days_per_interval

    dt_series = pd.to_datetime(dates.to_series(), errors='coerce')
    df = pd.DataFrame({'ids': dates.to_series().index}, index=dt_series)

    if start_date is not None:
        filter_before = pd.Timestamp(start_date)
        df = df.iloc[np.where(dt_series >= filter_before)]
        if filter_before not in df.index:
            # this will be stripped in _sample_group::_sampler
            # the purpose is to force Pandas to begin the window at this
            # time instead of the first observation (by making NaN the first
            # observation)
            df.loc[filter_before] = float('nan')

    grouped = df.groupby(pd.Grouper(freq=window_size, convention='start', 
                                    closed='left'), 
                         group_keys=False)
    filtered_df = grouped.apply(_sample_group(samples_per_interval, seed))

    df = df.dropna(axis=0)
    selection = pd.Series(False, index=dates.to_series().index)
    selection[filtered_df['ids']] = True

    md = qiime2.Metadata(dates.to_dataframe())
    return IDSelection(selection, md, 'subsample_longitudinal')
