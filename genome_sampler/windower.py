from collections import defaultdict

import pandas as pd
import biom


def _group_by_window(df, k, column='date', min_count=1):
    # pd.rolling and pd.Grouper don't fully account for this use case as best
    # as I can tell. Specifically, the goal here is a sliding window, such that
    # the indices within each window can be obtained rather than aggregated.
    index_name = df.index.name
    df = df.sort_values(column)
    df = df[~df[column].isnull()]

    if len(df) == 0:
        raise ValueError("DataFrame appears empty")

    # obtain the min / max observed dates
    min_dt = df[column].min()
    max_dt = df[column].max()

    # determine the last start day for the sliding window
    stop_dt = max_dt - pd.Timedelta(days=k - 1)

    # define a timedelta spanning the window
    k_timedelta = pd.Timedelta(days=k)

    # for each day
    windows = []
    for day in range((stop_dt - min_dt).days + 1):
        # determine the boundaries of the window
        window_left = (min_dt + pd.Timedelta(days=day))
        window_right = window_left + k_timedelta

        # gather all rows with a date of or greater than the left window
        # boundary, and "logical AND" them with the windows less than the right
        # boundary
        rows = (df[column] >= window_left) & (df[column] < window_right)

        if rows.sum() >= min_count:
            window = df[rows]
            windows.append({'starting_timepoint': window_left,
                            index_name: list(window.index)})

    return windows


def _create_window_df(window_metadata, outer_grouping, inner_grouping, dates):
    window_df = pd.DataFrame(window_metadata,
                             columns=['sample-id', outer_grouping,
                                      inner_grouping, dates])
    window_df['timepoint_gradient'] = window_df[dates].rank(method='dense')
    window_df[dates] = window_df[dates].apply(lambda x: x.strftime('%Y-%m-%d'))
    return window_df.set_index('sample-id')


def _make_id(outer_name, inner_name, timestamp):
    ts = timestamp.strftime('%Y-%m-%d')
    return f"{outer_name}:{inner_name}:{ts}".replace(' ', '_')


def _id_order(ids):
    order_map = {i: idx for idx, i in enumerate(sorted(ids))}
    ids_ordered = [k for k, _ in order_map.items()]
    return (order_map, ids_ordered)


def _create_table(rcv):
    # construct stable orderings for features and samples
    feature_order, feature_ids = _id_order({feature for feature, _ in rcv})
    sample_order, sample_ids = _id_order({sample for _, sample in rcv})

    # map into [[row_idx, col_idx, count]]
    rcv_flat = [[feature_order[i], sample_order[j], v]
                for (i, j), v in rcv.items()]

    return biom.Table(rcv_flat, feature_ids, sample_ids)


def sliding_window(metadata: pd.DataFrame,
                   dates: str,
                   outer_grouping: str,
                   inner_grouping: str,
                   window_size: int,
                   minimum_count: int,
                   feature_grouping: str = None) -> (biom.Table, pd.DataFrame):
    col_check = [dates, outer_grouping, inner_grouping]
    if feature_grouping is not None:
        col_check.append(feature_grouping)

    for k in col_check:
        if k not in metadata.columns:
            raise KeyError(f"{k} does not appear in the metadata")

    metadata = metadata.to_dataframe()
    metadata[dates] = pd.to_datetime(metadata[dates], errors='coerce')
    metadata = metadata[~metadata[dates].isnull()]

    if len(metadata) == 0:
        raise ValueError("No features left after removing invalid dates")

    if feature_grouping is not None:
        # index -> strain group
        feature_map = metadata[feature_grouping].to_dict()
    else:
        feature_map = {i: i for i in metadata.index}

    window_metadata = []
    rcv = defaultdict(int)
    index_name = metadata.index.name

    # e.g. for each country
    for outer_name, outer_group in metadata.groupby(outer_grouping):
        # e.g. for each city
        for inner_name, inner_group in outer_group.groupby(inner_grouping):
            for window in _group_by_window(inner_group, window_size, dates,
                                           minimum_count):
                window_id = _make_id(outer_name, inner_name,
                                     window['starting_timepoint'])

                # record the window metadata
                window_metadata.append([window_id, outer_name, inner_name,
                                        window['starting_timepoint']])

                # create row-column-value entries
                for feature in window[index_name]:
                    rcv[(feature_map[feature], window_id)] += 1

    window_df = _create_window_df(window_metadata, outer_grouping,
                                  inner_grouping, dates)
    table = _create_table(rcv)

    return (table, window_df)
