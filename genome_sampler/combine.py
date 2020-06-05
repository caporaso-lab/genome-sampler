import operator

import pandas as pd
import numpy as np

import qiime2

from genome_sampler.common import IDSelection


def _combine_df_error_if_not_equal(a, b):
    if set(a.index) != set(b.index):
        raise ValueError("Metadata id sets are not equal. Can't combine.")

    # if column isn't in a, return b
    if a.dropna().empty:
            return b

    # if column isn't in b, return a
    if b.dropna().empty:
            return a

    # if column is in both but aren't equal, error
    if not (a == b).all():
        raise ValueError("Can't combine inconsistent metadata.")

    # columns are equal, so it doesn't matter which we return
    return a


def combine_selections(selections: IDSelection):
    output_label = "combined_selections"
    if len(selections) == 1:
        return IDSelection(selections[0].inclusion,
                           selections[0].metadata,
                           label=output_label)

    inclusion = selections[0].inclusion
    inclusion_ids = set(inclusion.index)
    metadata = selections[0].metadata.to_dataframe()
    metadata_ids = set(metadata.index)
    for e in selections[1:]:
        if inclusion_ids != set(e.inclusion.index):
            raise ValueError("Inclusion id sets are not equal. Can't combine.")
        inclusion = inclusion.combine(e.inclusion, operator.or_)

        df = e.metadata.to_dataframe()
        if metadata_ids != set(df.index):
            raise ValueError("Metadata id sets are not equal. Can't combine.")
        metadata = metadata.combine(df, _combine_df_error_if_not_equal)

    return IDSelection(inclusion, qiime2.Metadata(metadata), output_label)
