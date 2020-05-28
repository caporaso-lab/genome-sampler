# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import numpy as np

import qiime2

from q2_types.feature_data import DNAIterator

def filter_seqs(sequences: DNAIterator,
                min_length: int=1,
                max_length: int=None,
                max_proportion_ambiguous: float=1.0) -> pd.Series:
    max_length = max_length or np.inf
    result = {}
    for sequence in sequences:
        sequence_length = len(sequence)
        ambiguous_proportion = sequence.degenerates().sum() / sequence_length
        too_long = sequence_length > max_length
        too_short = sequence_length < min_length
        too_ambiguous = ambiguous_proportion > max_proportion_ambiguous
        if not (too_short or too_long or too_ambiguous):
            result[sequence.metadata['id']] = sequence
    return pd.Series(result)
