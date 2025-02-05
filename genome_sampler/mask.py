# ----------------------------------------------------------------------------
# Copyright (c) 2020-2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio
import pandas as pd
import numpy as np


def _filter_mask_by_level(mask, level):
    if level == 'mask':
        return mask[mask['FILTER'] == 'mask']
    else:
        return mask


def _create_position_map(aln, refseq_id):
    try:
        non_gaps = ~aln.loc[refseq_id].gaps()
    except KeyError:
        raise KeyError('Reference sequence %s is not present in alignment.' %
                       refseq_id)
    return non_gaps.nonzero()[0]


# so: https://stackoverflow.com/a/32681075/579416
def _rle(inarray):
    """ run length encoding. Partial credit to R rle function.
        Multi datatype arrays catered for including non Numpy
        returns: tuple (runlengths, startpositions, values) """
    ia = np.asarray(inarray)                # force numpy
    n = len(ia)
    if n == 0:
        return (None, None, None)
    else:
        y = np.array(ia[1:] != ia[:-1])      # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)    # must include last element posi
        z = np.diff(np.append(-1, i))        # run lengths
        p = np.cumsum(np.append(0, z))[:-1]  # positions
        return (z, p, ia[i])


def _find_terminal_gaps(aligned_seq):
    run_len, pos, vals = _rle(aligned_seq.gaps())

    output = np.full(len(aligned_seq), False, dtype=bool)
    if vals[0]:
        output[:run_len[0]] = True

    if vals[-1]:
        output[pos[-1]:] = True

    return output


def _create_terminal_gap_mask(aln, mask):
    result = np.full(aln.shape[1], True, dtype=bool)
    for chrom in mask['CHROM'].unique():
        result &= _find_terminal_gaps(aln.loc[chrom])
    return result


def _create_mask(aln, mask):
    result = np.full(aln.shape[1], False, dtype=bool)
    for chrom in mask['CHROM'].unique():
        chrom_idx = mask['CHROM'] == chrom
        mask_idx = mask['POS'].loc[chrom_idx] - 1

        position_lookup = _create_position_map(aln, chrom)
        try:
            aligned_mask_idx = position_lookup[mask_idx]
        except IndexError:
            raise IndexError('Reference sequence position out of range '
                             'for sequence %s' % chrom)

        result[aligned_mask_idx] = True
    return result


def _apply_mask(aln, mask):
    # True values indicate positions to remove in input mask
    return aln[:, ~mask]


def mask(alignment: skbio.TabularMSA, mask: pd.DataFrame,
         level: str = 'mask', mask_terminal_gaps: bool = True
         ) -> skbio.TabularMSA:
    alignment.reassign_index(minter='id')

    mask = _filter_mask_by_level(mask, level)
    mask_vector = _create_mask(alignment, mask)
    if mask_terminal_gaps:
        terminal_gap_vector = _create_terminal_gap_mask(alignment, mask)
        mask_vector |= terminal_gap_vector
    masked_alignment = _apply_mask(alignment, mask_vector)
    return masked_alignment
