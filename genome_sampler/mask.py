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

def _refseq_to_aln_positions(aln, mask):
    # positions in mask will be one-indexed
    result = []
    position_maps = {}
    for e in mask.itertuples():
        refseq_id = e.CHROM
        refseq_position = e.POS
        if refseq_id not in position_maps:
            position_maps[refseq_id] = _create_position_map(aln, refseq_id)
        position_map = position_maps[refseq_id]
        try:
            result.append(position_map[refseq_position - 1])
        except IndexError:
            raise IndexError('Reference sequence position %d is out of range '
                             'for sequence %s' % (refseq_position, refseq_id))
    return result

def _compute_boolean_mask(aln, aln_positions_to_remove):
    result = np.ones(aln.shape[1])
    np.put(result, aln_positions_to_remove, [0])
    return result

def _apply_mask(aln, mask):
    return aln[:, mask]

def mask(alignment: skbio.TabularMSA, mask: pd.DataFrame,
         level: str='mask') -> skbio.TabularMSA:
    mask = _filter_mask_by_level(mask, level)
    aln_positions_to_remove = \
        _refseq_to_aln_positions(aln, refseq_id, refseq_positions)
    boolean_mask = _compute_mask(aln, aln_positions_to_remove)
    masked_alignment = _apply_mask(aln, boolean_mask)
    return masked_alignment
