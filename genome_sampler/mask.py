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


def _refseq_to_aln_positions(aln, mask, mask_gapped_ends):
    # positions in mask will be one-indexed
    # positions in result will be zero-indexed
    result = []
    position_maps = {}
    for e in mask.itertuples():
        refseq_id = e.CHROM
        refseq_position = e.POS

        if refseq_id not in position_maps:
            position_maps[refseq_id] = _create_position_map(aln, refseq_id)
        position_map = position_maps[refseq_id]

        try:
            aln_position = position_map[refseq_position - 1]
        except IndexError:
            raise IndexError('Reference sequence position %d is out of range '
                             'for sequence %s' % (refseq_position, refseq_id))
        result.append(aln_position)

        if mask_gapped_ends:
            # if the aligned reference sequences contains only gap characters
            # up to the position that is being masked, mask all alignment
            # positions preceeding this position (since they're outside the
            # range of the reference sequence)
            if aln.loc[refseq_id, :aln_position].gaps().all():
                result.extend(range(aln_position))
            # if the aligned reference sequences contains only gap characters
            # after the position that is being masked, mask all alignment
            # positions following this position (since they're outside the
            # range of the reference sequence)
            if aln.loc[refseq_id, aln_position + 1:].gaps().all():
                result.extend(range(aln_position+1, aln.shape[1]))

        result.sort()

    return result


def _compute_boolean_mask(aln, aln_positions_to_remove):
    result = np.full(aln.shape[1], True, dtype=bool)
    np.put(result, aln_positions_to_remove, [False])
    return result


def _apply_mask(aln, mask):
    return aln[:, mask]


def mask(alignment: skbio.TabularMSA, mask: pd.DataFrame,
         level: str = 'mask', mask_gapped_ends: bool = True
         ) -> skbio.TabularMSA:
    mask = _filter_mask_by_level(mask, level)
    aln_positions_to_remove = \
        _refseq_to_aln_positions(alignment, mask, mask_gapped_ends)
    boolean_mask = _compute_boolean_mask(
            alignment, aln_positions_to_remove)
    masked_alignment = _apply_mask(alignment, boolean_mask)
    return masked_alignment
