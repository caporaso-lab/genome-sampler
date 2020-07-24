import pandas as pd
import pandas.testing as pdt
import skbio
import numpy as np
import numpy.testing as npt

from qiime2.plugin.testing import TestPluginBase

from genome_sampler.plugin_setup import BrokenVCFFormat
from genome_sampler.mask import (
    _filter_mask_by_level, _create_position_map, _refseq_to_aln_positions,
    _compute_boolean_mask, _apply_mask, mask)


class MaskTests(TestPluginBase):
    package = 'genome_sampler.tests'

    def setUp(self):
        super().setUp()
        _, self.mask1 = self.transform_format(BrokenVCFFormat,
                                              pd.DataFrame,
                                              filename='mask1.tsv')
        _, self.mask2 = self.transform_format(BrokenVCFFormat,
                                              pd.DataFrame,
                                              filename='mask2.tsv')
        _, self.mask3 = self.transform_format(BrokenVCFFormat,
                                              pd.DataFrame,
                                              filename='mask3.tsv')

        seqs = [
            skbio.DNA('ACGT'),
            skbio.DNA('AG-T'),
            skbio.DNA('-C-T')]
        self.msa1 = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])

    def test_filter_mask_by_level_caution(self):
        obs = _filter_mask_by_level(self.mask1, "caution")
        pdt.assert_frame_equal(obs, self.mask1)
        self.assertTrue('caution' in list(obs['FILTER']))

    def test_filter_mask_by_level_mask(self):
        obs = _filter_mask_by_level(self.mask1, "mask")
        self.assertEqual(obs.shape[0], 16)
        self.assertFalse('caution' in list(obs['FILTER']))

    def test_create_position_map(self):
        obs = _create_position_map(self.msa1, 's1')
        exp = np.array([0, 1, 2, 3])
        npt.assert_array_equal(obs, exp)

        obs = _create_position_map(self.msa1, 's2')
        exp = np.array([0, 1, 3])
        npt.assert_array_equal(obs, exp)

        obs = _create_position_map(self.msa1, 's3')
        exp = np.array([1, 3])
        npt.assert_array_equal(obs, exp)

    def test_create_position_map_error(self):
        with self.assertRaisesRegex(KeyError, 'Reference sequence s4 is not'):
            _create_position_map(self.msa1, 's4')

    def test_refseq_to_aln_positions(self):
        obs = _refseq_to_aln_positions(self.msa1, self.mask2)
        self.assertEqual(obs, [1, 3])

    def test_refseq_to_aln_positions_error(self):
        with self.assertRaisesRegex(IndexError, 'sequence position 42.*e s1'):
            _refseq_to_aln_positions(self.msa1, self.mask3)

    def test_compute_boolean_mask(self):
        obs = _compute_boolean_mask(self.msa1, [])
        npt.assert_array_equal(obs, np.array([True, True, True, True]))

        obs = _compute_boolean_mask(self.msa1, [1, 3])
        npt.assert_array_equal(obs, np.array([True, False, True, False]))

        obs = _compute_boolean_mask(self.msa1, [0, 1, 2, 3])
        npt.assert_array_equal(obs, np.array([False, False, False, False]))

    def test_apply_mask(self):
        obs = _apply_mask(self.msa1, [True, True, True, True])
        self.assertEqual(obs, self.msa1)

        obs = _apply_mask(self.msa1, [True, False, False, False])
        seqs = [
            skbio.DNA('A'),
            skbio.DNA('A'),
            skbio.DNA('-')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

        obs = _apply_mask(self.msa1, [True, False, False, True])
        seqs = [
            skbio.DNA('AT'),
            skbio.DNA('AT'),
            skbio.DNA('-T')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

        obs = _apply_mask(self.msa1, [True, False, True, True])
        seqs = [
            skbio.DNA('AGT'),
            skbio.DNA('A-T'),
            skbio.DNA('--T')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

        obs = _apply_mask(self.msa1, [False, False, False, False])
        seqs = [
            skbio.DNA(''),
            skbio.DNA(''),
            skbio.DNA('')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)
