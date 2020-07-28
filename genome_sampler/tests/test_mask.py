import pandas as pd
import pandas.testing as pdt
import skbio
import numpy as np
import numpy.testing as npt

from qiime2.plugin.testing import TestPluginBase

from genome_sampler.plugin_setup import VCFLikeMaskFormat
from genome_sampler.mask import (
    _filter_mask_by_level, _create_position_map, _find_terminal_gaps,
    _create_terminal_gap_mask, _create_mask, _apply_mask, mask)


class MaskTests(TestPluginBase):
    package = 'genome_sampler.tests'

    def setUp(self):
        super().setUp()
        _, self.mask1 = self.transform_format(VCFLikeMaskFormat,
                                              pd.DataFrame,
                                              filename='mask1.tsv')
        _, self.mask2 = self.transform_format(VCFLikeMaskFormat,
                                              pd.DataFrame,
                                              filename='mask2.tsv')
        _, self.mask3 = self.transform_format(VCFLikeMaskFormat,
                                              pd.DataFrame,
                                              filename='mask3.tsv')
        _, self.mask4 = self.transform_format(VCFLikeMaskFormat,
                                              pd.DataFrame,
                                              filename='mask4.tsv')
        _, self.mask5 = self.transform_format(VCFLikeMaskFormat,
                                              pd.DataFrame,
                                              filename='mask5.tsv')
        _, self.mask6 = self.transform_format(VCFLikeMaskFormat,
                                              pd.DataFrame,
                                              filename='mask6.tsv')
        _, self.mask7 = self.transform_format(VCFLikeMaskFormat,
                                              pd.DataFrame,
                                              filename='mask7.tsv')

        seqs = [
            skbio.DNA('ACGT'),
            skbio.DNA('AG-T'),
            skbio.DNA('-C-T')]
        self.msa1 = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])

        seqs = [
            skbio.DNA('TCNTGNNNGGTGCCA-CC--AAA--'),
            skbio.DNA('TCNTGCTCGGTGCCA-CC--AAAT-'),
            skbio.DNA('TCNTGCTCGGTACCA-CC--AAA--'),
            skbio.DNA('-CN-GCTCGGTGCCA-CCGGAAACT'),
            skbio.DNA('TCNTGCTCGGTGCCA-CC--AAATT'),
            skbio.DNA('--NTGCTCGGTGCCA-CC--AAAT-')]
        self.msa2 = skbio.TabularMSA(
                seqs, index=['s1', 's2', 's3', 'S_4', 'seq5.555', 's11'])

    def test_find_terminal_gaps_5_prime(self):
        aseq = skbio.DNA('--ACGTAGTCGA-AGCT----GATCG')
        exp = np.asarray(([True] * 2) + ([False] * (len(aseq) - 2)))
        obs = _find_terminal_gaps(aseq)
        npt.assert_array_equal(obs, exp)

    def test_find_terminal_gaps_none(self):
        aseq = skbio.DNA('A---ACGTAGTCGA-AGCT----GATCG')
        exp = np.asarray([False] * len(aseq))
        obs = _find_terminal_gaps(aseq)
        npt.assert_array_equal(obs, exp)

    def test_find_terminal_gaps_3_prime(self):
        aseq = skbio.DNA('A---ACGTAGTCGA-AGCT----')
        exp = np.asarray(([False] * (len(aseq) - 4)) + ([True] * 4))
        obs = _find_terminal_gaps(aseq)
        npt.assert_array_equal(obs, exp)

    def test_find_terminal_gaps_both_ends(self):
        aseq = skbio.DNA('---ACGTA-----GT-CGA-AGCT----')
        exp = np.asarray(([True] * 3) + ([False] * (len(aseq) - 7))
                         + ([True] * 4))
        obs = _find_terminal_gaps(aseq)
        npt.assert_array_equal(obs, exp)

    def test_find_terminal_gaps_all(self):
        aseq = skbio.DNA('-------')
        exp = np.asarray([True] * len(aseq))
        obs = _find_terminal_gaps(aseq)
        npt.assert_array_equal(obs, exp)

    def test_filter_mask_by_level_caution(self):
        obs = _filter_mask_by_level(self.mask1, "caution")
        pdt.assert_frame_equal(obs, self.mask1)
        self.assertTrue('caution' in list(obs['FILTER']))

    def test_filter_mask_by_level_mask(self):
        obs = _filter_mask_by_level(self.mask1, "mask")
        self.assertEqual(obs.shape[0], 16)
        self.assertFalse('caution' in list(obs['FILTER']))

    def test_create_position_map_no_gaps(self):
        obs = _create_position_map(self.msa1, 's1')
        exp = np.array([0, 1, 2, 3])
        npt.assert_array_equal(obs, exp)

    def test_create_position_map_some_gaps(self):
        obs = _create_position_map(self.msa1, 's2')
        exp = np.array([0, 1, 3])
        npt.assert_array_equal(obs, exp)

        obs = _create_position_map(self.msa1, 's3')
        exp = np.array([1, 3])
        npt.assert_array_equal(obs, exp)

    def test_create_position_map_all_gaps(self):
        seqs = [
            skbio.DNA('ACGT'),
            skbio.DNA('AG-T'),
            skbio.DNA('----')]
        msa = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])

        obs = _create_position_map(msa, 's3')
        exp = np.array([])
        npt.assert_array_equal(obs, exp)

    def test_create_position_map_error(self):
        with self.assertRaisesRegex(KeyError, 'Reference sequence s4 is not'):
            _create_position_map(self.msa1, 's4')

    def test_create_terminal_gap_mask_one_chrome(self):
        obs = _create_terminal_gap_mask(self.msa1, self.mask5)
        npt.assert_array_equal(obs, [True, False, False, False])

        seqs = [
            skbio.DNA('ACG-'),
            skbio.DNA('AG-T'),
            skbio.DNA('-C-T')]
        msa = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        obs = _create_terminal_gap_mask(msa, self.mask5)
        npt.assert_array_equal(obs, [True, False, False, False])

    def test_create_terminal_gap_mask_two_chrome(self):
        obs = _create_terminal_gap_mask(self.msa1, self.mask2)
        npt.assert_array_equal(obs, [False, False, False, False])

        seqs = [
            skbio.DNA('-CGT'),
            skbio.DNA('AG-T'),
            skbio.DNA('-C-T')]
        msa = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        obs = _create_terminal_gap_mask(msa, self.mask2)
        npt.assert_array_equal(obs, [True, False, False, False])

        seqs = [
            skbio.DNA('-CG-'),
            skbio.DNA('AG-T'),
            skbio.DNA('-C--')]
        msa = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        obs = _create_terminal_gap_mask(msa, self.mask2)
        npt.assert_array_equal(obs, [True, False, False, True])

    def test_create_terminal_gap_mask_none(self):
        seqs = [
            skbio.DNA('AC--'),
            skbio.DNA('AG-T'),
            skbio.DNA('--TT')]
        msa = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        obs = _create_terminal_gap_mask(msa, self.mask2)
        npt.assert_array_equal(obs, [False, False, False, False])

        seqs = [
            skbio.DNA('AC--'),
            skbio.DNA('AG-T'),
            skbio.DNA('AATT')]
        msa = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        obs = _create_terminal_gap_mask(msa, self.mask5)
        npt.assert_array_equal(obs, [False, False, False, False])

    def test_create_terminal_gap_mask_all(self):
        seqs = [
            skbio.DNA('----'),
            skbio.DNA('AGAT'),
            skbio.DNA('----')]
        msa = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        obs = _create_terminal_gap_mask(msa, self.mask2)
        npt.assert_array_equal(obs, [True, True, True, True])

        seqs = [
            skbio.DNA('ACG-'),
            skbio.DNA('AG-T'),
            skbio.DNA('----')]
        msa = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        obs = _create_terminal_gap_mask(msa, self.mask5)
        npt.assert_array_equal(obs, [True, True, True, True])

    def test_create_mask_one_chrome(self):
        obs = _create_mask(self.msa1, self.mask5)
        npt.assert_array_equal(obs, [False, False, False, True])

    def test_create_mask_two_chrome(self):
        obs = _create_mask(self.msa1, self.mask2)
        npt.assert_array_equal(obs, [False, True, False, True])

    def test_create_mask_mask_all(self):
        obs = _create_mask(self.msa1, self.mask6)
        npt.assert_array_equal(obs, [True, True, True, True])

    def test_create_mask_mask_none(self):
        obs = _create_mask(self.msa1, self.mask7)
        npt.assert_array_equal(obs, [False, False, False, False])

    def test_create_mask_error(self):
        with self.assertRaisesRegex(IndexError, 'sequence position out .* s1'):
            _create_mask(self.msa1, self.mask3)

    def test_apply_mask_mask_none(self):
        obs = _apply_mask(self.msa1, np.array([False, False, False, False]))
        self.assertEqual(obs, self.msa1)

    def test_apply_mask_mask_some(self):
        obs = _apply_mask(self.msa1, np.array([False, True, True, True]))
        seqs = [
            skbio.DNA('A'),
            skbio.DNA('A'),
            skbio.DNA('-')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

        obs = _apply_mask(self.msa1, np.array([False, True, True, False]))
        seqs = [
            skbio.DNA('AT'),
            skbio.DNA('AT'),
            skbio.DNA('-T')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

        obs = _apply_mask(self.msa1, np.array([False, True, False, False]))
        seqs = [
            skbio.DNA('AGT'),
            skbio.DNA('A-T'),
            skbio.DNA('--T')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

    def test_apply_mask_mask_all(self):
        obs = _apply_mask(self.msa1, np.array([True, True, True, True]))
        seqs = [
            skbio.DNA(''),
            skbio.DNA(''),
            skbio.DNA('')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

    def test_mask2_wo_terminal_gap_mask(self):
        obs = mask(self.msa1, self.mask2, "mask", False)
        seqs = [
            skbio.DNA('ACG'),
            skbio.DNA('AG-'),
            skbio.DNA('-C-')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

        obs = mask(self.msa1, self.mask2, "caution", False)
        seqs = [
            skbio.DNA('AG'),
            skbio.DNA('A-'),
            skbio.DNA('--')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

    def test_mask2_w_terminal_gap_mask(self):
        obs = mask(self.msa1, self.mask2, "mask", True)
        seqs = [
            skbio.DNA('ACG'),
            skbio.DNA('AG-'),
            skbio.DNA('-C-')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

        obs = mask(self.msa1, self.mask2, "caution", True)
        seqs = [
            skbio.DNA('AG'),
            skbio.DNA('A-'),
            skbio.DNA('--')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

        obs = mask(self.msa1, self.mask5, "caution", True)
        seqs = [
            skbio.DNA('CG'),
            skbio.DNA('G-'),
            skbio.DNA('C-')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

        seqs = [
            skbio.DNA('ACGTA'),
            skbio.DNA('AG-TA'),
            skbio.DNA('-C-T-')]
        msa = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])

        obs = mask(msa, self.mask5, "caution", True)
        seqs = [
            skbio.DNA('CG'),
            skbio.DNA('G-'),
            skbio.DNA('C-')]
        exp = skbio.TabularMSA(seqs, index=['s1', 's2', 's3'])
        self.assertEqual(obs, exp)

    def test_mask4_wo_terminal_gap_mask(self):
        obs = mask(self.msa2, self.mask4, "mask", False)
        seqs = [
            skbio.DNA('TCNNNGGTGCCA-CC--A-'),
            skbio.DNA('TCCTCGGTGCCA-CC--A-'),
            skbio.DNA('TCCTCGGTACCA-CC--A-'),
            skbio.DNA('-CCTCGGTGCCA-CCGGAT'),
            skbio.DNA('TCCTCGGTGCCA-CC--AT'),
            skbio.DNA('--CTCGGTGCCA-CC--A-')]
        exp = skbio.TabularMSA(
                seqs, index=['s1', 's2', 's3', 'S_4', 'seq5.555', 's11'])
        self.assertEqual(obs, exp)

        obs = mask(self.msa2, self.mask4, "caution", False)
        seqs = [
            skbio.DNA('TCNNNGGCCA-CC--A-'),
            skbio.DNA('TCCTCGGCCA-CC--A-'),
            skbio.DNA('TCCTCGGCCA-CC--A-'),
            skbio.DNA('-CCTCGGCCA-CCGGAT'),
            skbio.DNA('TCCTCGGCCA-CC--AT'),
            skbio.DNA('--CTCGGCCA-CC--A-')]
        exp = skbio.TabularMSA(
                seqs, index=['s1', 's2', 's3', 'S_4', 'seq5.555', 's11'])
        self.assertEqual(obs, exp)

    def test_mask4_w_terminal_gap_mask(self):
        obs = mask(self.msa2, self.mask4, "mask", True)
        seqs = [
            skbio.DNA('NNNGGTGCCA-CC--A'),
            skbio.DNA('CTCGGTGCCA-CC--A'),
            skbio.DNA('CTCGGTACCA-CC--A'),
            skbio.DNA('CTCGGTGCCA-CCGGA'),
            skbio.DNA('CTCGGTGCCA-CC--A'),
            skbio.DNA('CTCGGTGCCA-CC--A')]
        exp = skbio.TabularMSA(
                seqs, index=['s1', 's2', 's3', 'S_4', 'seq5.555', 's11'])
        self.assertEqual(obs, exp)

        obs = mask(self.msa2, self.mask4, "caution", True)
        seqs = [
            skbio.DNA('NNNGGCCA-CC--A'),
            skbio.DNA('CTCGGCCA-CC--A'),
            skbio.DNA('CTCGGCCA-CC--A'),
            skbio.DNA('CTCGGCCA-CCGGA'),
            skbio.DNA('CTCGGCCA-CC--A'),
            skbio.DNA('CTCGGCCA-CC--A')]
        exp = skbio.TabularMSA(
                seqs, index=['s1', 's2', 's3', 'S_4', 'seq5.555', 's11'])
        self.assertEqual(obs, exp)
