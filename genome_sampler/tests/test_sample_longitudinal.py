# ----------------------------------------------------------------------------
# Copyright (c) 2020-2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import numpy as np

import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat

from genome_sampler.sample_longitudinal import sample_longitudinal


class TestSubsampleLongitudinal(TestPluginBase):
    package = 'genome_sampler.tests'

    _N_TEST_ITERATIONS = 50

    def setUp(self):
        super().setUp()

        s1 = pd.Series(['2019-12-31', '2020-01-09', '2020-01-10',
                        '2019-11-01', '2020-01-11', '2020-02-21',
                        '2020-02-21', '2020-02-21', '2020-03-15'],
                       index=[chr(x) for x in range(65, 74)])
        s1.index.name = 'id'
        s1.name = 'date-md'
        self.md1 = qiime2.CategoricalMetadataColumn(s1)

        s2 = pd.Series(['2020-01-02', '2019-11-01', '2020-02-21',
                        '2020-02-21', '2020-02-21', '2020-03-15',
                        '2020-01-03', '2020-01-04', '2020-01-05',
                        '2020-01-06', '2020-01-07', '2020-01-08',
                        '2020-01-09', '2020-01-10', '2020-01-11',
                        '2020-01-12', '2020-01-13', '2020-01-14',
                        '2020-01-15', '2020-01-16', '2020-01-17'],
                       index=[chr(x) for x in range(65, 86)])
        s2.index.name = 'id'
        s2.name = 'date-md'
        self.md2 = qiime2.CategoricalMetadataColumn(s2)

    def test_default(self):
        sel = sample_longitudinal(self.md1)

        self.assertEqual(sel.inclusion.sum(), 9)
        self.assertEqual(sel.metadata.get_column('date-md'), self.md1)
        self.assertEqual(sel.label, 'sample_longitudinal')

    def test_start_date_in_data(self):
        sel = sample_longitudinal(self.md1,
                                  start_date='2019-12-31')

        self.assertEqual(sel.inclusion.sum(), 8)
        self.assertEqual(sel.metadata.get_column('date-md'), self.md1)
        self.assertEqual(sel.label, 'sample_longitudinal')
        self.assertFalse(np.nan in list(sel.inclusion.index))

    def test_start_date_not_in_data(self):
        sel = sample_longitudinal(self.md1,
                                  start_date='2019-12-30')

        self.assertEqual(sel.inclusion.sum(), 8)
        self.assertEqual(sel.metadata.get_column('date-md'), self.md1)
        self.assertEqual(sel.label, 'sample_longitudinal')
        self.assertFalse(np.nan in list(sel.inclusion.index))

    def test_one_sample_per_interval(self):
        sel = sample_longitudinal(self.md1,
                                  samples_per_interval=1)

        self.assertEqual(sel.inclusion.sum(), 6)
        self.assertEqual(sel.metadata.get_column('date-md'), self.md1)
        self.assertEqual(sel.label, 'sample_longitudinal')

    def test_two_sample_per_interval(self):
        sel = sample_longitudinal(self.md1,
                                  samples_per_interval=2)

        self.assertEqual(sel.inclusion.sum(), 8)
        self.assertEqual(sel.metadata.get_column('date-md'), self.md1)
        self.assertEqual(sel.label, 'sample_longitudinal')

    def test_interval_bounds1(self):
        for _ in range(self._N_TEST_ITERATIONS):
            sel = sample_longitudinal(self.md2,
                                      samples_per_interval=1,
                                      start_date='2019-12-26')

            exp_int1_dates = ['2020-01-02', '2020-01-03', '2020-01-04',
                              '2020-01-05', '2020-01-06', '2020-01-07',
                              '2020-01-08']
            exp_int2_dates = ['2020-01-09', '2020-01-10', '2020-01-11',
                              '2020-01-12', '2020-01-13', '2020-01-14',
                              '2020-01-15']
            exp_int3_dates = ['2020-01-16', '2020-01-17']
            exp_int4_dates = ['2020-02-21']
            exp_int5_dates = ['2020-03-15']

            self.assertEqual(sel.inclusion.sum(), 5)
            self.assertEqual(sel.metadata.get_column('date-md'), self.md2)
            self.assertEqual(sel.label, 'sample_longitudinal')

            sampled_dates = set(self.md2.to_series()[sel.inclusion].values)
            self.assertEqual(len(sampled_dates & set(exp_int1_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int2_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int3_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int4_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int5_dates)), 1)

    def test_interval_bounds2(self):
        for _ in range(self._N_TEST_ITERATIONS):
            sel = sample_longitudinal(self.md2,
                                      samples_per_interval=1,
                                      start_date='2019-12-27')

            exp_int1_dates = ['2020-01-02']
            exp_int2_dates = ['2020-01-03', '2020-01-04', '2020-01-05',
                              '2020-01-06', '2020-01-07', '2020-01-08',
                              '2020-01-09']
            exp_int3_dates = ['2020-01-10', '2020-01-11', '2020-01-12',
                              '2020-01-13', '2020-01-14', '2020-01-15',
                              '2020-01-16']
            exp_int4_dates = ['2020-01-17']
            exp_int5_dates = ['2020-02-21']
            exp_int6_dates = ['2020-03-15']

            self.assertEqual(sel.inclusion.sum(), 6)
            self.assertEqual(sel.metadata.get_column('date-md'), self.md2)
            self.assertEqual(sel.label, 'sample_longitudinal')

            sampled_dates = set(self.md2.to_series()[sel.inclusion].values)
            self.assertEqual(len(sampled_dates & set(exp_int1_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int2_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int3_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int4_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int5_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int6_dates)), 1)

    def test_interval_bounds3(self):
        for _ in range(self._N_TEST_ITERATIONS):
            sel = sample_longitudinal(self.md2,
                                      samples_per_interval=1,
                                      start_date='2019-12-28')

            exp_int1_dates = ['2020-01-02', '2020-01-03']
            exp_int2_dates = ['2020-01-04', '2020-01-05',
                              '2020-01-06', '2020-01-07', '2020-01-08',
                              '2020-01-09', '2020-01-10']
            exp_int3_dates = ['2020-01-11', '2020-01-12',
                              '2020-01-13', '2020-01-14', '2020-01-15',
                              '2020-01-16', '2020-01-17']
            exp_int4_dates = ['2020-02-21']
            exp_int5_dates = ['2020-03-15']

            self.assertEqual(sel.inclusion.sum(), 5)
            self.assertEqual(sel.metadata.get_column('date-md'), self.md2)
            self.assertEqual(sel.label, 'sample_longitudinal')

            sampled_dates = set(self.md2.to_series()[sel.inclusion].values)
            self.assertEqual(len(sampled_dates & set(exp_int1_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int2_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int3_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int4_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int5_dates)), 1)

    def test_interval_size(self):
        for _ in range(self._N_TEST_ITERATIONS):
            sel = sample_longitudinal(self.md2,
                                      start_date='2019-12-19',
                                      samples_per_interval=1,
                                      days_per_interval=14)

            exp_int1_dates = ['2020-01-02', '2020-01-03', '2020-01-04',
                              '2020-01-05', '2020-01-06', '2020-01-07',
                              '2020-01-08', '2020-01-09', '2020-01-10',
                              '2020-01-11', '2020-01-12', '2020-01-13',
                              '2020-01-14', '2020-01-15']
            exp_int2_dates = ['2020-01-16', '2020-01-17']
            exp_int3_dates = ['2020-02-21']
            exp_int4_dates = ['2020-03-15']

            self.assertEqual(sel.inclusion.sum(), 4)
            self.assertEqual(sel.metadata.get_column('date-md'), self.md2)
            self.assertEqual(sel.label, 'sample_longitudinal')

            sampled_dates = set(self.md2.to_series()[sel.inclusion].values)
            self.assertEqual(len(sampled_dates & set(exp_int1_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int2_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int3_dates)), 1)
            self.assertEqual(len(sampled_dates & set(exp_int4_dates)), 1)

    def test_seed(self):
        sel1 = sample_longitudinal(self.md2,
                                   samples_per_interval=1,
                                   start_date='2019-12-26', seed=1)
        for _ in range(self._N_TEST_ITERATIONS):
            sel2 = sample_longitudinal(self.md2,
                                       samples_per_interval=1,
                                       start_date='2019-12-26', seed=1)
            self.assertEqual(list(sel1.inclusion.items()),
                             list(sel2.inclusion.items()))

    def test_seqs_restrict_metadata(self):
        context_seqs = self.get_data_path('context-seqs-4.fasta')
        context_seqs = DNAFASTAFormat(context_seqs, 'r')
        s = pd.Series(['2019-11-01', '2020-01-17'],
                      index=['B', 'U'])
        s.index.name = 'id'
        s.name = 'date-md'
        exp_md = qiime2.CategoricalMetadataColumn(s)

        for _ in range(self._N_TEST_ITERATIONS):
            sel = sample_longitudinal(self.md2, context_seqs)

            self.assertEqual(sel.inclusion.sum(), 2)
            self.assertTrue(sel.inclusion['B'])
            self.assertTrue(sel.inclusion['U'])
            self.assertEqual(sel.metadata.get_column('date-md'), exp_md)
            self.assertEqual(sel.label, 'sample_longitudinal')
