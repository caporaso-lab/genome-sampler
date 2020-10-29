import pandas as pd
import pandas.testing as pdt
import numpy as np

import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat

from genome_sampler.sample_neighbors import (
    sample_neighbors, _clusters_from_vsearch_out, _sample_cluster,
    _generate_weights)


class TestSubsampleNeighbors(TestPluginBase):
    package = 'genome_sampler.tests'
    _N_TEST_ITERATIONS = 250

    def setUp(self):
        super().setUp()
        focal_seqs1 = self.get_data_path('focal-seqs-1.fasta')
        self.focal_seqs1 = DNAFASTAFormat(focal_seqs1, 'r')

        context_seqs1 = self.get_data_path('context-seqs-1.fasta')
        self.context_seqs1 = DNAFASTAFormat(context_seqs1, 'r')

        context_md1 = self.get_data_path('context-metadata-1.tsv')
        self.context_md1 = qiime2.Metadata.load(context_md1)

        focal_seqs2 = self.get_data_path('focal-seqs-2.fasta')
        self.focal_seqs2 = DNAFASTAFormat(focal_seqs2, 'r')

        context_seqs2 = self.get_data_path('context-seqs-2.fasta')
        self.context_seqs2 = DNAFASTAFormat(context_seqs2, 'r')

        context_md2 = self.get_data_path('context-metadata-2.tsv')
        self.context_md2 = qiime2.Metadata.load(context_md2)

    def test_sample_neighbors_no_locale(self):
        sel = sample_neighbors(self.focal_seqs1,
                               self.context_seqs1,
                               percent_id=0.98,
                               samples_per_cluster=2)

        exp_inclusion = pd.Series([True, True, False, False, True, False],
                                  index=['c1', 'c2', 'c3', 'c4', 'c5', 'c6'],
                                  name='inclusion')
        exp_metadata = pd.DataFrame(index=['c1', 'c2', 'c3', 'c4', 'c5', 'c6'])
        exp_metadata.index.name = 'id'
        exp_metadata = qiime2.Metadata(exp_metadata)

        pdt.assert_series_equal(sel.inclusion, exp_inclusion)
        self.assertEqual(sel.metadata, exp_metadata)
        self.assertEqual(sel.label, 'sample_neighbors')

    def test_sample_neighbors_no_locale_alt_percent_id(self):
        sel = sample_neighbors(self.focal_seqs1,
                               self.context_seqs1,
                               percent_id=1.0,
                               samples_per_cluster=2)

        exp_inclusion = pd.Series([True, True, False, False, False, False],
                                  index=['c1', 'c2', 'c3', 'c4', 'c5', 'c6'],
                                  name='inclusion')
        exp_metadata = pd.DataFrame(index=['c1', 'c2', 'c3', 'c4', 'c5', 'c6'])
        exp_metadata.index.name = 'id'
        exp_metadata = qiime2.Metadata(exp_metadata)

        pdt.assert_series_equal(sel.inclusion, exp_inclusion)
        self.assertEqual(sel.metadata, exp_metadata)
        self.assertEqual(sel.label, 'sample_neighbors')

    def test_sample_neighbors_no_locale_alt_samples_per_cluster(self):
        sel = sample_neighbors(self.focal_seqs1,
                               self.context_seqs1,
                               percent_id=0.98,
                               samples_per_cluster=3)

        exp_inclusion = pd.Series([True, True, True, False, True, False],
                                  index=['c1', 'c2', 'c3', 'c4', 'c5', 'c6'],
                                  name='inclusion')
        exp_metadata = pd.DataFrame(index=['c1', 'c2', 'c3', 'c4', 'c5', 'c6'])
        exp_metadata.index.name = 'id'
        exp_metadata = qiime2.Metadata(exp_metadata)

        pdt.assert_series_equal(sel.inclusion, exp_inclusion)
        self.assertEqual(sel.metadata, exp_metadata)
        self.assertEqual(sel.label, 'sample_neighbors')

    def test_sample_neighbors_locale(self):
        count_obs_c2 = 0
        count_obs_c3 = 0
        count_obs_c4 = 0
        count_obs_c5 = 0

        exp_metadata = self.context_md1.to_dataframe()
        exp_metadata.index.name = 'id'
        exp_metadata = qiime2.Metadata(exp_metadata)

        for _ in range(self._N_TEST_ITERATIONS):
            sel = sample_neighbors(self.focal_seqs1,
                                   self.context_seqs1,
                                   percent_id=0.98,
                                   samples_per_cluster=2,
                                   locale=self.context_md1.get_column('x'))

            obs_sampled_context_seqs = sel.inclusion[sel.inclusion].keys()
            self.assertTrue('c1' in set(obs_sampled_context_seqs))
            self.assertEqual(sel.inclusion.sum(), 3)
            self.assertEqual(len(sel.inclusion), 6)

            self.assertEqual(sel.metadata, exp_metadata)
            self.assertEqual(sel.label, 'sample_neighbors')

            if 'c2' in obs_sampled_context_seqs:
                count_obs_c2 += 1
            if 'c3' in obs_sampled_context_seqs:
                count_obs_c3 += 1
            if 'c4' in obs_sampled_context_seqs:
                count_obs_c4 += 1
            if 'c5' in obs_sampled_context_seqs:
                count_obs_c5 += 1

        # since c2, c3, and c5 all have locale "def" and c4 has locale "hijk",
        # so we expect to see c4 more frequently than any of the other three
        self.assertTrue(count_obs_c4 > count_obs_c2)
        self.assertTrue(count_obs_c4 > count_obs_c3)
        self.assertTrue(count_obs_c4 > count_obs_c5)

    def test_sample_neighbors_locale_w_seed(self):
        exp_metadata = self.context_md1

        # since we're setting a random seed, the result we get the first
        # time is our expected every time
        exp_sel = sample_neighbors(self.focal_seqs1,
                                   self.context_seqs1,
                                   percent_id=0.98,
                                   samples_per_cluster=2,
                                   locale=self.context_md1.get_column('x'),
                                   seed=0)
        self.assertTrue(exp_sel.inclusion['c1'])
        self.assertEqual(exp_sel.inclusion.sum(), 3)
        self.assertEqual(exp_sel.metadata, exp_metadata)

        for _ in range(self._N_TEST_ITERATIONS):
            sel = sample_neighbors(self.focal_seqs1,
                                   self.context_seqs1,
                                   percent_id=0.98,
                                   samples_per_cluster=2,
                                   locale=self.context_md1.get_column('x'),
                                   seed=0)

            pdt.assert_series_equal(sel.inclusion,
                                    exp_sel.inclusion)

    def test_sample_neighbors_invalid_max_accepts(self):
        with self.assertRaisesRegex(ValueError, 'obtained per cluster'):
            sample_neighbors(self.focal_seqs1,
                             self.context_seqs1,
                             percent_id=0.98,
                             samples_per_cluster=11)

    def test_sample_neighbors_terminal_gaps_ignored(self):
        sel = sample_neighbors(self.focal_seqs2,
                               self.context_seqs2,
                               percent_id=1.0,
                               samples_per_cluster=2)

        exp_inclusion = pd.Series([True],
                                  index=['c1'],
                                  name='inclusion')
        exp_metadata = pd.DataFrame(index=['c1'])
        exp_metadata.index.name = 'id'
        exp_metadata = qiime2.Metadata(exp_metadata)

        pdt.assert_series_equal(sel.inclusion, exp_inclusion)
        self.assertEqual(sel.metadata, exp_metadata)
        self.assertEqual(sel.label, 'sample_neighbors')

    def test_sample_neighbors_metadata_superset(self):
        context_md = self.get_data_path('context-metadata-2-extra-ids.tsv')
        context_md = qiime2.Metadata.load(context_md)

        sel = sample_neighbors(self.focal_seqs2,
                               self.context_seqs2,
                               percent_id=1.0,
                               samples_per_cluster=2,
                               locale=context_md.get_column('x'))

        exp_inclusion = pd.Series([True],
                                  index=['c1'],
                                  name='inclusion')
        exp_metadata = context_md.filter_ids(['c1'])

        pdt.assert_series_equal(sel.inclusion, exp_inclusion)
        self.assertEqual(sel.metadata, exp_metadata)
        self.assertEqual(sel.label, 'sample_neighbors')

    def test_sample_neighbors_metadata_subset(self):
        context_md = self.get_data_path('context-metadata-1-missing-id.tsv')
        context_md = qiime2.Metadata.load(context_md)

        with self.assertRaisesRegex(ValueError, 'not present in the metadata'):
            sample_neighbors(self.focal_seqs1,
                             self.context_seqs1,
                             percent_id=0.98,
                             samples_per_cluster=1,
                             locale=context_md.get_column('x'))

    def test_clusters_from_vsearch_out_no_locale(self):
        vsearch_out = pd.DataFrame(
            [['a', 'c4', 5],
             ['a', 'c2', 0],
             ['a', 'c99', 1],
             ['bb', 'c2', 3],
             ['bb', 'c9', 0]],
            columns=['focal_id', 'context_id', 'n_mismatches']
        )
        exp_columns = ['context_id', 'n_mismatches', 'locale']
        exp = {'a': pd.DataFrame([['c4', 5, None],
                                  ['c2', 0, None],
                                  ['c99', 1, None]],
                                 columns=exp_columns),
               'bb': pd.DataFrame([['c2', 3, None],
                                   ['c9', 0, None]],
                                  columns=exp_columns)}
        obs = _clusters_from_vsearch_out(vsearch_out, locale=None)

        self.assertEqual(len(obs), 2)
        pdt.assert_frame_equal(obs['a'], exp['a'])
        pdt.assert_frame_equal(obs['bb'], exp['bb'])

    def test_clusters_from_vsearch_out_locale(self):
        vsearch_out = pd.DataFrame(
            [['a', 'c4', 5],
             ['a', 'c2', 0],
             ['a', 'c99', 1],
             ['bb', 'c2', 3],
             ['bb', 'c9', 0]],
            columns=['focal_id', 'context_id', 'n_mismatches']
        )
        locale = pd.Series(['ab', 'cd', 'ef', 'gh'],
                           index=['c2', 'c4', 'c9', 'c99'])
        exp_columns = ['context_id', 'n_mismatches', 'locale']
        exp = {'a': pd.DataFrame([['c4', 5, 'cd'],
                                  ['c2', 0, 'ab'],
                                  ['c99', 1, 'gh']],
                                 columns=exp_columns),
               'bb': pd.DataFrame([['c2', 3, 'ab'],
                                   ['c9', 0, 'ef']],
                                  columns=exp_columns)}
        obs = _clusters_from_vsearch_out(vsearch_out, locale=locale)

        self.assertEqual(len(obs), 2)
        pdt.assert_frame_equal(obs['a'], exp['a'])
        pdt.assert_frame_equal(obs['bb'], exp['bb'])

    def test_sample_cluster_empty(self):
        cluster = pd.DataFrame([])
        obs = _sample_cluster(cluster, 1, np.random.RandomState())
        exp = []
        self.assertEqual(obs, exp)

    def test_sample_cluster_select_all(self):
        columns = ['context_id', 'n_mismatches', 'locale']
        cluster = pd.DataFrame([['c4', 5, 'cd'],
                                ['c2', 0, 'ab'],
                                ['c99', 1, 'gh']],
                               columns=columns)
        obs = _sample_cluster(cluster, 3, np.random.RandomState())
        exp = ['c4', 'c2', 'c99']
        self.assertEqual(set(obs), set(exp))

        obs = _sample_cluster(cluster, 4, np.random.RandomState())
        exp = ['c4', 'c2', 'c99']
        self.assertEqual(set(obs), set(exp))

    def test_sample_cluster_no_locale_1_sample_per_cluster(self):
        columns = ['context_id', 'n_mismatches', 'locale']
        cluster = pd.DataFrame([['c4', 5, None],
                                ['c2', 0, None],
                                ['c99', 1, None],
                                ['c42', 2, None]],
                               columns=columns)

        obs = _sample_cluster(cluster, 1, np.random.RandomState())
        exp = ['c2']
        self.assertEqual(set(obs), set(exp))

    def test_sample_cluster_no_locale_2_samples_per_cluster(self):
        columns = ['context_id', 'n_mismatches', 'locale']
        cluster = pd.DataFrame([['c4', 5, None],
                                ['c2', 0, None],
                                ['c99', 1, None],
                                ['c42', 2, None]],
                               columns=columns)

        obs = _sample_cluster(cluster, 2, np.random.RandomState())
        exp = ['c2', 'c4']
        self.assertEqual(set(obs), set(exp))

    def test_sample_cluster_no_locale_3_samples_per_cluster(self):
        columns = ['context_id', 'n_mismatches', 'locale']
        cluster = pd.DataFrame([['c4', 5, None],
                                ['c2', 0, None],
                                ['c99', 1, None],
                                ['c42', 2, None]],
                               columns=columns)

        obs = _sample_cluster(cluster, 3, np.random.RandomState())
        exp = ['c2', 'c99', 'c4']
        self.assertEqual(set(obs), set(exp))

    def test_sample_cluster_single_locale_2_samples_per_cluster(self):
        columns = ['context_id', 'n_mismatches', 'locale']
        cluster = pd.DataFrame([['c4', 5, 'abc'],
                                ['c2', 0, 'abc'],
                                ['c99', 1, 'abc'],
                                ['c42', 2, 'abc']],
                               columns=columns)

        obs = _sample_cluster(cluster, 2, np.random.RandomState())
        exp = ['c2', 'c4']
        self.assertEqual(set(obs), set(exp))

    def test_sample_cluster_multiple_locales_2_samples_per_cluster(self):
        columns = ['context_id', 'n_mismatches', 'locale']
        cluster = pd.DataFrame([['c4', 5, 'abc'],
                                ['c2', 0, 'abc'],
                                ['c99', 1, 'def'],
                                ['c42', 2, 'abc']],
                               columns=columns)

        count_obs_c4 = 0
        count_obs_c2 = 0
        count_obs_c99 = 0
        count_obs_c42 = 0

        for _ in range(self._N_TEST_ITERATIONS):
            obs = _sample_cluster(cluster, 2, np.random.RandomState())
            self.assertEqual(len(obs), 2)
            if 'c4' in obs:
                count_obs_c4 += 1
            if 'c2' in obs:
                count_obs_c2 += 1
            if 'c99' in obs:
                count_obs_c99 += 1
            if 'c42' in obs:
                count_obs_c42 += 1

        # c4, c2, and c42 all have locale "abc" and c99 has locale "abc",
        # so we expect to see c99 more frequently than any of the other three
        self.assertTrue(count_obs_c99 > count_obs_c4)
        self.assertTrue(count_obs_c99 > count_obs_c2)
        self.assertTrue(count_obs_c99 > count_obs_c42)

    def test_sample_cluster_missing_locales(self):
        columns = ['context_id', 'n_mismatches', 'locale']
        cluster = pd.DataFrame([['c4', 5, 'abc'],
                                ['c2', 0, float('nan')],
                                ['c99', 1, float('nan')],
                                ['c42', 2, 'abc'],
                                ['c142', 2, 'abc'],
                                ['c242', 2, 'abc']],
                               columns=columns)

        count_obs_c4 = 0
        count_obs_c2 = 0
        count_obs_c99 = 0
        count_obs_c42 = 0
        count_obs_c142 = 0
        count_obs_c242 = 0

        for _ in range(self._N_TEST_ITERATIONS):
            obs = _sample_cluster(cluster, 3, np.random.RandomState())
            self.assertEqual(len(obs), 3)
            if 'c4' in obs:
                count_obs_c4 += 1
            if 'c2' in obs:
                count_obs_c2 += 1
            if 'c99' in obs:
                count_obs_c99 += 1
            if 'c42' in obs:
                count_obs_c42 += 1
            if 'c142' in obs:
                count_obs_c142 += 1
            if 'c242' in obs:
                count_obs_c242 += 1

        # c99 and c2 have unknown locale while all others have "abc" locale
        # so we expect to see c99 amd c2 more frequently
        self.assertTrue(count_obs_c99 > count_obs_c4)
        self.assertTrue(count_obs_c99 > count_obs_c42)
        self.assertTrue(count_obs_c99 > count_obs_c142)
        self.assertTrue(count_obs_c99 > count_obs_c242)
        self.assertTrue(count_obs_c2 > count_obs_c4)
        self.assertTrue(count_obs_c2 > count_obs_c42)
        self.assertTrue(count_obs_c2 > count_obs_c142)
        self.assertTrue(count_obs_c2 > count_obs_c242)

    def test_generate_weights_single_class(self):
        loc = pd.Series(['A'] * 4)

        results = _generate_weights(loc)

        self.assertEqual([1/4] * 4, list(results))

    def test_generate_weights_multiple_classes(self):
        loc = pd.Series(['A', 'B', 'A', 'B', 'C', 'D'])

        results = _generate_weights(loc)

        self.assertEqual([1/8, 1/8, 1/8, 1/8, 2/8, 2/8], list(results))

    def test_generate_weights_missing_data(self):
        loc = pd.Series(['A', 'B', 'A', 'B', float('nan'), float('nan')])

        results = _generate_weights(loc)

        self.assertEqual([1/8, 1/8, 1/8, 1/8, 2/8, 2/8], list(results))
