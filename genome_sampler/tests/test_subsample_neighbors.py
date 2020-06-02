import pandas as pd
import pandas.testing as pdt
import numpy as np

import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat

from genome_sampler.subsample_neighbors import (
    subsample_neighbors, _clusters_from_vsearch_out, _sample_cluster)


class TestSubsampleNeighbors(TestPluginBase):
    package = 'genome_sampler.tests'
    _N_TEST_ITERATIONS = 50

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

    def test_subsample_neighbors_no_locale(self):
        sel = subsample_neighbors(self.focal_seqs1,
                                  self.context_seqs1,
                                  self.context_md1,
                                  percent_id=0.98,
                                  samples_per_cluster=2)

        obs_sampled_context_seqs = sel.inclusion[sel.inclusion].keys()
        exp_sampled_context_seqs = ['c1', 'c2', 'c5']
        self.assertEqual(set(obs_sampled_context_seqs),
                         set(exp_sampled_context_seqs))
        self.assertEqual(len(sel.inclusion), len(self.context_md1.ids))
        self.assertEqual(sel.metadata, self.context_md1)
        self.assertEqual(sel.label, 'subsample_neighbors')

    def test_subsample_neighbors_no_locale_alt_percent_id(self):
        sel = subsample_neighbors(self.focal_seqs1,
                                  self.context_seqs1,
                                  self.context_md1,
                                  percent_id=1.0,
                                  samples_per_cluster=2)

        obs_sampled_context_seqs = sel.inclusion[sel.inclusion].keys()
        exp_sampled_context_seqs = ['c1', 'c2']
        self.assertEqual(set(obs_sampled_context_seqs),
                         set(exp_sampled_context_seqs))
        self.assertEqual(len(sel.inclusion), len(self.context_md1.ids))
        self.assertEqual(sel.metadata, self.context_md1)
        self.assertEqual(sel.label, 'subsample_neighbors')

    def test_subsample_neighbors_no_locale_alt_samples_per_cluster(self):
        sel = subsample_neighbors(self.focal_seqs1,
                                  self.context_seqs1,
                                  self.context_md1,
                                  percent_id=0.98,
                                  samples_per_cluster=3)

        obs_sampled_context_seqs = sel.inclusion[sel.inclusion].keys()
        exp_sampled_context_seqs = ['c1', 'c2', 'c3', 'c5']
        self.assertEqual(set(obs_sampled_context_seqs),
                         set(exp_sampled_context_seqs))
        self.assertEqual(len(sel.inclusion), len(self.context_md1.ids))
        self.assertEqual(sel.metadata, self.context_md1)
        self.assertEqual(sel.label, 'subsample_neighbors')

    def test_subsample_neighbors_locale(self):
        count_obs_c2 = 0
        count_obs_c3 = 0
        count_obs_c4 = 0
        count_obs_c5 = 0

        for _ in range(self._N_TEST_ITERATIONS):
            sel = subsample_neighbors(self.focal_seqs1,
                                      self.context_seqs1,
                                      self.context_md1,
                                      percent_id=0.98,
                                      samples_per_cluster=2,
                                      locale='x')

            obs_sampled_context_seqs = sel.inclusion[sel.inclusion].keys()
            self.assertTrue('c1' in set(obs_sampled_context_seqs))
            self.assertEqual(sel.inclusion.sum(), 3)
            self.assertEqual(len(sel.inclusion), len(self.context_md1.ids))
            self.assertEqual(sel.metadata, self.context_md1)
            self.assertEqual(sel.label, 'subsample_neighbors')
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

    def test_subsample_neighbors_locale_w_seed(self):
        # since we're setting a random seed, the result we get the first
        # time is our expected every time
        exp_sel = subsample_neighbors(self.focal_seqs1,
                                      self.context_seqs1,
                                      self.context_md1,
                                      percent_id=0.98,
                                      samples_per_cluster=2,
                                      locale='x',
                                      seed=0)
        exp_sampled_context_seqs = exp_sel.inclusion[exp_sel.inclusion].keys()
        self.assertTrue('c1' in set(exp_sampled_context_seqs))
        self.assertEqual(len(exp_sampled_context_seqs), 3)

        for _ in range(self._N_TEST_ITERATIONS):
            sel = subsample_neighbors(self.focal_seqs1,
                                      self.context_seqs1,
                                      self.context_md1,
                                      percent_id=0.98,
                                      samples_per_cluster=2,
                                      locale='x',
                                      seed=0)

            obs_sampled_context_seqs = sel.inclusion[sel.inclusion].keys()
            self.assertEqual(set(obs_sampled_context_seqs),
                             set(exp_sampled_context_seqs))

    def test_subsample_neighbors_invalid_max_accepts(self):
        with self.assertRaisesRegex(ValueError, 'obtained per cluster'):
            subsample_neighbors(self.focal_seqs1,
                                self.context_seqs1,
                                self.context_md1,
                                percent_id=0.98,
                                samples_per_cluster=11)

    def test_subsample_neighbors_invalid_locale(self):
        with self.assertRaisesRegex(ValueError, 'not a column'):
            subsample_neighbors(self.focal_seqs1,
                                self.context_seqs1,
                                self.context_md1,
                                percent_id=0.98,
                                samples_per_cluster=3,
                                locale='y')

    def test_subsample_neighbors_terminal_gaps_ignored(self):
        sel = subsample_neighbors(self.focal_seqs2,
                                  self.context_seqs2,
                                  self.context_md2,
                                  percent_id=1.0,
                                  samples_per_cluster=2)

        obs_sampled_context_seqs = sel.inclusion[sel.inclusion].keys()
        exp_sampled_context_seqs = ['c1']
        self.assertEqual(set(obs_sampled_context_seqs),
                         set(exp_sampled_context_seqs))
        self.assertEqual(len(sel.inclusion), len(self.context_md2.ids))
        self.assertEqual(sel.metadata, self.context_md2)
        self.assertEqual(sel.label, 'subsample_neighbors')

    def test_subsample_neighbors_metadata_superset(self):
        context_md2 = self.get_data_path('context-metadata-2-extra-ids.tsv')
        context_md2 = qiime2.Metadata.load(context_md2)

        sel = subsample_neighbors(self.focal_seqs2,
                                  self.context_seqs2,
                                  context_md2,
                                  percent_id=1.0,
                                  samples_per_cluster=2)

        obs_sampled_context_seqs = sel.inclusion[sel.inclusion].keys()
        exp_sampled_context_seqs = ['c1']
        self.assertEqual(set(obs_sampled_context_seqs),
                         set(exp_sampled_context_seqs))
        self.assertEqual(sel.metadata.ids, ['c1'])
        self.assertEqual(sel.label, 'subsample_neighbors')

    def test_subsample_neighbors_metadata_subset(self):
        context_md = self.get_data_path('context-metadata-1-missing-id.tsv')
        context_md = qiime2.Metadata.load(context_md)

        with self.assertRaisesRegex(ValueError, 'not contained in the index'):
            subsample_neighbors(self.focal_seqs1,
                                self.context_seqs1,
                                context_md,
                                percent_id=0.98,
                                samples_per_cluster=1)

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

        # since c4, c2, and c42 all have locale "abc" and c99 has locale "abc",
        # so we expect to see c99 more frequently than any of the other three
        self.assertTrue(count_obs_c99 > count_obs_c4)
        self.assertTrue(count_obs_c99 > count_obs_c2)
        self.assertTrue(count_obs_c99 > count_obs_c42)
