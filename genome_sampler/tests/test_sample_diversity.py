# ----------------------------------------------------------------------------
# Copyright (c) 2020-2025, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import pandas.testing as pdt

import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat

from genome_sampler.sample_diversity import sample_diversity


class TestSubsampleDiversity(TestPluginBase):
    package = 'genome_sampler.tests'

    def setUp(self):
        super().setUp()

    def test_sample_diversity(self):
        context_seqs1 = self.get_data_path('context-seqs-1.fasta')
        context_seqs1 = DNAFASTAFormat(context_seqs1, 'r')

        sel = sample_diversity(context_seqs1,
                               percent_id=0.98)

        exp_inclusion = pd.Series([True, True, False, False, False, True],
                                  index=['c1', 'c2', 'c3', 'c4', 'c5', 'c6'],
                                  name='inclusion')
        exp_metadata = pd.DataFrame(index=['c1', 'c2', 'c3', 'c4', 'c5', 'c6'])
        exp_metadata.index.name = 'id'
        exp_metadata = qiime2.Metadata(exp_metadata)

        pdt.assert_series_equal(sel.inclusion, exp_inclusion)
        self.assertEqual(sel.metadata, exp_metadata)
        self.assertEqual(sel.label, 'sample_diversity')

    def test_sample_diversity_alt_percent_id(self):
        context_seqs1 = self.get_data_path('context-seqs-1.fasta')
        context_seqs1 = DNAFASTAFormat(context_seqs1, 'r')

        sel = sample_diversity(context_seqs1,
                               percent_id=0.99)

        exp_inclusion = pd.Series([True, True, False, False, True, True],
                                  index=['c1', 'c2', 'c3', 'c4', 'c5', 'c6'],
                                  name='inclusion')
        exp_metadata = pd.DataFrame(index=['c1', 'c2', 'c3', 'c4', 'c5', 'c6'])
        exp_metadata.index.name = 'id'
        exp_metadata = qiime2.Metadata(exp_metadata)

        pdt.assert_series_equal(sel.inclusion, exp_inclusion)
        self.assertEqual(sel.metadata, exp_metadata)
        self.assertEqual(sel.label, 'sample_diversity')

    def test_sample_diversity_terminal_gaps_ignored(self):
        context_seqs3 = self.get_data_path('context-seqs-3.fasta')
        context_seqs3 = DNAFASTAFormat(context_seqs3, 'r')

        sel = sample_diversity(context_seqs3,
                               percent_id=1.0)

        exp_inclusion = pd.Series([True, False],
                                  index=['c1', 'c2'],
                                  name='inclusion')
        exp_metadata = pd.DataFrame(index=['c1', 'c2'])
        exp_metadata.index.name = 'id'
        exp_metadata = qiime2.Metadata(exp_metadata)

        pdt.assert_series_equal(sel.inclusion, exp_inclusion)
        self.assertEqual(sel.metadata, exp_metadata)
        self.assertEqual(sel.label, 'sample_diversity')
