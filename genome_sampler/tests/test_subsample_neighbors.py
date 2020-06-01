import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat

from genome_sampler.subsample_neighbors import subsample_neighbors


class TestSubsampleNeighbors(TestPluginBase):
    package = 'genome_sampler.tests'

    def setUp(self):
        super().setUp()
        focal_seqs1 = self.get_data_path('focal-seqs-1.fasta')
        self.focal_seqs1 = DNAFASTAFormat(focal_seqs1, 'r')

        context_seqs1 = self.get_data_path('context-seqs-1.fasta')
        self.context_seqs1 = DNAFASTAFormat(context_seqs1, 'r')

        context_md1 = self.get_data_path('context-metadata-1.tsv')
        self.context_md1 = qiime2.Metadata.load(context_md1)

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
