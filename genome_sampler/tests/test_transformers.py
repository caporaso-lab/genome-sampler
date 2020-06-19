import os.path

import skbio

from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNASequencesDirectoryFormat

from ..plugin_setup import GISAIDDNAFASTAFormat


class TestGISAIDDNAFASTAFormatTransformerTests(TestPluginBase):
    package = 'genome_sampler.tests'

    def test_gisaid_dna_fasta_format_to_dna_iterator(self):
        input, obs = self.transform_format(GISAIDDNAFASTAFormat,
                                           DNASequencesDirectoryFormat,
                                           filename='gisaid1.fasta')

        obs_fp = os.path.join(str(obs), 'dna-sequences.fasta')
        obs = skbio.io.read(obs_fp, format='fasta', constructor=skbio.DNA)
        obs = list(obs)

        self.assertEqual(len(obs), 3)

        self.assertEqual(obs[0].metadata['id'], 'USA/AZ-TGXXXX/2020')
        self.assertEqual(str(obs[0]), 'ACGTNTGCATNACANTGCTANNNNNNNN')

        self.assertEqual(obs[1].metadata['id'], 'pangolin/Asdf/JKL/2017')
        self.assertEqual(str(obs[1]), 'ACGTGACCANNNNNNNNNACGTCAGTACAGTACCANN')

        self.assertEqual(obs[2].metadata['id'], 'NC_045512.2')
        self.assertEqual(len(obs[2]), 29903)

    def test_transformer_character_replacements(self):
        input, obs = self.transform_format(GISAIDDNAFASTAFormat,
                                           DNASequencesDirectoryFormat,
                                           filename='gisaid2.fasta')

        obs_fp = os.path.join(str(obs), 'dna-sequences.fasta')
        obs = skbio.io.read(obs_fp, format='fasta', constructor=skbio.DNA)
        obs = list(obs)

        self.assertEqual(len(obs), 3)

        self.assertEqual(obs[0].metadata['id'], 's1')
        self.assertEqual(str(obs[0]), 'ACGTNTGCATNACA')

        self.assertEqual(obs[1].metadata['id'], 's2')
        self.assertEqual(str(obs[1]), 'ACGTNTGCATCA')

        self.assertEqual(obs[2].metadata['id'], 's3')
        self.assertEqual(str(obs[2]), 'ACGTNTGCATTACA')

    def test_transformer_sequence_exclusion(self):
        input, obs = self.transform_format(GISAIDDNAFASTAFormat,
                                           DNASequencesDirectoryFormat,
                                           filename='gisaid3.fasta')

        obs_fp = os.path.join(str(obs), 'dna-sequences.fasta')
        obs = skbio.io.read(obs_fp, format='fasta', constructor=skbio.DNA)
        obs = list(obs)

        self.assertEqual(len(obs), 3)

        self.assertEqual(obs[0].metadata['id'], 's1')
        self.assertEqual(str(obs[0]), 'ACGTNTGCATNACA')

        self.assertEqual(obs[1].metadata['id'], 's2')
        self.assertEqual(str(obs[1]), 'ACGTNTGCATCA')

        self.assertEqual(obs[2].metadata['id'], 's5')
        self.assertEqual(str(obs[2]), 'ACGTNTGCATTACA')

    def test_transformer_sequence_exclusion_last_record(self):
        input, obs = self.transform_format(GISAIDDNAFASTAFormat,
                                           DNASequencesDirectoryFormat,
                                           filename='gisaid4.fasta')

        obs_fp = os.path.join(str(obs), 'dna-sequences.fasta')
        obs = skbio.io.read(obs_fp, format='fasta', constructor=skbio.DNA)
        obs = list(obs)

        self.assertEqual(len(obs), 3)

        self.assertEqual(obs[0].metadata['id'], 's1')
        self.assertEqual(str(obs[0]), 'ACGTNTGCATNACA')

        self.assertEqual(obs[1].metadata['id'], 's2')
        self.assertEqual(str(obs[1]), 'ACGTNTGCATCA')

        self.assertEqual(obs[2].metadata['id'], 's5')
        self.assertEqual(str(obs[2]), 'ACGTNTGCATTACA')
