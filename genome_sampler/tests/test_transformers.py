# ----------------------------------------------------------------------------
# Copyright (c) 2020-2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path

import skbio
import pandas as pd

from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNASequencesDirectoryFormat

from ..plugin_setup import GISAIDDNAFASTAFormat, VCFMaskFormat


class GISAIDDNAFASTAFormatTransformerTests(TestPluginBase):
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


class VCFMaskFormatTransformerTests(TestPluginBase):
    package = 'genome_sampler.tests'

    def test_vcf_mask_to_df1(self):
        input, obs = self.transform_format(VCFMaskFormat,
                                           pd.DataFrame,
                                           filename='mask1.vcf')

        self.assertEqual(obs.shape, (20, 8))
        exp_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
                       'QUAL', 'FILTER', 'INFO']
        self.assertEqual(list(obs.columns), exp_columns)
        self.assertEqual(list(obs[['CHROM', 'POS', 'FILTER']].loc[0]),
                         ['MN908947.3', 1, 'mask'])
        self.assertEqual(list(obs[['CHROM', 'POS', 'FILTER']].loc[9]),
                         ['MN908947.3', 76, 'caution'])
        self.assertEqual(list(obs[['CHROM', 'POS', 'FILTER']].loc[19]),
                         ['MN908947.3', 29903, 'mask'])

    def test_vcf_mask_to_df2(self):
        input, obs = self.transform_format(VCFMaskFormat,
                                           pd.DataFrame,
                                           filename='mask2.vcf')

        self.assertEqual(obs.shape, (2, 8))
        exp_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
                       'QUAL', 'FILTER', 'INFO']
        self.assertEqual(list(obs.columns), exp_columns)
        self.assertEqual(list(obs[['CHROM', 'POS', 'FILTER']].loc[0]),
                         ['s3', 1, 'caution'])
        self.assertEqual(list(obs[['CHROM', 'POS', 'FILTER']].loc[1]),
                         ['s1', 4, 'mask'])

    def test_vcf_mask_validation_failure(self):
        filepath = self.get_data_path('gisaid1.fasta')
        format = VCFMaskFormat(filepath, mode='r')
        with self.assertRaisesRegex(ValidationError, "VCFMaskFormat"):
            format.validate()
