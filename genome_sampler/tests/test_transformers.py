# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
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
