import skbio
import pandas as pd
import pandas.testing as pdt

import qiime2
from qiime2.plugin.testing import TestPluginBase

from genome_sampler.label_seqs import label_seqs


class TestLabelSeqs(TestPluginBase):
    package = 'genome_sampler.tests'

    def setUp(self):
        super().setUp()

        df = pd.DataFrame({'COL1': ['A', 'B', 'C'], 'COL2': ['1', '2', '3']},
                          index=['id1', 'id2', 'id3'])
        df.index.name = 'id'
        self.md = qiime2.Metadata(df)
        self.seqs = pd.Series([skbio.DNA('ACGT', metadata={'id': 'id1'}),
                               skbio.DNA('ACGT', metadata={'id': 'id2'}),
                               skbio.DNA('ACGT', metadata={'id': 'id3'})],
                              index=['id1', 'id2', 'id3'])

    def test_label_works(self):
        exp_series = \
            pd.Series(['ACGT', 'ACGT', 'ACGT'],
                      index=['id1+A+1', 'id2+B+2', 'id3+C+3'])

        obs_series = label_seqs(self.seqs, self.md, ['COL1', 'COL2'], '+')

        pdt.assert_series_equal(obs_series, exp_series)
