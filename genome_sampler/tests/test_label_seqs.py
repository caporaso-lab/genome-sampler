import numpy as np
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

        df = pd.DataFrame({'COL1': ['A', np.nan, 'C'],
                           'COL2': ['1', '2', '3']},
                          index=['id1', 'id2', 'id3'])
        df.index.name = 'id'
        self.missing_md = qiime2.Metadata(df)

        self.seqs = pd.Series(['ACGT', 'ACGT', 'ACGT'],
                              index=['id1', 'id2', 'id3'])
        self.labeled_seqs = pd.Series(['ACGT', 'ACGT', 'ACGT'],
                                      index=['id1+A+1', 'id2+B+2', 'id3+C+3'])
        self.four_seqs = pd.Series(['ACGT', 'ACGT', 'ACGT', 'ACGT'],
                                   index=['id1', 'id2', 'id3', 'id4'])

    def test_label_works(self):
        obs_series = label_seqs(self.seqs, '+', self.md, ['COL1', 'COL2'])

        pdt.assert_series_equal(obs_series, self.labeled_seqs)

    def test_delabel_works(self):
        obs_series = label_seqs(self.labeled_seqs, '+')

        pdt.assert_series_equal(obs_series, self.seqs)

    def test_label_works_one_column(self):
        exp_series = pd.Series(['ACGT', 'ACGT', 'ACGT'],
                               index=['id1+A', 'id2+B', 'id3+C'])

        obs_series = label_seqs(self.seqs, '+', self.md, ['COL1'])

        pdt.assert_series_equal(obs_series, exp_series)

    def test_column_missing_value(self):
        exp_series = pd.Series(['ACGT', 'ACGT', 'ACGT'],
                               index=['id1+A+1', 'id2+missing+2', 'id3+C+3'])

        obs_series = label_seqs(self.seqs, '+', self.missing_md,
                                ['COL1', 'COL2'])

        pdt.assert_series_equal(obs_series, exp_series)

    def test_md_no_columns(self):
        with self.assertRaisesRegex(ValueError, 'Columns and metadata'):
            label_seqs(self.seqs, '+', self.md)

    def test_columns_no_md(self):
        with self.assertRaisesRegex(ValueError, 'Columns and metadata'):
            label_seqs(self.seqs, '+', columns=['COL1', 'COL2'])

    def test_md_missing_id(self):
        with self.assertRaisesRegex(ValueError, "The following.*'id4'"):
            label_seqs(self.four_seqs, '+', self.md, ['COL1', 'COL2'])

    def test_requested_col_not_present(self):
        with self.assertRaisesRegex(ValueError, "The column 'COL3' is not"):
            label_seqs(self.seqs, '+', self.md, ['COL3'])

    def test_more_than_ten_missing(self):
        seqs = ['ACGT'] * 11
        index = ['id%d' % d for d in range(11)]
        eleven_seqs = pd.Series(seqs, index=index)

        df = pd.DataFrame({'COL': [1]}, index=['x'])
        df.index.name = 'id'
        empty_md = qiime2.Metadata(df)

        with self.assertRaisesRegex(ValueError, 'omitted'):
            label_seqs(eleven_seqs, '+', empty_md, ['COL'])

    def test_column_missing_value_contains_delimiter(self):
        with self.assertRaisesRegex(ValueError, ':.*mis:sing.*not allowed'):
            label_seqs(self.seqs, ':', self.missing_md, ['COL1', 'COL2'],
                       'mis:sing')
