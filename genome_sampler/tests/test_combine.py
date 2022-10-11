import pandas as pd
import pandas.testing as pdt
import numpy as np

import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types._format import IDSelection

from genome_sampler.combine import combine_selections


class TestCombineSelections(TestPluginBase):
    package = 'genome_sampler.tests'

    def setUp(self):
        super().setUp()
        df1 = pd.DataFrame([['x'], ['y'], ['z'], ['a']],
                           index=['a', 'b', 'c', 'd'],
                           columns=['locale'])
        df1.index.name = 'id'
        self.md1 = qiime2.Metadata(df1)

        self.sel1 = IDSelection(pd.Series([True, False, False, False],
                                          index=['a', 'b', 'c', 'd'],
                                          name='inclusion'),
                                self.md1,
                                label='sel1')
        self.sel2 = IDSelection(pd.Series([False, True, False, False],
                                          index=['a', 'b', 'c', 'd'],
                                          name='inclusion'),
                                self.md1,
                                label='sel2')
        self.sel3 = IDSelection(pd.Series([False, False, True, False],
                                          index=['a', 'b', 'c', 'd'],
                                          name='inclusion'),
                                self.md1,
                                label='sel3')

    def test_combine_selections_single(self):
        sel = combine_selections([self.sel1])
        pdt.assert_series_equal(sel.inclusion, self.sel1.inclusion)
        self.assertEqual(sel.metadata, self.sel1.metadata)
        self.assertEqual(sel.label, 'combined_selections')

    def test_combine_selections(self):
        sel = combine_selections([self.sel1, self.sel2, self.sel3])

        exp_inclusion = pd.Series([True, True, True, False],
                                  index=['a', 'b', 'c', 'd'],
                                  name='inclusion')
        exp_metadata = self.md1

        pdt.assert_series_equal(sel.inclusion, exp_inclusion)
        self.assertEqual(sel.metadata, exp_metadata)
        self.assertEqual(sel.label, 'combined_selections')

    def test_combine_selections_alt_metadata(self):
        df = pd.DataFrame([[42, 88], [3, 88], [99, 88], [np.nan, 88]],
                          index=['a', 'b', 'c', 'd'],
                          columns=['value', 'time-travel-speed-mph'])
        df.index.name = 'id'
        alt_md = qiime2.Metadata(df)

        sel4 = IDSelection(self.sel3.inclusion, alt_md, 'abc')

        sel = combine_selections([self.sel1, self.sel2, sel4])

        exp_inclusion = pd.Series([True, True, True, False],
                                  index=['a', 'b', 'c', 'd'],
                                  name='inclusion')
        exp_df = pd.DataFrame([['x', 88, 42], ['y', 88, 3],
                               ['z', 88, 99], ['a', 88, np.nan]],
                              index=['a', 'b', 'c', 'd'],
                              columns=['locale',  'time-travel-speed-mph',
                                       'value'])
        exp_df.index.name = 'id'
        exp_md = qiime2.Metadata(exp_df)

        pdt.assert_series_equal(sel.inclusion, exp_inclusion)
        self.assertEqual(sel.metadata, exp_md)
        self.assertEqual(sel.label, 'combined_selections')

    def test_combine_selections_inconsistent_metadata(self):
        df = pd.DataFrame([['x'], ['y'], ['w'], ['a']],
                          index=['a', 'b', 'c', 'd'],
                          columns=['locale'])
        df.index.name = 'id'
        alt_md = qiime2.Metadata(df)

        sel4 = IDSelection(self.sel3.inclusion, alt_md, 'abc')

        with self.assertRaisesRegex(ValueError, 'inconsistent metadata'):
            combine_selections([self.sel1, self.sel2, sel4])

    def test_error_on_non_equal_inclusion_id_sets(self):
        bad_sel1 = IDSelection(pd.Series([False, False, True],
                                         index=['a', 'b', 'c'],
                                         name='inclusion'),
                               self.md1,
                               label='somthing')
        with self.assertRaisesRegex(ValueError, "id sets are not equal"):
            combine_selections([self.sel1, bad_sel1])

        with self.assertRaisesRegex(ValueError, "id sets are not equal"):
            combine_selections([bad_sel1, self.sel1])

        with self.assertRaisesRegex(ValueError, "id sets are not equal"):
            combine_selections([self.sel1, self.sel2, self.sel3, bad_sel1])

    def test_error_on_non_equal_metadata_id_sets(self):
        df = pd.DataFrame([['x'], ['y'], ['z']],
                          index=['a', 'b', 'c'],
                          columns=['locale'])
        df.index.name = 'id'
        bad_md1 = qiime2.Metadata(df)

        bad_sel1 = IDSelection(self.sel1.inclusion,
                               bad_md1,
                               label='somthing')

        with self.assertRaisesRegex(ValueError, "id sets are not equal"):
            combine_selections([self.sel1, bad_sel1])

        with self.assertRaisesRegex(ValueError, "id sets are not equal"):
            combine_selections([bad_sel1, self.sel1])

        with self.assertRaisesRegex(ValueError, "id sets are not equal"):
            combine_selections([self.sel1, self.sel2, self.sel3, bad_sel1])
