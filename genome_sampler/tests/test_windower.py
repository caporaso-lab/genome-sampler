import unittest

import qiime2
import pandas as pd
import numpy as np
import pandas.testing as pdt
import biom

from genome_sampler.windower import (sliding_window, _make_id, _id_order,
                                     _create_table, _create_window_df,
                                     _group_by_overlapping_windows,
                                     _group_by_nonoverlapping_windows)


class WindowTests(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame([
            ['a', '2020-07-01', 'X', 'Y', 'G1'],
            ['b', '2020-07-02', 'X', 'Y', 'G2'],
            ['c', '2020-07-03', 'X', 'Y', 'G3'],
            ['c2', '2020-07-03', 'X', 'Y', 'G3'],
            ['d', '2020-07-04', 'X', 'Z', 'G3'],
            ['e', '2020-07-07', 'X', 'Z', 'G3'],
            ['f', '2020-07-08', 'X', 'Z', 'G4'],
            ['g', '2020-07-10', 'X', 'Z', 'G4'],
            ['h', '2020-07-20', 'X', 'Z', 'G4'],
            ['i', '2020-07-30', 'W', 'Z', 'G4'],
            ['j', '2020-07-31', 'W', 'Z', 'G4'],
            ['k', '2020-08-01', 'X', 'Z', 'G5'],
            ['l', '2020-08-02', 'X', 'Z', 'G4'],
            ['m', '2020-08-03', 'X', 'Z', 'G5']],
            columns=['feature-id', 'date', 'location1', 'location2', 'strain'])
        self.df['date'] = pd.to_datetime(self.df['date'])
        self.df.set_index('feature-id', inplace=True)

    def test_group_by_nonoverlapping_windows_k1(self):
        exp = [{'starting_timepoint': r['date'],
                'feature-id': [idx]} for idx, r in self.df.iterrows()
               if not idx == 'c2']
        exp[2]['feature-id'].append('c2')

        obs = _group_by_nonoverlapping_windows(self.df, k=1)
        self.assertEqual(obs, exp)

    def test_group_by_nonoverlapping_windows_k2(self):
        exp = [{'starting_timepoint': pd.to_datetime('2020-07-01'),
                'feature-id': ['a', 'b']},
               {'starting_timepoint': pd.to_datetime('2020-07-03'),
                'feature-id': ['c', 'c2', 'd']},
               {'starting_timepoint': pd.to_datetime('2020-07-07'),
                'feature-id': ['e', 'f']},
               {'starting_timepoint': pd.to_datetime('2020-07-09'),
                'feature-id': ['g']},
               {'starting_timepoint': pd.to_datetime('2020-07-19'),
                'feature-id': ['h']},
               {'starting_timepoint': pd.to_datetime('2020-07-29'),
                'feature-id': ['i']},
               {'starting_timepoint': pd.to_datetime('2020-07-31'),
                'feature-id': ['j', 'k']},
               {'starting_timepoint': pd.to_datetime('2020-08-02'),
                'feature-id': ['l', 'm']}]

        obs = _group_by_nonoverlapping_windows(self.df, k=2)
        self.assertEqual(obs, exp)

    def test_group_by_overlapping_window_k1(self):
        exp = [{'starting_timepoint': r['date'],
                'feature-id': [idx]} for idx, r in self.df.iterrows()
               if not idx == 'c2']
        exp[2]['feature-id'].append('c2')

        obs = _group_by_overlapping_windows(self.df, k=1)
        self.assertEqual(obs, exp)

    def test_group_by_overlapping_window_k2(self):
        exp = [{'starting_timepoint': pd.to_datetime('2020-07-01'),
                'feature-id': ['a', 'b']},
               {'starting_timepoint': pd.to_datetime('2020-07-02'),
                'feature-id': ['b', 'c', 'c2']},
               {'starting_timepoint': pd.to_datetime('2020-07-03'),
                'feature-id': ['c', 'c2', 'd']},
               {'starting_timepoint': pd.to_datetime('2020-07-07'),
                'feature-id': ['e', 'f']},
               {'starting_timepoint': pd.to_datetime('2020-07-30'),
                'feature-id': ['i', 'j']},
               {'starting_timepoint': pd.to_datetime('2020-07-31'),
                'feature-id': ['j', 'k']},
               {'starting_timepoint': pd.to_datetime('2020-08-01'),
                'feature-id': ['k', 'l']},
               {'starting_timepoint': pd.to_datetime('2020-08-02'),
                'feature-id': ['l', 'm']}]
        obs = _group_by_overlapping_windows(self.df, k=2, min_count=2)
        self.assertEqual(obs, exp)

    def test_sliding_window_nostrain_k2(self):
        exp_md = pd.DataFrame([['X:Y:2020-07-01', 'X', 'Y', '2020-07-01', 1.],
                               ['X:Y:2020-07-02', 'X', 'Y', '2020-07-02', 2.],
                               ['X:Z:2020-07-07', 'X', 'Z', '2020-07-07', 3.],
                               ['W:Z:2020-07-30', 'W', 'Z', '2020-07-30', 4.],
                               ['X:Z:2020-08-01', 'X', 'Z', '2020-08-01', 5.],
                               ['X:Z:2020-08-02', 'X', 'Z', '2020-08-02', 6.]],
                              columns=['sample-id', 'location1', 'location2',
                                       'date', 'timepoint_gradient'])
        exp_md.set_index('sample-id', inplace=True)
        exp_table = biom.Table(np.array([[1, 0, 0, 0, 0, 0],  # a
                                         [1, 1, 0, 0, 0, 0],  # b
                                         [0, 1, 0, 0, 0, 0],  # c
                                         [0, 1, 0, 0, 0, 0],  # c2
                                         [0, 0, 1, 0, 0, 0],  # e
                                         [0, 0, 1, 0, 0, 0],  # f
                                         [0, 0, 0, 1, 0, 0],  # i
                                         [0, 0, 0, 1, 0, 0],  # j
                                         [0, 0, 0, 0, 1, 0],  # k
                                         [0, 0, 0, 0, 1, 1],  # l
                                         [0, 0, 0, 0, 0, 1]]),  # m
                               ['a', 'b', 'c', 'c2', 'e', 'f', 'i', 'j', 'k',
                                'l', 'm'],
                               list(exp_md.index))

        # qiime doesn't have datetime support right now
        sw_md = self.df.copy()
        sw_md['date'] = sw_md['date'].apply(lambda x: x.strftime('%Y-%m-%d'))

        obs_table, obs_md = sliding_window(qiime2.Metadata(sw_md), 'date',
                                           'location1', 'location2', 2, 2)
        obs_table = obs_table.sort_order(exp_table.ids())
        self.assertEqual(obs_table, exp_table)
        obs_md.sort_values('timepoint_gradient', inplace=True)
        pdt.assert_frame_equal(obs_md, exp_md)

    def test_make_id(self):
        exp = 'fo_o:bar:2020-02-01'
        obs = _make_id('fo o', 'bar', pd.to_datetime("2020-02-01 12:34"))
        self.assertEqual(obs, exp)

    def test_create_window_df(self):
        md = [['s1', 'foo', 'baz', pd.to_datetime('2020-01-10')],
              ['s2', 'bar', 'baz', pd.to_datetime('2020-01-08')],
              ['s3', 'bar', 'stuff', pd.to_datetime('2020-01-15')],
              ['s4', 'foo', 'baz', pd.to_datetime('2020-01-02')],
              ['s5', 'bar', 'baz', pd.to_datetime('2020-01-01')]]
        exp = pd.DataFrame(md, columns=['sample-id', 'c1', 'c2', 'dates'])
        exp.set_index('sample-id', inplace=True)
        exp['timepoint_gradient'] = np.array([4, 3, 5, 2, 1], dtype=float)
        exp['dates'] = ['2020-01-10', '2020-01-08', '2020-01-15',
                        '2020-01-02', '2020-01-01']
        obs = _create_window_df(md, 'c1', 'c2', 'dates')
        pdt.assert_frame_equal(obs, exp)

    def test_create_table(self):
        rcv = {('f1', 's1'): 10,
               ('f3', 's1'): 2,
               ('f2', 's2'): 1}
        exp = biom.Table(np.array([[10, 0], [0, 1], [2, 0]]),
                         ['f1', 'f2', 'f3'], ['s1', 's2'])
        obs = _create_table(rcv)
        self.assertEqual(obs, exp)

    def test_id_order(self):
        ids = set(list("abczyxm"))
        exp_map = {'a': 0, 'b': 1, 'c': 2, 'm': 3, 'x': 4, 'y': 5, 'z': 6}
        exp_ids = list('abcmxyz')

        obs_map, obs_ids = _id_order(ids)
        self.assertEqual(obs_map, exp_map)
        self.assertEqual(obs_ids, exp_ids)


if __name__ == '__main__':
    unittest.main()
