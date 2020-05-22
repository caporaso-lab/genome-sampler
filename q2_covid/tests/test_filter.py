# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import skbio
import pandas as pd

from qiime2.plugin.testing import TestPluginBase

from ..filter import filter_seqs


class FilterTests(TestPluginBase):
    package = 'q2_covid.tests'

    def test_no_filter(self):
        exp = pd.Series({'s1': skbio.DNA('ACGTTNGACA', metadata={'id': 's1'}),
                         's2': skbio.DNA('A', metadata={'id': 's2'}),
                         's3': skbio.DNA('NNNNNN', metadata={'id': 's3'})})
        obs = filter_seqs(exp)

        self.assertEqual(list(obs.index), list(exp.index))
        self.assertEqual(list(obs), list(exp))

    def test_too_short(self):
        inp = pd.Series({'s1': skbio.DNA('ACGTTGACA', metadata={'id': 's1'}),
                         's2': skbio.DNA('AA', metadata={'id': 's2'})})
        exp = pd.Series({'s1': skbio.DNA('ACGTTGACA', metadata={'id': 's1'})})
        obs = filter_seqs(inp, min_length=3)

        self.assertEqual(list(obs.index), list(exp.index))
        self.assertEqual(list(obs), list(exp))

    def test_too_long(self):
        inp = pd.Series({'s1': skbio.DNA('ACGTTGACA', metadata={'id': 's1'}),
                         's2': skbio.DNA('AA', metadata={'id': 's2'})})
        exp = pd.Series({'s2': skbio.DNA('AA', metadata={'id': 's2'})})
        obs = filter_seqs(inp, max_length=3)

        self.assertEqual(list(obs.index), list(exp.index))
        self.assertEqual(list(obs), list(exp))

    def test_too_ambiguous(self):
        inp = pd.Series({'s1': skbio.DNA('ACGTTGACANNNN', metadata={'id': 's1'}),
                         's2': skbio.DNA('AA', metadata={'id': 's2'})})
        exp = pd.Series({'s2': skbio.DNA('AA', metadata={'id': 's2'})})
        obs = filter_seqs(inp, max_ambiguous_fraction=.3)

        self.assertEqual(list(obs.index), list(exp.index))
        self.assertEqual(list(obs), list(exp))
    
    def test_too_long_and_too_ambiguous(self):
        inp = pd.Series({'s1': skbio.DNA('ACGTTGACANNNN', metadata={'id': 's1'}),
                         's2': skbio.DNA('AA', metadata={'id': 's2'})})
        exp = pd.Series({'s2': skbio.DNA('AA', metadata={'id': 's2'})})
        obs = filter_seqs(inp, max_ambiguous_fraction=.3, max_length=5)

        self.assertEqual(list(obs.index), list(exp.index))
        self.assertEqual(list(obs), list(exp))

    def test_too_short_and_too_ambiguous(self):
        inp = pd.Series({'s1': skbio.DNA('ACGTTGACA', metadata={'id': 's1'}),
                         's2': skbio.DNA('AAN', metadata={'id': 's2'})})
        exp = pd.Series({'s1': skbio.DNA('ACGTTGACA', metadata={'id': 's1'})})
        obs = filter_seqs(inp, max_ambiguous_fraction=.3, min_length=4)

        self.assertEqual(list(obs.index), list(exp.index))
        self.assertEqual(list(obs), list(exp))

    def test_empty_return(self):
        inp = pd.Series({'s1': skbio.DNA('ACGTTGACA', metadata={'id': 's1'}),
                         's2': skbio.DNA('AAN', metadata={'id': 's2'})})
        exp = pd.Series()
        obs = filter_seqs(inp, min_length=29000)

        self.assertEqual(list(obs.index), list(exp.index))
        self.assertEqual(list(obs), list(exp))
