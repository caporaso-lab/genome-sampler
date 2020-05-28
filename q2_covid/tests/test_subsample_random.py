import unittest

import pandas as pd
import qiime2

from q2_covid.subsample_random import subsample_random


class TestSubsampleRandom(unittest.TestCase):
    def setUp(self):
        df = pd.DataFrame({'foo': range(26)},
                          index=[chr(x) for x in range(65, 91)])
        df.index.name = 'id'
        self.md = qiime2.Metadata(df)

    def test_n_works(self):
        n = 15

        sel = subsample_random(self.md, n)

        self.assertEqual(sel.inclusion.sum(), n)
        self.assertEqual(sel.metadata, self.md)
        self.assertEqual(sel.label, 'subsample_random')

    def test_n_too_large(self):
        n = 27

        with self.assertRaisesRegex(ValueError, 'larger than'):
            subsample_random(self.md, n)

    def test_n_matches(self):
        n = self.md.id_count

        sel = subsample_random(self.md, n)

        self.assertEqual(sel.inclusion.sum(), n)

    def test_seed_works(self):
        n = 15
        seed = 123
        sel = subsample_random(self.md, n, seed=seed)

        exp = list(sel.inclusion.items())

        for _ in range(50):
            sel = subsample_random(self.md, n, seed=seed)
            self.assertEqual(list(sel.inclusion.items()), exp)



if __name__ == '__main__':
    unittest.main()
