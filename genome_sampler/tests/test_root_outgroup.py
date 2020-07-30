from skbio import TreeNode

from qiime2.plugin.testing import TestPluginBase

from genome_sampler.root_outgroup import root_outgroup


class TestRootOutgroup(TestPluginBase):
    package = 'genome_sampler.tests'

    def setUp(self):
        super().setUp()

        self.bifurcated_tree = TreeNode.read(['(((a,b)c,(d,e)f)g,i);'])
        self.trifurcated_tree = TreeNode.read(['(((a,b)c,(d,e)f)g,h,(j,k)i);'])

    def test_bifurcation_shallow_remove(self):
        obs = str(root_outgroup(self.bifurcated_tree, 'i'))
        exp = '(((a,b)c,(d,e)f)g);\n'

        self.assertEqual(obs, exp)

    def test_bifurcation_shallow_keep(self):
        obs = str(root_outgroup(self.bifurcated_tree, 'i', False))
        exp = '(i,(((a,b)c,(d,e)f)g))root;\n'

        self.assertEqual(obs, exp)

    def test_bifurcation_deep_remove(self):
        obs = str(root_outgroup(self.bifurcated_tree, 'a'))
        exp = '(b,((d,e)f,(i)g)c);\n'

        self.assertEqual(obs, exp)

    def test_bifurcation_deep_keep(self):
        obs = str(root_outgroup(self.bifurcated_tree, 'a', False))
        exp = '(a,(b,((d,e)f,(i)g)c))root;\n'

        self.assertEqual(obs, exp)

    def test_trifurcation_shallow_remove(self):
        obs = str(root_outgroup(self.trifurcated_tree, 'h'))
        exp = '(((a,b)c,(d,e)f)g,(j,k)i);\n'

        self.assertEqual(obs, exp)

    def test_trifurcation_shallow_keep(self):
        obs = str(root_outgroup(self.trifurcated_tree, 'h', False))
        exp = '(h,(((a,b)c,(d,e)f)g,(j,k)i))root;\n'

        self.assertEqual(obs, exp)

    def test_trifurcation_deep_remove(self):
        obs = str(root_outgroup(self.trifurcated_tree, 'a'))
        exp = '(b,((d,e)f,(h,(j,k)i)g)c);\n'

        self.assertEqual(obs, exp)

    def test_trifurcation_deep_keep(self):
        obs = str(root_outgroup(self.trifurcated_tree, 'a', False))
        exp = '(a,(b,((d,e)f,(h,(j,k)i)g)c))root;\n'

        self.assertEqual(obs, exp)
