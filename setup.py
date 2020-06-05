from setuptools import find_packages, setup

import versioneer

setup(
    name='genome-sampler',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    packages=find_packages(),
    author="Evan Bolyen",
    author_email="ebolyen@gmail.com",
    description="QIIME 2 plugin for sampling from collections of sequences.",
    url="https://github.com/caporaso-lab/genome-sampler",
    entry_points={
        'qiime2.plugins':
        ['genome-sampler=genome_sampler.plugin_setup:plugin']
    },
    package_data={
        'genome_sampler': ['citations.bib', 'assets/*'],
        'genome_sampler.tests': ['data/*'],
    },
    zip_safe=False,
)
