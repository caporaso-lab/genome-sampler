from setuptools import find_packages, setup

import versioneer

setup(
    name='q2-covid',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    packages=find_packages(),
    author="Evan Bolyen",
    author_email="ebolyen@gmail.com",
    description="TODO",
    url="https://github.com/caporaso-lab/q2-covid",
    entry_points={
        'qiime2.plugins':
        ['q2-covid=q2_covid.plugin_setup:plugin']
    },
    package_data={
        'q2_covid': ['citations.bib']
    },
    zip_safe=False,
)
