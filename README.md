# genome-sampler

![lint-build-test](https://github.com/caporaso-lab/genome-sampler/workflows/lint-build-test/badge.svg)

Tools for sampling viral genomes across time of isolation, location of isolation, and genome divergence.

This software is a [QIIME 2](https://qiime2.org) plugin implementing the subsampling workflow that was initially developed for the [Arizona COVID-19 Genomics Union analysis of SARS-CoV-2 genomes](https://www.medrxiv.org/content/10.1101/2020.05.08.20095935v1).

If you're interested in contributing to genome-sampler, please review the software project's [code of conduct](https://github.com/caporaso-lab/code-of-conduct/blob/master/code-of-conduct.md), which is adapted from the [Contributor Covenant](https://www.contributor-covenant.org), version 1.4.


## Installation (from source)
A conda package will be available in the near future, however for the moment only a source installation is supported.

First create a suitable conda environment:
```
conda create -y -n genome-sampler
conda activate genome-sampler
```

Next install dependencies:
```
conda install -c conda-forge -c bioconda -c qiime2 -c defaults \
  qiime2 q2cli q2templates q2-types q2-feature-table q2-metadata vsearch
```

Finally install from source:
```
pip install git+https://github.com/caporaso-lab/genome-sampler.git
```
