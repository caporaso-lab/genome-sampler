(install)=
# `genome-sampler` installation instructions

This document provides instructions for installing genome-sampler.

## Installation instructions

### Install Miniconda
[Miniconda](https://conda.io/miniconda.html) provides the conda environment
and package manager, and is the recommended way to install `genome-sampler`.
Follow the instructions for downloading and installing Miniconda. You may
choose either Miniconda2 or Miniconda3 (i.e. Miniconda Python 2 or 3).
`genome-sampler` will work with either version of Miniconda.

### Install `genome-sampler` from source
A conda package will be available in the near future. For the moment we only
provide a source installation.

First create a suitable conda environment:
```
conda create -y -n genome-sampler
conda activate genome-sampler
```

Next install dependencies:
```
conda install \
  -c conda-forge -c bioconda -c qiime2 -c defaults \
  qiime2 q2cli q2templates q2-types q2-feature-table q2-metadata vsearch snakemake
```

Finally install from source:
```
pip install git+https://github.com/caporaso-lab/genome-sampler.git
```
