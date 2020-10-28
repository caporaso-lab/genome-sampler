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

### Install `wget`

```bash
conda install wget
```

### Create a conda environment for `genome-sampler`

For linux installation environments, please run:

```bash
wget https://raw.githubusercontent.com/caporaso-lab/genome-sampler/master/conda-env-files/genome-sampler-py36-linux-conda.yml
conda env create -n genome-sampler-2020.8 --file genome-sampler-py36-linux-conda.yml
rm genome-sampler-py36-linux-conda.yml
```

For macOS installation environments, please run:

```bash
wget https://raw.githubusercontent.com/caporaso-lab/genome-sampler/master/conda-env-files/genome-sampler-py36-osx-conda.yml
conda env create -n genome-sampler-2020.8 --file genome-sampler-py36-osx-conda.yml
rm genome-sampler-py36-osx-conda.yml
```

### Activate the `genome-sampler` conda environment

```bash
conda activate genome-sampler-2020.8
```
