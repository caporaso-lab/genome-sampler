(install)=
# Installation

This document provides instructions for installing `genome-sampler`.

## Installation instructions

### Install Miniconda
[Miniconda](https://conda.io/miniconda.html) provides the conda environment and package manager, and is currently the only supported way to install `genome-sampler`.
Follow the instructions for downloading and installing Miniconda 3.

### Follow the installation instructions on the QIIME 2 Library

Refer to the [`genome-sampler` installation instructions](https://library.qiime2.org/plugin/caporaso-lab/genome-sampler) on the QIIME 2 Library.
We recommend naming your environment `genome-sampler` (by specifying `--name genome-sampler`) in the `conda env create` command.

### Activate the `genome-sampler` conda environment

```bash
conda activate genome-sampler
```

```{note}
If you used a different environment name when creating your conda environment, you'll need to specify that environment name here instead of `genome-sampler`.
```
