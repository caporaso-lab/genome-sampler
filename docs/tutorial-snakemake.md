(usage-tutorial-snakemake)=
# Snakemake tutorial

This document illustrates how to use `genome-sampler` on a small tutorial data set using Snakemake.
This makes genome-sampler very straight-forward to run, but more challenging to customize to your needs relative to using it in a step-by-step manner as illustrated in [](usage-tutorial).

## Download tutorial data
The tutorial data set used here is intended for educational purposes only.
If you're interested in using these sequences for other analyses, we recommend starting with sequence repositories such as GISAID or NCBI Genbank, which have much more recent versions.

Download the tutorial sequences and corresponding metadata using the following commands:

```
mkdir tutorial-data/

wget -O tutorial-data/context-metadata.tsv https://raw.githubusercontent.com/caporaso-lab/genome-sampler/r2020.8/snakemake/tutorial-data/context-metadata.tsv
wget -O tutorial-data/focal-metadata.tsv https://raw.githubusercontent.com/caporaso-lab/genome-sampler/r2020.8/snakemake/tutorial-data/focal-metadata.tsv
wget -O tutorial-data/context-seqs.fasta https://raw.githubusercontent.com/caporaso-lab/genome-sampler/r2020.8/snakemake/tutorial-data/context-seqs.fasta
wget -O tutorial-data/focal-seqs.fasta https://raw.githubusercontent.com/caporaso-lab/genome-sampler/r2020.8/snakemake/tutorial-data/focal-seqs.fasta
```

## Using `genome-sampler` (Snakemake workflow)

The full `genome-sampler` workflow can be run using [Snakemake](https://snakemake.readthedocs.io/en/stable/) [@snakemake-bioi].
Begin by [installing Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Download the Snakemake and associated config file using `curl` as follows:

```
wget -O Snakefile https://raw.githubusercontent.com/caporaso-lab/genome-sampler/r2020.8/snakemake/Snakefile
wget -O config.yaml https://raw.githubusercontent.com/caporaso-lab/genome-sampler/r2020.8/snakemake/config.yaml
```

Place the resulting Snakemake file in the same folder as the sequence and metadata files that you downloaded above.
Then run:

```
snakemake
```

When this workflow completes, there will be two primary output files that you'll use.
`sequences.fasta` will contain your subsampled context sequences and your focal sequences.
You should use this file for downstream analyses, such as alignment and phylogenetic analyses.
`selection-summary.qzv` will provide a summary of the sampling run.
You can view this file using [QIIME 2 View](https://view.qiime2.org).

```{note}
After loading `selection-summary.qzv` with [QIIME 2 View](https://view.qiime2.org), click the _Provenance_ tab to see full details on the workflow that was executed.
This allows you to review exactly what steps were performed, what parameters were used for each step, and what versions of `genome-sampler` and its dependencies were installed when you ran `genome-sampler`.
Keep this file for your records so you can refer to it if any of this information is needed in the future (for example, when publishing your findings).
```

You should now be able to move on to analysis of your own data.
See [](#adapting-tutorial) to learn about what changes you might want to make to your `Snakefile` before running your own analysis.

(adapting-tutorial)=
## Adapting Snakemake workflow for application to your own data
The Snakemake workflow presented above is a good starting point for your own analyses.
There are typically a few things to do to adapt `Snakefile` for your own data.
These changes will be made to the `config.yaml` file which should be found alongside the `Snakefile`.

1. Modify input filepaths as needed.
   The input filepaths listed correspond to the names of the files provided for the tutorial.
   You can either name your files using those names, or update the input filepath values.
2. Modify output filepaths if you'd like these to be different from the ones used in the tutorial.
3. Modify longitudinal, neighbor, and diversity sampling parameters as desired.
   If you end up experimenting with different values for these parameters, which we encourage, we would love to hear about your findings.
   Be aware that increasing the `*_percent_id` parameters will increase the runtime of your analysis, and decreasing those values will decrease the runtime of your analysis.

🐍