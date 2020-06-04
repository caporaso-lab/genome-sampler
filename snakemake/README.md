# Snakemake genome-sampler workflow

This workflow will run the entire genome-sampler workflow using Snakemake and QIIME 2.

## 1. Download example data

Download metadata as tsv from the two tabs in [this spreadsheet](https://docs.google.com/spreadsheets/d/18IyZK6gvwcqKrl2U1FnucrC71Q5VSy_qFTz4ktffNK4/edit#gid=0). Name the file downloaded from the context-metadata tab context-metadata.tsv and the file downloaded from the focal-metadata tab focal-metadata.tsv.

Download the two fasta files [here](https://www.dropbox.com/sh/tkb0c4snk5zodj8/AABLCykSiEe5zqv8gTeOSegna?dl=0).

## 2. Edit the Snakefile [Optional]

The `Snakefile` in this directory contains all commands necessary to run the workflow. It also contains the configuration information used to customize this workflow on your local system, instead of relying on a separate configurations file. No edits are necessary to process the example data provided above, but if processing other files the filepaths and parameters should be edited to customize your workflow.

## 3. Run the workflow

From within the `/genome-sampler/snakemake/` directory, run the following command:
```
snakemake
```

And watch the magic happen.

The fully processed and subsampled sequences will be output to `sequences.qza`. All other intermediate files will be saved to the same directory.

🐍
