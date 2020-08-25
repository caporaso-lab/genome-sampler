(usage-tutorial)=
# `genome-sampler` usage tutorial

This document provides illustrates how to use `genome-sampler` on a small
tutorial data set.

## Download tutorial data
The tutorial data set used here is intended for educational purposes only. If
you're interested in using these sequences for other analyses, we recommend
starting with sequence repositories such as GISAID or NCBI Genbank, which may
have updated versions.

Download metadata as tsv from the two tabs in [this
spreadsheet](https://docs.google.com/spreadsheets/d/18IyZK6gvwcqKrl2U1FnucrC71Q5VSy_qFTz4ktffNK4/edit#gid=0).
Name the file downloaded from the `context-metadata` tab
`context-metadata.tsv` and the file downloaded from the `focal-metadata` tab
`focal-metadata.tsv`.

Download the two fasta files from [this Dropbox
folder](https://www.dropbox.com/sh/tkb0c4snk5zodj8/AABLCykSiEe5zqv8gTeOSegna?dl=0).

## Using `genome-sampler` (Snakemake workflow)

The full `genome-sampler` workflow can be run using Snakemake. If you'd like
to get started quickly and use default parameters, start here. If you'd like
more control over your analysis or want to work through the steps
individually, move on to the next section.

Download the Snakemake file using `curl` as follows:

```
curl -sL https://raw.githubusercontent.com/caporaso-lab/genome-sampler/master/snakemake/Snakefile --output Snakefile
```

Place the resulting Snakemake file in the same folder as the sequence and
metadata files that you downloaded above. Then run:

```
snakemake
```

When this workflow completes, there will be two primary output files that
you'll use. `sequences.fasta` will contain your subsampled context sequences
and your focal sequences. You should use this file for downstream analyses,
such as alignment and phylogenetic analyses. `selection-summary.qzv` will
provide a summary of the sampling run. You can view this file using [QIIME 2
View](https://view.qiime2.org). Be sure to click the _Provenance_ tab after
loading the file on that page - that will provide full details on the
workflow that was executed.

You should now be able to move on to analysis of your own data. If you'd like
to modify parameters of the workflow, you can do so by opening the Snakemake
file in a text editor and editing the values in the `CONFIGS` section.

🐍

## Using `genome-sampler` (step-by-step instructions)

You'll begin the workflow by importing the fasta files into QIIME 2
Artifacts. QIIME 2 Artifacts are structured zip files which will contain the
data (still in fasta format), but also some QIIME 2-specific information that
(among other things) automates your method recording so you don't have to do
that.

```
qiime tools import \
  --input-path context-seqs.fasta \
  --output-path context-seqs.qza \
  --input-format GISAIDDNAFASTAFormat \
  --type FeatureData[Sequence]

qiime tools import \
  --input-path focal-seqs.fasta \
  --output-path focal-seqs.qza \
  --input-format GISAIDDNAFASTAFormat \
  --type FeatureData[Sequence]
```

If you're obtaining context sequences from a public repository, you may
encounter context sequences that don't have associated metadata records. Steps
in this workflow that require context sequence metadata would fail as a
result. At this stage we therefore filter any context sequences that don't
have metadata records.

```
qiime feature-table filter-seqs \
  --i-data context-seqs.qza \
  --m-metadata-file context-metadata.tsv \
  --o-filtered-data context-seqs-w-metadata.qza
```

Next, we'll apply a quality filter to the sequence data. Technically this is
an optional step for both the context and focal sequences, but in practice if
you obtain your context sequences from a public repository you should
consider the context sequence filtering to be essential. You can use this to
remove sequences based on their length (if they're too short or too long),
and based on the proportion of ambiguous (e.g., `N`) nucleotide characters
that they contain. Here we'll apply this both to the context sequences and
the focal sequences. We'll remove sequences with a proportion of 1% or more
ambiguous characters.

```
qiime genome-sampler filter-seqs \
  --i-sequences context-seqs-w-metadata.qza \
  --p-max-proportion-ambiguous 0.01 \
  --o-filtered-sequences filtered-context-seqs.qza

qiime genome-sampler filter-seqs \
  --i-sequences focal-seqs.qza \
  --p-max-proportion-ambiguous 0.01 \
  --o-filtered-sequences filtered-focal-seqs.qza
```

If you'd like to see what other options you can control during this filtering
step, you can run:

```
qiime genome-sampler filter-seqs --help
```

The `--help` parameter can be provided to any of the commands that are used
in this tutorial.

At any time, you could get some summary information about your sequences
using the following command:

```
qiime feature-table tabulate-seqs \
  --i-data filtered-focal-seqs.qza \
  --o-visualization filtered-focal-seqs.qza
```

That command will create a QIIME 2 visualization that you can view using
[QIIME 2 View](https://view.qiime2.org). Try viewing that file and finding
information such as the number of sequences present in this file and the
median length of the sequences. (Your data is not uploaded to a server when
you visit [QIIME 2 View](https://view.qiime2.org), so you don't need to be
conerned about exposing sensitive research data.)

We're now ready to start sampling our data, and we'll do this in three steps.
These subsampling steps are independent, so can be run in any order.

First, we'll sample across time. By default, this will select seven genomes
per seven day period. The file that gets generated as a result here isn't
something that is very useful to view directly, but we'll use it in a few
minutes to see how many sequences were retained at this step.

```
qiime genome-sampler sample-longitudinal \
  --i-context-seqs filtered-context-seqs.qza \
  --m-dates-file context-metadata.tsv \
  --m-dates-column date \
  --o-selection date-selection.qza
```

Next, we'll sample across viral diversity. This will cluster sequences at a
percent identity threshold of 99.95%, and select the centroid sequence from
each cluster to include in downstream analyses.

```
qiime genome-sampler sample-diversity \
  --i-context-seqs filtered-context-seqs.qza \
  --p-percent-id 0.9995 \
  --o-selection diversity-selection.qza
```

Last, we'll sample near neighbors of the focal sequences from the context
sequences. Of the near-neighbor sequences that are identified for each focal
sequence (up to 10 by default), 3 sequences will be selected at random such
that each `location` (defined in the `context-metadata.tsv` file) has an
equal probability of being selected for inclusion in downstream analysis.

```
qiime genome-sampler sample-neighbors \
  --i-focal-seqs filtered-focal-seqs.qza \
  --i-context-seqs filtered-context-seqs.qza \
  --m-locale-file context-metadata.tsv \
  --m-locale-column location \
  --p-percent-id 0.9999 \
  --p-samples-per-cluster 3 \
  --o-selection neighbor-selection.qza
```

Now, we'll combine the results of the three sampling approaches and generate
a summary of the full selection process.

```
qiime genome-sampler combine-selections \
  --i-selections date-selection.qza diversity-selection.qza neighbor-selection.qza \
  --o-combined-selection master-selection.qza

qiime genome-sampler summarize-selections \
  --i-selections date-selection.qza diversity-selection.qza neighbor-selection.qza \
  --o-visualization selection-summary.qzv
```

The `selection-summary.qzv` file provides a summary of the sampling run. You
can view this file using [QIIME 2 View](https://view.qiime2.org). Be sure to
click the _Provenance_ tab after loading the file on that page - that will
provide full details on the workflow that was executed. How many sequences
were retained by each sampling step?

We're now ready to start compiling our final data set. To do this, we'll use
the `master-selection.qza` file to select specific context sequences from the
`filtered-context-seqs.qza` file that was generated earlier.

```
qiime feature-table filter-seqs \
  --i-data filtered-context-seqs.qza \
  --m-metadata-file master-selection.qza \
  --o-filtered-data subsampled-context-seqs.qza
```

Next, merge the resulting subsampled context sequences with the focal
sequences to create the final QIIME 2 artifact. Then, export that to generate
a fasta file that you can use in downstream analysis.

```
qiime feature-table merge-seqs \
  --i-data subsampled-context-seqs.qza filtered-focal-seqs.qza \
  --o-merged-data sequences.qza

qiime tools export \
  --input-path sequences.qza \
  --output-path sequences.fasta \
  --output-format DNAFASTAFormat
```

## Optional: alignment and phylogenetic reconstruction with QIIME 2

QIIME 2 contains tools for preliminary sequence alignment and
phylogenetic tree generation in the
[q2-alignment](https://docs.qiime2.org/2020.2/plugins/available/alignment/)
and
[q2-phylogeny](https://docs.qiime2.org/2020.2/plugins/available/phylogeny/)
plugins (which are not installed by default with genome-sampler). Building
high quality trees is not the focus of QIIME 2, but it can help you to
generate and visualize initial alignments and trees.

These steps would take the `sequences.qza` file as input, so you could
postpone or skip the export step that you ran above.

If you'd like to learn more about this, refer to {ref}`downstream`.
