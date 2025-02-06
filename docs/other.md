(other)=
# Other topics

```{note}
This document hasn't yet been updated to use the Usage API as part of the February 2025 updates.
As a result, only instructions for command line usage have been prepared and the document isn't yet doc-tested (so some commands may become out of date).
Reach out on the [issue tracker](https://github.com/caporaso-lab/genome-sampler/issues) if you notice any issues.
```

(sampling-focal-sequences)=
## Sampling focal sequences

`genome-sampler` can be used for sampling focal sequences in addition to context sequences.
In fact, all of the exact same sampling approaches can be applied in exactly the same way.
The only additional file you'll need to do this (relative to the files used in [](#usage-tutorial)) is a focal sequence metadata file.
This can mirror the format of the context sequence metadata file.
For more details on the format of these files, see [Metadata file format](https://use.qiime2.org/en/latest/references/metadata.html).

(parallel)=
## Running in parallel

`genome-sampler` can be run in parallel to speed it up.
This is done in different ways depending on whether you're running the steps individually or through Snakemake.

If you're using Snakemake, you need only pass additional cores (via `snakemake --cores N`).
The parallelizable steps will automatically use the provided resources.

If you're running the steps individually you can pass the `--p-n-threads` option to several of the commands.
For example, `sample-diversity` is the slowest step in the workflow.
You can provide the `--p-n-threads` parameter to run it in parallel:

```
qiime genome-sampler sample-diversity \
 --i-context-seqs filtered-context-seqs.qza \
 --p-percent-id 0.9995 \
 --o-selection diversity-selection.qza \
 --p-n-threads N
```

When running this command, you should set `N` to be the number of available processors or cores on a single node of your system.
For example, on a cluster that has nodes with 28 cores, you could run:

```
qiime genome-sampler sample-diversity \
 --i-context-seqs filtered-context-seqs.qza \
 --p-percent-id 0.9995 \
 --o-selection diversity-selection.qza \
 --p-n-threads 28
```

This would use all of the resources on a single node of that cluster.

(alignment-mask)=
## Pre-computed alignment mask

We provide support for masking alignments with the pre-computed mask provided [here](https://github.com/W-L/ProblematicSites_SARS-CoV2).
Like the [](downstream) presented here, this is experimental for working with viral genomes.
Let us know how it works for you. It would be very useful to hear about how your final results (e.g., a phylogenetic tree) computed using this mask compares to one created with your own custom mask.

(removing-sequences)=
## Removing sequences present in both focal and context sequence collections

As many researchers are rapidly submitting their sequence data to public repositories, sometimes your focal sequences may be present in context sequences that you obtain from a public repository.
You can use QIIME 2 to remove those sequences from your context sequence collection.

First, visualize your focal sequence and context sequence collections with QIIME 2 to confirm that you know what QIIME 2 considers to be the sequence
identifiers.
You can do this by running the following two commands:

```
qiime feature-table tabulate-seqs \
 --i-data focal-seqs.qza \
 --o-visualization focal-seqs.qzv

qiime feature-table tabulate-seqs \
 --i-data context-seqs.qza \
 --o-visualization context-seqs.qzv
```

Then, load the two visualizations with [QIIME 2 View](https://view.qiime2.org).
The values in the _Feature ID_ column are the sequence ids. If the ids of the sequences that you want to remove are the same across your two data sets, you can follow the steps under [](#shared-identifiers).
If the identifiers you want to remove are not the same across your two data sets, or you're not sure, follow the steps under [](#differing-identifiers).

(shared-identifiers)=
### Shared identifiers

If your sequence identifiers are identical across your focal and context sequence collections, you can remove context sequences that are also focal sequences as follows:

```
qiime feature-table filter-seqs\
 --i-data context-seqs.qza \
 --m-metadata-file focal-metadata.tsv \
 --p-exclude-ids \
 --o-filtered-data filtered-context-seqs.qza
```

All of your downstream analysis steps should then use the resulting `filtered-context-seqs.qza` file.

Alternatively, if you want to remove these sequences from your focal sequence collection, you could do this as follows:

```
qiime feature-table filter-seqs\
 --i-data focal-seqs.qza \
 --m-metadata-file context-metadata.tsv \
 --p-exclude-ids \
 --o-filtered-data filtered-focal-seqs.qza
```

(differing-identifiers)=
### Differing identifiers
If your sequence identifiers differ across your focal and context sequences (or you're not sure if they're the same or if they differ), you can perform a metadata-based filtering of your sequences.
The easiest way to approach this is to add a new column to your context sequence metadata file (`context-metadata.tsv` in the tutorial), with a column header like `is-focal-sequence`.
You can do this with your favorite spreadsheet editor.
For all of the sequences you would like to remove, add the value `TRUE` to this column.
For all of the sequences you would like to retain, add the value `FALSE` to this column.
Then, save the file and run the following command:

```
qiime feature-table filter-seqs \
 --i-data context-seqs.qza \
 --m-metadata-file context-metadata.tsv \
 --p-where "[is_focal_sequence]='FALSE'" \
 --o-filtered-data filtered-context-seqs.qza
```

This will create a new file, `filtered-context-seqs.qza`, which only contains the sequences that you indicated were not focal sequences in the context metadata.
All of your downstream analysis steps should then use the resulting `filtered-context-seqs.qza` file.

Alternatively, if you want to remove sequences from your focal sequence collection instead of from your context sequence collection, you can edit your focal sequence metadata file, and then use that to filter in the same way that context sequences were filtered here.

```{note}
If your metadata file already contains a column that differentiates the sequences you want to remove from the sequences you want to retain, it's not necessary to add a new column.
Instead, you can design a value for the ``--p-where`` parameter that will allow the filter to operate on data from that column.
```
