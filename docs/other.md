(other)=
# Other topics

(sampling-focal-sequences)=
## Sampling focal sequences

`genome-sampler` can be used for sampling focal sequences in addition to
context sequences. In fact, all of the exact same sampling approaches can
be applied in exactly the same way. The only additional file you'll need to do
this (relative to the files used in {ref}`usage-tutorial`) is a focal sequence
metadata file. This can mirror the format of the context sequence metadata
file. For more details on the format of these files, see
[Metadata in QIIME 2](https://docs.qiime2.org/2020.6/tutorials/metadata/).

(parallel)=
## Running in parallel

`genome-sampler` can be run in parallel to speed it up. This is done in
different ways depending on whether you're running the steps individually or
through Snakemake.

If you're using Snakemake, you need to edit `Snakefile` to set the `N_THREADS`
value to the number of threads you'd like genome-sampler to use.

If you're running the steps individually you can pass the `--p-n-threads`
option to several of the commands. For example, `sample-diversity` is the
slowest step in the workflow. You can provide the `--p-n-threads` parameter to
run it in parallel:

```
qiime genome-sampler sample-diversity \
 --i-context-seqs filtered-context-seqs.qza \
 --p-percent-id 0.9995 \
 --o-selection diversity-selection.qza \
 --p-n-threads N
```

When running this command, you should set `N` to be the number of available
processors or cores on a single node of your system. For example, on a
cluster that has nodes with 28 cores, you could run:

```
qiime genome-sampler sample-diversity \
 --i-context-seqs filtered-context-seqs.qza \
 --p-percent-id 0.9995 \
 --o-selection diversity-selection.qza \
 --p-n-threads 28
```

This would use all of the resources on a single node of that cluster. In the
future QIIME 2 will have support for splitting workflows like this across
multiple cluster nodes, but we do not have this support now.

(adapting-tutorial)=
## Adapting Snakemake workflow for application to your own data
The Snakemake workflow presented in the tutorial is a good starting point for
your own analyses. There are typically a few things to do to adapt
`Snakefile` for your own data. The changes that you'll make to the `Snakefile`
are all in the section at the top that is annotated as ``CONFIGS``.

1. Modify input filepaths as needed. The input filepaths listed correspond to
the names of the files provided for the tutorial. You can either name your
files using those names, or update the input filepath values.
2. Modify output filepaths if you'd like these to be different from the ones
used in the tutorial.
3. Modify the `N_THREADS` value to match the number of threads that are
accessible to you for your analysis. The default is `1`. This will work, but
may take a _very_ long time to run. See {ref}`parallel` for more information
on this topic.
4. Modify longitudinal, neighbor, and diversity sampling parameters as
desired. If you end up experimenting with different values for these
parameters, which we encourage, we would love to hear about your findings. Be
aware that increasing the `*_PERCENT_ID` parameters will increase the runtime
of your analysis, and decreasing those values will decrease the runtime of your
analysis.

(alignment-mask)=
## Pre-computed alignment mask

We provide support for masking alignments with the pre-computed mask
provided [here](https://github.com/W-L/ProblematicSites_SARS-CoV2). Like the
rest of our {ref}`alignment and phylogenetics workflow <downstream>`,
this is experimental. Let us know how it works for you. It would be very useful
to hear about how your final results (e.g., a phylogenetic tree) computed
using this mask compares to one created with your own custom mask.
