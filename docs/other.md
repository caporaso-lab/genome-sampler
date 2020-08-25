(other)=
# Other topics

(sampling-focal-sequences)=
## Sampling focal sequences

xyz

## Running in parallel

`genome-sampler` can be run in parallel to speed it up. This is done in
different ways depending on whether you're running the steps individually or
through Snakemake.

If you're using Snakemake, you need to edit `Snakefile` and set the `N_THREADS`
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
 --p-n-threads n
```

When running this command, you should set n to be the number of available
processors or cores on a single node of your system. For example, on a
cluster that has nodes with 28 cores, I would run:

```
qiime genome-sampler sample-diversity \
 --i-context-seqs filtered-context-seqs.qza \
 --p-percent-id 0.9995 \
 --o-selection diversity-selection.qza \
 --p-n-threads 28
```

This would use all of the resources on a single node of that cluster. In the
future QIIME 2 will have support for splitting workflows like this across
multiple cluster nodes, but we do not have this support at this time.

## Pre-computed alignment mask

We provide support for masking alignments with the pre-computed mask
provided [here](https://github.com/W-L/ProblematicSites_SARS-CoV2). Let us
know how it works for you. It would be very useful to hear about how your
final results (e.g., a phylogenetic tree) computed using this mask compares to
one created with your own custom mask.
