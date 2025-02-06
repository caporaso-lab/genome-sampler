(downstream)=
# Downstream alignment and phylogenetics

```{note}
This document hasn't yet been updated to use the Usage API as part of the February 2025 updates.
As a result, only instructions for command line usage have been prepared and the document isn't yet doc-tested (so some commands may become out of date).
Reach out on the [issue tracker](https://github.com/caporaso-lab/genome-sampler/issues) if you notice any issues.
```

QIIME 2 contains tools for sequence alignment, alignment filtering, and phylogenetic reconstruction in the q2-alignment and q2-phylogeny plugins.
These are not installed with genome-sampler by default, but can be accessed by installing the [QIIME 2 amplicon distribution](https://docs.qiime2.org).

QIIME 2 often wraps other widely used tools in QIIME 2 plugins rather than implement them directly.
For example, under-the-hood, `genome-sampler`'s `sample-diversity` action is using `vsearch` [@vsearch].
In this document we'll build an alignment with MAFFT [@mafft7], apply a pre-computed alignment position mask, and then build a phylogenetic tree with IQTree 2 [@iqtree2].

```{warning}
The QIIME 2 used here are not designed for high-accuracy phylogenetic analysis.
If you use the workflow illustrated here, we recommend that you consider the resulting trees to be experimental.
In other words, this workflow can give you a quick idea of what is happening with your data, but you may want to use other tools for your final trees.
```

```{note}
This tutorial assumes that you have already run [](#usage-tutorial).
If you haven't, you'll be missing files that are needed for this tutorial.
```

## Obtain reference sequence and alignment mask

```
wget -O alignment-mask.qza https://raw.githubusercontent.com/caporaso-lab/genome-sampler/r2020.8/snakemake/tutorial-data/alignment-mask.qza
wget -O sarscov2-reference-genome.qza https://raw.githubusercontent.com/caporaso-lab/genome-sampler/r2020.8/snakemake/tutorial-data/sarscov2-reference-genome.qza
```

## Align sequences and build a tree

First, we'll add the SARS-CoV-2 reference sequence to the sequence collection obtained in the tutorial.
Notice that we're working with the `.qza` file that was created in that tutorial, not the `.fasta` file that we exported.

```
qiime feature-table merge-seqs \
  --i-data sequences.qza \
  --i-data sarscov2-reference-genome.qza \
  --o-merged-data sequences-w-ref.qza
```

Next, we'll perform sequence alignment using MAFFT.

```
qiime alignment mafft \
  --i-sequences sequences-w-ref.qza \
  --o-alignment aligned-sequences-w-ref.qza
```

After aligning the sequences, a "mask" can be applied to filter positions from the alignment that are likely to be uninformative.
At present, we're experimenting with a [alignment-mask](#pre-computed alignment mask).

```
qiime genome-sampler mask \
  --i-alignment aligned-sequences-w-ref.qza \
  --i-mask alignment-mask.qza \
  --o-masked-alignment masked-aligned-sequences-w-ref.qza
```

Finally, we build a tree from the resulting alignment.
This will generate an unrooted phylogenetic tree.

```
qiime phylogeny iqtree \
  --i-alignment masked-aligned-sequences-w-ref.qza \
  --o-tree unrooted-tree.qza
```

This `.qza` file can be viewed directly with [iTOL](https://itol.embl.de/) to get a quick look.

All of the `.qza` files that were generated in this example can be exported using `qiime tools export`.
Exporting of sequence or alignment files will provide you with fasta files by default, and exporting of the phylogenetic tree will provide you with a newick file by default.
See the [QIIME 2 exporting documentation](https://docs.qiime2.org/2024.10/tutorials/exporting/) for more details.
