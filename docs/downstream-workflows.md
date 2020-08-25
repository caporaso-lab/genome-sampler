(downstream)=
# Downstream alignment and phylogenetics with QIIME 2

QIIME 2 contains tools for sequence alignment, alignment filtering, and
phylogenetic reconstruction in the
[q2-alignment](https://docs.qiime2.org/2020.2/plugins/available/alignment/)
and
[q2-phylogeny](https://docs.qiime2.org/2020.2/plugins/available/phylogeny/)
plugins (which are not installed by default with genome-sampler). This document
illustrates how these steps can be applied after working through
{ref}`usage-tutorial`.

```{warning}
QIIME 2 is not designed for high-accuracy phylogenetic analysis. If you use the
workflow illustrated here, we recommend that you consider the resulting trees
to be experimental. In other words, this workflow can give you a quick idea of
what is happening with your data, but you may want to use other tools for your
final trees.
```

QIIME 2 often wraps other widely used tools in QIIME 2 plugins rather than
implement them directly. For example, under-the-hood, `genome-sampler`'s
`sample-diversity` action is using `vsearch` {cite}`vsearch-peerj`. In this
document we'll build an alignment with MAFFT {cite}`mafft7`,
apply a pre-computed alignment position mask, and then build a phylogenetic
tree with IQTree 2 {cite}`iqtree2`.

## Installing other plugins

TODO: Fill this in at release time...

## Obtain reference sequence and alignment mask

Download the two `.qza` files from
[this Dropbox folder](https://www.dropbox.com/sh/tkb0c4snk5zodj8/AABLCykSiEe5zqv8gTeOSegna?dl=0).

## Align sequences and build a tree

First, we'll add the SARS-CoV-2 reference sequence to the sequence collection
obtained in the tutorial. Notice that we're working with the `.qza` file that
was created in that tutorial, not the `.fasta` file that we exported.

```
qiime feature-table merge-seqs \
  --i-data sequences.qza \
  --i-data sarscov2-reference-genome.qza \
  --o-merged-data sequences-w-ref.qza
```

Perform sequence alignment using MAFFT.

```
qiime alignment mafft \
  --i-sequences sequences-w-ref.qza \
  --o-alignment aligned-sequences-w-ref.qza
```

After aligning the sequences, a "mask" can be applied to filter positions from
the alignment that are likely to be uninformative. At present, we're
experimenting with a {ref}`pre-computed alignment mask <alignment-mask>`.

```
qiime genome-sampler mask \
  --i-alignment aligned-sequences-w-ref.qza \
  --i-mask alignment-mask.qza \
  --o-masked-alignment masked-aligned-sequences-w-ref.qza
```

Finally, we build a tree from the resulting alignment. This `.qza` file can
be viewed directly with [iTOL](https://itol.embl.de/) to get a quick look.
We'll also soon be adding support for viewing this with
[Empress](https://github.com/biocore/empress).

```
qiime phylogeny iqtree \
  --i-alignment masked-aligned-sequences-w-ref.qza \
  --o-tree unrooted-tree.qza
```

All of the `.qza` files that were generated in this example can be exported
using `qiime tools export`. Exporting of sequence or alignment files will
provide you with fasta files by default, and exporting of the phylogenetic
tree will provide you with a newick file by default. See the
[QIIME 2 exporting documentation](https://docs.qiime2.org/2020.8/tutorials/exporting/)
for more details.
