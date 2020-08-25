(downstream)=
# Downstream alignment and phylogenetics with QIIME 2

**Optional**: QIIME 2 contains tools for sequence alignment and
phylogenetic reconstruction in the
[q2-alignment](https://docs.qiime2.org/2020.2/plugins/available/alignment/)
and
[q2-phylogeny](https://docs.qiime2.org/2020.2/plugins/available/phylogeny/)
plugins (which are not installed by default with genome-sampler).

In this document we'll build an alignment with MAFFT {cite}`Katoh2013-lz`,
apply a pre-computed alignment position mask, and then build a phylogenetic
tree with IQTree 2 {cite}`Minh2020-mz`.

```
qiime alignment mafft \
  --i-sequences sequences.qza \
  --o-aligned-sequences aligned-sequences.qza
```

```
qiime phylogeny fasttree \
  --i-alignment aligned-sequences.qza \
  --o-tree tree.qza
```
