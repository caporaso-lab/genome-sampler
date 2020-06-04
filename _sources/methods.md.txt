# Supplementary Bioinformatics Methods

**NOTE: this is a work in progress**

## Environment setup

```
conda env create -n az-covid-1 --file environment.yml
```

## Sequence data preparation

```
$ md5sum Global_allsequences_2020-04-16.fasta
98834fc11a1de1f2ce265c4f36b8531b  Global_allsequences_2020-04-16.fasta

$ md5sum nextstrain-metadata.tsv
1751d2b16e5fa9d62b7f8f8e17bd68dc  nextstrain-metadata.tsv
```

Clean up GISAID sequences
```
$ az-covid-1/gisaid-cleanup.py Global_allsequences_2020-04-16.fasta gisaid-unaligned.fasta
Count read: 9230
Count discarded (contains spaces): 408
Count discarded (>10% N): 163
Count written: 8659
```

Sample GISAID sequences over time
```
$ az-covid-1/subsample-time.py --rate 7,7 --hist time-hist.pdf nextstrain-metadata.tsv time-sampled.tsv
Subsampled 9841 down to 106 samples.

$ az-covid-1/filter-seqs.py --tsv time-sampled.tsv gisaid-unaligned.fasta time-sampled.fasta
For 106 ids to keep, 97 were found out of 8659 total samples
```

Create a single file containing the TGEN/AZ sequences, the Genbank reference
sequence (NC_045512.2), and the time-sampled GISAID sequences to create our
“reference set” of sequences we know we want to include.
```
$ cat tgen-unaligned.fasta reference-genome.fasta time-sampled.fasta > reference-centroids.fasta
```

Search GISAID sequences against the “reference set”, to identify close
matches for all reference sequences. This ensures that for all sequences that
we plan to include, we have their nearest neighbors. This protects against
artificial monophylies. Store the unmatched GISAID sequences from this step
for de novo clustering.
```
$ vsearch --threads 28 --usearch_global gisaid-unaligned.fasta --id 0.9999 --db reference-centroids.fasta --uc reference-search-results-0.9999.uc --qmask none --output_no_hits --notmatched unmatched-0.9999.fasta
```

Cluster the unmatched GISAID sequences de novo at multiple different percent
identity thresholds. This allows us to sample the diversity of the SARS-CoV-2
genomes. We will identify the output that we want to use based on the number
of sequences obtained from that search. Note that the 2nd through 5th de novo
clustering commands are clustering the output of the 1st de novo clustering
command to reduce runtime of each.
```
$ vsearch --threads 28 --cluster_fast unmatched-0.9999.fasta --id 0.9999 --qmask none --xsize --uc unmatched-0.9999.uc --centroids denovo-centroids-0.9999.fasta

$ vsearch --threads 28 --cluster_fast denovo-centroids-0.9999.fasta --id 0.9995 --qmask none --xsize --uc unmatched-0.9995.uc --centroids denovo-centroids-0.9995.fasta

$ vsearch --threads 28 --cluster_fast denovo-centroids-0.9999.fasta --id 0.9990 --qmask none --xsize --uc unmatched-0.9990.uc --centroids denovo-centroids-0.9990.fasta

$ vsearch --threads 28 --cluster_fast denovo-centroids-0.9999.fasta --id 0.9950 --qmask none --xsize --uc unmatched-0.9950.uc --centroids denovo-centroids-0.9950.fasta

$ vsearch --threads 28 --cluster_fast denovo-centroids-0.9999.fasta --id 0.9900 --qmask none --xsize --uc unmatched-0.9900.uc --centroids denovo-centroids-0.9900.fasta
```

We will sample the remaining SARS-CoV-2 diversity at 99.9% identity threshold
as 48 sequences were obtained at that threshold, which is about the number we
expect will be reasonable to use with BEAST.

```
$ grep -c '^>' denovo-centroids-0.99*.fasta
denovo-centroids-0.9900.fasta:9
denovo-centroids-0.9950.fasta:11
denovo-centroids-0.9990.fasta:48
denovo-centroids-0.9995.fasta:106
denovo-centroids-0.9999.fasta:1869
```

For each TGEN/AZ sequence, select up to its three closest matches based on
the global alignment search.
```
$ az-covid-1/sample-clusters.py --query gisaid-unaligned.fasta --uc reference-search-results-0.9999.uc --target tgen-unaligned.fasta --n 3 --tsv nextstrain-metadata.tsv tgen-unaligned-sampled-clusters.tsv

$ az-covid-1/filter-seqs.py --tsv tgen-unaligned-sampled-clusters.tsv gisaid-unaligned.fasta mrca-centroids-siblings.fasta
```

Count the number of sequences in each of our sequence collections. Then merge
them to create the sequence collection we plan to use for this analysis.
```
$ grep -c "^>" reference-centroids.fasta mrca-centroids-siblings.fasta denovo-centroids-0.9990.fasta
reference-centroids.fasta:174
mrca-centroids-siblings.fasta:244
denovo-centroids-0.9990.fasta:48
```

```
$ cat reference-centroids.fasta mrca-centroids-siblings.fasta denovo-centroids-0.9990.fasta > unaligned-preprint-seqs.fasta
```

```
$ grep -c "^>" unaligned-preprint-seqs.fasta
466
```

```
$ md5sum unaligned-preprint-seqs.fasta*
6506af95df6c4c3b8671489418b6773f  unaligned-preprint-seqs.fasta
e801322e946c1b81782a2d6b962be15e  unaligned-preprint-seqs.fasta.gz
```

Brendan identified sequences to exclude based on low quality alignments or
missing metadata. He also requested some specific sequences be added. These
additions and subtractions result in the data used for our analyses.
```
$ az-covid-1/clean-seqs.py --remove brendans-ids-to-drop.tsv --fix-date True unaligned-preprint-seqs.fasta unaligned-cleaned-preprint-seqs-2020-04-28.fasta
54 sequences were duplicated, 37 were removed, and 71 were renamed out of 466 total.

$ grep -c '^>' unaligned-cleaned-preprint-seqs-2020-04-28.fasta
375

$ grep -A 1 EPI_ISL_407215 gisaid-unaligned.fasta >> unaligned-cleaned-preprint-seqs-2020-04-28.fasta

$ grep -A 1 EPI_ISL_406223 gisaid-unaligned.fasta >> unaligned-cleaned-preprint-seqs-2020-04-28.fasta

$ grep -A 1 EPI_ISL_424671 gisaid-unaligned.fasta >> unaligned-cleaned-preprint-seqs-2020-04-28.fasta

$ grep -A 1 EPI_ISL_424668 gisaid-unaligned.fasta >> unaligned-cleaned-preprint-seqs-2020-04-28.fasta

$ grep -c '^>' unaligned-cleaned-preprint-seqs-2020-04-28.fasta
379
```

```
$ gzip -c unaligned-cleaned-preprint-seqs-2020-04-28.fasta > unaligned-cleaned-preprint-seqs-2020-04-28.fasta.gz

$ md5sum unaligned-cleaned-preprint-seqs-2020-04-28.fasta*
07bce6efa67969113c46dc91597e367b unaligned-cleaned-preprint-seqs-2020-04-28.fasta
8bef74159d699415023b829efbd6670d unaligned-cleaned-preprint-seqs-2020-04-28.fasta.gz
```

An additional five genomes sequenced at TGEN were added following this
workflow as these were assembled following the execution of this workflow.
These sequences are ids:
TGEN-CoV-AZ-WMTS-TG268002|Coconino|2020-03-16
TGEN-CoV-AZ-WMTS-TG271878|UNKNOWN|2020-03-20
TGEN-CoV-AZ-WMTS-TG268099|Coconino|2020-03-17
TGEN-CoV-AZ-WMTS-TG271862|Cochise|2020-03-22
TGEN-CoV-AZ-WMTS-TG271435|Coconino|2020-03-2
