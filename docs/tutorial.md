(usage-tutorial)=
# Step-by-step tutorial

This document illustrates how to use `genome-sampler` on a small tutorial data set using step-by-step instructions.
This gives you complete control over the analysis, but is more complex to run relative to using the genome-sampler Snakemake workflow as illustrated in [](usage-tutorial-snakemake).

## Download tutorial data
The tutorial data set used here is intended for educational purposes only.
If you're interested in using these sequences for other analyses, we recommend starting with sequence repositories such as GISAID or NCBI Genbank, which have much more recent versions.

Download the tutorial sequences and corresponding metadata using the following commands:

```{describe-usage}
:scope: tutorial

context_metadata = use.init_metadata_from_url(
   'context-metadata',
   'https://raw.githubusercontent.com/caporaso-lab/genome-sampler/r2020.8/snakemake/tutorial-data/context-metadata.tsv')

focal_metadata = use.init_metadata_from_url(
   'focal-metadata',
   'https://raw.githubusercontent.com/caporaso-lab/genome-sampler/r2020.8/snakemake/tutorial-data/focal-metadata.tsv')
```

```{describe-usage}

def fasta_factory(url):
    import urllib.request
    import tempfile
    from io import TextIOWrapper

    from genome_sampler.common import GISAIDDNAFASTAFormat
    from qiime2.plugin.util import transform

    data = urllib.request.urlopen(url)
    ff = GISAIDDNAFASTAFormat()
    with ff.open() as fh:
      fh.write(TextIOWrapper(data).read())
    return ff

def context_fasta_factory():
    url = 'https://raw.githubusercontent.com/caporaso-lab/genome-sampler/r2020.8/snakemake/tutorial-data/context-seqs.fasta'
    return fasta_factory(url)

context_seqs_raw = use.init_format('context-seqs-raw', context_fasta_factory, ext='fasta')

def focal_fasta_factory():
    url = 'https://raw.githubusercontent.com/caporaso-lab/genome-sampler/r2020.8/snakemake/tutorial-data/focal-seqs.fasta'
    return fasta_factory(url)

focal_seqs_raw = use.init_format('focal-seqs-raw', focal_fasta_factory, ext='fasta')
```

## Importing

You'll begin the workflow by importing the fasta files into QIIME 2 Artifacts.
QIIME 2 Artifacts are structured zip files which will contain the data (still in fasta format), but also some QIIME 2-specific information that (among other things) automates your method recording so you don't have to do that.

```{describe-usage}
context_seqs = use.import_from_format('context_seqs', 'FeatureData[Sequence]', context_seqs_raw, 'GISAIDDNAFASTAFormat')
focal_seqs = use.import_from_format('focal_seqs', 'FeatureData[Sequence]', focal_seqs_raw, 'GISAIDDNAFASTAFormat')
```

If you're obtaining context sequences from a public repository, you may encounter context sequences that don't have associated metadata records.
Steps in this workflow that require context sequence metadata would fail as a result.
At this stage we therefore filter any context sequences that don't have metadata records.

```{describe-usage}
context_seqs_w_metadata, = use.action(
  use.UsageAction(plugin_id='feature_table',
                  action_id='filter_seqs'),
  use.UsageInputs(data=context_seqs, metadata=context_metadata),
  use.UsageOutputNames(filtered_data='context_seqs_w_metadata')
)
```

## Quality filtering

Next, we'll apply a quality filter to the sequence data.
Technically this is an optional step for both the context and focal sequences, but in practice if you obtain your context sequences from a public repository you should consider the context sequence filtering to be essential.
You can use this to remove sequences based on their length (if they're too short or too long), and based on the proportion of ambiguous (e.g., `N`) nucleotide characters that they contain.
Here we'll apply this both to the context sequences and the focal sequences.
We'll remove sequences with a proportion of 1% or more ambiguous characters.

```{describe-usage}
filtered_context_seqs, = use.action(
  use.UsageAction(plugin_id='genome_sampler',
                  action_id='filter_seqs'),
  use.UsageInputs(sequences=context_seqs_w_metadata,
                  max_proportion_ambiguous=0.01),
  use.UsageOutputNames(filtered_sequences='filtered_context_seqs')
)

filtered_focal_seqs, = use.action(
  use.UsageAction(plugin_id='genome_sampler',
                  action_id='filter_seqs'),
  use.UsageInputs(sequences=focal_seqs,
                  max_proportion_ambiguous=0.01),
  use.UsageOutputNames(filtered_sequences='filtered_focal_seqs')
)
```

```{tip}
The `--help` parameter can be provided to any of the commands that are used in this tutorial.
[See here](https://use.qiime2.org/en/latest/tutorials/intro.html#exploring-the-available-functionality) for additional detail.
```

At any time, you could get some summary information about your sequences using the following command:

```{describe-usage}
filtered_focal_seqs_summary, = use.action(
  use.UsageAction(plugin_id='feature_table',
                  action_id='tabulate_seqs'),
  use.UsageInputs(data=filtered_focal_seqs),
  use.UsageOutputNames(visualization='filtered_focal_seqs_summary')
)
```

That command will create a [QIIME 2 visualization](https://use.qiime2.org/en/latest/back-matter/glossary.html#term-visualization), which most frequently would be viewed using [QIIME 2 View](https://view.qiime2.org).
Try viewing that file and finding information such as the number of sequences present in this file and the median length of the sequences.
(Your data is not uploaded to a server when you visit [QIIME 2 View](https://view.qiime2.org), so you don't need to be concerned about exposing sensitive research data.)
For additional information on how to view QIIME 2 visualizations, see [here](https://use.qiime2.org/en/latest/how-to-guides/view-visualizations.html).

```{note}
If some of your focal sequences are present in your context sequence collection (for example because you submitted your sequences to GISAID before downloading GISAID), you should remove those sequences from either your focal or context sequence collection.
See [](#removing-sequences) for instructions on how to do that.
```

## Sampling genomes

We're now ready to start sampling our data, and we'll do this in three steps.
These subsampling steps are independent, so can be run in any order.

### Sampling across time ⏱️

First, we'll sample across time.
By default, this will select seven genomes per seven day period.
The file that gets generated as a result here isn't something that is very useful to view directly, but we'll use it in a few minutes to see how many sequences were retained at this step.

```{describe-usage}
dates = use.get_metadata_column('date', 'date', context_metadata)

date_selection, = use.action(
  use.UsageAction(plugin_id='genome_sampler',
                  action_id='sample_longitudinal'),
  use.UsageInputs(context_seqs=filtered_context_seqs,
                  dates=dates),
  use.UsageOutputNames(selection='date_selection')
)
```

### Sampling across biological diversity 🌲

Next, we'll sample across viral diversity.
This will cluster sequences at a percent identity threshold of 99.95%, and select the centroid sequence from each cluster to include in downstream analyses.

```{describe-usage}
diversity_selection, = use.action(
  use.UsageAction(plugin_id='genome_sampler',
                  action_id='sample_diversity'),
  use.UsageInputs(context_seqs=filtered_context_seqs,
                  percent_id=0.9995),
  use.UsageOutputNames(selection='diversity_selection')
)
```

### Sampling across location of isolation 🗺️

Last, we'll sample near neighbors of the focal sequences from the context sequences.
Of the near-neighbor sequences that are identified for each focal sequence (up to 10 by default), 3 sequences will be selected at random such that each `location` (defined in the `context-metadata.tsv` file) has an equal probability of being selected for inclusion in downstream analysis.

```{describe-usage}
locale = use.get_metadata_column('location', 'location', context_metadata)

neighbor_selection, = use.action(
  use.UsageAction(plugin_id='genome_sampler',
                  action_id='sample_neighbors'),
  use.UsageInputs(focal_seqs=filtered_focal_seqs,
                  context_seqs=filtered_context_seqs,
                  locale=locale,
                  percent_id=0.9999,
                  samples_per_cluster=3),
  use.UsageOutputNames(selection='neighbor_selection')
)
```

### Combine and summarize samples

Now, we'll combine the results of the three sampling approaches:

```{describe-usage}
master_selection, = use.action(
  use.UsageAction(plugin_id='genome_sampler',
                  action_id='combine_selections'),
  use.UsageInputs(selections=[date_selection, diversity_selection, neighbor_selection]),
  use.UsageOutputNames(combined_selection='master_selection')
)
```

And we'll generate a summary of the full selection process:

```{describe-usage}
selection_summary, = use.action(
  use.UsageAction(plugin_id='genome_sampler',
                  action_id='summarize_selections'),
  use.UsageInputs(selections=[date_selection, diversity_selection, neighbor_selection]),
  use.UsageOutputNames(visualization='selection_summary')
)
```

The `selection-summary.qzv` file provides a summary of the sampling run.
You can view this file using [QIIME 2 View](https://view.qiime2.org).
Be sure to click the _Provenance_ tab after loading the file on that page - that will provide full details on the workflow that was executed.
How many sequences were retained by each sampling step?

## Preparing sampled data for downstream applications

We're now ready to start compiling our final data set.
To do this, we'll use the `master-selection.qza` file to select specific context sequences from the `filtered-context-seqs.qza` file that was generated earlier.

```{describe-usage}
subsampled_context_seqs, = use.action(
  use.UsageAction(plugin_id='feature_table',
                  action_id='filter_seqs'),
  use.UsageInputs(data=filtered_context_seqs,
                  metadata=use.view_as_metadata('_', master_selection)),
  use.UsageOutputNames(filtered_data='subsampled_context_seqs')
)
```

Next, merge the resulting subsampled context sequences with the focal sequences to create the final QIIME 2 artifact.

```{describe-usage}
sequences, = use.action(
  use.UsageAction(plugin_id='feature_table',
                  action_id='merge_seqs'),
  use.UsageInputs(data=[subsampled_context_seqs, filtered_focal_seqs]),
  use.UsageOutputNames(merged_data='sequences')
)
```

## Exporting data

````{note}
At this point, you're done and you can export a fasta file that you can use in downstream analysis if you plan to leave QIIME 2.
Using the command line interface, this can be done as follows:

```
qiime tools export \
  --input-path sequences.qza \
  --output-path sequences.fasta \
  --output-format DNAFASTAFormat
```

This step needs to be updated so examples can be provided for different interfaces.

````

## Optional: alignment and phylogenetic reconstruction with QIIME 2

QIIME 2 contains tools for preliminary sequence alignment and phylogenetic tree generation in the [q2-alignment](https://docs.qiime2.org/2020.2/plugins/available/alignment/) and [q2-phylogeny](https://docs.qiime2.org/2020.2/plugins/available/phylogeny/) plugins.
QIIME 2 is not designed for high-accuracy phylogenetic analysis, but it can help you to generate and visualize initial alignments and trees.

These steps would take the `sequences.qza` file as input, so if you're interested in using these you could postpone or skip the export step that you ran above.

If you'd like to learn more about this, refer to [](#downstream).
