import skbio
import pandas as pd

import qiime2
from qiime2.plugin import (
    Plugin,
    Metadata,
    Int,
    Range,
    Str,
    MetadataColumn,
    Categorical,
    Float,
    List,
)
from q2_types.feature_data import (
    FeatureData,
    DNAFASTAFormat,
    DNASequencesDirectoryFormat,
    Sequence,
)

import genome_sampler
from genome_sampler.common import (
    IDSelectionDirFmt,
    IDSelection,
    Selection,
    IDMetadataFormat,
    UNIXListFormat,
    GISAIDDNAFASTAFormat,
)
from genome_sampler.sample_random import sample_random
from genome_sampler.sample_longitudinal import sample_longitudinal
from genome_sampler.sample_neighbors import sample_neighbors
from genome_sampler.sample_diversity import sample_diversity
from genome_sampler.filter import filter_seqs
from genome_sampler.combine import combine_selections
from genome_sampler.summarize import summarize_selection

plugin = Plugin(
    name='genome-sampler',
    website='https://caporasolab.us/genome-sampler',
    package='genome_sampler',
    version=genome_sampler.__version__,
    description='Tools for sampling from collections of genomes.',
    short_description='Genome sampler.'
)

plugin.register_formats(IDSelectionDirFmt)
plugin.register_formats(GISAIDDNAFASTAFormat)
plugin.register_semantic_types(Selection)
plugin.register_semantic_type_to_format(FeatureData[Selection],
                                        artifact_format=IDSelectionDirFmt)


@plugin.register_transformer
def _1(obj: IDSelection) -> IDSelectionDirFmt:
    result = IDSelectionDirFmt()

    inclusion = obj.inclusion
    assert not inclusion.index.has_duplicates

    include = inclusion.index[inclusion]
    exclude = inclusion.index[~inclusion]
    with open(result.included.path_maker(), 'w') as fh:
        fh.write('\n'.join(include))
    with open(result.excluded.path_maker(), 'w') as fh:
        fh.write('\n'.join(exclude))

    obj.metadata.save(result.metadata.path_maker())

    with open(result.label.path_maker(), 'w') as fh:
        fh.write(obj.label)

    return result


@plugin.register_transformer
def _2(fmt: IDSelectionDirFmt) -> qiime2.Metadata:
    md = fmt.metadata.view(IDMetadataFormat).to_metadata()
    return md.filter_ids(fmt.included.view(UNIXListFormat).to_list())


@plugin.register_transformer
def _3(fmt: IDSelectionDirFmt) -> IDSelection:
    md = fmt.metadata.view(IDMetadataFormat).to_metadata()
    inclusion = pd.Series(False, index=md.to_dataframe().index)
    included = fmt.included.view(UNIXListFormat).to_list()
    inclusion[included] = True
    with fmt.label.view(UNIXListFormat).open() as fh:
        label = fh.read().strip()
    return IDSelection(inclusion, md, label)


def _read_gisaid_dna_fasta(path):
    def _cleanup_gen():
        with open(path) as input_f:
            for line in input_f:
                if line.startswith('>'):
                    yield line
                else:
                    # Due to a bug in skbio 0.5.5, the lowercase option can't
                    # be used with skbio.io.read for reading DNA sequences.
                    # Convert sequences to uppercase here.
                    line = line.upper()
                    # Spaces and gap characters can appear in unaligned GISAID
                    # sequence records, so we strip those.
                    line = line.replace('-', '')
                    line = line.replace('.', '')
                    line = line.replace(' ', '')
                    yield line

    result = skbio.io.read(_cleanup_gen(), verify=False,
                           format='fasta', constructor=skbio.DNA)
    return result


@plugin.register_transformer
def _4(fmt: GISAIDDNAFASTAFormat) -> DNASequencesDirectoryFormat:
    data = _read_gisaid_dna_fasta(str(fmt))
    df = DNASequencesDirectoryFormat()
    ff = DNAFASTAFormat()

    with ff.open() as file:
        skbio.io.write(data, format='fasta', into=file)

    df.file.write_data(ff, DNAFASTAFormat)
    return df


plugin.methods.register_function(
    function=sample_random,
    inputs={},
    parameters={
        'ids': Metadata,
        'n': Int % Range(1, None),
        'seed': Int % Range(0, None)
    },
    outputs=[('selection', FeatureData[Selection])],
    parameter_descriptions={
        'ids': 'IDs to subsample from.',
        'n': 'Number of IDs to sample.',
        'seed': 'Random seed to use to initialize random number generator.'
    },
    input_descriptions={},
    output_descriptions={
        'selection': 'The selected IDs.'
    },
    name='Sample IDs at random.',
    description=('Randomly sample IDs without replacement. This is useful '
                 'for evaluation purposes (e.g., to test whether another '
                 'sampling approach is more useful than random sampling).')
)


plugin.methods.register_function(
    function=sample_longitudinal,
    inputs={},
    parameters={
        'dates': MetadataColumn[Categorical],
        'start_date': Str,
        'samples_per_interval': Int % Range(1, None),
        'days_per_interval': Int % Range(1, None),
        'seed': Int % Range(0, None),
    },
    outputs=[('selection', FeatureData[Selection])],
    input_descriptions={},
    parameter_descriptions={
        'dates': 'Dates to sample from.',
        'start_date': 'Start date of first interval. Dates before this date '
                      ' will be excluded. The start date plus the '
                      '`days_per_interval` defines the bounds of the '
                      'sampling intervals. If not provided, this will '
                      'default to the first date in metadata.',
        'samples_per_interval': 'The number of random dates to select in '
                                'each interval.',
        'days_per_interval': 'The length of each interval in days.',
        'seed': 'Seed used for random number generators.',
    },
    output_descriptions={
        'selection': 'The selected ids (i.e., the subsampled dates).'
    },
    name='Sample dates uniformly across time',
    description='Sample dates at random without replacement '
                'from each user-defined interval. Dates should be provided '
                'in ISO-8601 format (see '
                'https://en.wikipedia.org/wiki/ISO_8601) both in metadata '
                'and for `start_date`.',
)


plugin.methods.register_function(
    function=sample_neighbors,
    inputs={'focal_seqs': FeatureData[Sequence],
            'context_seqs': FeatureData[Sequence]},
    parameters={
        'percent_id': Float % Range(0, 1, inclusive_end=True),
        'samples_per_cluster': Int % Range(1, None),
        'locale': MetadataColumn[Categorical],
        'max_accepts': Int % Range(1, None),
        'n_threads': Int % Range(1, None),
        'seed': Int % Range(0, None)
    },
    outputs=[('selection', FeatureData[Selection])],
    input_descriptions={
        'focal_seqs': 'The focal sequences.',
        'context_seqs': 'The context sequences to be sampled from.'},
    parameter_descriptions={
        'percent_id': ('The percent identity threshold for searching. If a '
                       'context sequence matches a focal sequence at greater '
                       'than or equal to this percent identity, the context '
                       'sequence will be considered a neighbor of the focal '
                       'sequence.'),
        'samples_per_cluster': ('The number of context sequences to sample '
                                'per cluster, where clusters are the up-to '
                                '`max_accepts` context sequences that match '
                                'at `percent_id` to a given focal sequence.'),
        'locale': ('The metadata column that contains locale '
                   'data. If provided, sampling will be performed across '
                   'locales. (While this was designed for locale sampling, '
                   'any categorical metadata column could be provided.)'),
        'max_accepts': ('The maximum number of context sequences that match '
                        'a focal sequence at `percent_id` or higher that '
                        'will be identified. Up to `samples_per_cluster` of '
                        'these will be sampled.'),
        'n_threads': 'The number of threads to use for processing.',
        'seed': 'Seed used for random number generators.',
    },
    output_descriptions={
        'selection': 'The selected ids (i.e., the subsampled neighbors).'
    },
    name=('Sample context sequences that are near-neighbors of focal '
          'sequences.'),
    description=('Sample context sequences that are near-neighbors of focal '
                 'sequences, including sampling over locales if provided. '
                 'This is useful for avoiding apparent monophylies of '
                 'focal sequences.'),
)


plugin.methods.register_function(
    function=sample_diversity,
    inputs={'context_seqs': FeatureData[Sequence]},
    parameters={
        'percent_id': Float % Range(0, 1, inclusive_end=True),
        'max_accepts': Int % Range(1, None),
        'n_threads': Int % Range(1, None)
    },
    outputs=[('selection', FeatureData[Selection])],
    input_descriptions={
        'context_seqs': 'The context sequences to be sampled from.'},
    parameter_descriptions={
        'percent_id': ('The percent identity threshold for clustering. '
                       'Context sequences will be dereplicated such that no '
                       'pair of retained sequences will have a percent '
                       'identity to one another that is this high.'),
        'max_accepts': ('The maximum number of context sequences that will '
                        'be queried before a new cluster is defned.'),
        'n_threads': 'The number of threads to use for processing.',
    },
    output_descriptions={
        'selection': ('The selected ids (i.e., the diversity-sampled context '
                      'sequences).')
    },
    name='Identify a divergent collection of context sequences.',
    description=('Sample context sequences to a collection of divergent '
                 'sequences. This is useful for retaining the diversity of '
                 'the context sequences in a smaller data set, and for '
                 'downsampling abundant lineages.'),
)


plugin.methods.register_function(
    function=filter_seqs,
    inputs={'sequences': FeatureData[Sequence]},
    parameters={
        'min_length': Int % Range(1, None),
        'max_length': Int % Range(1, None),
        'max_proportion_ambiguous': Float % Range(0, 1, inclusive_end=True)
    },
    outputs=[('filtered_sequences', FeatureData[Sequence])],
    parameter_descriptions={
        'min_length': 'The minimum length of a sequence that will allow'
                      ' it be retained.',
        'max_length': 'The maximum length of a sequence that will allow'
                      ' the sequence to be retained.',
        'max_proportion_ambiguous':
            'The maximum proportion of sequence characters that can be '
            'ambiguous (e.g., N) that will allow it be retained.',
    },
    input_descriptions={'sequences': 'The sequences to be filtered.'},
    output_descriptions={
        'filtered_sequences': 'The sequences retained after filtering.'
    },
    name='Filter sequences.',
    description='Filter sequences based on their length and ambiguity.',
)


plugin.methods.register_function(
    function=combine_selections,
    inputs={'selections': List[FeatureData[Selection]]},
    parameters={},
    outputs=[('combined_selection', FeatureData[Selection])],
    parameter_descriptions={},
    input_descriptions={'selections': 'The id selections to be combined.'},
    output_descriptions={'combined_selection': 'The combined id selection.'},
    name='Combine id selections.',
    description='Combine list of id selections into single id selection.'
)

plugin.visualizers.register_function(
    function=summarize_selection,
    inputs={'selections': List[FeatureData[Selection]]},
    parameters={},
    input_descriptions={'selections': 'Selections to summarize.'},
    parameter_descriptions={},
    name='Summarize one or more selections.',
    description='Provide basic summary statistics on the number of IDs'
                ' selected by the provided selections.'
)
