import qiime2
from qiime2.plugin import (Plugin, Metadata, Int, Range, Str, MetadataColumn,
                           Categorical)

from q2_types.feature_data import FeatureData

import q2_covid
from q2_covid.common import (IDSelectionDirFmt, IDSelection, Selection,
                             IDMetadataFormat, UNIXListFormat)
from q2_covid.subsample_random import subsample_random
from q2_covid.subsample_longitudinal import subsample_longitudinal

plugin = Plugin(
    name='covid',
    website='https://github.com/caporaso-lab/q2-covid',
    package='q2_covid',
    version=q2_covid.__version__,
    description='Tools for genomic epidemiology focused on the SARS-CoV-2'
                ' virus.',
    short_description='Tools for genomic epidemiology.'
)

plugin.register_formats(IDSelectionDirFmt)
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


plugin.methods.register_function(
    function=subsample_random,
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
    name='Randomly sample IDs',
    description='Randomly sample IDs without replacement.'
)


plugin.methods.register_function(
    function=subsample_longitudinal,
    inputs={},
    parameters={
        'dates': MetadataColumn[Categorical],
        'start_date': Str,
        'samples_per_interval': Int % Range(1, None),
        'days_per_interval': Int % Range(1, None),
        'seed': Int % Range(0, None)
    },
    outputs=[('selection', FeatureData[Selection])],
    input_descriptions={},
    parameter_descriptions={
        'dates': 'Dates to sample from.',
        'start_date': ('Start date of first interval. Dates before this date '
                       ' will be excluded. The start date plus the '
                       '`days_per_interval` defines the bounds of the '
                       'sampling intervals. If not provided, this will '
                       'default to the first date in metadata.'),
        'samples_per_interval': ('The number of random dates to select in each '
                                 'interval.'),
        'days_per_interval': ('The length of each interval in days.'),
        'seed': ('Seed used for random number generators.')
    },
    output_descriptions={
        'selection': 'The subsampled dates.'
    },
    name='Subsample dates across time',
    description=('Sample dates at random without replacement '
                 'from each user-defined interval.')
)
