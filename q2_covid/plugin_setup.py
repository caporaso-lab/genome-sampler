import qiime2
from qiime2.plugin import Plugin, Metadata, Int, Range

from q2_types.feature_data import FeatureData

import q2_covid
from q2_covid.common import (IDSelectionDirFmt, IDSelection, Selection)
from q2_covid.subsample_random import subsample_random

plugin = Plugin(
    name='covid',
    website='https://github.com/caporaso-lab/q2-covid',
    package='q2_covid',
    version=q2_covid.__version__,
    description='TODO',
    short_description='TODO'
)


plugin.register_semantic_types(Selection)
plugin.register_semantic_type_to_format(FeatureData[Selection],
                                        artifact_format=IDSelectionDirFmt)

@plugin.register_transformer
def _1(obj: IDSelection) -> IDSelectionDirFmt:
    result = IDSelectionDirFmt()

    inclusion = obj.inclusion
    include = inclusion.index[inclusion]
    exclude = inclusion.index[~inclusion]
    with open(result.included.path_maker(), 'w') as fh:
        fh.write('\n'.join(include))
    with open(result.excluded.path_maker(), 'w') as fh:
        fh.write('\n'.join(exclude))

    md = obj.metadata.copy()
    md.index.name = 'id'
    qiime2.Metadata(md).save(result.metadata.path_maker())

    with open(result.label.path_maker(), 'w') as fh:
        fh.write(obj.label)

    return result


plugin.methods.register_function(
    function=subsample_random,
    inputs={},
    parameters={
        'ids': Metadata,
        'n': Int % Range(1, None),
        'seed': Int
    },
    outputs=[('selection', FeatureData[Selection])],
    parameter_descriptions={
        'ids': 'IDs to subsample from.',
        'n': 'How many IDs to sample.',
        'seed': 'Random seed to use to initialize random number generator'
    },
    input_descriptions={},
    output_descriptions={
        'selection': 'The selected IDs'
    },
    name='Randomly sample IDs',
    description=''
)