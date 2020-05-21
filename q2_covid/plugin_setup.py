import qiime2
from qiime2.plugin import Plugin, Metadata, Int, Range

import skbio
import tempfile

from q2_types.feature_data import (
        FeatureData, DNAIterator, DNAFASTAFormat, 
        DNASequencesDirectoryFormat)

import q2_covid
from q2_covid.common import (IDSelectionDirFmt, IDSelection, Selection,
                             IDMetadataFormat, UNIXListFormat, 
                             GISAIDDNAFASTAFormat)
from q2_covid.subsample_random import subsample_random

plugin = Plugin(
    name='covid',
    website='https://github.com/caporaso-lab/q2-covid',
    package='q2_covid',
    version=q2_covid.__version__,
    description='Tools for sampling collections of genomes.',
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
def _3(fmt: GISAIDDNAFASTAFormat) -> DNASequencesDirectoryFormat:
    data = _read_gisaid_dna_fasta(str(fmt))
    df = DNASequencesDirectoryFormat()
    ff = DNAFASTAFormat()

    with ff.open() as file:
        skbio.io.write(data, format='fasta', into=file)

    df.file.write_data(ff, DNAFASTAFormat)
    return df


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
