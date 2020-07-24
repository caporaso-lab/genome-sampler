import subprocess

import pandas as pd
import skbio

import qiime2
import qiime2.plugin.model as model
from qiime2.plugin import SemanticType
from q2_types.feature_data import FeatureData
from qiime2.plugin import ValidationError


class UNIXListFormat(model.TextFileFormat):
    def _validate_(self, level):
        # any file with lines will be valid
        return True

    def to_list(self):
        with self.open() as fh:
            return [s.strip() for s in fh]


class BrokenVCFFormat(model.TextFileFormat):
    def _validate_(self, level):
        # any file with lines will be valid
        return True

    def to_list(self):
        with self.open() as fh:
            return [s.strip() for s in fh if not s.startswith('##')]


class IDMetadataFormat(model.TextFileFormat):
    def _validate_(self, level):
        try:
            self.to_metadata()
        except qiime2.metadata.MetadataFileError as md_exc:
            raise model.ValidationError(md_exc) from md_exc

    def to_metadata(self):
        return qiime2.Metadata.load(str(self))


class IDSelectionDirFmt(model.DirectoryFormat):
    included = model.File('included.txt', format=UNIXListFormat)
    excluded = model.File('excluded.txt', format=UNIXListFormat)
    metadata = model.File('metadata.tsv', format=IDMetadataFormat)
    label = model.File('label.txt', format=UNIXListFormat)


class IDSelection:
    def __init__(self, inclusion: pd.Series, metadata: qiime2.Metadata,
                 label: str):
        self.inclusion = inclusion
        self.metadata = metadata
        self.label = label


Selection = SemanticType('Selection', variant_of=FeatureData.field['type'])


# modified from DNAFastaFormat in q2_types to allow lowercase characters
# https://github.com/qiime2/q2-types/blob/058ee0e40e38edaa02b1aad034df37456aeb4ddf/q2_types/feature_data/_format.py#L146
class GISAIDDNAFASTAFormat(model.TextFileFormat):
    def _validate_lines(self, max_lines):
        last_line_was_ID = False
        ids = {}

        with open(str(self), 'rb') as fh:
            try:
                first = fh.read(6)
                if first[:3] == b'\xEF\xBB\xBF':
                    first = first[3:]
                # Empty files should validate
                if first.strip() == b'':
                    return
                if first[0] != ord(b'>'):
                    raise ValidationError("First line of file is not a valid "
                                          "description. Descriptions must "
                                          "start with '>'")
                fh.seek(0)
                for line_number, line in enumerate(fh, 1):
                    line = line.strip()
                    if line_number >= max_lines:
                        return
                    line = line.decode('utf-8-sig')
                    if line.startswith('>'):
                        if last_line_was_ID:
                            raise ValidationError('Multiple consecutive '
                                                  'descriptions starting on '
                                                  f'line {line_number-1!r}')
                        line = line.split()
                        if line[0] == '>':
                            if len(line) == 1:
                                raise ValidationError(
                                    f'Description on line {line_number} is '
                                    'missing an ID.')
                            else:
                                raise ValidationError(
                                    f'ID on line {line_number} starts with a '
                                    'space. IDs may not start with spaces')
                        if line[0] in ids:
                            raise ValidationError(
                                f'ID on line {line_number} is a duplicate of '
                                f'another ID on line {ids[line[0]]}.')
                        ids[line[0]] = line_number
                        last_line_was_ID = True
                    else:
                        last_line_was_ID = False
            except UnicodeDecodeError as e:
                raise ValidationError(f'utf-8 cannot decode byte on line '
                                      f'{line_number}') from e

    def _validate_(self, max_lines):
        level_map = {'min': 100, 'max': float('inf')}
        self._validate_lines(level_map[max_lines])


def run_command(cmd, verbose=True):
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command being run is below. This command cannot "
              "be manually re-run as it will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)


def ids_from_fasta(fasta):
    seqs = skbio.io.read(fasta, format='fasta', constructor=skbio.DNA)
    return [s.metadata['id'] for s in seqs]
