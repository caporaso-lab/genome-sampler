import pandas as pd

import qiime2
import qiime2.plugin.model as model
from qiime2.plugin import SemanticType

from q2_types.feature_data import FeatureData


class UNIXListFormat(model.TextFileFormat):
    def _validate_(self, level):
        # any file with lines will be valid
        return True

    def to_list(self):
        with self.open() as fh:
            return [s.strip() for s in fh]


class IDMetadata(model.TextFileFormat):
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
    metadata = model.File('metadata.tsv', format=IDMetadata)
    label = model.File('label.txt', format=UNIXListFormat)


class IDSelection:
    def __init__(self, inclusion: pd.Series, metadata: qiime2.Metadata,
                 label: str):
        self.inclusion = inclusion
        self.metadata = metadata
        self.label = label


Selection = SemanticType('Selection', variant_of=FeatureData.field['type'])