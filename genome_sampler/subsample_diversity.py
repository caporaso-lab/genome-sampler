import tempfile

import pandas as pd

import qiime2
from q2_types.feature_data import DNAFASTAFormat

from genome_sampler.common import IDSelection, run_command, ids_from_fasta


# According to the vsearch 2.14.2 documentation, percent_id is defined as:
#  (matching columns) / (alignment length - terminal gaps)
def subsample_diversity(context_seqs: DNAFASTAFormat,
                        percent_id: float,
                        max_accepts: int = 10,
                        n_threads: int = 1) -> IDSelection:

    context_ids = ids_from_fasta(str(context_seqs))
    inclusion = pd.Series(False, index=context_ids, name='inclusion')
    metadata = pd.DataFrame(index=pd.Index(inclusion.index))
    metadata.index.name = 'id'

    with tempfile.NamedTemporaryFile() as uc_out_f:
        command = ['vsearch',
                   '--threads', str(n_threads),
                   '--cluster_fast', str(context_seqs),
                   '--id', str(percent_id),
                   '--uc', uc_out_f.name,
                   '--qmask', 'none',
                   '--maxaccepts', str(max_accepts),
                   ]
        run_command(command)

        uc = pd.read_csv(
            uc_out_f.name, sep='\t', na_values='*',
            names=['type', 'cluster_id', 'length', 'perc_id', 'strand',
                   'BLANK1', 'BLANK2', 'cigar', 'query', 'target'])

    # the S lines define the cluster centroids
    context_seqs_to_keep = uc[uc['type'] == 'S'].index
    inclusion[context_seqs_to_keep] = True

    return IDSelection(inclusion,
                       qiime2.Metadata(metadata),
                       "subsample_diversity")
