import tempfile

import pandas as pd
import numpy as np

import qiime2
from q2_types.feature_data import DNAFASTAFormat

from genome_sampler.common import IDSelection, run_command

def _clusters_from_uc(uc):
    hits = uc[uc['type'] == 'H']
    clusters = {}
    for r in hits.itertuples():
        if r.target not in clusters:
            clusters[r.target] = []
        clusters[r.target].append(r)
    return clusters

# According to the vsearch 2.14.2 documentation, percent_id is defined as:
#  (matching columns) / (alignment length - terminal gaps)
def subsample_neighbors(context_seqs: DNAFASTAFormat,
                        ids: qiime2.Metadata,
                        percent_id: float,
                        max_accepts: int = 10,
                        n_threads: int = 1) -> IDSelection:

    df = ids.to_dataframe()
    inclusion = pd.Series(False, index=df.index)

# vsearch --threads 28 --cluster_fast denovo-centroids-0.9999.fasta --id 0.9800 --qmask none --xsize --uc unmatched-0.9800.uc --centroids denovo-centroids-0.9800.fasta

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

        vsearch_out = pd.read_csv(
            open(vsearch_out_f.name), sep='\t', na_values='*',
            names=['focal_id', 'context_id', 'n_mismatches'])

    uc = pd.read_csv(
        uc, sep='\t', na_values='*',
        names=['type', 'cluster_id', 'length', 'perc_id', 'strand', 'BLANK1',
               'BLANK2', 'cigar', 'query', 'target'])

        clusters = _clusters_from_vsearch_out(vsearch_out, locale)
        context_seqs_to_keep = \
            _sample_clusters(clusters, samples_per_cluster, seed=seed)
        inclusion[context_seqs_to_keep] = True

    return IDSelection(inclusion, ids, "subsample_neighbors")
