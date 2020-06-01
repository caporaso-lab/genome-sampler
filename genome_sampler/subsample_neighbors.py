import tempfile

import pandas as pd
import numpy as np

import qiime2
from q2_types.feature_data import DNAFASTAFormat

from genome_sampler.common import IDSelection, run_command


def _clusters_from_vsearch_out(vsearch_out, locale):
    clusters = {}
    if locale is None:
        def f(context_id):
            return None
    else:
        def f(context_id):
            return locale[context_id]
    for r in vsearch_out.itertuples():
        if r.focal_id not in clusters:
            clusters[r.focal_id] = []
        clusters[r.focal_id].append((r.context_id,
                                     r.n_mismatches,
                                     f(r.context_id)))
    return clusters


def _sort_on_n_mismatches(element):
    return element[1]


def _sample_cluster(sorted_cluster, samples_per_cluster):

    len_sorted_cluster = len(sorted_cluster)
    if len_sorted_cluster == 0:
        return []

    indices = map(int,
                  np.linspace(0, len_sorted_cluster-1, samples_per_cluster))
    return [sorted_cluster[i][0] for i in indices]


def _sample_clusters(clusters, samples_per_cluster):
    result = []

    for cluster in clusters.values():
        cluster.sort(key=_sort_on_n_mismatches)
        sampled_cluster = _sample_cluster(cluster, samples_per_cluster)
        result += sampled_cluster
    return set(result)


# def make_snp_tables(clusters, seqs, df):
#     snps = {}
#     for target, hits in clusters.items():
#         idx = [get_id(h.query) for h in hits]
#         table = pd.DataFrame(
#             {'snp': [snp_count(h, seqs) for h in hits],
#              'location': [location(get_id(h.query), df) for h in hits],
#              'date': pd.to_datetime(df['date'][df.index.isin(idx)],
#                                     errors='coerce')},
#             index=idx)
#         table = table.sort_values(by=['snp', 'date'])
#         snps[get_id(target)] = table
#     return snps

# According to the vsearch 2.14.2 documentation, percent_id is defined as:
#  (matching columns) / (alignment length - terminal gaps)
def subsample_neighbors(focal_seqs: DNAFASTAFormat,
                        context_seqs: DNAFASTAFormat,
                        ids: qiime2.Metadata,
                        percent_id: float,
                        samples_per_cluster: int,
                        locale: str = None,
                        max_accepts: int = 10,
                        n_threads: int = 1) -> IDSelection:

    df = ids.to_dataframe()
    inclusion = pd.Series(False, index=df.index)

    with tempfile.NamedTemporaryFile() as vsearch_out_f:
        command = ['vsearch',
                   '--threads', str(n_threads),
                   '--usearch_global', str(focal_seqs),
                   '--id', str(percent_id),
                   '--db', str(context_seqs),
                   '--userout', vsearch_out_f.name,
                   '--qmask', 'none',
                   '--maxaccepts', str(max_accepts),
                   '--uc_allhits',
                   '--userfields', 'query+target+mism']
        run_command(command, verbose=False)

        vsearch_out = pd.read_csv(
            open(vsearch_out_f.name), sep='\t', na_values='*',
            names=['focal_id', 'context_id', 'n_mismatches'])

        clusters = _clusters_from_vsearch_out(vsearch_out, locale)
        context_seqs_to_keep = _sample_clusters(clusters, samples_per_cluster)
        inclusion[context_seqs_to_keep] = True

    return IDSelection(inclusion, ids, "subsample_neighbors")
