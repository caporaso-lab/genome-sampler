import tempfile

import pandas as pd
import numpy as np

import qiime2
from q2_types.feature_data import DNAFASTAFormat

from genome_sampler.common import IDSelection, run_command, ids_from_fasta


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
    columns = ['context_id', 'n_mismatches', 'locale']
    return {focal_id: pd.DataFrame(hits, columns=columns)
            for focal_id, hits in clusters.items()}


def _sample_cluster(cluster, samples_per_cluster, random_state):
    len_cluster = cluster.shape[0]
    if len_cluster == 0:
        return []

    if len_cluster <= samples_per_cluster:
        return list(cluster['context_id'])

    unique_locales = cluster['locale'].unique()
    if len(unique_locales) == 1:
        # if only one locale, select context sequences spanning the range
        # of number of mismatches
        sorted_cluster = cluster.sort_values(by='n_mismatches')
        indices = map(int, np.linspace(0, len_cluster-1, samples_per_cluster))
        return [sorted_cluster.iloc[i]['context_id'] for i in indices]

    # convert locations into number of obs, take the inverse, then scale by
    # number of unique locations so it all sums to 1. Infrequent location obs
    # will have a high weight, frequent location obs will have an individually
    # lower weight, but collectively all locations have the same weight when
    # all obs are summed within a location class.
    loc = cluster['locale']
    weights = ((1 / loc.map(loc.value_counts())) / len(loc.unique()))
    # TODO: if multiple samples are being selected from the same locale, it
    # would be good to sample across the number of mismatches to the focal
    # sequence, as we're doing where there is only a single locale. To do that,
    # I think we'd need to get the count of how many of each locale we want
    # at random (like the rarefaction algorithm) and then use the approach
    # above (which would then generalize to one or more locales).
    result = cluster.sample(samples_per_cluster,
                            weights=weights,
                            random_state=random_state)
    return list(result['context_id'])


def _sample_clusters(clusters, samples_per_cluster, seed):
    result = []
    random_state = np.random.RandomState(seed=seed)

    for cluster in clusters.values():
        sampled_cluster = _sample_cluster(
                cluster, samples_per_cluster, random_state)
        result += sampled_cluster
    return set(result)


# According to the vsearch 2.14.2 documentation, percent_id is defined as:
#  (matching columns) / (alignment length - terminal gaps)
def subsample_neighbors(focal_seqs: DNAFASTAFormat,
                        context_seqs: DNAFASTAFormat,
                        ids: qiime2.Metadata,
                        percent_id: float,
                        samples_per_cluster: int,
                        locale: str = None,
                        max_accepts: int = 10,
                        n_threads: int = 1,
                        seed: int = None) -> IDSelection:

    if max_accepts < samples_per_cluster:
        raise ValueError('max_accepts (%d) must be greater than or equal to '
                         'samples_per_cluster (%d), since it is determines '
                         'the largest number of samples that could be '
                         'obtained per cluster.' %
                         (max_accepts, samples_per_cluster))

    context_ids = ids_from_fasta(str(context_seqs))

    inclusion = pd.Series(False, index=context_ids)
    ids = ids.filter_ids(inclusion.index)
    if locale is not None:
        locale = ids.get_column(locale).to_series()

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
        run_command(command)

        vsearch_out = pd.read_csv(
            open(vsearch_out_f.name), sep='\t', na_values='*',
            names=['focal_id', 'context_id', 'n_mismatches'])

        clusters = _clusters_from_vsearch_out(vsearch_out, locale)
        context_seqs_to_keep = \
            _sample_clusters(clusters, samples_per_cluster, seed=seed)
        inclusion[context_seqs_to_keep] = True

    return IDSelection(inclusion, ids, "subsample_neighbors")
