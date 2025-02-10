# ----------------------------------------------------------------------------
# Copyright (c) 2020-2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources

import pandas as pd
import q2templates


from genome_sampler.common import IDSelection

SUMMARY_TEMPLATE = pkg_resources.resource_filename(
    'genome_sampler', 'assets/summarize/index.html')


def _build_summary_table(selections):
    rows = []
    for selection in selections:
        rows.append({'action': selection.label,
                     'included': selection.inclusion.sum(),
                     'total': len(selection.inclusion)})
    table = pd.DataFrame(rows, columns=['action', 'included', 'total'])
    return table


def summarize_selections(output_dir: str, selections: IDSelection):
    table = _build_summary_table(selections)
    html_table = q2templates.df_to_html(table, index=False)

    table_fn = 'table.tsv'
    # Not using qiime2.Metadata b/c we don't have a meaningful ID col
    table.to_csv(
        os.path.join(output_dir, table_fn), sep='\t', encoding='utf-8')

    context = {
        'table': html_table,
        'table_fn': table_fn,
    }

    q2templates.render(SUMMARY_TEMPLATE, output_dir, context=context)
