import os
import pkg_resources

import pandas as pd
import q2templates


from genome_sampler.common import IDSelection

SUMMARY_TEMPLATE = pkg_resources.resource_filename(
    'genome_sampler', 'assets/summarize/index.html')


def _build_summary_table(selections):
    table = pd.DataFrame({}, columns=['action', 'included', 'total'])
    for selection in sorted(selections, key=lambda x: x.label):
        table = table.append({'action': selection.label,
                              'included': selection.inclusion.sum(),
                              'total': len(selection.inclusion)},
                             ignore_index=True)

    return table


def summarize_selection(output_dir: str, selections: IDSelection):
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
