import pkg_resources

import pandas as pd
import q2templates


from genome_sampler.common import IDSelection

SUMMARY_TEMPLATE = pkg_resources.resource_filename(
    'genome_sampler', 'assets/index.html')


def _build_summary_table(selections):
    table = pd.DataFrame({}, columns=['action', 'included', 'total'])
    for selection in selections:
        table = table.append({'action': selection.label,
                              'included': selection.inclusion.sum(),
                              'total': len(selection.inclusion)},
                             ignore_index=True)

    return table


def summarize_selection(output_dir: str, selections: IDSelection):
    context = {}

    table = _build_summary_table(selections)
    html_table = q2templates.df_to_html(table, index=False)
    context['table'] = html_table

    q2templates.render(SUMMARY_TEMPLATE, output_dir, context=context)
