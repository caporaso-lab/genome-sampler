import os
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


def _make_longitudinal_plot(dates):
    # TODO: return an altair plot for this.
    pass

def _merge_time_to_dates(selections):
    to_merge = []
    for selection in selections:
        if selection.label == 'sample_longitudinal':
            to_merge.append(selection)

    if not to_merge:
        return None



def summarize_selection(output_dir: str, selections: IDSelection):
    context = {}

    table = _build_summary_table(selections)
    html_table = q2templates.df_to_html(table, index=False)
    context['table'] = html_table

    dates = _merge_time_to_dates(selections)
    if dates is not None:
        longitudinal_plot = _make_longitudinal_plot(dates)
        longitudinal_plot.save(os.path.join(output_dir, 'longitudinal.html'))
        context['longitudinal_plot'] = True

    q2templates.render(SUMMARY_TEMPLATE, output_dir, context=context)
