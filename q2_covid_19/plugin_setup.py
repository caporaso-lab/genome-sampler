from qiime2.plugin import Plugin

import q2_covid_19

plugin = Plugin(
    name='covid-19',
    website='https://github.com/caporaso-lab/q2-covid-19',
    package='q2_covid_19',
    version=q2_covid_19.__version__,
    description='TODO',
    short_description='TODO'
)
