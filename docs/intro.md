# `genome-sampler` documentation

Thanks for your interest in `genome-sampler`!

This page is main entry point into the `genome-sampler` documentation. After
reading this documentation, you'll know how to install and use
`genome-sampler` and where to turn if you need help.

(getting-started)=
## Getting started with `genome-sampler`
Installation instructions for `genome-sampler` are available in our
{ref}`install`.

To learn how to use `genome-sampler` you can work through our
{ref}`usage-tutorial` after installing.


(help)=
## Getting help, contributing, and our community code of conduct

If you need technical support, please post a question to the
[QIIME 2 Forum](https://forum.qiime2.org).

```{note}
You'll get help more quickly on the [QIIME 2 Forum](https://forum.qiime2.org)
than if you email the developers directly, since many people are monitoring
the forum for support requests, while your message might get lost in an
individual developers email inbox for days or weeks if they're busy with
other projects.
```

We are very interested in contributions to
`genome-sampler` from the community. Please get in touch via the
[GitHub issue tracker](https://github.com/caporaso-lab/genome-sampler/issues)
or the [QIIME 2 Forum](https://forum.qiime2.org) if you’re interested in
contributing.

If `genome-sampler` is missing a feature that would be helpful for your work,
please post to the
[GitHub issue tracker](https://github.com/caporaso-lab/genome-sampler/issues).

Before getting in touch, please review the software project's
[code of
conduct](https://github.com/caporaso-lab/code-of-conduct/blob/master/code-of-conduct.md),
which is adapted from the
[Contributor
Covenant](https://www.contributor-covenant.org), version 1.4.


(license)=
## Licensing and source code
`genome-sampler` is open-source and free for all use. Software and unit tests
are available in our
[GitHub repository](https://github.com/caporaso-lab/genome-sampler) under the
BSD 3-clause license.

(about)=
## About `genome-sampler`
`genome-sampler` is a [QIIME 2](https://qiime2.org) plugin. QIIME 2 offers
useful features for bioinformatics software, including that it ensures
reproducibility of analyses and it is interface agnostic, meaning that the
same functionality can be accessed through a command line interface, a Python 3
API, and various graphical interfaces that are currently in different stages
of development. If you're interested in learning more about QIIME 2 and how it
can help with your bioinformatics software, read the
[QIIME 2 paper](https://pubmed.ncbi.nlm.nih.gov/31341288/)
{cite}`qiime2-nbt` and then the
[QIIME 2 developer documentation](https://dev.qiime2.org).

`genome-sampler`'s documentation is written using
[Myst](https://myst-parser.readthedocs.io/en/latest) and rendered using
[Jupyter Book](https://jupyterbook.org/).

`genome-sampler` is primarily developed at the
[Pathogen and Microbiome Institute](http://pmi.nau.edu) at
[Northern Arizona University](http://www.nau.edu).

(citing)=
## Citing `genome-sampler`
If you use `genome-sampler` in published work, we ask that you
cite [our paper](https://f1000research.com/articles/9-657)
{cite}`genome-sampler-f1000`.

The primary workflow implemented in `genome-sampler` also makes extensive use
of `vsearch` {cite}`vsearch-peerj`, so you should also cite that in your
published work.

If you use other components of QIIME 2 (as discussed in
{ref}`downstream-workflows`) you may end up using other tools that need to be
cited. If you load a `.qza` or `.qzv` file with
[QIIME 2 View](https://view.qiime2.org), you can obtain a list of the papers
you should cite under the _Details_ tab in Bibtex format. That file can be
loaded into most citation managers, such as
[Paperpile](https://paperpile.com/app) or EndNote.
