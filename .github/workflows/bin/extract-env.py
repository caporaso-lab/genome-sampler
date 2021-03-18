#!/usr/bin/env python

import os
import json
import argparse

import yaml


def get_conda_prefix():
    conda_prefix = os.getenv('CONDA_PREFIX')
    if conda_prefix is None:
        raise Exception("Not in a conda environment.")

    return conda_prefix


class CondaMeta:
    def __init__(self, prefix):
        self.prefix = prefix
        self.meta = os.path.join(self.prefix, 'conda-meta')
        self._cache = {}

        self.meta_lookup = {}
        for filename in os.listdir(self.meta):
            if filename.endswith('.json'):
                name = filename.rsplit('-', 2)[0]
                self.meta_lookup[name] = os.path.join(self.meta, filename)

    def __getitem__(self, package):
        if package not in self._cache:
            try:
                with open(self.meta_lookup[package]) as fh:
                    self._cache[package] = json.load(fh)
            except KeyError:
                raise Exception(
                    "Package %r not found in current environment." % package)

        return self._cache[package]

    def iter_primary_deps(self, package):
        for dep in self[package]['depends']:
            dep = dep.split(' ')[0]
            # https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-virtual.html
            if dep in ('__cuda', '__osx', '__glibc', '__unix', '__win'):
                continue
            else:
                yield dep

    def iter_deps(self, package, *, include_self=True, _seen=None):
        if _seen is None:
            _seen = set()

        if include_self:
            yield package

        for dependency in self.iter_primary_deps(package):
            if dependency in _seen:
                continue
            else:
                _seen.add(dependency)
                yield from self.iter_deps(dependency, _seen=_seen)

    def get_version(self, package):
        return self[package]['version']


def write_env_file(env, fh):
    import qiime2

    channels = {
        'channels': [
            'https://packages.qiime2.org/qiime2/latest/tested/',
            'qiime2/label/r%s' % qiime2.__release__,
            'conda-forge',
            'bioconda',
            'defaults']}
    yaml.dump(channels, fh, default_flow_style=False)

    dependencies = {'dependencies': ['='.join(r) for r in sorted(env.items())]}
    yaml.dump(dependencies, fh, default_flow_style=False)


def main(packages):
    prefix = get_conda_prefix()
    meta = CondaMeta(prefix)

    env = {}
    for package in packages:
        for d in meta.iter_deps(package, include_self=True):
            env[d] = meta.get_version(d)

    return env


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract environment from an existing conda environment')
    parser.add_argument('packages', metavar='PACKAGE', type=str, nargs='+',
                        help='Packages to extract.')
    parser.add_argument('--env-file', dest='environment_file', action='store',
                        default=None, help='Path to write an env file to.')
    parser.add_argument('--json', dest='json', action='store_true',
                        help='Write as json')

    args = parser.parse_args()

    env = main(args.packages)

    if args.environment_file is not None:
        with open(args.environment_file, 'w') as fh:
            write_env_file(env, fh)
    elif args.json:
        print(json.dumps(env))
    else:
        for key, value in env.items():
            print("%s\t%s" % (key, value))
