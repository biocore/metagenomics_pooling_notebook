#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, metapool development team.
#
# Distributed under the terms of the Modified BSD License.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages
from glob import glob
from os.path import dirname
import versioneer

# NOTE: if you change this section,
# ALSO change the environment.yml file
base = [
    'biom-format >= 2.1.16',
    'matplotlib >= 3.9.2',
    'numpy >= 2.0.2',
    'openpyxl >= 3.1.5',
    'pandas >= 3.0',
    # seems seqtk cannot be installed through the setup mechanism?
    # 'seqtk >= 1.4',
    'seaborn >= 0.13.2',
    'scikit-learn >= 1.5.2',
    # below can't be installed by conda
    'sample-sheet >= 0.13.0',
    # git-only; pinned to the commit Qiita-MIINT consumes, so the pre-flight
    # this notebook writes matches the version that reads it
    'run-preflight @ git+https://github.com/the-miint/kl-run-preflight'
    '@03ee9e3b86132b97d42298380a43d9e52eaaad23']

test = [
    'flake8 >= 7.1.1',
    'pytest >= 7.0.0',
    'pytest-cov >= 4.0.0',
    'papermill >= 2.6.0',
    'pep8 >= 1.7.1']

coverage = [
    'coverage >= 7.6.8',
    'coveralls']

notebook = [
    'jupyter >= 1.1.1',
    'notebook >= 6.5.7',
    'watermark >= 2.5.0',
    'ipyfilechooser']

all_deps = base + test + coverage + notebook

# collect notebook data files
notebooks_fp = []
for fp in glob('notebooks/*.ipynb'):
    notebooks_fp.append((dirname(fp), [fp]))

setup(
    name='metapool',
    python_requires='>=3.11',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    package_data={
        'metapool': ['data/*.tsv', 'data/*.xlsx', 'tests/data/*.csv']},
    include_package_data=True,
    # adding all the notebooks fps
    data_files=notebooks_fp,
    install_requires=base,
    extras_require={'test': test,
                    'coverage': coverage,
                    'all': all_deps},
    entry_points={
        'console_scripts': [
            'seqpro=metapool.scripts.seqpro:format_preparation_files',
            ('seqpro_mf=metapool.scripts.seqpro_mf:'
                'format_preparation_files_mf')],
    })
