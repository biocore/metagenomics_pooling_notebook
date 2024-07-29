#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, Metapool development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup
from glob import glob
from os.path import dirname

import versioneer


classifiers = [
    'Development Status :: 5 - Production/Stable',
    'License :: OSI Approved :: MIT License',
    'Environment :: Console',
    'Topic :: Software Development :: Libraries :: Application Frameworks',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Operating System :: Unix',
    'Operating System :: POSIX',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows']


description = 'Metagenomic pooling Jupyter notebook helper'

with open('README.md') as f:
    long_description = f.read()

keywords = 'microbiome wetlab bioinformatics',

base = ['numpy', 'pandas', 'matplotlib >= 2.0', 'seaborn >= 0.7.1', 'click',
        'sample_sheet', 'openpyxl', 'qiita_client @ https://github.com/'
        'qiita-spots/qiita_client/archive/master.zip', 'scikit-learn',
        'biom-format']
test = ["nose", "pep8", "flake8"]
coverage = ['coverage']
notebook = ['jupyter', 'notebook', 'jupyter_contrib_nbextensions', 'watermark']
all_deps = base + test + coverage + notebook

notebooks_fp = []
for fp in glob('notebooks/*.ipynb'):
    notebooks_fp.append((dirname(fp), [fp]))

setup(name='metapool',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      license='MIT',
      description=description,
      long_description=long_description,
      keywords=keywords,
      classifiers=classifiers,
      author="Metapool development team.",
      maintainer="Metapool development team.",
      url='https://github.com/biocore/metagenomics_pooling_notebook',
      test_suite='nose.collector',
      packages=find_packages(),
      package_data={
          'metapool': ['data/*.tsv', 'data/*.xlsx', 'tests/data/*.csv']},
      # addingt all the notebooks fps
      data_files=notebooks_fp,
      install_requires=base,
      extras_require={'test': test,
                      'coverage': coverage,
                      'all': all_deps}
      )
