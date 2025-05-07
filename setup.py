#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, metapool development team.
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
    'Development Status :: 2 - Pre-Alpha',
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

keywords = 'microbiome wetlab bioinformatics'

base = ['qiita_client @ https://github.com/qiita-spots/qiita_client/archive/'
        'master.zip', 'biom-format >= 2.1.16', 'click >= 8.1.7',
        'matplotlib >= 3.9.2', 'numpy >= 2.0.2', 'openpyxl >= 3.1.5',
        'pandas >= 2.2.3', 'sample-sheet >= 0.13.0', 'scikit-learn >= 1.5.2',
        'seaborn >= 0.13.2']

test = ['nose >= 1.3.7', 'pep8 >= 1.7.1', 'flake8 >= 7.1.1']

coverage = ['coverage >= 7.6.8']

notebook = ['jupyter >= 1.1.1', 'jupyter_contrib_nbextensions >= 0.7.0',
            'notebook >= 6.5.7', 'watermark >= 2.5.0']

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
      author="Jon Sanders",
      maintainer="Amanda Birmingham",
      url='https://github.com/biocore/kl_metapool',
      test_suite='nose.collector',
      packages=find_packages(),
      package_data={
          'metapool': ['data/*.tsv', 'data/*.xlsx', 'tests/data/*.csv']},
      # adding all the notebooks fps
      data_files=notebooks_fp,
      install_requires=base,
      extras_require={'test': test,
                      'coverage': coverage,
                      'all': all_deps},
      entry_points={
          'console_scripts': [
              'seqpro=metapool.scripts.seqpro:format_preparation_files',
              ('seqpro_mf=metapool.scripts.seqpro_mf:format_preparation_'
               'files_mf'),
          ],

      })
