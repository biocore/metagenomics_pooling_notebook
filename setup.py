#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, calour development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re
import ast
from setuptools import find_packages, setup


# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('metapool/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

classifiers = [
    'Development Status :: 2 - Pre-Alpha',
    'License :: OSI Approved :: MIT License',
    'Environment :: Console',
    'Topic :: Software Development :: Libraries :: Application Frameworks',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Operating System :: Unix',
    'Operating System :: POSIX',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows']


description = 'Metagenomics pooling Jupyter notebook helper'

with open('README.md') as f:
    long_description = f.read()

keywords = 'microbiome wetlab bioinformatics',

setup(name='metapool',
      version=version,
      license='MIT',
      description=description,
      long_description=long_description,
      keywords=keywords,
      classifiers=classifiers,
      author="Jon Sanders",
      maintainer="Jon Sanders",
      url='https://github.com/tanaes/metagenomics_pooling_notebook',
      test_suite='nose.collector',
      packages=find_packages(),
      package_data={
        'metapool': ['tests/data/*.csv']},
      install_requires=[
          'numpy',
          'pandas',
          'matplotlib >= 2.0',
          'jupyter',
          'notebook',
          'jupyter_contrib_nbextensions',
          'seaborn >= 0.7.1',
          'sample_sheet'],
      extras_require={'test': ["nose", "pep8", "flake8"],
                      'coverage': ["coverage"]},
      entry_points={
          'console_scripts': [
              'calour=calour.cli:cmd',
          ]})
