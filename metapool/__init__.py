#!/usr/bin/env python

from .prep import preparations_for_run, parse_prep
from .sample_sheet import (sample_sheet_to_dataframe, KLSampleSheet,
                           validate_and_scrub_sample_sheet, make_sample_sheet)
from .plate import validate_plate_metadata
from .amplipool import assign_emp_index
from .igm import IGMManifest

__credits__ = ("https://github.com/biocore/metagenomics_pooling_notebook/"
               "graphs/contributors")
__all__ = ['preparations_for_run', 'parse_prep',
           'sample_sheet_to_dataframe', 'KLSampleSheet',
           'validate_and_scrub_sample_sheet', 'make_sample_sheet',
           'parse_sample_sheet',
           'validate_plate_metadata',
           'assign_emp_index',
           'IGMManifest']


from . import _version
__version__ = _version.get_versions()['version']
