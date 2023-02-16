#!/usr/bin/env python
from .prep import (preparations_for_run, parse_prep, generate_qiita_prep_file,
                   preparations_for_run_mapping_file)
from .sample_sheet import (sample_sheet_to_dataframe, KLSampleSheet,
                           validate_and_scrub_sample_sheet, make_sample_sheet,
                           quiet_validate_and_scrub_sample_sheet)
from .plate import (validate_plate_metadata, requires_dilution, dilute_gDNA,
                    autopool, find_threshold)
from .amplipool import assign_emp_index
from .igm import IGMManifest
from .count import run_counts
from .metapool import extract_stats_metadata, sum_lanes

__credits__ = ("https://github.com/biocore/metagenomics_pooling_notebook/"
               "graphs/contributors")
__all__ = ['preparations_for_run', 'parse_prep', 'generate_qiita_prep_file',
           'sample_sheet_to_dataframe', 'KLSampleSheet',
           'validate_and_scrub_sample_sheet', 'make_sample_sheet',
           'quiet_validate_and_scrub_sample_sheet', 'validate_plate_metadata',
           'assign_emp_index', 'IGMManifest', 'run_counts',
           'requires_dilution', 'dilute_gDNA', 'autopool', 'find_threshold',
           'extract_stats_metadata', 'sum_lanes',
           'preparations_for_run_mapping_file']

from . import _version

__version__ = _version.get_versions()['version']
