#!/usr/bin/env python
from .prep import (preparations_for_run, parse_prep, generate_qiita_prep_file,
                   preparations_for_run_mapping_file, remove_qiita_id,
                   demux_pre_prep, pre_prep_needs_demuxing)
from .sample_sheet import (SAMPLE_PROJECT_KEY, SAMPLES_KEY, ORIG_NAME_KEY,
                           DATA_SAMPLE_ID_KEY, PROJECT_SHORT_NAME_KEY,
                           QIITA_ID_KEY, CONTAINS_REPLICATES_KEY,
                           sample_sheet_to_dataframe, make_sample_sheet,
                           AmpliconSampleSheet, MetagenomicSampleSheetv90,
                           MetagenomicSampleSheetv100,
                           MetagenomicSampleSheetv101, AbsQuantSampleSheetv10,
                           MetatranscriptomicSampleSheetv0, demux_sample_sheet,
                           sheet_needs_demuxing, KLSampleSheet,
                           load_sample_sheet, MetatranscriptomicSampleSheetv10)
from .plate import (validate_plate_metadata, requires_dilution, dilute_gDNA,
                    autopool, find_threshold)
from .amplipool import assign_emp_index
from .igm import IGMManifest
from .count import run_counts
from .metapool import (extract_stats_metadata, sum_lanes,
                       compress_plates, add_controls)
from .controls import QIITA_SAMPLE_NAME_KEY, is_blank


__credits__ = ("https://github.com/biocore/metagenomics_pooling_notebook/"
               "graphs/contributors")

__all__ = ['IGMManifest', 'add_controls', 'assign_emp_index', 'autopool',
           'compress_plates', 'demux_pre_prep', 'demux_sample_sheet',
           'dilute_gDNA', 'extract_stats_metadata', 'find_threshold',
           'generate_qiita_prep_file', 'make_sample_sheet', 'parse_prep',
           'pre_prep_needs_demuxing', 'preparations_for_run',
           'preparations_for_run_mapping_file', 'remove_qiita_id',
           'requires_dilution', 'run_counts', 'sample_sheet_to_dataframe',
           'sheet_needs_demuxing', 'sum_lanes', 'validate_plate_metadata',
           'MetagenomicSampleSheetv90', 'MetagenomicSampleSheetv100',
           'MetagenomicSampleSheetv101', 'AmpliconSampleSheet',
           'MetatranscriptomicSampleSheetv0',
           'MetatranscriptomicSampleSheetv10', 'AbsQuantSampleSheetv10',
           # KLSampleSheet is needed for instance() calls.
           'KLSampleSheet', 'load_sample_sheet', 'is_blank',
           'QIITA_SAMPLE_NAME_KEY', 'SAMPLE_PROJECT_KEY', 'SAMPLES_KEY',
           'ORIG_NAME_KEY', 'DATA_SAMPLE_ID_KEY', 'PROJECT_SHORT_NAME_KEY',
           'QIITA_ID_KEY', 'CONTAINS_REPLICATES_KEY']

from . import _version

__version__ = _version.get_versions()['version']
