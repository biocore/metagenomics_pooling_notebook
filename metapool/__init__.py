#!/usr/bin/env python
from .prep import (preparations_for_run, parse_prep, generate_qiita_prep_file,
                   preparations_for_run_mapping_file, remove_qiita_id,
                   demux_pre_prep, pre_prep_needs_demuxing)
from .mp_strings import (SAMPLE_NAME_KEY, QIITA_ID_KEY, SAMPLES_DETAILS_KEY,
                         ORIG_NAME_KEY, PRIMARY_STUDY_KEY,
                         SECONDARY_STUDIES_KEY, PROJECT_SHORT_NAME_KEY,
                         PROJECT_FULL_NAME_KEY, CONTAINS_REPLICATES_KEY,
                         parse_project_name, get_short_name_and_id,
                         get_qiita_id_from_project_name)
from .sample_sheet import (SS_SAMPLE_ID_KEY,
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

from .controls import is_blank
from .metapool import (extract_stats_metadata, sum_lanes, compress_plates,
                       add_controls, generate_override_cycles_value,
                       TUBECODE_KEY)


__credits__ = ("https://github.com/biocore/metagenomics_pooling_notebook/"
               "graphs/contributors")

__all__ = ['IGMManifest', 'add_controls', 'assign_emp_index', 'autopool',
           'compress_plates', 'demux_pre_prep', 'demux_sample_sheet',
           'dilute_gDNA', 'extract_stats_metadata', 'find_threshold',
           'generate_qiita_prep_file', 'make_sample_sheet', 'parse_prep',
           'pre_prep_needs_demuxing', 'preparations_for_run',
           'preparations_for_run_mapping_file', 'remove_qiita_id',
           'parse_project_name', 'get_short_name_and_id',
           'get_qiita_id_from_project_name',
           'requires_dilution', 'run_counts', 'sample_sheet_to_dataframe',
           'sheet_needs_demuxing', 'sum_lanes', 'validate_plate_metadata',
           'MetagenomicSampleSheetv90', 'MetagenomicSampleSheetv100',
           'MetagenomicSampleSheetv101', 'AmpliconSampleSheet',
           'MetatranscriptomicSampleSheetv0', 'generate_override_cycles_value',
           'MetatranscriptomicSampleSheetv10', 'AbsQuantSampleSheetv10',
           # KLSampleSheet is needed for instance() calls.
           'KLSampleSheet', 'load_sample_sheet', 'is_blank',
           'SAMPLE_NAME_KEY', 'PRIMARY_STUDY_KEY', 'SECONDARY_STUDIES_KEY',
           'SAMPLES_DETAILS_KEY', 'SS_SAMPLE_ID_KEY', 'ORIG_NAME_KEY',
           'TUBECODE_KEY', 'QIITA_ID_KEY', 'PROJECT_SHORT_NAME_KEY',
           'PROJECT_FULL_NAME_KEY', 'CONTAINS_REPLICATES_KEY']

from . import _version

__version__ = _version.get_versions()['version']
