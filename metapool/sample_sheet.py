import csv
import collections
import warnings
from datetime import datetime
from itertools import chain, repeat, islice
from json import loads as json_loads
import sample_sheet
import pandas as pd
from types import MappingProxyType
from metapool.mp_strings import parse_project_name, \
    get_qiita_id_from_project_name, \
    SAMPLES_DETAILS_KEY, SAMPLE_NAME_KEY, SAMPLE_PROJECT_KEY, \
    CONTAINS_REPLICATES_KEY, ORIG_NAME_KEY, EXPT_DESIGN_DESC_KEY, \
    PM_PROJECT_NAME_KEY, PM_PROJECT_PLATE_KEY, PM_BLANK_KEY, QIITA_ID_KEY, \
    PROJECT_FULL_NAME_KEY
from metapool.metapool import (bcl_scrub_name, sequencer_i5_index,
                               REVCOMP_SEQUENCERS, TUBECODE_KEY)
from metapool.plate import ErrorMessage, WarningMessage, PlateReplication
from metapool.controls import SAMPLE_CONTEXT_COLS, \
    get_all_projects_in_context, is_blank, get_controls_details_from_context, \
    get_delimited_controls_details_from_compressed_plate, \
    make_manual_control_details, denormalize_controls_details

_BIOINFORMATICS_KEY = 'Bioinformatics'
_CONTACT_KEY = 'Contact'
_SAMPLE_CONTEXT_KEY = 'SampleContext'
_HEADER_KEY = 'Header'
_READS_KEY = 'Reads'
_SETTINGS_KEY = 'Settings'
_DATA_KEY = 'Data'
_ASSAY_KEY = 'Assay'
_SS_SAMPLE_PROJECT_KEY = 'Sample_Project'
_SS_QIITA_ID_KEY = 'QiitaID'
_SS_SAMPLE_NAME_KEY = 'Sample_Name'
_EXPERIMENT_NAME_KEY = 'Experiment Name'
_EMAIL_KEY = 'Email'
_HUMAN_FILTERING_KEY = 'HumanFiltering'
_SHEET_TYPE_KEY = 'SheetType'
_SHEET_VERSION_KEY = 'SheetVersion'

STANDARD_METAG_SHEET_TYPE = 'standard_metag'
STANDARD_METAT_SHEET_TYPE = 'standard_metat'
ABSQUANT_SHEET_TYPE = 'abs_quant_metag'
_AMPLICON = 'TruSeq HT'
_DUMMY_SHEET_TYPE = 'dummy_amp'
_PLACEHOLDER_EXPT_NAME = 'RKL_experiment'

# NB: These two string literals are assigned to (totally independent) public
# constants in the qp_klp Step module, and then used throughout qp_klp.  I
# don't know if those literals and these come from the same place (in
# which case they should share a single definition), or whether they just
# happen coincidentally to be the same, in which case they probably shouldn't.
_METAGENOMIC = 'Metagenomic'
_METATRANSCRIPTOMIC = 'Metatranscriptomic'

SS_SAMPLE_ID_KEY = 'Sample_ID'


# MappingProxyType is used to make the dictionary immutable.
# NOTE that key order is (a) important and (b) preserved; as of Python 3.7,
# insertion order is preserved in dictionaries (cite:
# https://mail.python.org/pipermail/python-dev/2017-December/151283.html ).
_BASE_BIOINFORMATICS_COLS = MappingProxyType(
    {_SS_SAMPLE_PROJECT_KEY: str,
     _SS_QIITA_ID_KEY: str,
     'BarcodesAreRC': bool,
     'ForwardAdapter': str,
     'ReverseAdapter': str,
     _HUMAN_FILTERING_KEY: bool,
     'library_construction_protocol': str,
     EXPT_DESIGN_DESC_KEY: str
     })
_BIOINFORMATICS_COLS_W_REP_SUPPORT = MappingProxyType(
    _BASE_BIOINFORMATICS_COLS | {CONTAINS_REPLICATES_KEY: bool})
_CONTACT_COLS = MappingProxyType({
    _SS_SAMPLE_PROJECT_KEY: str,
    _EMAIL_KEY: str})

# Note that there doesn't appear to be a difference between 95, 99, and 100
# beyond the value observed in 'Well_description' column. The real
# difference is between standard_metag and abs_quant_metag.
_BASE_DATA_COLUMNS = (SS_SAMPLE_ID_KEY, _SS_SAMPLE_NAME_KEY, 'Sample_Plate',
                      'well_id_384', 'I7_Index_ID', 'index', 'I5_Index_ID',
                      'index2', _SS_SAMPLE_PROJECT_KEY, 'Well_description')

_BASE_CARRIED_PREP_COLUMNS = (EXPT_DESIGN_DESC_KEY, 'i5_index_id',
                              'i7_index_id', 'index', 'index2',
                              'library_construction_protocol',
                              SAMPLE_NAME_KEY, 'sample_plate',
                              'sample_project', 'well_description',
                              'well_id_384')

_ELUTION_VOL_KEY = 'vol_extracted_elution_ul'

_BASE_METAG_REMAPPER = MappingProxyType({
            'sample sheet Sample_ID': SS_SAMPLE_ID_KEY,
            'Sample': _SS_SAMPLE_NAME_KEY,
            PM_PROJECT_PLATE_KEY: 'Sample_Plate',
            'Well': 'well_id_384',
            'i7 name': 'I7_Index_ID',
            'i7 sequence': 'index',
            'i5 name': 'I5_Index_ID',
            'i5 sequence': 'index2',
            PM_PROJECT_NAME_KEY: _SS_SAMPLE_PROJECT_KEY,
        })


class KLSampleSheet(sample_sheet.SampleSheet):
    _ASSAYS = frozenset({_AMPLICON, _METAGENOMIC, _METATRANSCRIPTOMIC})
    _BIOINFORMATICS_AND_CONTACT = MappingProxyType({
        _BIOINFORMATICS_KEY: None,
        _CONTACT_KEY: None
    })

    _KL_ADDTL_DF_SECTIONS = MappingProxyType({
        _BIOINFORMATICS_KEY: _BASE_BIOINFORMATICS_COLS,
        _CONTACT_KEY: _CONTACT_COLS,
    })

    _HEADER = MappingProxyType({
        'IEMFileVersion': '4',
        _SHEET_TYPE_KEY: None,
        _SHEET_VERSION_KEY: None,
        'Investigator Name': 'Knight',
        _EXPERIMENT_NAME_KEY: _PLACEHOLDER_EXPT_NAME,
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        _ASSAY_KEY: None,
        'Description': '',
        'Chemistry': 'Default',
    })

    _READS = MappingProxyType({
        'Read1': 151,
        'Read2': 151
    })

    _SETTINGS = MappingProxyType({
        'ReverseComplement': '0',
        'MaskShortReads': '1',
        'OverrideCycles': 'Y151;I8N2;I8N2;Y151'
    })

    _ALL_METADATA = MappingProxyType({
        **_HEADER, **_SETTINGS, **_READS, **_BIOINFORMATICS_AND_CONTACT})

    # If modifying, see issue #233
    sections = (_HEADER_KEY, _READS_KEY, _SETTINGS_KEY, _DATA_KEY,
                _BIOINFORMATICS_KEY, _CONTACT_KEY)

    _ORDERED_BY_DATA_COLUMNS = False

    # NB: Inside `make_sample_sheet`, the 'Well_description' column is
    # (over)written by concatenating project plate, sample, and well.
    # So it is required in the sense that it has to be present in the sample
    # sheet, but it is not required in the sense that it has to be provided
    # by a user when they create a sample sheet through this module.
    _data_columns = (SS_SAMPLE_ID_KEY, _SS_SAMPLE_NAME_KEY, 'Sample_Plate',
                     'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID',
                     'index2', _SS_SAMPLE_PROJECT_KEY, 'Well_description')

    _column_alts = MappingProxyType({'well_description': 'Well_description',
                                     'description': 'Well_description',
                                     'Description': 'Well_description',
                                     'sample_plate': 'Sample_Plate'})

    _CARRIED_PREP_COLUMNS = (EXPT_DESIGN_DESC_KEY, 'i5_index_id',
                             'i7_index_id', 'index', 'index2',
                             'library_construction_protocol',
                             SAMPLE_NAME_KEY,
                             'sample_plate', 'sample_project',
                             'well_description', 'Sample_Well', 'Lane')

    _GENERATED_PREP_COLUMNS = ('center_name', 'center_project_name',
                               'instrument_model', 'lane', 'platform',
                               'run_center', 'run_date', 'run_prefix', 'runid',
                               'sequencing_meth')

    @property
    def CARRIED_PREP_COLUMNS(self):
        return list(self._CARRIED_PREP_COLUMNS)

    @property
    def GENERATED_PREP_COLUMNS(self):
        return list(self._GENERATED_PREP_COLUMNS)

    def __new__(cls, path=None, *args, **kwargs):
        """
            Override so that base class cannot be instantiated.
        """
        if cls is KLSampleSheet:
            raise TypeError(
                f"only children of '{cls.__name__}' may be instantiated")

        instance = super(KLSampleSheet, cls).__new__(cls, *args, **kwargs)
        return instance

    def __init__(self, path=None, *args, **kwargs):
        """Knight Lab's SampleSheet subclass

        Includes a number of parsing and writing changes to allow for the
        inclusion of the Bioinformatics and Contact sections.

        The majority of the code in the _parse and write methods are copied
        from release 0.11.0 https://tinyurl.com/clintval-sample-sheet

        Parameters
        ----------
        path: str, optional
            File path to the sample sheet to load.

        Attributes
        ----------
        Header: sample_sheet.Section
            Header section of the sample sheet.
        Reads: sample_sheet.Section
            Reads section of the sample sheet
        Settings: sample_sheet.Section
            Settings section of the sample sheet.
        Bioinformatics: pd.DataFrame
            Bioinformatics section of the sample sheet.
        Contact: pd.DataFrame
            Contact section of the sample sheet.
        path: str
            File path where the data was parsed from.
        """
        # don't pass the path argument to avoid the superclass from parsing
        # the data.
        super().__init__()
        self._remapper = None

        self.Bioinformatics = None
        self.Contact = None
        self.path = path

        if self.path:
            self._parse(self.path)

            # Convert the boolean parameters from strings to booleans.
            # Ignore any messages returned because we are not validating,
            # just converting datatypes.
            self._normalize_df_sections_booleans()

    def _parse(self, path):
        section_name = ''
        section_header = None

        def _is_empty(csv_line):
            return not ''.join(csv_line).strip()

        with open(path, encoding=self._encoding) as handle:
            lines = list(csv.reader(handle, skipinitialspace=True))

            # Comments at the top of the file are no longer supported. Only
            # handle comments if they are contiguous and in the first block of
            # lines. Otherwise if comments are found halfway through the file
            # that will result in an error.
            #
            # Empty lines are ignored
            show_warning = False
            while len(lines) and (_is_empty(lines[0]) or
                                  lines[0][0].startswith('# ')):
                lines.pop(0)
                show_warning = True
            if show_warning:
                message = (f"Comments at the beginning of the sample sheet "
                           f"are no longer supported. This information will "
                           f"be ignored. Please use the {_CONTACT_KEY} "
                           f"section instead")
                warnings.warn(message)

            for i, line in enumerate(lines):
                # Skip to next line if this line is empty to support formats of
                # sample sheets with multiple newlines as section seperators.
                #
                #   https://github.com/clintval/sample-sheet/issues/46
                #
                if _is_empty(line):
                    continue

                # Raise exception if we encounter invalid characters.
                if any(
                    character not in sample_sheet.VALID_ASCII
                    for character in set(''.join(line))
                ):
                    raise ValueError(
                        f'Sample sheet contains invalid characters on line '
                        f'{i + 1}: {"".join(line)}'
                    )

                header_match = self._section_header_re.match(line[0])

                # If we enter a section save its name and continue to next
                # line.
                if header_match:
                    section_name, *_ = header_match.groups()
                    if (
                        section_name not in self._sections
                        and section_name not in type(self).sections
                    ):
                        self.add_section(section_name)

                    if section_header is not None:
                        section_header = None
                    continue

                # [Reads] - vertical list of integers.
                if section_name == _READS_KEY:
                    self.Reads.append(int(line[0]))
                    continue

                # [Data] - delimited data with the first line a header.
                elif section_name == _DATA_KEY:
                    if section_header is not None:
                        self.add_sample(
                            sample_sheet.Sample(dict(zip(section_header,
                                                         line))))
                    elif any(key == '' for key in line):
                        raise ValueError(
                            f'Header for [{_DATA_KEY}] section is not allowed '
                            f'to have empty fields: {line}'
                        )
                    else:
                        section_header = self._process_section_header(line)
                    continue

                elif section_name in self._KL_ADDTL_DF_SECTIONS:
                    if getattr(self, section_name) is not None:
                        # vals beyond the header are empty values so don't add
                        # them
                        line = line[:len(getattr(self, section_name).columns)]
                        df = getattr(self, section_name)
                        df.loc[len(df)] = line
                    else:
                        # CSV rows are padded to include commas for the longest
                        # line in the file, so we remove them to avoid creating
                        # empty columns
                        col_names = [value for value in line if value != '']
                        setattr(self, section_name,
                                pd.DataFrame(columns=col_names))
                    continue

                # [<Other>]
                else:
                    key, value, *_ = line
                    section = getattr(self, section_name)
                    section[key] = value
                    continue

    def _extend_remapper(self, addtl_remapper):
        curr_remapper = getattr(self, '_remapper', {})
        return MappingProxyType(curr_remapper | addtl_remapper)

    def set_override_cycles(self, value):
        # assume that any value including None is valid.
        # None should be silently converted to empty string as a truly
        # empty value would be better than None or na.
        if value is None:
            value = ''

        self.Settings['OverrideCycles'] = value

    def _process_section_header(self, columns):
        for i in range(0, len(columns)):
            if columns[i] in KLSampleSheet._column_alts:
                # overwrite existing alternate name w/internal representation.
                columns[i] = KLSampleSheet._column_alts[columns[i]]
        return columns

    def write(self, handle, blank_lines=1) -> None:
        """Write to a file-like object.

        Parameters
        ----------
        handle: file_like
            Object with a write method.
        blank_lines: int
            Number of blank lines to write between sections
        """
        if not isinstance(blank_lines, int) or blank_lines <= 0:
            raise ValueError('Number of blank lines must be a positive int.')

        def _get_section_len(section_name):
            section_df = getattr(self, section_name)
            if section_df is None:
                return 0
            return len(section_df.columns)

        section_lens = list(map(_get_section_len, self._KL_ADDTL_DF_SECTIONS))
        # NB: the added [2] ensures the list always contains at least the value
        # 2, which is the default number of columns if all others are zero.
        all_lens = [len(self.all_sample_keys)] + section_lens + [2]
        csv_width = max(all_lens)
        writer = csv.writer(handle)

        def pad_iterable(iterable, size, padding=''):
            return list(islice(chain(iterable, repeat(padding)), size))

        def write_blank_lines(a_writer, n=blank_lines, width=None):
            if width is None:
                width = csv_width
            for i in range(n):
                a_writer.writerow(pad_iterable([], width))

        for title in self.sections:
            writer.writerow(pad_iterable([f'[{title}]'], csv_width))

            # Data is not a section in this class
            if title != _DATA_KEY:
                section = getattr(self, title)

            if title == _READS_KEY:
                for read in self.Reads:
                    writer.writerow(pad_iterable([read], csv_width))
            elif title == _DATA_KEY:
                if self._ORDERED_BY_DATA_COLUMNS:
                    # order according to the expected column order.  If there
                    # are columns other than the expected ones, put them at the
                    # end in alphabetical order.
                    extra_cols = set(self.all_sample_keys) - \
                                 set(self._get_expected_data_columns())
                    expected_cols = set(self.all_sample_keys) - extra_cols
                    col_order = \
                        [x for x in self._get_expected_data_columns()
                         if x in expected_cols] + sorted(list(extra_cols))
                else:
                    # legacy behavior
                    col_order = self.all_sample_keys

                writer.writerow(pad_iterable(col_order, csv_width))

                for sample in self.samples:
                    line = [getattr(sample, k) for k in col_order]
                    writer.writerow(pad_iterable(line, csv_width))

            elif title in self._KL_ADDTL_DF_SECTIONS:
                if section is not None:
                    # these sections are represented as DataFrame objects
                    writer.writerow(pad_iterable(section.columns.tolist(),
                                                 csv_width))

                    for _, row in section.iterrows():
                        writer.writerow(pad_iterable(row.values.tolist(),
                                                     csv_width))

            else:
                for key, value in section.items():
                    writer.writerow(pad_iterable([key, value], csv_width))
            write_blank_lines(writer)

    def merge(self, sheets):
        """Merge the Data section of multiple sample sheets

        For the Date field in the Header section, we only keep the date of the
        current sheet.

        Parameters
        ----------
        sheets: list of KLSampleSheet
            The sample sheets to merge into `self`.

        Raises
        ------
        ValueError
            If the Header, Settings or Reads section is different between
            merged sheets.
        """
        for number, sheet in enumerate(sheets):
            for section in [_HEADER_KEY, _SETTINGS_KEY, _READS_KEY]:
                this, that = getattr(self, section), getattr(sheet, section)

                # For the Header section we'll ignore the Date field since that
                # is likely to be different but shouldn't be a condition to
                # prevent merging two sheets.
                if section == _HEADER_KEY:
                    if this is not None:
                        this = {k: v for k, v in this.items() if k != 'Date'}
                    if that is not None:
                        that = {k: v for k, v in that.items() if k != 'Date'}

                if this != that:
                    raise ValueError(('The %s section is different for sample '
                                     'sheet %d ') % (section, 1 + number))

            for sample in sheet.samples:
                self.add_sample(sample)

            # these sections are data frames
            for section in self._KL_ADDTL_DF_SECTIONS:
                this, that = getattr(self, section), getattr(sheet, section)

                # if both frames are not None then we concatenate the rows.
                if this is not None and that is not None:
                    for _, row in that.iterrows():
                        this.loc[len(this)] = row.copy()

                    this.drop_duplicates(keep='first', ignore_index=True,
                                         inplace=True)

                # if self's frame is None then assign a copy
                elif this is None and that is not None:
                    setattr(self, section, that.copy())
                # means that either self's is the only frame that's not None,
                # so we don't need to merge anything OR that both frames are
                # None so we have nothing to merge.
                else:
                    pass

    def _remap_table(self, table, strict):
        result = table.copy(deep=True)

        if strict:
            # All columns not defined in _remapper will be filtered result.
            result = table[self._remapper.keys()].copy()
            result.rename(self._remapper, axis=1, inplace=True)
        else:
            # if a column named 'index' is present in table, assume it is a
            # numeric index and not a sequence of bases, which is required in
            # the output. Assume the column that will become 'index' is
            # defined in _remapper.
            if 'index' in set(result.columns):
                result.drop(columns=['index'], inplace=True)

            _remapper = KLSampleSheet._column_alts | self._remapper
            result.rename(_remapper, axis=1, inplace=True)

            # result may contain additional columns that aren't allowed in the
            # [Data] section of a sample-sheet e.g.: 'Extraction Kit Lot'.
            # There may also be required columns that aren't defined in result.

            # once all columns have been renamed to their preferred names, we
            # must determine the proper set of column names for this sample-
            # sheet. For legacy classes this is simply the list of columns
            # defined in each sample-sheet version. For newer classes, this is
            # defined at run-time and requires examining the metadata that
            # will define the [Data] section.
            required_columns = self._get_expected_data_columns(table=result)
            subset = list(set(required_columns) & set(result.columns))
            result = result[subset]

        return result

    def _add_data_to_sheet(self, table, sequencer, lanes, assay, strict=True):
        if self._remapper is None:
            raise ValueError("sample-sheet does not contain a valid Assay"
                             " type.")

        # Well_description column is now defined here as the concatenation
        # of the following columns. If the column existed previously it will
        # be overwritten, otherwise it will be created here.
        well_description = table[PM_PROJECT_PLATE_KEY].astype(str) + "." + \
            table['Sample'].astype(str) + "." + table['Well'].astype(str)

        table = self._remap_table(table, strict)

        table['Well_description'] = well_description

        for column in self._get_expected_data_columns():
            if column not in table.columns:
                warnings.warn('The column %s in the sample sheet is empty' %
                              column)
                table[column] = ''

        if assay != _AMPLICON:
            table['index2'] = sequencer_i5_index(sequencer, table['index2'])

            self.Bioinformatics['BarcodesAreRC'] = str(
                sequencer in REVCOMP_SEQUENCERS)

        for lane in lanes:
            for sample in table.to_dict(orient='records'):
                sample['Lane'] = lane
                self.add_sample(sample_sheet.Sample(sample))

        return table

    def _add_metadata_to_sheet(self, metadata, sequencer):
        # set the default to avoid index errors if only one of the two is
        # provided.
        self.Reads = [self._READS['Read1'],
                      self._READS['Read2']]

        for metadata_key in self._ALL_METADATA:
            if metadata_key in self._READS:
                if metadata_key == 'Read1':
                    self.Reads[0] = metadata.get(metadata_key, self.Reads[0])
                else:
                    self.Reads[1] = metadata.get(metadata_key, self.Reads[1])

            elif metadata_key in self._SETTINGS:
                self.Settings[metadata_key] = \
                    metadata.get(metadata_key, self._SETTINGS[metadata_key])

            elif metadata_key in self._HEADER:
                if metadata_key == 'Date':
                    # we set the default value here and not in the global
                    # _HEADER dictionary to make sure the date is when the
                    # metadata is written not when the module is imported.
                    self.Header[metadata_key] = metadata.get(
                        metadata_key, datetime.today().strftime('%Y-%m-%d'))
                else:
                    self.Header[metadata_key] = \
                        metadata.get(metadata_key, self._HEADER[metadata_key])

            elif metadata_key in self._KL_ADDTL_DF_SECTIONS:
                # order the df columns to match the expected order
                temp_df = pd.DataFrame(metadata[metadata_key])
                col_order = \
                    list(self._KL_ADDTL_DF_SECTIONS[metadata_key].keys())
                if temp_df.empty:
                    temp_df = pd.DataFrame(columns=col_order)
                else:
                    if set(col_order) != set(temp_df.columns):
                        raise ValueError(
                            f"Columns in {metadata_key} section do not match "
                            f"expected columns: {col_order}")
                    temp_df = temp_df[col_order]
                setattr(self, metadata_key, temp_df)

        # Per MacKenzie's request for 16S don't include Investigator Name and
        # Experiment Name
        if metadata[_ASSAY_KEY] == _AMPLICON:
            if 'Investigator Name' in self.Header:
                del self.Header['Investigator Name']
            if _EXPERIMENT_NAME_KEY in self.Header:
                del self.Header[_EXPERIMENT_NAME_KEY]

            # these are only relevant for metagenomics because they are used in
            # bclconvert
            del self.Settings['MaskShortReads']
            del self.Settings['OverrideCycles']

        # 'MaskShortReads' and 'OverrideCycles' are not relevant for iSeq runs,
        # and can cause issues downstream.

        # Note: 'iseq' should remain at the tail of this list, since it
        # is a substring of the others.
        # NB: If modifying this list, see issue ##234!
        sequencer_types = ['novaseq', 'hiseq', 'miseq', 'miniseq', 'iseq']
        type_found = None
        for sequencer_type in sequencer_types:
            if sequencer_type in sequencer.lower():
                type_found = sequencer_type
                break

        if type_found is None:
            # if even the 'iSeq' substring could not be found, this is an
            # unlikely and unexpected value for sequencer.
            raise ValueError(f"{sequencer} isn't a known sequencer")
        elif type_found == 'iseq':
            #   Verify the settings exist before deleting them.
            if 'MaskShortReads' in self.Settings:
                del self.Settings['MaskShortReads']
            if 'OverrideCycles' in self.Settings:
                del self.Settings['OverrideCycles']

        return self

    def get_lane_number(self):
        lanes = []

        for sample in self.samples:
            lanes.append(sample.Lane)

        lanes = list(set(lanes))

        if len(lanes) > 1:
            raise ValueError("This sample-sheet contains more than one lane")

        return int(lanes[0])

    def _get_expected_data_columns(self, table=None):
        # this base (general) implementation of this method does nothing w/
        # the table parameter. It is present only for compatibility with child
        # methods.
        return self._data_columns

    def validate_and_scrub_sample_sheet(self, echo_msgs=True):
        """Validate the sample sheet and scrub invalid characters

        The character scrubbing is only applied to the Sample_Project and the
        Sample_ID columns. The main difference between this function and
        quiet_validate_and_scrub_sample_sheet is that this function will
        *always* print errors and warnings to standard output.

        Returns
        -------
        Boolean
            True if sample-sheet is valid or if only warnings were reported,
            False if one or more errors were reported.
        """
        msgs = self.quiet_validate_and_scrub_sample_sheet()

        # display Errors and Warnings directly to stdout:
        # this function is used in both Jupyter notebooks where msgs should
        # be displayed, and by other functions that simply want True or False.
        if echo_msgs:
            [msg.echo() for msg in msgs]

        # in addition to displaying any messages, return False if any Errors
        # were found, or True if there were just Warnings or no messages at
        # all.
        if not any([isinstance(m, ErrorMessage) for m in msgs]):
            return True
        else:
            return False

    def quiet_validate_and_scrub_sample_sheet(self):
        """Quietly validate the sample sheet and scrub invalid characters

        The character scrubbing is only applied to the Sample_Project and the
        Sample_ID columns.

        Returns
        -------
        list
            List of error or warning messages.
        """
        msgs = []

        # we print an error return None and exit when this happens otherwise
        # we won't be able to run other checks
        for column in self._get_expected_data_columns():
            if column not in self.all_sample_keys:
                msgs.append(ErrorMessage(f'The {column} column in the '
                                         f'{_DATA_KEY} section is missing'))

        # All children of sample_sheet.SampleSheet will have the following four
        # sections defined: ['Header', 'Reads', 'Settings', 'Data']. All
        # children of KLSampleSheet will have their columns defined in
        # child.sections. We will test only for the difference between these
        # two sets.
        for section in set(type(self).sections).difference({_HEADER_KEY,
                                                            _READS_KEY,
                                                            _SETTINGS_KEY,
                                                            _DATA_KEY}):
            if getattr(self, section) is None:
                msgs.append(ErrorMessage(f'The {section} section cannot be '
                                         'empty'))

        # For cases where a child of KLSampleSheet() is instantiated w/out a
        # filepath, the four base sections will be defined but will be empty
        # unless they are manually populated.
        for attribute in type(self)._HEADER:
            if attribute not in self.Header:
                msgs.append(ErrorMessage(f"'{attribute}' is not declared in "
                                         f"{_HEADER_KEY} section"))

        # Manually populated entries can be arbitrary. Ensure a minimal degree
        # of type consistency.
        expected_assay_type = type(self)._HEADER[_ASSAY_KEY]
        if _ASSAY_KEY in self.Header:
            if self.Header[_ASSAY_KEY] != expected_assay_type:
                msgs.append(ErrorMessage(f"'{_ASSAY_KEY}' value is not "
                                         f"'{expected_assay_type}'"))

        # For sheets that were created by loading in a sample-sheet file,
        # confirm that the SheetType in the file is what is expected from
        # the child class. This helps w/trial-and-error loads that use
        # validation to load a random sample-sheet into the correct class.
        expected_sheet_type = type(self)._HEADER[_SHEET_TYPE_KEY]
        if self.Header[_SHEET_TYPE_KEY] != expected_sheet_type:
            msgs.append(ErrorMessage(f"'{_SHEET_TYPE_KEY}' value is not "
                                     f"'{expected_sheet_type}'"))

        expected_sheet_version = int(type(self)._HEADER[_SHEET_VERSION_KEY])

        # sanitize sample-sheet SheetVersion before attempting to convert to
        # int() type. Remove any additional enclosing quotes.
        sheet_version = list(self.Header[_SHEET_VERSION_KEY])
        sheet_version = [c for c in sheet_version if c not in ['"', "'"]]
        try:
            sheet_version = int(''.join(sheet_version))
        except ValueError:
            msgs.append(ErrorMessage(f"'{self.Header[_SHEET_VERSION_KEY]}' "
                                     f"does not look like a valid value"))

        if sheet_version != expected_sheet_version:
            msgs.append(ErrorMessage(f"'{_SHEET_VERSION_KEY}' value is not "
                                     f"'{expected_sheet_version}'"))

        # if any errors are found up to this point then we can't continue with
        # the validation process.
        if msgs:
            return msgs

        # we track the updated projects as a dictionary so we can propagate
        # these changes to the Bioinformatics and Contact sections.
        # I think it's not necessary to update the _SAMPLE_CONTEXT_KEY section
        # because it uses qiita study ids, not project names, and qiita
        # study ids are not scrubbed.
        updated_samples, updated_projects = [], {}
        for sample in self.samples:
            new_sample = bcl_scrub_name(sample.Sample_ID)
            new_project = bcl_scrub_name(sample.Sample_Project)

            if new_sample != sample.Sample_ID:
                updated_samples.append(sample.Sample_ID)
                sample.Sample_ID = new_sample
            if new_project != sample.Sample_Project:
                updated_projects[sample.Sample_Project] = new_project
                sample[_SS_SAMPLE_PROJECT_KEY] = new_project

        if updated_samples:
            msgs.append(WarningMessage('The following sample names were '
                                       'scrubbed for bcl2fastq compatibility'
                                       ':\n%s' % ', '.join(updated_samples)))
        if updated_projects:
            msgs.append(WarningMessage(
                f"The following project names were scrubbed for bcl2fastq "
                f"compatibility. If the same invalid characters are "
                f"also found in the {_BIOINFORMATICS_KEY} and {_CONTACT_KEY} "
                f"sections, those will be automatically scrubbed too:\n"
                f"{', '.join(sorted(updated_projects))}"))

            # make the changes to prevent useless errors where the scrubbed
            # names fail to match between sections.
            # new pandas won't let you set value inplace on a slice
            project_remapper = {'Sample_Project': updated_projects}
            self.Contact.replace(project_remapper, inplace=True)
            self.Bioinformatics.replace(project_remapper, inplace=True)

        pairs = collections.Counter([(s.Lane, s.Sample_Project)
                                     for s in self.samples])
        # warn users when there's missing lane values
        empty_projects = [project for lane, project in pairs
                          if str(lane).strip() == '']
        if empty_projects:
            msgs.append(
                ErrorMessage(
                    'The following projects are missing a Lane value: '
                    '%s' % ', '.join(sorted(empty_projects))))

        data_project_names = {s.Sample_Project for s in self.samples}

        # project identifiers are digit groups at the end of the project name
        # preceded by an underscore, as in: CaporasoIllumina_550
        bad_projects = []
        for project_name in data_project_names:
            try:
                # don't actually need the return value here :)
                parse_project_name(project_name)
            except ValueError:
                bad_projects.append(project_name)

        if bad_projects:
            msgs.append(
                ErrorMessage(
                    f"The following project names in the "
                    f"{_SS_SAMPLE_PROJECT_KEY} column are missing a Qiita "
                    f"study identifier: "
                    f"{', '.join(sorted(bad_projects))}"))

        # check that the bioinformatics and data sections have the exact
        # same list of projects in them
        bfx_project_names = set(self.Bioinformatics[_SS_SAMPLE_PROJECT_KEY])
        not_shared = data_project_names ^ bfx_project_names
        if not_shared:
            msgs.append(
                ErrorMessage(
                    f"The following projects need to be in the {_DATA_KEY} "
                    f"and {_BIOINFORMATICS_KEY} sections: "
                    f"{', '.join(sorted(not_shared))}"))

        # contact and sample context don't have to have every project in the
        # bioinformatics section, but they can't have any project that ISN'T
        # in the bioinformatics section!
        # NB:below logic works even if there ISN'T a sample context section
        contact_project_names = set(self.Contact[_SS_SAMPLE_PROJECT_KEY])
        sample_context_project_ids = get_all_projects_in_context(
            getattr(self, _SAMPLE_CONTEXT_KEY, None))

        # for each section, the below lists the expected superset and the
        # expected subset.  Note that contacts compares project names, while
        # sample context compares qiita study ids.
        subset_sections = {
            _CONTACT_KEY: (bfx_project_names, contact_project_names),
            _SAMPLE_CONTEXT_KEY: (set(self.Bioinformatics[_SS_QIITA_ID_KEY]),
                                  set(sample_context_project_ids))}
        for curr_section, curr_sets in subset_sections.items():
            curr_missing_projects = curr_sets[1] - curr_sets[0]
            if len(curr_missing_projects) > 0:
                msgs.append(ErrorMessage((
                    f"The following projects were only found in the "
                    f"{curr_section} section: "
                    f"{', '.join(sorted(curr_missing_projects))}. "
                    f"Projects need to be listed in the {_DATA_KEY} and "
                    f"{_BIOINFORMATICS_KEY} section in order to be included "
                    f"in the {curr_section} section.")))
            # end if there are missing projects for this section
        # next section to check

        # silently convert boolean values to either True or False and generate
        # messages for all unrecognizable values.
        msgs += self._normalize_df_sections_booleans()

        # return all collected Messages, even if it's an empty list.
        return msgs

    def sample_is_a_blank(self, sample_name):
        """Check if a sample is a blank.

        Uses sample context section if available; if not, uses sample name
        format.

        Parameters
        ----------
        sample_name: str
            The name of the sample to check.

        Returns
        -------
        bool
            True if the sample is a blank, False otherwise.
        """

        # 1st ensure the input sample name occurs in the list of samples
        samples_details = self._get_samples_details()
        if sample_name not in samples_details:
            raise ValueError(f"Sample '{sample_name}' not found in the "
                             f"{_DATA_KEY} section")
        # endif sample not found

        sample_context = getattr(self, _SAMPLE_CONTEXT_KEY, None)
        # NB: this can be run with a null sample context, so no worries :)
        return is_blank(sample_name, sample_context)

    def get_controls_details(self):
        """Get the control details from the sample sheet."""
        controls_details = {}
        sample_context = getattr(self, _SAMPLE_CONTEXT_KEY, None)
        if sample_context is not None:
            controls_details = \
                get_controls_details_from_context(sample_context)
        else:
            # get the project-name-to-qiita-id mapping from bioinformatics.
            # NB: the logic in this section assumes that all sample project to
            # qiita id mappings are unique in the project-to-qiita_id direction
            # (e.g., MIDAS_10317:10317 and TMI_10317:10317 is acceptable, but
            # MIDAS_10317:10317 and MIDAS_10317:10318 is not).
            proj_name_to_qiita_id_df = self.Bioinformatics[
                [_SS_SAMPLE_PROJECT_KEY, _SS_QIITA_ID_KEY]]

            # get the control info out of the Data section
            samples_details = self._get_samples_details()
            for curr_sample_name, curr_details in samples_details.items():
                if self.sample_is_a_blank(curr_sample_name):
                    curr_name_mask = \
                        proj_name_to_qiita_id_df[_SS_SAMPLE_PROJECT_KEY] == \
                        curr_details[SAMPLE_PROJECT_KEY]
                    curr_qiita_id = proj_name_to_qiita_id_df.loc[
                        curr_name_mask, _SS_QIITA_ID_KEY].tolist()[0]
                    new_control_details = make_manual_control_details(
                        curr_sample_name, curr_qiita_id)
                    controls_details[curr_sample_name] = new_control_details
                # endif this is a blank
            # next sample
        # endif is/isn't a sample context section

        return controls_details

    def get_denormalized_controls_list(self):
        # get info on each control, showing its primary and secondary studies
        normalized_details = self.get_controls_details()

        # create info on each individual control+study combination
        # (so a control that has a primary study and two secondary studies
        # will have three entries in the denormalized_details dictionary).
        # There are no longer primary and secondary study id keys, just a
        # single qiita id key for each control. Result is a list not a dict
        # since a control name may now have >1 record.
        denormalized_details = denormalize_controls_details(normalized_details)

        # get info on all projects in the sample sheet
        project_details = self.get_projects_details()

        def _get_matching_full_project_name(qiita_id):
            matching = \
                [x[PROJECT_FULL_NAME_KEY] for x in project_details.values()
                 if x[QIITA_ID_KEY] == qiita_id]
            if len(matching) != 1:
                raise ValueError(f"Expected 1 matching project for qiita id "
                                 f"{qiita_id}, found {len(matching)}")
            return matching[0]

        # add full project name into each denormalized entry, using qiita id
        for curr_control_details in denormalized_details:
            curr_qiita_id = curr_control_details[QIITA_ID_KEY]
            curr_full_project_name = \
                _get_matching_full_project_name(curr_qiita_id)
            curr_control_details[PROJECT_FULL_NAME_KEY] = \
                curr_full_project_name
        # next denormalized control

        return denormalized_details

    def get_projects_details(self):
        # parse bioinformatics section and data section to generate a
        # durable list of project identifiers and associated sample identifiers
        results = {}

        bioinformatics = self.Bioinformatics
        for curr_project_record in bioinformatics.to_dict(orient='records'):
            curr_full_project_name = \
                curr_project_record[_SS_SAMPLE_PROJECT_KEY]
            curr_proj_dict = parse_project_name(curr_full_project_name)
            curr_proj_dict[SAMPLES_DETAILS_KEY] = {}

            if CONTAINS_REPLICATES_KEY in curr_project_record:
                curr_proj_dict[CONTAINS_REPLICATES_KEY] = \
                    curr_project_record[CONTAINS_REPLICATES_KEY]

            results[curr_full_project_name] = curr_proj_dict
        # next project

        samples_details = self._get_samples_details()
        for curr_sample_name, curr_sample_details in samples_details.items():
            curr_sample_project = curr_sample_details[SAMPLE_PROJECT_KEY]
            results[curr_sample_project][SAMPLES_DETAILS_KEY][
                curr_sample_name] = curr_sample_details
        # next sample

        return results

    def _get_samples_details(self):
        samples_dict = {}

        # Since the project-name is stored in an internal variable
        # in a third-party library, convert samples to JSON using the exposed
        # method first.
        samples = json_loads(self.to_json())[_DATA_KEY]
        for curr_sample in samples:
            curr_sample_name = curr_sample[_SS_SAMPLE_NAME_KEY]
            curr_sample_dict = {
                SAMPLE_NAME_KEY: curr_sample_name,
                SAMPLE_PROJECT_KEY: curr_sample[_SS_SAMPLE_PROJECT_KEY],
                SS_SAMPLE_ID_KEY: curr_sample[SS_SAMPLE_ID_KEY]}

            if ORIG_NAME_KEY in curr_sample:
                curr_sample_dict[ORIG_NAME_KEY] = curr_sample[ORIG_NAME_KEY]

            samples_dict[curr_sample_name] = curr_sample_dict
        # next sample

        return samples_dict

    def _normalize_df_sections_booleans(self):
        msgs = []
        for curr_section_name in self._KL_ADDTL_DF_SECTIONS:
            if hasattr(self, curr_section_name) and \
                    getattr(self, curr_section_name) is not None:
                msgs += self._normalize_section_booleans(curr_section_name)
        return msgs

    def _normalize_section_booleans(self, section_name):
        msgs = []

        def func(x):
            if type(x) is bool:
                # column type is already correct.
                return x
            elif type(x) is str:
                # strings should be converted to bool if possible.
                if x.strip().lower() == 'true':
                    return True
                elif x.strip().lower() == 'false':
                    return False

            # if value isn't recognizably True or False, leave it
            # unchanged and leave a message for the user.
            msgs.append(f"'{x}' is not 'True' or 'False'")
            return x

        section = getattr(self, section_name)
        col_info = self._KL_ADDTL_DF_SECTIONS[section_name]
        for col_name, col_type in col_info.items():
            if col_type is bool:
                if col_name in section:
                    section[col_name] = section[col_name].apply(func)
        setattr(self, section_name, section)

        return msgs

    def _validate_sample_sheet_metadata(self, metadata):
        msgs = []

        # Note: this method is used by all sample sheets, and not all
        # sample sheets include SampleContext, so it is not listed here.
        for req in [_ASSAY_KEY, _BIOINFORMATICS_KEY, _CONTACT_KEY]:
            if req not in metadata:
                msgs.append(ErrorMessage('%s is a required attribute' % req))

        # if both sections are found, then check that all the columns in all
        # extra dataframe sections are present; note that checks for the
        # contents (as opposed to mere presence) are done in the sample sheet
        # validation routine
        if _BIOINFORMATICS_KEY in metadata and _CONTACT_KEY in metadata:
            for section, cols_info in self._KL_ADDTL_DF_SECTIONS.items():
                columns = frozenset(cols_info.keys())

                for i, project in enumerate(metadata[section]):
                    if set(project.keys()) != columns:
                        message = (('In the %s section Project #%d does not '
                                    'have exactly these keys %s') %
                                   (section, i + 1, ', '.join(sorted(columns)))
                                   )
                        msgs.append(ErrorMessage(message))
                    if section == _BIOINFORMATICS_KEY:
                        if (project['library_construction_protocol'] is None or
                                project[
                                    'library_construction_protocol'] == ''):
                            message = (('In the %s section Project #%d does '
                                        'not have library_construction_'
                                        'protocol specified') %
                                       (section, i + 1))
                            msgs.append(ErrorMessage(message))
                        if (project[EXPT_DESIGN_DESC_KEY] is None or
                                project[
                                    EXPT_DESIGN_DESC_KEY] == ''):
                            message = (('In the %s section Project #%d does '
                                        'not have experiment_design_'
                                        'description specified') %
                                       (section, i + 1))
                            msgs.append(ErrorMessage(message))
        if metadata.get(_ASSAY_KEY) is not None and metadata[_ASSAY_KEY] \
                not in self._ASSAYS:
            msgs.append(ErrorMessage(f"{metadata[_ASSAY_KEY]} is not a "
                                     f"supported {_ASSAY_KEY}"))

        keys = set(metadata.keys())
        if not keys.issubset(self._ALL_METADATA):
            extra = sorted(keys - set(self._ALL_METADATA))
            msgs.append(
                ErrorMessage('These metadata keys are not supported: %s'
                             % ', '.join(extra)))

        return msgs


class KLSampleSheetWithSampleContext(KLSampleSheet):
    _KL_ADDTL_DF_SECTIONS = MappingProxyType({
        _BIOINFORMATICS_KEY: _BIOINFORMATICS_COLS_W_REP_SUPPORT,
        _CONTACT_KEY: _CONTACT_COLS,
        _SAMPLE_CONTEXT_KEY: SAMPLE_CONTEXT_COLS
    })

    _ALL_METADATA = MappingProxyType(
        KLSampleSheet._ALL_METADATA | {_SAMPLE_CONTEXT_KEY: None})

    sections = (_HEADER_KEY, _READS_KEY, _SETTINGS_KEY, _DATA_KEY,
                _BIOINFORMATICS_KEY, _CONTACT_KEY, _SAMPLE_CONTEXT_KEY)

    _ORDERED_BY_DATA_COLUMNS = True

    def __new__(cls, path=None, *args, **kwargs):
        """
            Override so that base class cannot be instantiated.
        """
        if cls is KLSampleSheetWithSampleContext:
            raise TypeError(
                f"only children of '{cls.__name__}' may be instantiated")

        instance = super(KLSampleSheetWithSampleContext, cls).__new__(
            cls, *args, **kwargs)
        return instance

    def __init__(self, path=None):
        """Knight Lab's SampleSheet subclass that includes SampleContext

        Expands Knight Lab SampleSheet to include a new (required) section
        called 'SampleContext'. This section is used to store information
        about the blanks and other controls used in the sequencing run.

        Parameters
        ----------
        path: str, optional
            File path to the sample sheet to load.
        """

        # NB: it matters that this come *before* the __init__ call;
        # __init__ calls _parse, which will automatically populate the
        # SampleContext section if it is present in the file--but only if
        # it is defined here first.
        self.SampleContext = None
        super().__init__(path=path)
        self._remapper = _BASE_METAG_REMAPPER
        self._data_columns = _BASE_DATA_COLUMNS
        self._CARRIED_PREP_COLUMNS = _BASE_CARRIED_PREP_COLUMNS


class AbsQuantMixin(object):
    _ABS_SYNDNA_INPUT_MASS_KEY = 'mass_syndna_input_ng'
    _ABS_GDNA_CONC_KEY = 'extracted_gdna_concentration_ng_ul'
    _ABS_SYNDNA_POOL_NUM_KEY = 'syndna_pool_number'
    _ABSQUANT_SPECIFIC_COLUMNS = (
        _ABS_SYNDNA_INPUT_MASS_KEY, _ABS_GDNA_CONC_KEY,
        _ELUTION_VOL_KEY, _ABS_SYNDNA_POOL_NUM_KEY)

    _ABSQUANT_REMAPPER = MappingProxyType(
        _BASE_METAG_REMAPPER | {
            _ABS_SYNDNA_POOL_NUM_KEY: _ABS_SYNDNA_POOL_NUM_KEY,
            _ABS_SYNDNA_INPUT_MASS_KEY: _ABS_SYNDNA_INPUT_MASS_KEY,
            _ABS_GDNA_CONC_KEY: _ABS_GDNA_CONC_KEY,
            _ELUTION_VOL_KEY: _ELUTION_VOL_KEY
        })

    def __init__(self, path=None):
        super().__init__(path=path)
        self._remapper = self._extend_remapper(self._ABSQUANT_REMAPPER)
        self._data_columns = \
            self._data_columns + self._ABSQUANT_SPECIFIC_COLUMNS
        self._CARRIED_PREP_COLUMNS = \
            self._CARRIED_PREP_COLUMNS + self._ABSQUANT_SPECIFIC_COLUMNS


# NB: Must be mixed in to something that inherits from KLSampleSheetWithContext
# since it relies on having self._remapper set in the super's __init__, and on
# the `samples` and `data_columns` attributes
class KatharoseqMixin(object):
    _KATH_RACK_ID_KEY = 'Kathseq_RackID'
    _KATHARO_PREFIX = 'katharo'

    # columns present in a pre-prep file (amplicon) that included katharoseq
    # controls. Presumably we will need these same columns in a sample-sheet.
    _optional_katharoseq_columns = (_KATH_RACK_ID_KEY, TUBECODE_KEY,
                                    'katharo_description',
                                    'number_of_cells',
                                    'platemap_generation_date',
                                    'project_abbreviation',
                                    _ELUTION_VOL_KEY, 'well_id_96')

    @staticmethod
    def _is_katharo_name(sample_name):
        return sample_name.lower().startswith(KatharoseqMixin._KATHARO_PREFIX)

    def __init__(self, path=None):
        super().__init__(path=path)
        self._remapper = self._extend_remapper(
            {self._KATH_RACK_ID_KEY: self._KATH_RACK_ID_KEY})

    def contains_katharoseq_samples(self):
        # when creating samples manually, as opposed to loading a sample-sheet
        # from file, whether or not a sample-sheet contains katharoseq
        # controls can change from add_sample() to add_sample() and won't be
        # determined when MetagenomicSampleSheetv101() is created w/out a
        # file. Hence, perform this check on demand() as opposed to once at
        # init().
        for sample in self.samples:
            # assume any sample-name beginning with 'katharo' in any form of
            # case is a katharoseq sample.
            if self._is_katharo_name(sample.Sample_Name):
                return True

        return False

    def _table_contains_katharoseq_samples(self, table):
        # for instances when a MetagenomicSampleSheetv101() object contains
        # no samples, and the samples will be added in a single method call.
        # this helper method will return True only if a katharo-control
        # sample is found. Note criteria for this method should be kept
        # consistent w/the above method (contains_katharoseq_samples).
        is_katharos = table[_SS_SAMPLE_NAME_KEY].apply(self._is_katharo_name)
        return is_katharos.any()

    def _get_expected_data_columns(self, table=None):
        if table is None:
            # if [Data] section contains katharoseq samples, add the expected
            # additional katharoseq columns to the official list of expected
            # columns before validation or other processing begins.
            if self.contains_katharoseq_samples():
                return self._data_columns + self._optional_katharoseq_columns

        else:
            # assume that there are no samples added to this object yet. This
            # means that self.contains_katharoseq_samples() will always return
            # False. Assume table contains a list of samples that may or may
            # not contain katharoseq controls.
            if self._table_contains_katharoseq_samples(table):
                return self._data_columns + self._optional_katharoseq_columns

        return self._data_columns


class AmpliconSampleSheet(KLSampleSheet):
    _HEADER = {
        'IEMFileVersion': '4',
        _SHEET_TYPE_KEY: _DUMMY_SHEET_TYPE,
        _SHEET_VERSION_KEY: '0',
        # Per MacKenzie's request, these are removed if found.
        # 'Investigator Name': 'Knight',
        # 'Experiment Name': 'RKL_experiment',
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        _ASSAY_KEY: _AMPLICON,
        'Description': '',
        'Chemistry': 'Default',
    }

    _CARRIED_PREP_COLUMNS = (EXPT_DESIGN_DESC_KEY, 'i5_index_id',
                             'i7_index_id', 'index', 'index2',
                             'library_construction_protocol',
                             SAMPLE_NAME_KEY,
                             'sample_plate', 'sample_project',
                             'well_description', 'Sample_Well')

    def __init__(self, path=None):
        super().__init__(path)
        self._remapper = {
            'sample sheet Sample_ID': SS_SAMPLE_ID_KEY,
            'Sample': _SS_SAMPLE_NAME_KEY,
            PM_PROJECT_PLATE_KEY: 'Sample_Plate',
            'Well': 'Sample_Well',
            'Name': 'I7_Index_ID',
            'Golay Barcode': 'index',
            PM_PROJECT_NAME_KEY: _SS_SAMPLE_PROJECT_KEY,
        }


class MetagenomicSampleSheetv102(
        KatharoseqMixin, KLSampleSheetWithSampleContext):
    # Adds support for optional KATHAROSEQ columns in [Data] section.

    _HEADER = KLSampleSheet._HEADER.copy()
    _HEADER[_SHEET_TYPE_KEY] = STANDARD_METAG_SHEET_TYPE
    _HEADER[_SHEET_VERSION_KEY] = '102'
    _HEADER[_ASSAY_KEY] = _METAGENOMIC


class MetagenomicSampleSheetv101(KLSampleSheetWithSampleContext):
    # A copy of MetagenomicSampleSheetv100 but inherits from
    # KLSampleSheetWithSampleContext. This is the first version of the
    # metagenomic sample sheet that includes the SampleContext section.
    _HEADER = KLSampleSheet._HEADER.copy()
    _HEADER[_SHEET_TYPE_KEY] = STANDARD_METAG_SHEET_TYPE
    _HEADER[_SHEET_VERSION_KEY] = '101'
    _HEADER[_ASSAY_KEY] = _METAGENOMIC


class MetagenomicSampleSheetv100(KLSampleSheet):
    _HEADER = {
        'IEMFileVersion': '4',
        _SHEET_TYPE_KEY: STANDARD_METAG_SHEET_TYPE,
        _SHEET_VERSION_KEY: '100',
        'Investigator Name': 'Knight',
        _EXPERIMENT_NAME_KEY: _PLACEHOLDER_EXPT_NAME,
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        _ASSAY_KEY: _METAGENOMIC,
        'Description': '',
        'Chemistry': 'Default',
    }

    # Note that there doesn't appear to be a difference between 95, 99, and 100
    # beyond the value observed in 'Well_description' column. The real
    # difference is between standard_metag and abs_quant_metag.
    _data_columns = _BASE_DATA_COLUMNS

    _KL_ADDTL_DF_SECTIONS = MappingProxyType({
        _BIOINFORMATICS_KEY: _BIOINFORMATICS_COLS_W_REP_SUPPORT,
        _CONTACT_KEY: _CONTACT_COLS,
    })

    _CARRIED_PREP_COLUMNS = (EXPT_DESIGN_DESC_KEY, 'i5_index_id',
                             'i7_index_id', 'index', 'index2',
                             'library_construction_protocol',
                             SAMPLE_NAME_KEY,
                             'sample_plate', 'sample_project',
                             'well_description', 'well_id_384')

    def __init__(self, path=None):
        super().__init__(path=path)
        self._remapper = _BASE_METAG_REMAPPER


class MetagenomicSampleSheetv90(KLSampleSheet):
    """
    MetagenomicSampleSheetv90 is meant to be a class to handle legacy
    Metagenomic type sample-sheets, since KLSampleSheet() itself can't be
    instantiated anymore. What makes it unique is that it specifies a version
    number and defines the classic values for self._remapper.
    """
    _HEADER = {
        'IEMFileVersion': '4',
        _SHEET_TYPE_KEY: STANDARD_METAG_SHEET_TYPE,
        _SHEET_VERSION_KEY: '90',
        'Investigator Name': 'Knight',
        _EXPERIMENT_NAME_KEY: _PLACEHOLDER_EXPT_NAME,
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        _ASSAY_KEY: _METAGENOMIC,
        'Description': '',
        'Chemistry': 'Default',
    }

    # data_columns are the same as base KLSampleSheet so they will not be
    # overridden here. _BIOINFORMATICS_COLUMNS as well.
    _CARRIED_PREP_COLUMNS = (EXPT_DESIGN_DESC_KEY, 'i5_index_id',
                             'i7_index_id', 'index', 'index2',
                             'library_construction_protocol',
                             SAMPLE_NAME_KEY,
                             'sample_plate', 'sample_project',
                             'well_description', 'Sample_Well')

    def __init__(self, path=None):
        super().__init__(path=path)
        self._remapper = {
            'sample sheet Sample_ID': SS_SAMPLE_ID_KEY,
            'Sample': _SS_SAMPLE_NAME_KEY,
            PM_PROJECT_PLATE_KEY: 'Sample_Plate',
            'Well': 'Sample_Well',
            'i7 name': 'I7_Index_ID',
            'i7 sequence': 'index',
            'i5 name': 'I5_Index_ID',
            'i5 sequence': 'index2',
            PM_PROJECT_NAME_KEY: _SS_SAMPLE_PROJECT_KEY
        }


class AbsQuantSampleSheetv10(KLSampleSheet):
    _HEADER = {
        'IEMFileVersion': '4',
        _SHEET_TYPE_KEY: ABSQUANT_SHEET_TYPE,
        _SHEET_VERSION_KEY: '10',
        'Investigator Name': 'Knight',
        _EXPERIMENT_NAME_KEY: _PLACEHOLDER_EXPT_NAME,
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        _ASSAY_KEY: _METAGENOMIC,
        'Description': '',
        'Chemistry': 'Default',
    }

    _data_columns = \
        _BASE_DATA_COLUMNS + AbsQuantMixin._ABSQUANT_SPECIFIC_COLUMNS

    _KL_ADDTL_DF_SECTIONS = MappingProxyType({
        _BIOINFORMATICS_KEY: _BIOINFORMATICS_COLS_W_REP_SUPPORT,
        _CONTACT_KEY: _CONTACT_COLS,
    })

    _CARRIED_PREP_COLUMNS = \
        _BASE_CARRIED_PREP_COLUMNS + AbsQuantMixin._ABSQUANT_SPECIFIC_COLUMNS

    def __init__(self, path=None):
        super().__init__(path=path)
        self._remapper = AbsQuantMixin._ABSQUANT_REMAPPER


class AbsQuantSampleSheetv11(AbsQuantMixin, KLSampleSheetWithSampleContext):
    _HEADER = AbsQuantSampleSheetv10._HEADER.copy()
    _HEADER[_SHEET_TYPE_KEY] = ABSQUANT_SHEET_TYPE
    _HEADER[_SHEET_VERSION_KEY] = '11'


class MetatranscriptomicSampleSheetv0(KLSampleSheet):
    _HEADER = {
        'IEMFileVersion': '4',
        _SHEET_TYPE_KEY: STANDARD_METAG_SHEET_TYPE,
        _SHEET_VERSION_KEY: '0',
        'Investigator Name': 'Knight',
        _EXPERIMENT_NAME_KEY: _PLACEHOLDER_EXPT_NAME,
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        _ASSAY_KEY: _METATRANSCRIPTOMIC,
        'Description': '',
        'Chemistry': 'Default',
    }

    _data_columns = _BASE_DATA_COLUMNS

    _KL_ADDTL_DF_SECTIONS = MappingProxyType({
        _BIOINFORMATICS_KEY: _BIOINFORMATICS_COLS_W_REP_SUPPORT,
        _CONTACT_KEY: _CONTACT_COLS,
    })

    _CARRIED_PREP_COLUMNS = _BASE_CARRIED_PREP_COLUMNS

    def __init__(self, path=None):
        super().__init__(path=path)
        self._remapper = _BASE_METAG_REMAPPER


class MetatranscriptomicSampleSheetv10(KLSampleSheet):
    _HEADER = {
        'IEMFileVersion': '4',
        _SHEET_TYPE_KEY: STANDARD_METAT_SHEET_TYPE,
        _SHEET_VERSION_KEY: '10',
        'Investigator Name': 'Knight',
        _EXPERIMENT_NAME_KEY: _PLACEHOLDER_EXPT_NAME,
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        _ASSAY_KEY: _METATRANSCRIPTOMIC,
        'Description': '',
        'Chemistry': 'Default',
    }

    # MaskShortReads and OverrideCycles are present
    # "Well_description" column contains concatenated information
    # (Sample_Plate + Sample_Name + well_id_384) vs. just the sample_name
    # in previous iterations.

    _data_columns = (SS_SAMPLE_ID_KEY, _SS_SAMPLE_NAME_KEY, 'Sample_Plate',
                     'well_id_384', 'I7_Index_ID', 'index', 'I5_Index_ID',
                     'index2', _SS_SAMPLE_PROJECT_KEY,
                     'total_rna_concentration_ng_ul',
                     _ELUTION_VOL_KEY, 'Well_description')

    _CARRIED_PREP_COLUMNS = _BASE_CARRIED_PREP_COLUMNS + (
                            'total_rna_concentration_ng_ul',
                            _ELUTION_VOL_KEY)

    def __init__(self, path=None):
        super().__init__(path=path)
        self._remapper = _BASE_METAG_REMAPPER | {
                'Sample RNA Concentration': 'total_rna_concentration_ng_ul',
                _ELUTION_VOL_KEY: _ELUTION_VOL_KEY
            }


def _parse_header(fp):
    df = pd.read_csv(fp, dtype="str", sep=",", header=None,
                     names=range(100))

    # pandas will have trouble reading a sample-sheet if the csv has a
    # variable number of columns. This occurs in legacy sheets when a user
    # has introduced one too many ',' characters in a line.
    #
    # the solution is to fix the number of initial columns at a high enough
    # value to include all columns, name them with integers, and later
    # truncate all columns that are entirely empty.
    df.dropna(how='all', axis=1, inplace=True)

    # remove all whitespace rows (drop all rows that are entirely empty)
    df.dropna(how='all', axis=0, inplace=True)

    # remove all comments rows, whether they are at the top of the file (no
    # longer supported, technically), or not.
    comment_rows = df.index[df[0].str.startswith("#")].tolist()
    df = df.drop(index=comment_rows)
    # reset the index to make it easier to post-process.
    df.reset_index(inplace=True, drop=True)

    # for simplicity's sake, assume the first row marks the [Header]
    # column and raise an Error if not. By convention it should be, once
    # legacy comments and whitespace rows are removed.
    if df[0][0] != '[Header]':
        raise ValueError("Top section is not [Header]")

    # identify the beginning of the following section and remove everything
    # from the start of the first section on down. Remove the now redundant
    # [Header] from the top row as well.
    next_section_start = df.index[df[0].str.startswith("[")].tolist()[1]
    df = df.iloc[1:next_section_start]

    # lastly, trim off the additional all-empty columns that are now
    # present after the removal of the other sections.
    df.dropna(how='all', axis=1, inplace=True)

    # set the index to the attributes column of the sample-sheet, replacing
    # the numeric index which now isn't needed. The dataframe will now just
    # contain the index column and a single column named 1.
    df.set_index(0, inplace=True)

    # return the value of key '1'. This will return an immediately
    # recognizable dictionary of key/value pairs.
    results = df.to_dict()[1]

    # conversion to dict causes SheetVersion to be wrapped in single ticks.
    # e.g.: "'100'". These should be removed if present.
    if 'SheetVersion' in results:
        results['SheetVersion'] = results['SheetVersion'].replace("'", "")

    return results


def load_sample_sheet(sample_sheet_path):
    types = [AmpliconSampleSheet,
             MetagenomicSampleSheetv102, MetagenomicSampleSheetv101,
             MetagenomicSampleSheetv100, MetagenomicSampleSheetv90,
             AbsQuantSampleSheetv10, MetatranscriptomicSampleSheetv0,
             MetatranscriptomicSampleSheetv10]

    header = _parse_header(sample_sheet_path)

    required_attributes = ['Assay', 'SheetType', 'SheetVersion']
    missing_attributes = []
    for attribute in required_attributes:
        if attribute not in header:
            missing_attributes.append(f"'{attribute}'")

    if len(missing_attributes) != 0:
        raise ValueError("The following fields must be defined in [Header]: "
                         " %s" % ", ".join(missing_attributes))

    sheet = None

    for type in types:
        m = True
        for attribute in required_attributes:
            if type._HEADER[attribute] != header[attribute]:
                m = False
                break

        if m:
            # header matches all the attributes for the type.
            sheet = type(sample_sheet_path)
            break

    # return a SampleSheet() object if the metadata in the file was
    # successfully matched to a sample-sheet type. this allows the user to
    # call validate_and_scrub() or quiet_validate_and_scrub() on the sample-
    # sheet to determine its correctness or receive warnings and errors.
    if sheet is not None:
        return sheet

    raise ValueError(f"'{sample_sheet_path}' does not appear to be a valid "
                     "sample-sheet.")


def _create_sample_sheet(sheet_type, sheet_version, assay_type):
    def _make_version_err_msg(sheet_type, sheet_version):
        return f"'{sheet_version}' is an unrecognized SheetVersion for " \
               f"'{sheet_type}'"

    def _make_assay_err_msg(assay_type):
        return f"'{assay_type}' is an unrecognized Assay type"

    if sheet_type == STANDARD_METAG_SHEET_TYPE:
        if assay_type == _METAGENOMIC:
            if sheet_version == '102':
                sheet = MetagenomicSampleSheetv102()
            elif sheet_version == '101':
                sheet = MetagenomicSampleSheetv101()
            elif sheet_version == '90':
                sheet = MetagenomicSampleSheetv90()
            elif sheet_version in ['95', '99', '100']:
                # 95, 99, and v100 are functionally the same type.
                sheet = MetagenomicSampleSheetv100()
            else:
                raise ValueError(
                    _make_version_err_msg(sheet_type, sheet_version))
        elif assay_type == _METATRANSCRIPTOMIC:
            sheet = MetatranscriptomicSampleSheetv0()
        else:
            raise ValueError(_make_assay_err_msg(assay_type))
    elif sheet_type == STANDARD_METAT_SHEET_TYPE:
        if assay_type == _METATRANSCRIPTOMIC:
            if sheet_version == '0':
                sheet = MetatranscriptomicSampleSheetv0()
            elif sheet_version == '10':
                sheet = MetatranscriptomicSampleSheetv10()
            else:
                raise ValueError(
                    _make_version_err_msg(sheet_type, sheet_version))
        else:
            raise ValueError(_make_assay_err_msg(assay_type))
    elif sheet_type == ABSQUANT_SHEET_TYPE:
        if assay_type != _METAGENOMIC:
            raise ValueError(_make_assay_err_msg(assay_type))

        if sheet_version == '11':
            sheet = AbsQuantSampleSheetv11()
        elif sheet_version == '10':
            sheet = AbsQuantSampleSheetv10()
        else:
            raise ValueError(_make_version_err_msg(sheet_type, sheet_version))
    elif sheet_type == _DUMMY_SHEET_TYPE:
        sheet = AmpliconSampleSheet()
    else:
        raise ValueError("'%s' is an unrecognized SheetType" % sheet_type)

    return sheet


def make_sample_sheet(metadata, table, sequencer, lanes, strict=None):
    """Write a valid sample sheet

    Parameters
    ----------
    metadata: dict
        Metadata describing the sample sheet with the following fields.
        If a value is omitted from this dictionary the values in square
        brackets are used as defaults.

        - Bioinformatics: List of dictionaries describing each project's
          attributes: Sample_Project, QiitaID, BarcodesAreRC, ForwardAdapter,
          ReverseAdapter, HumanFiltering, library_construction_protocol,
          experiment_design_description
        - Contact: List of dictionaries describing the e-mails to send to
          external stakeholders: Sample_Project, Email

        - IEMFileVersion: Illumina's Experiment Manager version [4]
        - Investigator Name: [Knight]
        - Experiment Name: [RKL_experiment]
        - Date: Date when the sheet is prepared [Today's date]
        - Workflow: how the sample sheet should be used [GenerateFASTQ]
        - Application: sample sheet's application [FASTQ Only]
        - Assay: assay type for the sequencing run. No default value will be
          set, this is required.
        - Description: additional information []
        - Chemistry: chemistry's description [Default]
        - read1: Length of forward read [151]
        - read2: Length of forward read [151]
        - ReverseComplement: If the reads in the FASTQ files should be reverse
          complemented by bcl2fastq [0]
    table: pd.DataFrame
        The Plate's data with one column per variable: sample name ('sample
        sheet Sample_ID'), forward and reverse barcodes ('i5 sequence', 'i7
        sequence'), forward and reverse barcode names ('i5 name', 'i7
        name'), description ('Sample'), well identifier ('Well'), project
        plate (PM_PROJECT_PLATE_KEY), project name (PM_PROJECT_NAME_KEY),
        and synthetic DNA pool number ('syndna_pool_number').
    sequencer: string
        A string representing the sequencer used.
    lanes: list of integers
        A list of integers representing the lanes used.
    strict: boolean
        If True, a subset of columns based on Assay type will define the
        columns in the [Data] section of the sample-sheet. Otherwise all
        columns in table will pass through into the sample-sheet. Either way
        some columns will be renamed as needed by Assay type.

    Returns
    -------
    samplesheet.SampleSheet
        SampleSheet object containing both the metadata and table

    Raises
    ------
    SampleSheetError
        If one of the required columns is missing.
    ValueError
        If the newly-created sample-sheet fails validation.
    """
    required_attributes = [_SHEET_TYPE_KEY, _SHEET_VERSION_KEY, _ASSAY_KEY]

    for attribute in required_attributes:
        if attribute not in metadata:
            raise ValueError("'%s' is not defined in metadata" % attribute)

    sheet_type = metadata[_SHEET_TYPE_KEY]
    sheet_version = metadata[_SHEET_VERSION_KEY]
    assay_type = metadata[_ASSAY_KEY]

    sheet = _create_sample_sheet(sheet_type, sheet_version, assay_type)

    messages = sheet._validate_sample_sheet_metadata(metadata)

    if len(messages) == 0:
        # if the user did not *explicitly* set the strict value
        if strict is None:
            # NB: the below is duck-typing.  It isn't checking whether the
            # sheet's data actually contains any katharoseq samples, but rather
            # whether the sheet *itself* can check whether it
            # contains any katharoseq samples.  If it has this ability, it
            # needs to go through the strict=False handling; see issue #236.
            strict = getattr(
                sheet, 'contains_katharoseq_samples', None) is None

        sheet._add_metadata_to_sheet(metadata, sequencer)
        sheet._add_data_to_sheet(table, sequencer, lanes, metadata[_ASSAY_KEY],
                                 strict)

        # now that we have a SampleSheet() object, validate it for any
        # additional errors that may have been present in the data and/or
        # metadata.
        messages = sheet.quiet_validate_and_scrub_sample_sheet()

        if not any([isinstance(m, ErrorMessage) for m in messages]):
            # No error messages equals success.
            # Echo any warning messages.
            for warning_msg in messages:
                warning_msg.echo()
            return sheet

    # Continue legacy behavior of echoing ErrorMessages and WarningMessages.
    msgs = []
    for message in messages:
        msgs.append(str(message))
        message.echo()

    # Introduce an exception raised for API calls that aren't reporting echo()
    # to the user. Specifically, calls from other modules rather than
    # notebooks. These legacy calls may or may not be testing the returned
    # value for None.
    raise ValueError("\n".join(msgs))


def sample_sheet_to_dataframe(sheet):
    """Converts the [Data] section of a sample sheet into a DataFrame

    Parameters
    ----------
    sheet: sample_sheet.KLSampleSheet
        Object from where to extract the data.

    Returns
    -------
    pd.DataFrame
        DataFrame object with the sample information.
    """

    # Get the columns names for the first sample so we have them in a list and
    # we can retrieve data in the same order on every iteration
    columns = sheet.all_sample_keys

    data = []
    for sample in sheet.samples:
        data.append([sample[column] for column in columns])

    out = pd.DataFrame(data=data, columns=[c.lower() for c in columns])
    out = out.merge(sheet.Bioinformatics[[_SS_SAMPLE_PROJECT_KEY,
                                          'library_construction_protocol',
                                          EXPT_DESIGN_DESC_KEY]],
                    left_on='sample_project', right_on=_SS_SAMPLE_PROJECT_KEY)
    out.drop(columns=_SS_SAMPLE_PROJECT_KEY, inplace=True)

    # it is 'sample_well' and not 'Sample_Well' because of c.lower() above.
    if 'sample_well' in out.columns:
        out.sort_values(by='sample_well', inplace=True)
    elif 'well_id_384' in out.columns:
        out.sort_values(by='well_id_384', inplace=True)
    else:
        raise ValueError("'Sample_Well' and 'well_id_384' columns are not "
                         "present")

    return out.set_index('sample_id')


def sheet_needs_demuxing(sheet):
    """Returns True if sample-sheet needs to be demultiplexed.

    Parameters
    ----------
    sheet: sample_sheet.KLSampleSheet
        Object from where to extract the data.

    Returns
    -------
    bool
        True if sample-sheet needs to be demultiplexed.
    """
    if CONTAINS_REPLICATES_KEY in sheet.Bioinformatics.columns:
        return _get_contains_replicates_value(sheet)

    # legacy sample-sheet does not handle replicates or no replicates were
    # found.
    return False


def _get_contains_replicates_value(sheet):
    contains_replicates = sheet.Bioinformatics[
        CONTAINS_REPLICATES_KEY].unique().tolist()

    # by convention, all projects in the sample-sheet are either going
    # to be True or False. If some projects are True while others are
    # False, we should raise an Error.
    if len(contains_replicates) > 1:
        raise ValueError(f"All projects in {_BIOINFORMATICS_KEY} section "
                         f"must either contain replicates or not.")

    # return either True or False, depending on the values found in
    # Bioinformatics section.
    return list(contains_replicates)[0]


def _demux_sample_sheet(sheet):
    """
    internal function that performs the actual demuxing.
    :param sheet: A valid KLSampleSheet confirmed to have replicates
    :return: a list of DataFrames.
    """
    df = sample_sheet_to_dataframe(sheet)

    # modify df to remove 'library_construction_protocol' and
    # 'experiment_design_description' columns that we don't want for the
    # [Data] section of this sample-sheet.

    df = df.drop(columns=['library_construction_protocol',
                          EXPT_DESIGN_DESC_KEY])

    # use PlateReplication object to convert each sample's 384 well location
    # into a 96-well location + quadrant. Since replication is performed at
    # the plate-level, this will identify which replicates belong in which
    # new sample-sheet.
    plate = PlateReplication(None)

    df['quad'] = df.apply(lambda row: plate.get_96_well_location_and_quadrant(
        row.destination_well_384)[0], axis=1)

    res = []

    for quad in sorted(df['quad'].unique()):
        # for each unique quadrant found, create a new dataframe that's a
        # subset containing only members of that quadrant. Delete the temporary
        # 'quad' column afterwards and reset the index to an integer value
        # starting at zero; the current-index will revert to a column named
        # 'sample_id'. Return the list of new dataframes.
        res.append(df[df['quad'] == quad].drop(['quad'], axis=1))

    return res


def demux_sample_sheet(sheet):
    """Given a sample-sheet w/samples that are plate-replicates, generate new
       sample-sheets for each unique plate-replicate.

        Parameters
        ----------
        sheet: sample_sheet.KLSampleSheet
            Object from where to extract the data.

        Returns
        -------
        list of sheets
    """
    if CONTAINS_REPLICATES_KEY not in sheet.Bioinformatics:
        raise ValueError("sample-sheet does not contain replicates")

    contains_repl_value = _get_contains_replicates_value(sheet)

    # contains_repl_value is of type 'np.bool_' rather than 'bool'. Hence,
    # the syntax below reflects what appears to be common practice for such
    # types.
    if not contains_repl_value:
        raise ValueError(f"No projects in {_BIOINFORMATICS_KEY} section "
                         f"contain replicates")

    demuxed_sheets = []

    # create new sample-sheets, one for each set of replicates. Since
    # replication is performed at the plate level (e.g.: plate 2 is a
    # replicate of plate 1, BLANKS and all), we can split replicates
    # according to their destination quadrant number.
    for df in _demux_sample_sheet(sheet):
        new_sheet = _create_sample_sheet(sheet.Header[_SHEET_TYPE_KEY],
                                         sheet.Header[_SHEET_VERSION_KEY],
                                         sheet.Header[_ASSAY_KEY])
        new_sheet.Header = sheet.Header
        new_sheet.Reads = sheet.Reads
        new_sheet.Settings = sheet.Settings

        # Add the SampleContext section to the new sheet. This is per-sample.
        if _SAMPLE_CONTEXT_KEY in sheet.sections:
            new_context_df = _get_demuxed_sample_context(sheet, df)
            new_sheet.SampleContext = new_context_df
            ctx_projects = \
                _get_sample_context_project_names(sheet, new_context_df)
        else:
            ctx_projects = set()

        projects = set(df.sample_project) | ctx_projects

        # Generate a list of projects associated with each set of samples.
        # Construct bioinformatics and contact sections for each set so that
        # projects not referenced in the sample-set are not included in the
        # Bioinformatics and Contact sections.
        # NB: Don't handle SampleContext section here bc it is sample-, not
        # project-specific, so needs to happen after we set the samples below.
        new_sheet.Bioinformatics = sheet.Bioinformatics.loc[
            sheet.Bioinformatics[_SS_SAMPLE_PROJECT_KEY].isin(projects)].drop(
            [CONTAINS_REPLICATES_KEY], axis=1).reset_index(drop=True)
        new_sheet.Contact = sheet.Contact.loc[
            sheet.Contact[_SS_SAMPLE_PROJECT_KEY].isin(projects)].reset_index(
            drop=True)

        # Add the SampleContext section to the new sheet. This is per-sample.
        if _SAMPLE_CONTEXT_KEY in sheet.sections:
            new_context_df = _get_demuxed_sample_context(sheet, df)
            new_sheet.SampleContext = new_context_df

        # for our purposes here, we want to reindex df so that the index
        # becomes Sample_ID and a new numeric index is created before
        # turning it into a dict. In other situations it remains beneficial
        # for _demux_sample_sheet to return a dataframe with sample_id as
        # the index, such as seqpro.
        df[SS_SAMPLE_ID_KEY] = df.index

        # remove the existing sample_name column that includes appended
        # well-ids. Replace further down w/orig_name column.
        df = df.drop(SAMPLE_NAME_KEY, axis=1)

        df.rename(columns={ORIG_NAME_KEY: _SS_SAMPLE_NAME_KEY,
                           'i7_index_id': 'I7_Index_ID',
                           'i5_index_id': 'I5_Index_ID',
                           'sample_project': _SS_SAMPLE_PROJECT_KEY},
                  inplace=True)
        for sample in df.to_dict(orient='records'):
            new_sheet.add_sample(sample_sheet.Sample(sample))

        demuxed_sheets.append(new_sheet)

    return demuxed_sheets


def _get_demuxed_sample_context(sheet, df):
    # The SampleContext section is a per-sample table, so we want to
    # leave out any samples that are in the ur-SampleContext but have
    # sample_name values that don't match one of the samples that will go in
    # sheet (note we have to match on the sample name that includes the well
    # id, not the well-id-stripped value that becomes the sample name in the
    # new sheet). Then we replace the sample_name column in the new
    # SampleContext section with the well-id-stripped sample_name so it matches
    # what goes in the revised Data section.
    relevant_samples_mask = \
        sheet.SampleContext[SAMPLE_NAME_KEY].isin(
            df[SAMPLE_NAME_KEY])
    temp_context_df = sheet.SampleContext.loc[
        relevant_samples_mask].reset_index(drop=True)
    expanded_temp_context_df = pd.merge(
        temp_context_df, df[[SAMPLE_NAME_KEY, ORIG_NAME_KEY]],
        how="left", on=SAMPLE_NAME_KEY)
    expanded_temp_context_df[SAMPLE_NAME_KEY] = \
        expanded_temp_context_df[ORIG_NAME_KEY]
    expanded_temp_context_df.drop(columns=ORIG_NAME_KEY, inplace=True)
    return expanded_temp_context_df


def _get_sample_context_project_names(sheet, external_context=None):
    ctx_projects = set()
    sample_context = external_context
    if external_context is None:
        if hasattr(sheet, _SAMPLE_CONTEXT_KEY):
            sample_context = getattr(sheet, _SAMPLE_CONTEXT_KEY)

    # The sample context section contains qiita study *ids*, not project
    # names. We need to match these to their corresponding project names in the
    # Bioinformatics section to get project name values useful for comparisons
    # with other parts of the sample sheet.
    ctx_project_ids = get_all_projects_in_context(sample_context)
    if ctx_project_ids is not None:
        bioinformatics = getattr(sheet, _BIOINFORMATICS_KEY)
        ctx_projects_mask = \
            bioinformatics[_SS_QIITA_ID_KEY].isin(ctx_project_ids)
        ctx_projects = \
            set(bioinformatics.loc[ctx_projects_mask, _SS_SAMPLE_PROJECT_KEY])
    # end if there are any projects in the sample context section

    return ctx_projects


def make_sections_dict(plate_df, studies_info, expt_name, expt_type,
                       expt_version, bioinfo_section_base):
    sections_dict = {
        _EXPERIMENT_NAME_KEY: expt_name,
        _SHEET_TYPE_KEY: expt_type,
        _SHEET_VERSION_KEY: expt_version,
        _ASSAY_KEY: _METAGENOMIC
    }

    bioinfo_dicts = []
    contacts_dicts = []
    for curr_study_dict in studies_info:
        # bioinformatics dict
        curr_bioinfo_dict = bioinfo_section_base.copy()
        curr_proj_name = curr_study_dict[PM_PROJECT_NAME_KEY]
        curr_qiita_id = get_qiita_id_from_project_name(curr_proj_name)
        curr_bioinfo_adds = {
            _SS_SAMPLE_PROJECT_KEY: curr_proj_name,
            _SS_QIITA_ID_KEY: curr_qiita_id,
            _HUMAN_FILTERING_KEY: curr_study_dict[_HUMAN_FILTERING_KEY],
            EXPT_DESIGN_DESC_KEY: curr_study_dict[EXPT_DESIGN_DESC_KEY],
            # Per Charlie 20240715: it is accurate that 'contains_replicates'
            # is added individually to each study dict but has the same value
            # in each because it is really a 384-well-plate-level value.
            # Someday maybe it can move into the header.
            CONTAINS_REPLICATES_KEY: plate_df[CONTAINS_REPLICATES_KEY].all(),
        }
        curr_bioinfo_dict.update(curr_bioinfo_adds)
        bioinfo_dicts.append(curr_bioinfo_dict)

        curr_contact_dict = {
            _SS_SAMPLE_PROJECT_KEY: curr_proj_name,
            _EMAIL_KEY: curr_study_dict[_EMAIL_KEY]
        }
        contacts_dicts.append(curr_contact_dict)
    # next study in the run

    sections_dict[_BIOINFORMATICS_KEY] = bioinfo_dicts
    sections_dict[_CONTACT_KEY] = contacts_dicts

    # rows that have True in the Blanks column are blanks :)
    sections_dict[_SAMPLE_CONTEXT_KEY] = \
        get_delimited_controls_details_from_compressed_plate(
            plate_df, blanks_mask=plate_df[PM_BLANK_KEY])

    return sections_dict
