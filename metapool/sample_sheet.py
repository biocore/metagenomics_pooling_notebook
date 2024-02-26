import re
import csv
import collections
import warnings
from datetime import datetime
from itertools import chain, repeat, islice
import sample_sheet
import pandas as pd
from metapool.metapool import (bcl_scrub_name, sequencer_i5_index,
                               REVCOMP_SEQUENCERS)
from metapool.plate import ErrorMessage, WarningMessage, PlateReplication


_AMPLICON = 'TruSeq HT'
_METAGENOMIC = 'Metagenomic'
_METATRANSCRIPTOMIC = 'Metatranscriptomic'
_STANDARD_METAG_SHEET_TYPE = 'standard_metag'
_STANDARD_METAT_SHEET_TYPE = 'standard_metat'
_DUMMY_SHEET_TYPE = 'dummy_amp'
_ABSQUANT_SHEET_TYPE = 'abs_quant_metag'
SHEET_TYPES = (_STANDARD_METAG_SHEET_TYPE, _ABSQUANT_SHEET_TYPE,
               _STANDARD_METAT_SHEET_TYPE)


class KLSampleSheet(sample_sheet.SampleSheet):
    _ASSAYS = frozenset({_AMPLICON, _METAGENOMIC, _METATRANSCRIPTOMIC})
    _BIOINFORMATICS_AND_CONTACT = {
        'Bioinformatics': None,
        'Contact': None
    }

    _CONTACT_COLUMNS = frozenset({
        'Sample_Project', 'Email'
    })

    _BIOINFORMATICS_COLUMNS = frozenset({
        'Sample_Project', 'QiitaID', 'BarcodesAreRC', 'ForwardAdapter',
        'ReverseAdapter', 'HumanFiltering', 'library_construction_protocol',
        'experiment_design_description'
    })

    _HEADER = {
        'IEMFileVersion': '4',
        'SheetType': None,
        'SheetVersion': None,
        'Investigator Name': 'Knight',
        'Experiment Name': 'RKL_experiment',
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        'Assay': None,
        'Description': '',
        'Chemistry': 'Default',
    }

    _READS = {
        'Read1': 151,
        'Read2': 151
    }

    _SETTINGS = {
        'ReverseComplement': '0',
        'MaskShortReads': '1',
        'OverrideCycles': 'Y151;I8N2;I8N2;Y151'
    }

    _ALL_METADATA = {**_HEADER, **_SETTINGS, **_READS,
                     **_BIOINFORMATICS_AND_CONTACT}

    sections = ('Header', 'Reads', 'Settings', 'Data', 'Bioinformatics',
                'Contact')

    data_columns = ('Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well',
                    'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                    'Sample_Project', 'Well_description')

    column_alts = {'well_description': 'Well_description',
                   'description': 'Well_description',
                   'Description': 'Well_description',
                   'sample_plate': 'Sample_Plate'}

    CARRIED_PREP_COLUMNS = ['experiment_design_description', 'i5_index_id',
                            'i7_index_id', 'index', 'index2',
                            'library_construction_protocol', 'sample_name',
                            'sample_plate', 'sample_project',
                            'well_description', 'Sample_Well', 'Lane']

    GENERATED_PREP_COLUMNS = ['center_name', 'center_project_name',
                              'instrument_model', 'lane', 'platform',
                              'run_center', 'run_date', 'run_prefix', 'runid',
                              'sequencing_meth']

    def __new__(cls, path=None, *args, **kwargs):
        """
            Override new() so that base class cannot be instantiated.
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
        self.remapper = None

        self.Bioinformatics = None
        self.Contact = None
        self.path = path

        if self.path:
            self._parse(self.path)

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
                message = ('Comments at the beginning of the sample sheet '
                           'are no longer supported. This information will '
                           'be ignored. Please use the Contact section '
                           'instead')
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
                if section_name == 'Reads':
                    self.Reads.append(int(line[0]))
                    continue

                # [Data] - delimited data with the first line a header.
                elif section_name == 'Data':
                    if section_header is not None:
                        self.add_sample(
                            sample_sheet.Sample(dict(zip(section_header,
                                                         line))))
                    elif any(key == '' for key in line):
                        raise ValueError(
                            f'Header for [Data] section is not allowed to '
                            f'have empty fields: {line}'
                        )
                    else:
                        section_header = self._process_section_header(line)
                    continue

                elif section_name in {'Bioinformatics', 'Contact'}:
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
                        line = [value for value in line if value != '']

                        setattr(self, section_name, pd.DataFrame(columns=line))
                    continue

                # [<Other>]
                else:
                    key, value, *_ = line
                    section = getattr(self, section_name)
                    section[key] = value
                    continue

    def _process_section_header(self, columns):
        for i in range(0, len(columns)):
            if columns[i] in KLSampleSheet.column_alts:
                # overwrite existing alternate name w/internal representation.
                columns[i] = KLSampleSheet.column_alts[columns[i]]
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

        writer = csv.writer(handle)
        csv_width = max([
            len(self.all_sample_keys),
            len(self.Bioinformatics.columns)
            if self.Bioinformatics is not None else 0,
            len(self.Contact.columns) if self.Contact is not None else 0,
            2])

        # custom Illumina sections will go between header reads
        section_order = (['Header', 'Reads', 'Settings', 'Data',
                          'Bioinformatics', 'Contact'])

        def pad_iterable(iterable, size, padding=''):
            return list(islice(chain(iterable, repeat(padding)), size))

        def write_blank_lines(writer, n=blank_lines, width=csv_width):
            for i in range(n):
                writer.writerow(pad_iterable([], width))

        for title in section_order:
            writer.writerow(pad_iterable([f'[{title}]'], csv_width))

            # Data is not a section in this class
            if title != 'Data':
                section = getattr(self, title)

            if title == 'Reads':
                for read in self.Reads:
                    writer.writerow(pad_iterable([read], csv_width))
            elif title == 'Data':
                writer.writerow(pad_iterable(self.all_sample_keys, csv_width))

                for sample in self.samples:
                    line = [getattr(sample, k) for k in self.all_sample_keys]
                    writer.writerow(pad_iterable(line, csv_width))

            elif title == 'Bioinformatics' or title == 'Contact':
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
            for section in ['Header', 'Settings', 'Reads']:
                this, that = getattr(self, section), getattr(sheet, section)

                # For the Header section we'll ignore the Date field since that
                # is likely to be different but shouldn't be a condition to
                # prevent merging two sheets.
                if section == 'Header':
                    if this is not None:
                        this = {k: v for k, v in this.items() if k != 'Date'}
                    if that is not None:
                        that = {k: v for k, v in that.items() if k != 'Date'}

                if this != that:
                    raise ValueError(('The %s section is different for sample '
                                     'sheet %d ') % (section, 1 + number))

            for sample in sheet.samples:
                self.add_sample(sample)

            # these two sections are data frames
            for section in ['Bioinformatics', 'Contact']:
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

    def _remap_table(self, table, strict=True):
        if self.remapper is None:
            raise ValueError("sample-sheet does not contain a valid Assay"
                             " type.")

        # Well_description column is now defined here as the concatenation
        # of the following columns. If the column existed previously it will
        # be overwritten, otherwise it will be created here. Alternate versions
        # of the column name have already been resolved at this point.

        # Note that the amplicon notebook currently generates the same values
        # for this column. If the functionality in the notebook changes, the
        # output will continue to be redfined with the current values here.
        well_description = table['Project Plate'].astype(str) + "." + table[
            'Sample'].astype(str) + "." + table['Well'].astype(str)

        if strict:
            # legacy operation. All columns not defined in remapper will be
            # filtered out.
            out = table[self.remapper.keys()].copy()
            out.rename(self.remapper, axis=1, inplace=True)
        else:
            out = table.copy(deep=True)

            # if a column named 'index' is present in table, assume it is a
            # numeric index and not a sequence of bases, which is required in
            # the output. Assume the column that will become 'index' is
            # defined in remapper.
            if 'index' in set(out.columns):
                out.drop(columns=['index'], inplace=True)

            # if an alternate form of a column name defined in
            # _KL_SAMPLE_SHEET_COLUMN_ALTS is found in table, assume it should
            # be renamed to its proper form and be included in the output e.g.:
            # 'sample_plate' -> 'Sample_Plate'.

            # assume keys in _KL_SAMPLE_SHEET_COLUMN_ALTS do not overlap w/
            # remapper (they currently do not). Define the full set of
            # potential columns to rename in table.

            # new syntax in 3.9 allows us to merge two dicts together w/OR.
            remapper = KLSampleSheet.column_alts | self.remapper
            out.rename(remapper, axis=1, inplace=True)

            # out may contain additional columns that aren't allowed in the
            # [Data] section of a sample-sheet e.g.: 'Extraction Kit Lot'.
            # There may also be required columns that aren't defined in out.
            subset = list(
                set(self.data_columns) & set(
                    out.columns))

            out = out[subset]

        # append the new 'Well_description' column, now that alternates have
        # been removed and non-essential columns have been dropped.
        out['Well_description'] = well_description

        for column in self.data_columns:
            if column not in out.columns:
                warnings.warn('The column %s in the sample sheet is empty' %
                              column)
                out[column] = ''

        return out

    def _add_data_to_sheet(self, table, sequencer, lanes, assay, strict=True):
        table = self._remap_table(table, strict)
        if assay != _AMPLICON:
            table['index2'] = sequencer_i5_index(sequencer, table['index2'])

            self.Bioinformatics['BarcodesAreRC'] = str(
                sequencer in REVCOMP_SEQUENCERS)

        for lane in lanes:
            for sample in table.to_dict(orient='records'):
                sample['Lane'] = lane
                self.add_sample(sample_sheet.Sample(sample))

    def _add_metadata_to_sheet(self, metadata, sequencer):
        # set the default to avoid index errors if only one of the two is
        # provided.
        self.Reads = [self._READS['Read1'],
                      self._READS['Read2']]

        for key in self._ALL_METADATA:
            if key in self._READS:
                if key == 'Read1':
                    self.Reads[0] = metadata.get(key, self.Reads[0])
                else:
                    self.Reads[1] = metadata.get(key, self.Reads[1])

            elif key in self._SETTINGS:
                self.Settings[key] = metadata.get(key,
                                                  self._SETTINGS[key])

            elif key in self._HEADER:
                if key == 'Date':
                    # we set the default value here and not in the global
                    # _HEADER dictionary to make sure the date is when the
                    # metadata is written not when the module is imported.
                    self.Header[key] = metadata.get(
                        key, datetime.today().strftime('%Y-%m-%d'))
                else:
                    self.Header[key] = metadata.get(key,
                                                    self._HEADER[key])

            elif key in self._BIOINFORMATICS_AND_CONTACT:
                setattr(self, key, pd.DataFrame(metadata[key]))

        # Per MacKenzie's request for 16S don't include Investigator Name and
        # Experiment Name
        if metadata['Assay'] == _AMPLICON:
            if 'Investigator Name' in self.Header:
                del self.Header['Investigator Name']
            if 'Experiment Name' in self.Header:
                del self.Header['Experiment Name']

            # these are only relevant for metagenomics because they are used in
            # bclconvert
            del self.Settings['MaskShortReads']
            del self.Settings['OverrideCycles']

        # 'MaskShortReads' and 'OverrideCycles' are not relevant for iSeq runs,
        # and can cause issues downstream.

        # Note: 'iseq' should remain at the tail of this list, since it
        # is a substring of the others.
        sequencer_types = ['novaseq', 'hiseq', 'miseq', 'miniseq', 'iseq']
        type_found = None
        for sequencer_type in sequencer_types:
            if sequencer_type in sequencer.lower():
                type_found = sequencer_type
                break

        if type_found is None:
            # if even the 'iSeq' substring could not be found, this is an
            # unlikely and unexpected value for sequencer.
            raise ErrorMessage(f"{sequencer} isn't a known sequencer")
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
        for column in self.data_columns:
            if column not in self.all_sample_keys:
                msgs.append(ErrorMessage(f'The {column} column in the Data '
                                         'section is missing'))

        # All children of sample_sheet.SampleSheet will have the following four
        # sections defined: ['Header', 'Reads', 'Settings', 'Data']. All
        # children of KLSampleSheet will have their columns defined in
        # child.sections. We will test only for the difference between these
        # two sets.
        for section in set(type(self).sections).difference({'Header',
                                                            'Reads',
                                                            'Settings',
                                                            'Data'}):
            if getattr(self, section) is None:
                msgs.append(ErrorMessage(f'The {section} section cannot be '
                                         'empty'))

        # For cases where a child of KLSampleSheet() is instantiated w/out a
        # filepath, the four base sections will be defined but will be empty
        # unless they are manually populated.
        for attribute in type(self)._HEADER:
            if attribute not in self.Header:
                msgs.append(ErrorMessage(f"'{attribute}' is not declared in "
                                         "Header section"))

        # Manually populated entries can be arbitrary. Ensure a minimal degree
        # of type consistency.
        expected_assay_type = type(self)._HEADER['Assay']
        if 'Assay' in self.Header:
            if self.Header['Assay'] != expected_assay_type:
                msgs.append(ErrorMessage("'Assay' value is not "
                                         f"'{expected_assay_type}'"))

        # For sheets that were created by loading in a sample-sheet file,
        # confirm that the SheetType in the file is what is expected from
        # the child class. This helps w/trial-and-error loads that use
        # validation to load a random sample-sheet into the correct class.
        expected_sheet_type = type(self)._HEADER['SheetType']
        if self.Header['SheetType'] != expected_sheet_type:
            msgs.append(ErrorMessage("'SheetType' value is not "
                                     f"'{expected_sheet_type}'"))

        # if any errors are found up to this point then we can't continue with
        # the validation process.
        if msgs:
            return msgs

        # we track the updated projects as a dictionary so we can propagate
        # these changes to the Bioinformatics and Contact sections
        updated_samples, updated_projects = [], {}
        for sample in self.samples:
            new_sample = bcl_scrub_name(sample.Sample_ID)
            new_project = bcl_scrub_name(sample.Sample_Project)

            if new_sample != sample.Sample_ID:
                updated_samples.append(sample.Sample_ID)
                sample.Sample_ID = new_sample
            if new_project != sample.Sample_Project:
                updated_projects[sample.Sample_Project] = new_project
                sample['Sample_Project'] = new_project

        if updated_samples:
            msgs.append(WarningMessage('The following sample names were '
                                       'scrubbed for bcl2fastq compatibility'
                                       ':\n%s' % ', '.join(updated_samples)))
        if updated_projects:
            msgs.append(WarningMessage('The following project names were '
                                       'scrubbed for bcl2fastq compatibility. '
                                       'If the same invalid characters are '
                                       'also found in the Bioinformatics and '
                                       'Contacts sections those will be '
                                       'automatically scrubbed too:\n%s' %
                                       ', '.join(sorted(updated_projects))))

            # make the changes to prevent useless errors where the scurbbed
            # names fail to match between sections.
            self.Contact.Sample_Project.replace(updated_projects,
                                                inplace=True)
            self.Bioinformatics.Sample_Project.replace(updated_projects,
                                                       inplace=True)

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

        projects = {s.Sample_Project for s in self.samples}

        # project identifiers are digit groups at the end of the project name
        # preceded by an underscore CaporasoIllumina_550
        qiita_id_re = re.compile(r'(.+)_(\d+)$')
        bad_projects = []
        for project_name in projects:
            if re.search(qiita_id_re, project_name) is None:
                bad_projects.append(project_name)
        if bad_projects:
            msgs.append(
                ErrorMessage(
                    'The following project names in the Sample_Project '
                    'column are missing a Qiita study '
                    'identifier: %s' % ', '.join(sorted(bad_projects))))

        # check Sample_project values match across sections
        bfx = set(self.Bioinformatics['Sample_Project'])
        contact = set(self.Contact['Sample_Project'])
        not_shared = projects ^ bfx
        if not_shared:
            msgs.append(
                ErrorMessage(
                    'The following projects need to be in the Data and '
                    'Bioinformatics sections %s' %
                    ', '.join(sorted(not_shared))))
        elif not contact.issubset(projects):
            msgs.append(
                ErrorMessage(('The following projects were only found in the '
                              'Contact section: %s projects need to be listed '
                              'in the Data and Bioinformatics section in order'
                              ' to be included in the Contact section.') %
                             ', '.join(sorted(contact - projects))))

        # return all collected Messages, even if it's an empty list.
        return msgs

    def _validate_sample_sheet_metadata(self, metadata):
        msgs = []

        for req in ['Assay', 'Bioinformatics', 'Contact']:
            if req not in metadata:
                msgs.append(ErrorMessage('%s is a required attribute' % req))

        # if both sections are found, then check that all the columns are
        # present, checks for the contents are done in the sample sheet
        # validation routine
        if 'Bioinformatics' in metadata and 'Contact' in metadata:
            for section in ['Bioinformatics', 'Contact']:
                if section == 'Bioinformatics':
                    columns = self._BIOINFORMATICS_COLUMNS
                else:
                    columns = self._CONTACT_COLUMNS

                for i, project in enumerate(metadata[section]):
                    if set(project.keys()) != columns:
                        message = (('In the %s section Project #%d does not '
                                    'have exactly these keys %s') %
                                   (section, i + 1, ', '.join(sorted(columns)))
                                   )
                        msgs.append(ErrorMessage(message))
                    if section == 'Bioinformatics':
                        if (project['library_construction_protocol'] is None or
                                project[
                                    'library_construction_protocol'] == ''):
                            message = (('In the %s section Project #%d does '
                                        'not have library_construction_'
                                        'protocol specified') %
                                       (section, i + 1))
                            msgs.append(ErrorMessage(message))
                        if (project['experiment_design_description'] is None or
                                project[
                                    'experiment_design_description'] == ''):
                            message = (('In the %s section Project #%d does '
                                        'not have experiment_design_'
                                        'description specified') %
                                       (section, i + 1))
                            msgs.append(ErrorMessage(message))
        if metadata.get('Assay') is not None and metadata['Assay'] \
                not in self._ASSAYS:
            msgs.append(ErrorMessage('%s is not a supported Assay' %
                                     metadata['Assay']))

        keys = set(metadata.keys())
        if not keys.issubset(self._ALL_METADATA):
            extra = sorted(keys - set(self._ALL_METADATA))
            msgs.append(
                ErrorMessage('These metadata keys are not supported: %s'
                             % ', '.join(extra)))

        return msgs


class AmpliconSampleSheet(KLSampleSheet):
    _HEADER = {
        'IEMFileVersion': '4',
        'SheetType': _DUMMY_SHEET_TYPE,
        'SheetVersion': '0',
        # Per MacKenzie's request, these are removed if found.
        # 'Investigator Name': 'Knight',
        # 'Experiment Name': 'RKL_experiment',
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        'Assay': _AMPLICON,
        'Description': '',
        'Chemistry': 'Default',
    }

    CARRIED_PREP_COLUMNS = ['experiment_design_description', 'i5_index_id',
                            'i7_index_id', 'index', 'index2',
                            'library_construction_protocol', 'sample_name',
                            'sample_plate', 'sample_project',
                            'well_description', 'Sample_Well']

    def __init__(self, path=None):
        super().__init__(path)
        self.remapper = {
            'sample sheet Sample_ID': 'Sample_ID',
            'Sample': 'Sample_Name',
            'Project Plate': 'Sample_Plate',
            'Well': 'Sample_Well',
            'Name': 'I7_Index_ID',
            'Golay Barcode': 'index',
            'Project Name': 'Sample_Project',
        }


class MetagenomicSampleSheetv100(KLSampleSheet):
    _HEADER = {
        'IEMFileVersion': '4',
        'SheetType': _STANDARD_METAG_SHEET_TYPE,
        'SheetVersion': '100',
        'Investigator Name': 'Knight',
        'Experiment Name': 'RKL_experiment',
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        'Assay': _METAGENOMIC,
        'Description': '',
        'Chemistry': 'Default',
    }

    # Note that there doesn't appear to be a difference between 95, 99, and 100
    # beyond the value observed in 'Well_description' column. The real
    # difference is between standard_metag and abs_quant_metag.

    # Marks change from 'Metagenomics' to 'Metagenomic' - encapsulate this
    # change: TODO.

    # Note: Remove syndna_pool_number as that was part of the purpose of
    # making this change. Also, it's always going to be empty or worse have
    # a value that won't be checked.
    data_columns = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                    'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                    'Sample_Project', 'Well_description']

    # For now, assume only MetagenomicSampleSheetv100 (and v95, v99) contains
    # 'contains_replicates' column. Assume AbsQuantSampleSheetv10 doesn't.

    _BIOINFORMATICS_COLUMNS = {'Sample_Project', 'QiitaID', 'BarcodesAreRC',
                               'ForwardAdapter', 'ReverseAdapter',
                               'HumanFiltering', 'contains_replicates',
                               'library_construction_protocol',
                               'experiment_design_description'}

    CARRIED_PREP_COLUMNS = ['experiment_design_description', 'i5_index_id',
                            'i7_index_id', 'index', 'index2',
                            'library_construction_protocol', 'sample_name',
                            'sample_plate', 'sample_project',
                            'well_description', 'well_id_384']

    def __init__(self, path=None):
        super().__init__(path=path)
        self.remapper = {
            'sample sheet Sample_ID': 'Sample_ID',
            'Sample': 'Sample_Name',
            'Project Plate': 'Sample_Plate',
            'Well': 'well_id_384',
            'i7 name': 'I7_Index_ID',
            'i7 sequence': 'index',
            'i5 name': 'I5_Index_ID',
            'i5 sequence': 'index2',
            'Project Name': 'Sample_Project',
        }


class MetagenomicSampleSheetv90(KLSampleSheet):
    '''
    MetagenomicSampleSheetv90 is meant to be a class to handle legacy
    Metagenomic type sample-sheets, since KLSampleSheet() itself can't be
    instantiated anymore. What makes it unique is that it specifies a version
    number and defines the classic values for self.remapper.
    '''
    _HEADER = {
        'IEMFileVersion': '4',
        'SheetType': _STANDARD_METAG_SHEET_TYPE,
        'SheetVersion': '90',
        'Investigator Name': 'Knight',
        'Experiment Name': 'RKL_experiment',
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        'Assay': _METAGENOMIC,
        'Description': '',
        'Chemistry': 'Default',
    }

    # data_columns are the same as base KLSampleSheet so they will not be
    # overridden here. _BIOINFORMATICS_COLUMNS as well.
    CARRIED_PREP_COLUMNS = ['experiment_design_description', 'i5_index_id',
                            'i7_index_id', 'index', 'index2',
                            'library_construction_protocol', 'sample_name',
                            'sample_plate', 'sample_project',
                            'well_description', 'Sample_Well']

    def __init__(self, path=None):
        super().__init__(path=path)
        self.remapper = {
            'sample sheet Sample_ID': 'Sample_ID',
            'Sample': 'Sample_Name',
            'Project Plate': 'Sample_Plate',
            'Well': 'Sample_Well',
            'i7 name': 'I7_Index_ID',
            'i7 sequence': 'index',
            'i5 name': 'I5_Index_ID',
            'i5 sequence': 'index2',
            'Project Name': 'Sample_Project'
        }


class AbsQuantSampleSheetv10(KLSampleSheet):
    _HEADER = {
        'IEMFileVersion': '4',
        'SheetType': _ABSQUANT_SHEET_TYPE,
        'SheetVersion': '10',
        'Investigator Name': 'Knight',
        'Experiment Name': 'RKL_experiment',
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        'Assay': _METAGENOMIC,
        'Description': '',
        'Chemistry': 'Default',
    }

    data_columns = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                    'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                    'Sample_Project', 'mass_syndna_input_ng',
                    'extracted_gdna_concentration_ng_ul',
                    'vol_extracted_elution_ul', 'syndna_pool_number',
                    'Well_description']

    _BIOINFORMATICS_COLUMNS = frozenset({
        'Sample_Project', 'QiitaID', 'BarcodesAreRC', 'ForwardAdapter',
        'ReverseAdapter', 'HumanFiltering', 'library_construction_protocol',
        'experiment_design_description', 'contains_replicates'
    })

    CARRIED_PREP_COLUMNS = ['experiment_design_description',
                            'extracted_gdna_concentration_ng_ul',
                            'i5_index_id', 'i7_index_id', 'index', 'index2',
                            'library_construction_protocol',
                            'mass_syndna_input_ng', 'sample_name',
                            'sample_plate', 'sample_project',
                            'syndna_pool_number', 'vol_extracted_elution_ul',
                            'well_description', 'well_id_384']

    def __init__(self, path=None):
        super().__init__(path=path)
        self.remapper = {
            'sample sheet Sample_ID': 'Sample_ID',
            'Sample': 'Sample_Name',
            'Project Plate': 'Sample_Plate',
            'Well': 'well_id_384',
            'i7 name': 'I7_Index_ID',
            'i7 sequence': 'index',
            'i5 name': 'I5_Index_ID',
            'i5 sequence': 'index2',
            'Project Name': 'Sample_Project',
            'syndna_pool_number': 'syndna_pool_number',
            'mass_syndna_input_ng': 'mass_syndna_input_ng',
            'extracted_gdna_concentration_ng_ul':
                'extracted_gdna_concentration_ng_ul',
            'vol_extracted_elution_ul':
                'vol_extracted_elution_ul'
        }


class MetatranscriptomicSampleSheetv0(KLSampleSheet):
    _HEADER = {
        'IEMFileVersion': '4',
        'SheetType': _STANDARD_METAG_SHEET_TYPE,
        'SheetVersion': '0',
        'Investigator Name': 'Knight',
        'Experiment Name': 'RKL_experiment',
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        'Assay': _METATRANSCRIPTOMIC,
        'Description': '',
        'Chemistry': 'Default',
    }

    data_columns = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                    'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                    'Sample_Project', 'Well_description']

    _BIOINFORMATICS_COLUMNS = frozenset({'Sample_Project', 'QiitaID',
                                         'BarcodesAreRC', 'ForwardAdapter',
                                         'ReverseAdapter', 'HumanFiltering',
                                         'contains_replicates',
                                         'library_construction_protocol',
                                         'experiment_design_description'})

    CARRIED_PREP_COLUMNS = ['experiment_design_description', 'i5_index_id',
                            'i7_index_id', 'index', 'index2',
                            'library_construction_protocol', 'sample_name',
                            'sample_plate', 'sample_project',
                            'well_description', 'well_id_384']

    def __init__(self, path=None):
        super().__init__(path=path)
        self.remapper = {
            'sample sheet Sample_ID': 'Sample_ID',
            'Sample': 'Sample_Name',
            'Project Plate': 'Sample_Plate',
            'Well': 'well_id_384',
            'i7 name': 'I7_Index_ID',
            'i7 sequence': 'index',
            'i5 name': 'I5_Index_ID',
            'i5 sequence': 'index2',
            'Project Name': 'Sample_Project',
        }


class MetatranscriptomicSampleSheetv10(KLSampleSheet):
    _HEADER = {
        'IEMFileVersion': '4',
        'SheetType': _STANDARD_METAT_SHEET_TYPE,
        'SheetVersion': '10',
        'Investigator Name': 'Knight',
        'Experiment Name': 'RKL_experiment',
        'Date': None,
        'Workflow': 'GenerateFASTQ',
        'Application': 'FASTQ Only',
        'Assay': _METATRANSCRIPTOMIC,
        'Description': '',
        'Chemistry': 'Default',
    }

    # MaskShortReads and OverrideCycles are present
    # "Well_description" column contains concatenated information
    # (Sample_Plate + Sample_Name + well_id_384) vs. just the sample_name
    # in previous iterations.

    data_columns = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                    'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                    'Sample_Project', 'total_rna_concentration_ng_ul',
                    'vol_extracted_elution_ul', 'Well_description']

    _BIOINFORMATICS_COLUMNS = frozenset({'Sample_Project', 'QiitaID',
                                         'BarcodesAreRC', 'ForwardAdapter',
                                         'ReverseAdapter', 'HumanFiltering',
                                         'contains_replicates',
                                         'library_construction_protocol',
                                         'experiment_design_description',
                                         'contains_replicates'})

    CARRIED_PREP_COLUMNS = ['experiment_design_description', 'i5_index_id',
                            'i7_index_id', 'index', 'index2',
                            'library_construction_protocol', 'sample_name',
                            'sample_plate', 'sample_project',
                            'well_description', 'well_id_384']

    def __init__(self, path=None):
        super().__init__(path=path)
        self.remapper = {
            'sample sheet Sample_ID': 'Sample_ID',
            'Sample': 'Sample_Name',
            'Project Plate': 'Sample_Plate',
            'Well': 'well_id_384',
            'i7 name': 'I7_Index_ID',
            'i7 sequence': 'index',
            'i5 name': 'I5_Index_ID',
            'i5 sequence': 'index2',
            'Project Name': 'Sample_Project',
        }


def load_sample_sheet(sample_sheet_path):
    # Load the sample-sheet using various KLSampleSheet children and return
    # the first instance that produces a valid sample-sheet. We assume that
    # because of specific SheetType and SheetVersion values, no one sample-
    # sheet can match more than one KLSampleSheet child.

    sheet = AbsQuantSampleSheetv10(sample_sheet_path)
    if sheet.validate_and_scrub_sample_sheet(echo_msgs=False):
        return sheet

    sheet = AmpliconSampleSheet(sample_sheet_path)
    if sheet.validate_and_scrub_sample_sheet(echo_msgs=False):
        return sheet

    sheet = MetagenomicSampleSheetv100(sample_sheet_path)
    if sheet.validate_and_scrub_sample_sheet(echo_msgs=False):
        return sheet

    sheet = MetagenomicSampleSheetv90(sample_sheet_path)
    if sheet.validate_and_scrub_sample_sheet(echo_msgs=False):
        return sheet

    sheet = MetatranscriptomicSampleSheetv10(sample_sheet_path)

    if sheet.validate_and_scrub_sample_sheet(echo_msgs=False):
        return sheet

    sheet = MetatranscriptomicSampleSheetv0(sample_sheet_path)

    if sheet.validate_and_scrub_sample_sheet(echo_msgs=False):
        return sheet

    sheet = MetatranscriptomicSampleSheetv0(sample_sheet_path)
    if sheet.validate_and_scrub_sample_sheet(echo_msgs=False):
        return sheet

    raise ValueError(f"'{sample_sheet_path}' does not appear to be a valid "
                     "sample-sheet.")


def _create_sample_sheet(sheet_type, sheet_version, assay_type):
    if sheet_type == _STANDARD_METAG_SHEET_TYPE:
        if assay_type == _METAGENOMIC:
            if sheet_version == '90':
                sheet = MetagenomicSampleSheetv90()
            elif sheet_version in ['95', '99', '100']:
                # 95, 99, and v100 are functionally the same type.
                sheet = MetagenomicSampleSheetv100()
            else:
                raise ValueError(f"'{sheet_version}' is an unrecognized Sheet"
                                 f"Version for '{sheet_type}'")
        elif assay_type == _METATRANSCRIPTOMIC:
            sheet = MetatranscriptomicSampleSheetv0()
        else:
            raise ValueError("'%s' is an unrecognized Assay type" % assay_type)
    elif sheet_type == _STANDARD_METAT_SHEET_TYPE:
        if assay_type == _METATRANSCRIPTOMIC:
            if sheet_version == '0':
                sheet = MetatranscriptomicSampleSheetv0()
            elif sheet_version == '10':
                sheet = MetatranscriptomicSampleSheetv10()
            else:
                raise ValueError(f"'{sheet_version}' is an unrecognized Sheet"
                                 f"Version for '{sheet_type}'")
        else:
            raise ValueError("'%s' is an unrecognized Assay type" % assay_type)
    elif sheet_type == _ABSQUANT_SHEET_TYPE:
        sheet = AbsQuantSampleSheetv10()
    elif sheet_type == _DUMMY_SHEET_TYPE:
        sheet = AmpliconSampleSheet()
    else:
        raise ValueError("'%s' is an unrecognized SheetType" % sheet_type)

    return sheet


def make_sample_sheet(metadata, table, sequencer, lanes, strict=True):
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
        plate ('Project Plate'), project name ('Project Name'), and synthetic
        DNA pool number ('syndna_pool_number').
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
    required_attributes = ['SheetType', 'SheetVersion', 'Assay']

    for attribute in required_attributes:
        if attribute not in metadata:
            raise ValueError("'%s' is not defined in metadata" % attribute)

    sheet_type = metadata['SheetType']
    sheet_version = metadata['SheetVersion']
    assay_type = metadata['Assay']

    sheet = _create_sample_sheet(sheet_type, sheet_version, assay_type)

    messages = sheet._validate_sample_sheet_metadata(metadata)

    if len(messages) == 0:
        sheet._add_metadata_to_sheet(metadata, sequencer)
        sheet._add_data_to_sheet(table, sequencer, lanes, metadata['Assay'],
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
    out = out.merge(sheet.Bioinformatics[['Sample_Project',
                                          'library_construction_protocol',
                                          'experiment_design_description']],
                    left_on='sample_project', right_on='Sample_Project')
    out.drop(columns='Sample_Project', inplace=True)
    out.sort_values(by='well_id_384', inplace=True)
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
    if 'contains_replicates' in sheet.Bioinformatics.columns:
        contains_replicates = sheet.Bioinformatics.contains_replicates.apply(
            lambda x: x.lower() == 'true').unique()
        if len(contains_replicates) > 1:
            raise ValueError("all projects in Bioinformatics section must "
                             "either contain replicates or not.")

        # return either True or False, depending on the values found in
        # Bioinformatics section.
        return list(contains_replicates)[0]

    # legacy sample-sheet does not handle replicates or no replicates were
    # found.
    return False


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
                          'experiment_design_description'])

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
    if 'contains_replicates' not in sheet.Bioinformatics:
        raise ValueError("sample-sheet does not contain replicates")

    # by convention, all projects in the sample-sheet are either going
    # to be True or False. If some projects are True while others are
    # False, we should raise an Error.
    contains_replicates = sheet.Bioinformatics.contains_replicates.apply(
        lambda x: x.lower() == 'true').unique()
    if len(contains_replicates) > 1:
        raise ValueError("all projects in Bioinformatics section must "
                         "either contain replicates or not.")

    # contains_replicates[0] is of type 'np.bool_' rather than 'bool'. Hence,
    # the syntax below reflects what appears to be common practice for such
    # types.
    if not contains_replicates[0]:
        raise ValueError("all projects in Bioinformatics section do not "
                         "contain replicates")

    demuxed_sheets = []

    # create new sample-sheets, one for each set of replicates. Since
    # replication is performed at the plate level (e.g.: plate 2 is a
    # replicate of plate 1, BLANKS and all), we can split replicates
    # according to their destination quadrant number.
    for df in _demux_sample_sheet(sheet):
        new_sheet = _create_sample_sheet(sheet.Header['SheetType'],
                                         sheet.Header['SheetVersion'],
                                         sheet.Header['Assay'])
        new_sheet.Header = sheet.Header
        new_sheet.Reads = sheet.Reads
        new_sheet.Settings = sheet.Settings
        projects = set(df.sample_project)

        # Generate a list of projects associated with each set of samples.
        # Construct bioinformatics and contact sections for each set so that
        # projects not referenced in the sample-set are not included in the
        # Bioinformatics and Contact sections.

        new_sheet.Bioinformatics = sheet.Bioinformatics.loc[
            sheet.Bioinformatics['Sample_Project'].isin(projects)].drop(
            ['contains_replicates'], axis=1).reset_index(drop=True)
        new_sheet.Contact = sheet.Contact.loc[
            sheet.Contact['Sample_Project'].isin(projects)].reset_index(
            drop=True)

        # for our purposes here, we want to reindex df so that the index
        # becomes Sample_ID and a new numeric index is created before
        # turning it into a dict. In other situations it remains beneficial
        # for _demux_sample_sheet to return a dataframe with sample_id as
        # the index, such as seqpro.
        df['Sample_ID'] = df.index
        df.rename(columns={'sample_name': 'Sample_Name',
                           'i7_index_id': 'I7_Index_ID',
                           'i5_index_id': 'I5_Index_ID',
                           'sample_project': 'Sample_Project'}, inplace=True)
        for sample in df.to_dict(orient='records'):
            new_sheet.add_sample(sample_sheet.Sample(sample))

        demuxed_sheets.append(new_sheet)

    return demuxed_sheets
