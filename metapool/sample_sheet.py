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


_KL_SAMPLE_SHEET_SECTIONS = [
    'Header', 'Reads', 'Settings', 'Data', 'Bioinformatics', 'Contact'
]

_KL_SAMPLE_SHEET_DATA_COLUMNS = [
    'Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384', 'I7_Index_ID',
    'index', 'I5_Index_ID', 'index2', 'Sample_Project', 'syndna_pool_number',
    'Well_description'
]

_KL_SAMPLE_SHEET_COLUMN_ALTS = {'well_description': 'Well_description',
                                'description': 'Well_description',
                                'Description': 'Well_description',
                                'sample_plate': 'Sample_Plate'}

_KL_AMPLICON_REMAPPER = {
    'sample sheet Sample_ID': 'Sample_ID',
    'Sample': 'Sample_Name',
    'Project Plate': 'Sample_Plate',
    'Well': 'well_id_384',
    'Name': 'I7_Index_ID',
    'Golay Barcode': 'index',
    'Project Name': 'Sample_Project',
}

_KL_METAGENOMIC_REMAPPER = {
    'sample sheet Sample_ID': 'Sample_ID',
    'Sample': 'Sample_Name',
    'Project Plate': 'Sample_Plate',
    'Well': 'well_id_384',
    'i7 name': 'I7_Index_ID',
    'i7 sequence': 'index',
    'i5 name': 'I5_Index_ID',
    'i5 sequence': 'index2',
    'Project Name': 'Sample_Project',
    'synDNA pool number': 'syndna_pool_number'
}

_KL_METATRANSCRIPTOMIC_REMAPPER = {
    'sample sheet Sample_ID': 'Sample_ID',
    'Sample': 'Sample_Name',
    'Project Plate': 'Sample_Plate',
    'Well': 'well_id_384',
    'i7 name': 'I7_Index_ID',
    'i7 sequence': 'index',
    'i5 name': 'I5_Index_ID',
    'i5 sequence': 'index2',
    'Project Name': 'Sample_Project',
    'synDNA pool number': 'syndna_pool_number'
}

_AMPLICON = 'TruSeq HT'
_METAGENOMIC = 'Metagenomic'
_METATRANSCRIPTOMIC = 'Metatranscriptomic'
_ASSAYS = {_AMPLICON, _METAGENOMIC, _METATRANSCRIPTOMIC}

_READS = {
    'Read1': 151,
    'Read2': 151
}

_SETTINGS = {
    'ReverseComplement': '0',

    # these are needed ever since we moved from bcl2fastq -> bclconvert
    'MaskShortReads': '1',
    'OverrideCycles': 'Y151;I8N2;I8N2;Y151'
}


_HEADER = {
    'IEMFileVersion': '4',
    'Investigator Name': 'Knight',
    'Experiment Name': 'RKL_experiment',
    'Date': None,
    'Workflow': 'GenerateFASTQ',
    'Application': 'FASTQ Only',
    'Assay': None,
    'Description': '',
    'Chemistry': 'Default',
}

_BIOINFORMATICS_COLUMNS = {
    'Sample_Project', 'QiitaID', 'BarcodesAreRC', 'ForwardAdapter',
    'ReverseAdapter', 'HumanFiltering',
    'library_construction_protocol', 'experiment_design_description'
}

_CONTACT_COLUMNS = {
    'Sample_Project', 'Email'
}

_BIOINFORMATICS_AND_CONTACT = {
    'Bioinformatics': None,
    'Contact': None
}

_ALL_METADATA = {**_HEADER, **_SETTINGS, **_READS,
                 **_BIOINFORMATICS_AND_CONTACT}


class KLSampleSheet(sample_sheet.SampleSheet):
    def __init__(self, path=None):
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
        # the data
        super().__init__()

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
                        and section_name not in _KL_SAMPLE_SHEET_SECTIONS
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
            if columns[i] in _KL_SAMPLE_SHEET_COLUMN_ALTS:
                # overwrite existing alternate name w/internal representation.
                columns[i] = _KL_SAMPLE_SHEET_COLUMN_ALTS[columns[i]]
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


def _validate_sample_sheet_metadata(metadata):
    msgs = []

    for req in ['Assay', 'Bioinformatics', 'Contact']:
        if req not in metadata:
            msgs.append(ErrorMessage('%s is a required attribute' % req))

    # if both sections are found, then check that all the columns are present,
    # checks for the contents are done in the sample sheet validation routine
    if 'Bioinformatics' in metadata and 'Contact' in metadata:
        for section in ['Bioinformatics', 'Contact']:
            if section == 'Bioinformatics':
                columns = _BIOINFORMATICS_COLUMNS
            else:
                columns = _CONTACT_COLUMNS

            for i, project in enumerate(metadata[section]):
                if set(project.keys()) != columns:
                    message = (('In the %s section Project #%d does not have '
                               'exactly these keys %s') %
                               (section, i+1, ', '.join(sorted(columns))))
                    msgs.append(ErrorMessage(message))
                if section == 'Bioinformatics':
                    if (project['library_construction_protocol'] is None or
                            project['library_construction_protocol'] == ''):
                        message = (('In the %s section Project #%d does not '
                                    'have library_construction_protocol '
                                    'specified') % (section, i+1))
                        msgs.append(ErrorMessage(message))
                    if (project['experiment_design_description'] is None or
                            project['experiment_design_description'] == ''):
                        message = (('In the %s section Project #%d does not '
                                    'have experiment_design_description '
                                    'specified') % (section, i+1))
                        msgs.append(ErrorMessage(message))
    if metadata.get('Assay') is not None and metadata['Assay'] not in _ASSAYS:
        msgs.append(ErrorMessage('%s is not a supported Assay' %
                                 metadata['Assay']))

    keys = set(metadata.keys())
    if not keys.issubset(_ALL_METADATA):
        extra = sorted(keys - set(_ALL_METADATA))
        msgs.append(ErrorMessage('These metadata keys are not supported: %s'
                                 % ', '.join(extra)))

    return msgs


def _add_metadata_to_sheet(metadata, sheet, sequencer):
    # set the default to avoid index errors if only one of the two is provided
    sheet.Reads = [_READS['Read1'], _READS['Read2']]

    for key in _ALL_METADATA:
        if key in _READS:
            if key == 'Read1':
                sheet.Reads[0] = metadata.get(key, sheet.Reads[0])
            else:
                sheet.Reads[1] = metadata.get(key, sheet.Reads[1])

        elif key in _SETTINGS:
            sheet.Settings[key] = metadata.get(key, _SETTINGS[key])

        elif key in _HEADER:
            if key == 'Date':
                # we set the default value here and not in the global _HEADER
                # dictionary to make sure the date is when the metadata is
                # written not when the module is imported
                sheet.Header[key] = metadata.get(
                    key, datetime.today().strftime('%Y-%m-%d'))
            else:
                sheet.Header[key] = metadata.get(key, _HEADER[key])

        elif key in _BIOINFORMATICS_AND_CONTACT:
            setattr(sheet, key, pd.DataFrame(metadata[key]))

    # Per MacKenzie's request for 16S don't include Investigator Name and
    # Experiment Name
    if metadata['Assay'] == _AMPLICON:
        del sheet.Header['Investigator Name']
        del sheet.Header['Experiment Name']

        # these are only relevant for metagenomics because they are used in
        # bclconvert
        del sheet.Settings['MaskShortReads']
        del sheet.Settings['OverrideCycles']

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
        if 'MaskShortReads' in sheet.Settings:
            del sheet.Settings['MaskShortReads']
        if 'OverrideCycles' in sheet.Settings:
            del sheet.Settings['OverrideCycles']

    return sheet


def _remap_table(table, assay, strict=True):
    if assay == _AMPLICON:
        remapper = _KL_AMPLICON_REMAPPER
    elif assay == _METAGENOMIC:
        remapper = _KL_METAGENOMIC_REMAPPER
    elif assay == _METATRANSCRIPTOMIC:
        remapper = _KL_METATRANSCRIPTOMIC_REMAPPER
    else:
        raise ValueError(f"'{assay} is not a valid Assay type.")

    # Well_description column is now defined here as the concatenation
    # of the following columns. If the column existed previously it will
    # be overwritten, otherwise it will be created here. Alternate versions
    # of the column name have already been resolved at this point.

    # Note that the amplicon notebook currently generates the same values for
    # this column. If the functionality in the notebook changes, the output
    # will continue to be redfined with the current values here.
    well_description = table['Project Plate'].astype(str) + "." + \
        table['Sample'].astype(str) + "." + table['Well'].astype(str)

    if strict:
        # legacy operation. All columns not defined in remapper will be
        # filtered out.
        out = table[remapper.keys()].copy()
        out.rename(remapper, axis=1, inplace=True)
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
        # remapper (they currently do not). Define the full set of potential
        # columns to rename in table.

        # new syntax in 3.9 allows us to merge two dicts together w/OR.
        remapper = _KL_SAMPLE_SHEET_COLUMN_ALTS | remapper
        out.rename(remapper, axis=1, inplace=True)

        # out may contain additional columns that aren't allowed in the [Data]
        # section of a sample-sheet e.g.: 'Extraction Kit Lot'. There may also
        # be required columns that aren't defined in out.
        subset = list(set(_KL_SAMPLE_SHEET_DATA_COLUMNS) & set(out.columns))
        out = out[subset]

    # append the new 'Well_description' column, now that alternates have been
    # removed and non-essential columns have been dropped.
    out['Well_description'] = well_description

    for column in _KL_SAMPLE_SHEET_DATA_COLUMNS:
        if column not in out.columns:
            warnings.warn('The column %s in the sample sheet is empty' %
                          column)
            out[column] = ''

    return out


def _add_data_to_sheet(table, sheet, sequencer, lanes, assay, strict=True):
    table = _remap_table(table, assay, strict)

    if assay != _AMPLICON:
        table['index2'] = sequencer_i5_index(sequencer, table['index2'])

        sheet.Bioinformatics['BarcodesAreRC'] = str(
            sequencer in REVCOMP_SEQUENCERS)

    for lane in lanes:
        for sample in table.to_dict(orient='records'):
            sample['Lane'] = lane
            sheet.add_sample(sample_sheet.Sample(sample))

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
        DNA pool number ('synDNA pool number').
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
        SampleSheet object containing both the metadata and table.

    Raises
    ------
    SampleSheetError
        If one of the required columns is missing.
    """
    messages = _validate_sample_sheet_metadata(metadata)

    if len(messages) == 0:
        sheet = KLSampleSheet()
        sheet = _add_metadata_to_sheet(metadata, sheet, sequencer)
        sheet = _add_data_to_sheet(table, sheet, sequencer, lanes,
                                   metadata['Assay'], strict)
        return sheet
    else:
        for message in messages:
            message.echo()


def quiet_validate_and_scrub_sample_sheet(sheet):
    """Quietly validate the sample sheet and scrub invalid characters

    The character scrubbing is only applied to the Sample_Project and the
    Sample_ID columns.

    Parameters
    ----------
    sheet: sample_sheet.KLSampleSheet
        The sample sheet object to validate and scrub.

    Returns
    -------
    list
        List of error or warning messages.
    sample_sheet.SampleSheet or None
        Corrected and validated sample sheet if no errors are found. Otherwise
        None is returned.
    """
    msgs = []

    # we print an error return None and exit when this happens otherwise we
    # won't be able to run some of the other checks
    for column in _KL_SAMPLE_SHEET_DATA_COLUMNS:
        if column not in sheet.all_sample_keys:
            msgs.append(
                ErrorMessage('The %s column in the Data section is missing' %
                             column))
    for section in ['Bioinformatics', 'Contact']:
        if getattr(sheet, section) is None:
            msgs.append(ErrorMessage('The %s section cannot be empty' %
                                     section))

    # if any errors are found up to this point then we can't continue with the
    # validation
    if msgs:
        return msgs, None

    # we track the updated projects as a dictionary so we can propagate these
    # changes to the Bioinformatics and Contact sections
    updated_samples, updated_projects = [], {}
    for sample in sheet.samples:
        new_sample = bcl_scrub_name(sample.Sample_ID)
        new_project = bcl_scrub_name(sample.Sample_Project)

        if new_sample != sample.Sample_ID:
            updated_samples.append(sample.Sample_ID)
            sample.Sample_ID = new_sample
        if new_project != sample.Sample_Project:
            updated_projects[sample.Sample_Project] = new_project
            sample['Sample_Project'] = new_project

    if updated_samples:
        msgs.append(
            WarningMessage('The following sample names were scrubbed for'
                           ' bcl2fastq compatibility:\n%s' %
                           ', '.join(updated_samples)))
    if updated_projects:
        msgs.append(
            WarningMessage('The following project names were scrubbed for'
                           ' bcl2fastq compatibility. If the same invalid '
                           'characters are also found in the Bioinformatics '
                           'and Contacts sections those will be automatically '
                           'scrubbed too:\n%s' %
                           ', '.join(sorted(updated_projects))))

        # make the changes to prevent useless errors where the scurbbed names
        # fail to match between sections.
        sheet.Contact.Sample_Project.replace(updated_projects, inplace=True)
        sheet.Bioinformatics.Sample_Project.replace(updated_projects,
                                                    inplace=True)

    pairs = collections.Counter([(s.Lane, s.Sample_Project)
                                 for s in sheet.samples])
    # warn users when there's missing lane values
    empty_projects = [project for lane, project in pairs
                      if str(lane).strip() == '']
    if empty_projects:
        msgs.append(
            ErrorMessage('The following projects are missing a Lane value: '
                         '%s' % ', '.join(sorted(empty_projects))))

    projects = {s.Sample_Project for s in sheet.samples}

    # project identifiers are digit groups at the end of the project name
    # preceded by an underscore CaporasoIllumina_550
    qiita_id_re = re.compile(r'(.+)_(\d+)$')
    bad_projects = []
    for project_name in projects:
        if re.search(qiita_id_re, project_name) is None:
            bad_projects.append(project_name)
    if bad_projects:
        msgs.append(
            ErrorMessage('The following project names in the Sample_Project '
                         'column are missing a Qiita study '
                         'identifier: %s' % ', '.join(sorted(bad_projects))))

    # check Sample_project values match across sections
    bfx = set(sheet.Bioinformatics['Sample_Project'])
    contact = set(sheet.Contact['Sample_Project'])
    not_shared = projects ^ bfx
    if not_shared:
        msgs.append(
            ErrorMessage('The following projects need to be in the Data and '
                         'Bioinformatics sections %s' %
                         ', '.join(sorted(not_shared))))
    elif not contact.issubset(projects):
        msgs.append(
            ErrorMessage(('The following projects were only found in the '
                          'Contact section: %s projects need to be listed in '
                          'the Data and Bioinformatics section in order to '
                          'be included in the Contact section.') %
                         ', '.join(sorted(contact - projects))))

    # if there are no error messages then return the sheet
    if not any([isinstance(m, ErrorMessage) for m in msgs]):
        return msgs, sheet
    else:
        return msgs, None


def validate_and_scrub_sample_sheet(sheet):
    """Validate the sample sheet and scrub invalid characters

    The character scrubbing is only applied to the Sample_Project and the
    Sample_ID columns. The main difference between this function and
    quiet_validate_and_scrub_sample_sheet is that this function will *always*
    print errors and warnings to standard output.

    Parameters
    ----------
    sheet: sample_sheet.KLSampleSheet
        The sample sheet object to validate and scrub.

    Returns
    -------
    sample_sheet.SampleSheet
        Corrected and validated sample sheet if no errors are found.
    """
    msgs, sheet = quiet_validate_and_scrub_sample_sheet(sheet)

    [msg.echo() for msg in msgs]

    if sheet is not None:
        return sheet


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
        new_sheet = KLSampleSheet()
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
        for sample in df.to_dict(orient='records'):
            new_sheet.add_sample(sample_sheet.Sample(sample))

        demuxed_sheets.append(new_sheet)

    return demuxed_sheets
