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

_KL_SAMPLE_SHEET_SECTIONS = [
    'Header', 'Reads', 'Settings', 'Data', 'Bioinformatics', 'Contact'
]

_KL_SAMPLE_SHEET_DATA_COLUMNS = [
    'Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well', 'I7_Index_ID',
    'index', 'I5_Index_ID', 'index2', 'Sample_Project', 'Well_description'
]

_KL_AMPLICON_REMAPPER = {
    'sample sheet Sample_ID': 'Sample_ID',
    'Sample': 'Sample_Name',
    'Project Plate': 'Sample_Plate',
    'Well': 'Sample_Well',
    'Name': 'I7_Index_ID',
    'Golay Barcode': 'index',
    'Project Name': 'Sample_Project',
}

_KL_METAGENOMICS_REMAPPER = {
    'sample sheet Sample_ID': 'Sample_ID',
    'Sample': 'Sample_Name',
    'Project Plate': 'Sample_Plate',
    'Well': 'Sample_Well',
    'i7 name': 'I7_Index_ID',
    'i7 sequence': 'index',
    'i5 name': 'I5_Index_ID',
    'i5 sequence': 'index2',
    'Project Name': 'Sample_Project',
}

_AMPLICON = 'Amplicon'
_METAGENOMICS = 'Metagenomics'
_ASSAYS = {_AMPLICON, _METAGENOMICS}

_READS = {
    'Read1': 151,
    'Read2': 151
}

_SETTINGS = {
    'ReverseComplement': 0
}

_HEADER = {
    'IEMFileVersion': 4,
    'Investigator Name': 'Knight',
    'Experiment Name': 'RKL_experiment',
    'Date': None,
    'Workflow': 'GenerateFASTQ',
    'Application': 'FASTQ Only',
    'Assay': None,
    'Description': '',
    'Chemistry': 'Default',
}

_BIOINFORMATICS_AND_CONTACT = {
    'Bioinformatics': None,
    'Contact': None
}


class KLSampleSheet(sample_sheet.SampleSheet):
    def __init__(self, path=None):
        """Knight Lab's SampleSheet subclass

        Includes a number of parsing and writing changes to allow for the
        inclusion of the Bioinformatics and Contact sections.

        The majority of the code in the _parse and write methods are copied
        from release 0.11.0 https://github.com/clintval/sample-sheet/blob/\
c3df258541a384a5058f8aa46b343ff032d8e247/sample_sheet/__init__.py

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

        with open(path, encoding=self._encoding) as handle:
            lines = list(csv.reader(handle, skipinitialspace=True))

        for i, line in enumerate(lines):
            # Skip to next line if this line is empty to support formats of
            # sample sheets with multiple newlines as section seperators.
            #
            #   https://github.com/clintval/sample-sheet/issues/46
            #
            if not ''.join(line).strip():
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

            # If we enter a section save it's name and continue to next line.
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
                        sample_sheet.Sample(dict(zip(section_header, line))))
                elif any(key == '' for key in line):
                    raise ValueError(
                        f'Header for [Data] section is not allowed to '
                        f'have empty fields: {line}'
                    )
                else:
                    section_header = line
                continue

            elif section_name in {'Bioinformatics', 'Contact'}:
                # CSV rows are padded to include commas for the longest line in
                # the file, so we remove them to avoid creating empty columns
                line = [value for value in line if value != '']

                if getattr(self, section_name) is not None:
                    df = getattr(self, section_name)
                    df.loc[len(df)] = line
                else:
                    setattr(self, section_name, pd.DataFrame(columns=line))
                continue

            # [<Other>]
            else:
                key, value, *_ = line
                section = getattr(self, section_name)
                section[key] = value
                continue

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
        csv_width = max([len(self.all_sample_keys),
                         len(self.Bioinformatics.columns),
                         len(self.Contact.columns), 2])

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

        Parameters
        ----------
        sheets: list of KLSampleSheet
            The sample sheets to merge into `self`.

        Returns
        -------
        KLSampleSheet
            The combined sample sheet.
        """
        for sheet in sheets:
            for sample in sheet.samples:
                self.add_sample(sample)

            # these two sections are data frames
            for section in ['Bioinformatics', 'Contact']:
                if (getattr(self, section) is not None and
                   getattr(sheet, section) is not None):
                    section = getattr(self, section)

                    for _, row in getattr(sheet, section).iterrows():
                        section.loc[len(section)] = row


def _validate_sample_sheet_metadata(metadata, table):
    # check there's a subset of keys
    # check assay is not None in the metadata
    return metadata


def _add_metadata_to_sheet(metadata, sheet):
    # set the default to avoid index errors if only one of the two is provided
    sheet.Reads = [_READS['Read1'], _READS['Read2']]

    # combine all the dictionaries with the default values
    combined = {**_HEADER, **_SETTINGS, **_READS,
                **_BIOINFORMATICS_AND_CONTACT}

    for key in combined:
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

    return sheet


def _remap_table(table, assay):
    if assay == _AMPLICON:
        remapper = _KL_AMPLICON_REMAPPER
    elif assay == _METAGENOMICS:
        remapper = _KL_METAGENOMICS_REMAPPER

    # make a copy because we are going to modify the data
    table = table[remapper.keys()].copy()

    table.rename(remapper, axis=1, inplace=True)

    for column in _KL_SAMPLE_SHEET_DATA_COLUMNS:
        if column not in table.columns:
            warnings.warn('The column %s in the sample sheet is empty' %
                          column)
            table[column] = ''

    return table


def _add_data_to_sheet(table, sheet, sequencer, lanes, assay):
    table = _remap_table(table, assay)

    # for amplicon we don't have reverse barcodes
    if assay == _METAGENOMICS:
        table['index2'] = sequencer_i5_index(sequencer, table['index2'])

    sheet.Bioinformatics['BarcodesAreRC'] = sequencer in REVCOMP_SEQUENCERS

    for lane in lanes:
        for sample in table.to_dict(orient='records'):
            sample['Lane'] = lane
            sheet.add_sample(sample_sheet.Sample(sample))

    return sheet


def make_sample_sheet(metadata, table, sequencer, lanes):
    """Write a valid sample sheet

    Parameters
    ----------
    metadata: dict
        Metadata describing the sample sheet with the following fileds.
        If a value is omitted from this dictionary the values in square
        brackets are used as defaults.

        - Bioinformatics: List of dictionaries describing each project's
          attributes: Sample_Project, QiitaID, BarcodesAreRC, ForwardAdapter,
          ReverseAdapter, HumanFiltering
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
        plate ('Project Plate'), and project name ('Project Name').

    Returns
    -------
    samplesheet.SampleSheet
        SampleSheet object containing both the metadata and table.

    Raises
    ------
    SampleSheetError
        If one of the required columns is missing.
    """
    metadata = _validate_sample_sheet_metadata(metadata, table)

    sheet = KLSampleSheet()

    sheet = _add_metadata_to_sheet(metadata, sheet)

    # make sure that we can use the lanes attribute as expected
    sheet = _add_data_to_sheet(table, sheet, sequencer, lanes,
                               metadata['Assay'])

    return sheet


def validate_sample_sheet(sheet):
    """Validate and correct some aspects of the sample sheet

    Parameters
    ----------
    sheet: sample_sheet.SampleSheet
        The sample sheet container as parsed from `parse_sample_sheet`.

    Returns
    -------
    sample_sheet.SampleSheet
        Corrected and validated sample sheet.

    Raises
    ------
    ValueError
        If the Sample_ID column or the Sample_project columns are missing.
        If duplicated elements are found in the Sample_ID column after name
        scrubbing.
    UserWarning
        When sample identifiers are scrubbed for bcl2fastq compatibility.
        If project names don't include a Qiita study identifier.
    """

    if sheet.samples[0].Sample_project is None:
        raise ValueError('The Sample_project column in the Data section is '
                         'missing')

    corrected_names = []
    for sample in sheet.samples:
        new = bcl_scrub_name(sample.Sample_ID)
        if new != sample.Sample_ID:
            corrected_names.append(sample.Sample_ID)
            sample.Sample_ID = new

    if corrected_names:
        warnings.warn('The following sample names were scrubbed for bcl2fastq'
                      ' compatibility:\n%s' % ', '.join(corrected_names))

    # based on Illumina's recommendation Sample_ID is the required column to
    # name samples
    samples = collections.Counter([s.Sample_ID for s in sheet.samples])

    # check all values are unique
    duplicates = {k: v for k, v in samples.items() if v > 1}

    if duplicates:
        names = '\n'.join(['%s: %d' % (k, v) for k, v in duplicates.items()])
        message = ('The following names in the Sample_ID column are listed '
                   'multiple times:\n%s' % names)

        # make it clear if the sample collisions are caused by the scrubbing
        if corrected_names:
            message = ('After scrubbing samples for bcl2fastq compatibility '
                       'there are repeated identifiers. ' + message)

        raise ValueError(message)
    lanes = collections.Counter([s.Lane for s in sheet.samples])
    # warn users when there's missing lane values
    if any([lane.strip() == '' for lane in lanes]):
        warnings.warn('The following projects are missing a Lane value')

    projects = collections.Counter([s.Sample_project for s in sheet.samples])

    # project identifiers are digit groups at the end of the project name
    # preceded by an underscore CaporasoIllumina_550
    qiita_id_re = re.compile(r'(.+)_(\d+)$')
    bad_projects = []
    for project_name in projects:
        if re.search(qiita_id_re, project_name) is None:
            bad_projects.append(project_name)
    if bad_projects:
        warnings.warn('The following project names in the Sample_project '
                      'column are missing a Qiita study '
                      'identifier: %s' % ', '.join(sorted(bad_projects)))

    # check sample projects match with the stuff present in the table

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

    return out.set_index('sample_id')
