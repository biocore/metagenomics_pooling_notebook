import re
import csv
import tempfile
import collections
import warnings

from datetime import datetime
from itertools import chain, repeat, islice

import sample_sheet
import pandas as pd

from metapool.metapool import bcl_scrub_name

_KL_SECTIONS = ['Header', 'Reads', 'Settings', 'Data', 'Bioinformatics',
                'Contact']


class KLSampleSheet(sample_sheet.SampleSheet):
    """Knight Lab's SampleSheet subclass

    Includes changes to the write method so the [Bioinformatics] and [Contact]
    sections can be found at the end of the document and follow the same
    standard as the [Data] section.
    """
    def __init__(self, path):
        super().__init__()

        # print('something is happening')

        self.Bioinformatics = None
        self.Contact = None
        self.path = path

        if self.path:
            self._parse(self.path)

        # TODO: Maybe check that bioinformatics and contact aren't none

    def _parse(self, path):
        # print('i am parsing')
        section_name = ''
        section_header = None

        # print('who')
        with open(path, encoding=self._encoding) as handle:
            lines = list(csv.reader(handle, skipinitialspace=True))

        # print('almost')
        for i, line in enumerate(lines):
            # print(i)
            # print('this is self', self)
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
                    and section_name not in _KL_SECTIONS
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
                line = [value for value in line if value != '']
                if getattr(self, section_name) is not None:
                    # add row to the dataframe
                    df = getattr(self, section_name)
                    df.loc[len(df)] = line

                    # we need a dropna or something like that
                    # elif any(key == '' for key in line):
                    #     raise ValueError(
                    #         f'Header for [{section_name}] section allowed '
                    #         f'to have empty fields: {line}'
                    #     )
                else:
                    setattr(self, section_name, pd.DataFrame(columns=line))
                continue

            # [<Other>]
            else:
                key, value, *_ = line
                section = getattr(self, section_name)
                section[key] = value
                continue

    # copied from release 0.11.0
    # https://github.com/clintval/sample-sheet/blob/
    # c3df258541a384a5058f8aa46b343ff032d8e247/sample_sheet/__init__.py
    def write(self, handle, blank_lines=1) -> None:
        """Write to a file-like object.

        Parameters
        ----------

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
    # Information for the Bioinformatics and Contact section goes here
    'sequencer': None,
    'lanes': None,
    'Sample_Projects': None
}


def _validate_sample_sheet_metadata(metadata, table):
    # check there's an exact match of keys
    # check assay is provided
    # check sample projects match with the stuff present in the table
    # check lanes are correct
    # check the sequencer is known

    return metadata


def _add_metadata_to_sheet(metadata, sheet):

    # combined dictionaries are useful for setting defaults and updating
    # whenever needed
    for key, value in _HEADER.update(_READS).update(_SETTINGS).items():

        if key in _READS:
            if sheet.Reads is None:
                sheet.Reads = [_READS['Read1'], _READS['Read2']]

                sheet.Reads[0] = metadata.get('Read1', sheet.Reads[0])
                sheet.Reads[1] = metadata.get('Read2', sheet.Reads[1])

        elif key in _SETTINGS:
            sheet.Settings[key]: metadata.get(key, _SETTINGS[key])

        elif key in _HEADER:
            if key == 'Date':
                sheet.Header[key] = metadata.get(
                    key, datetime.today().strftime('%Y-%m-%d'))
            else:
                sheet.Header[key]: metadata.get(key, _HEADER[key])

        elif key in _BIOINFORMATICS_AND_CONTACT:
            pass

    return sheet


def make_sample_sheet(metadata, table):
    """Write a valid sample sheet

    Parameters
    ----------
    metadata: dict
        Metadata describing the sample sheet with the following fileds.
        When a field is omitted the values in square brackets are used
        instead.

        - sequencer: ??
        - lanes: ??

        - e-mails: comma-separated list of e-mail addresses to notifyx
          [jdereus@ucsd.edu, ackermag@ucsd.edu, mmbryant@ucsd.edu]
        - PI: Principal investigator's email [robknight@ucsd.edu]

        - IEMFileVersion: Illumina's Experiment Manager version [4]
        - Investigator Name: [Knight]
        - Experiment Name: [RKL_experiment]
        - Date: Date when the sheet is prepared [Today's date]
        - Workflow: how the sample sheet should be used [GenerateFASTQ]
        - Application: sample sheet's application [FASTQ Only]
        - Assay: assay type for the sequencing run [Metagenomics]
        - Description: additional information []
        - Chemistry: chemistry's description [Default]
        - read1: Length of forward read [151]
        - read2: Length of forward read [151]
        - ReverseComplement: If the barcodes should be reverse complemented
          by bcl2fastq [0]
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
    # sheet = add_data_to_sheet(table, sheet)

    return sheet


def write_sample_sheet(sheet, path):
    """Writes a sample sheet container object to a path

    Parameters
    ----------
    sheet: sample_sheet.SampleSheet
        The sheet to be saved.
    path: str
        Path where the sample sheet will be saved to
    """
    if not isinstance(sheet, sample_sheet.SampleSheet):
        raise ValueError('The input sample sheet should be a SampleSheet '
                         'instance')

    with open(path, 'w') as f:
        # the SampleSheet class does not support writing comments
        f.write(sheet.comments)

        sheet.write(f)


def parse_sample_sheet(file_like):
    """Parse a sample sheet with comments

    Comments are not defined as part of Illumina's sample sheet specification.
    Hence this function explicitly handles comments by stripping them and
    parsing the rest of the contents using the `sample_sheet` package.

    Parameters
    ----------
    file_like: file path or open filehandle
        Object pointing to the sample sheet.

    Returns
    -------
    sample_sheet.SampleSheet
        An object with the sample sheet contents and comments added to a custom
        attribute `.comments`.
    """

    # read comments first
    comments = []
    with open(file_like) as f, tempfile.NamedTemporaryFile(mode='w+') as temp:
        for line in f:
            if line.startswith('# '):
                comments.append(line)
            else:
                temp.write(line)

        # return to the beginning of the file so contents can be parsed
        temp.seek(0)

        # important to parse before leaving the context manager
        sheet = sample_sheet.SampleSheet(temp.name)

        # save the comments as a custom attribute
        sheet.comments = ''.join(comments)

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
        If the sample sheet object has not comments attribute.
        If duplicated elements are found in the Sample_ID column after name
        scrubbing.
    UserWarning
        When sample identifiers are scrubbed for bcl2fastq compatibility.
        If project names don't include a Qiita study identifier.
    """

    if not hasattr(sheet, 'comments') or sheet.comments == '':
        raise ValueError('This sample sheet does not include comments, these '
                         'are required to notify interested parties upon '
                         'completed processing of the data')

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

    return sheet


def sample_sheet_to_dataframe(sheet):
    """Converts the sample section of a sample sheet into a DataFrame

    Parameters
    ----------
    sheet: sample_sheet.SampleSheet
        Object from where to extract the data.

    Returns
    -------
    pd.DataFrame
        DataFrame object with the sample information.
    """

    # Get the columns names for the first sample so we have them in a list and
    # we can retrieve data in the same order on every iteration
    columns = list(sheet.samples[0].keys())

    data = []
    for sample in sheet.samples:
        data.append([sample[column] for column in columns])

    out = pd.DataFrame(data=data, columns=[c.lower() for c in columns])

    return out.set_index('sample_id')


# The old stuff
def ss_temp():
    """Sample sheet template
    """
    s = ('{comments}[Header]\n'
         'IEMFileVersion{sep}{IEMFileVersion}\n'
         'Investigator Name{sep}{Investigator Name}\n'
         'Experiment Name{sep}{Experiment Name}\n'
         'Date{sep}{Date}\n'
         'Workflow{sep}{Workflow}\n'
         'Application{sep}{Application}\n'
         'Assay{sep}{Assay}\n'
         'Description{sep}{Description}\n'
         'Chemistry{sep}{Chemistry}\n\n'
         '[Reads]\n'
         '{read1}\n'
         '{read2}\n\n'
         '[Settings]\n'
         'ReverseComplement{sep}{ReverseComplement}\n\n'
         '[Data]\n'
         '{data}')

    return s


def format_sheet_comments(PI=None, contacts=None, other=None, sep=','):

    comments = ''

    if PI is not None:
        comments += 'PI{0}{1}\n'.format(
            sep,
            sep.join('{0}{1}{2}'.format(x, sep, PI[x]) for x in PI.keys()))

    if contacts is not None:
        comments += 'Contact{0}{1}\n{0}{2}\n'.format(
            sep,
            sep.join(x for x in sorted(contacts.keys())),
            sep.join(contacts[x] for x in sorted(contacts.keys())))

    if other is not None:
        comments += '%s\n' % other

    return(comments)


def format_sample_sheet(sample_sheet_dict, sep=',', template=ss_temp()):
    """Formats Illumina-compatible sample sheet.

    Parameters
    ----------
    sample_sheet_dict : 2-level dict
        dict with 1st level headers 'Comments', 'Header', 'Reads', 'Settings',
        and 'Data'.

    Returns
    -------
    sample_sheet : str
        the sample sheet string
    """
    if sample_sheet_dict['comments']:
        sample_sheet_dict['comments'] = re.sub(
            '^', '# ', sample_sheet_dict['comments'].rstrip(),
            flags=re.MULTILINE) + '\n'

    sample_sheet = template.format(**sample_sheet_dict, **{'sep': sep})

    return(sample_sheet)


def format_sample_data(sample_ids, i7_name, i7_seq, i5_name, i5_seq,
                       sample_plate, sample_proj, wells=None,
                       description=None, lanes=[1], sep=','):
    """
    Creates the [Data] component of the Illumina sample sheet from plate
    information

    Parameters
    ----------
    sample_sheet_dict : 2-level dict
        dict with 1st level headers 'Comments', 'Header', 'Reads', 'Settings',
        and 'Data'.

    Returns
    -------
    data : str
        the sample sheet string
    """
    data = ''

    if (len(sample_ids) != len(i7_name) != len(i7_seq) != len(i5_name) !=
       len(i5_seq)):
        raise ValueError('Sample information lengths are not all equal')

    if wells is None:
        wells = [''] * len(sample_ids)
    if description is None:
        description = [''] * len(sample_ids)
    if isinstance(sample_plate, str):
        sample_plate = [sample_plate] * len(sample_ids)
    if isinstance(sample_proj, str):
        sample_proj = [sample_proj] * len(sample_ids)

    header = ','.join(['Lane', 'Sample_ID', 'Sample_Name', 'Sample_Plate',
                       'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID',
                       'index2', 'Sample_Project', 'Description'])

    data += header

    sample_plate = list(sample_plate)
    wells = list(wells)
    i7_name = list(i7_name)
    i7_seq = list(i7_seq)
    i5_name = list(i5_name)
    i5_seq = list(i5_seq)
    sample_proj = list(sample_proj)
    description = list(description)

    for lane in lanes:
        for i, sample in enumerate(sample_ids):
            line = sep.join([str(lane),
                             sample,
                             sample,
                             sample_plate[i],
                             wells[i],
                             i7_name[i],
                             i7_seq[i],
                             i5_name[i],
                             i5_seq[i],
                             sample_proj[i],
                             description[i]])
            data += '\n' + line

    return(data)
