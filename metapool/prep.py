import re
import os
import gzip
import warnings
import pandas as pd

from glob import glob
from datetime import datetime


PREP_COLUMNS = ['sample_name', 'experiment_design_description',
                'library_construction_protocol', 'platform', 'run_center',
                'run_date', 'run_prefix', 'sequencing_meth', 'center_name',
                'center_project_name', 'instrument_model', 'runid',
                'sample_plate', 'sample_well', 'i7_index_id', 'index',
                'i5_index_id', 'index2', 'lane', 'sample_project',
                'well_description']

# put together by Gail, based on the instruments we know of
INSTRUMENT_LOOKUP = pd.DataFrame({
    'A00953': {'machine prefix': 'A', 'Vocab': 'Illumina NovaSeq',
               'Machine type': 'NovaSeq', 'run_center': 'IGM'},
    'A00169': {'machine prefix': 'A', 'Vocab': 'Illumina NovaSeq',
               'Machine type': 'NovaSeq', 'run_center': 'LJI'},
    'M05314': {'machine prefix': 'M', 'Vocab': 'Illumina MiSeq',
               'Machine type': 'MiSeq', 'run_center': 'KLM'},
    'K00180': {'machine prefix': 'K', 'Vocab': 'Illumina HiSeq 4000',
               'Machine type': 'HiSeq', 'run_center': 'IGM'},
    'D00611': {'machine prefix': 'D', 'Vocab': 'Illumina HiSeq 2500',
               'Machine type': 'HiSeq/RR', 'run_center': 'IGM'},
    'MN01225': {'machine prefix': 'MN', 'Vocab': 'Illumina MiniSeq',
                'Machine type': 'MiniSeq', 'run_center': 'CMI'}}).T


def parse_illumina_run_id(run_id):
    """Parse a run identifier

    Parameters
    ----------
    run_id: str
        The name of a run

    Returns
    -------
    str:
        When the run happened (YYYY-MM-DD)
    str:
        Instrument code
    """

    # of the form
    # YYMMDD_machinename_XXXX_FC
    tokens = run_id.split('_', 1)

    if len(tokens) != 2:
        raise ValueError('Unrecognized run identifier format "%s". The '
                         'expected format is YYMMDD_machinename_XXXX_FC.' %
                         run_id)

    # convert illumina's format to qiita's format
    run_date = datetime.strptime(tokens[0], '%y%m%d').strftime('%Y-%m-%d')

    return run_date, tokens[1]


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

    # used to get the right order every time
    columns = list(sheet.samples[0].keys())

    data = []
    for sample in sheet.samples:
        data.append([sample[column] for column in columns])

    out = pd.DataFrame(data=data, columns=[c.lower() for c in columns])

    return out.set_index('sample_id')


def is_nonempty_gz_file(name):
    """Taken from https://stackoverflow.com/a/37878550/379593"""
    with gzip.open(name, 'rb') as f:
        try:
            file_content = f.read(1)
            return len(file_content) > 0
        except Exception:
            return False


def remove_qiita_id(project_name):
    # project identifiers are digit groups at the end of the project name
    # preceded by an underscore CaporasoIllumina_550
    qiita_id_re = re.compile(r'(.+)_(\d+)$')

    # no matches
    matches = re.search(qiita_id_re, project_name)
    if matches is None:
        return project_name
    else:
        # group 1 is the project name
        return matches[1]


def get_run_prefix(run_path, project, sample, lane):
    """For a sample find the run prefix

    Parameters
    ----------
    run_path: str
        Base path for the run
    project: str
        Name of the project
    sample: str
        Sample name
    lane: str
        Lane number

    Returns
    -------
    str
        The run prefix of the sequence file in the lane, only if the sequence
        file is not empty.
    """
    base = os.path.join(run_path, project)
    path = base

    # check for the existence of the qc or human filtering folders
    qc = os.path.join(base, 'atropos_qc')
    if os.path.exists(qc) and len(_file_list(qc)):
        path = qc

        human_filtering = os.path.join(base, 'filtered_sequences')
        if (os.path.exists(human_filtering) and
           len(_file_list(human_filtering))):
            path = human_filtering

    results = glob(os.path.join(path, '%s_S*_L*%s_R*.fastq.gz' %
                                (sample, lane)))

    # at this stage there should only be two files forward and reverse
    if len(results) == 2:
        forward, reverse = sorted(results)

        if is_nonempty_gz_file(forward) and is_nonempty_gz_file(reverse):
            return os.path.basename(forward.split('_R1_')[0])
        else:
            return None

    else:
        return None


def _file_list(path):
    return [f for f in os.listdir(path) if not os.path.isdir(os.path.join(path, f))]


def get_machine_code(instrument_model):
    """Get the machine code for an instrument's code

    Parameters
    ----------
    instrument_model: str
        An instrument's model of the form A999999 or AA999999

    Returns
    -------
    """
    # the machine code represents the first 1 to 2 letters of the
    # instrument model
    machine_code = re.compile(r'^([a-zA-Z]{1,2})')
    matches = re.search(machine_code, instrument_model)
    if matches is None:
        raise ValueError('Cannot find a machine code. This instrument '
                         'model is malformed %s. The machine code is a '
                         'one or two character prefix.' % instrument_model)
    else:
        return matches[0]


def get_model_and_center(instrument_code):
    """Determine instrument model and center based on a lookup

    Parameters
    ----------
    instrument_code: str
        Instrument code from a run identifier.

    Returns
    -------
    str
        Instrument model.
    str
        Run center based on the machine's id.
    """
    run_center = "UCSDMI"
    instrument_model = instrument_code.split('_')[0]

    if instrument_model in INSTRUMENT_LOOKUP.index:
        instrument_model = INSTRUMENT_LOOKUP.loc[instrument_model, 'Vocab']
        run_center = INSTRUMENT_LOOKUP.loc[instrument_model, 'run_center']
    else:
        instrument_prefix = get_machine_code(instrument_model)

        if instrument_prefix not in INSTRUMENT_LOOKUP['machine prefix']:
            ValueError('Unrecognized machine prefix %s' % instrument_prefix)

        instrument_model = INSTRUMENT_LOOKUP[
            INSTRUMENT_LOOKUP['machine prefix'] == instrument_prefix
        ]['Vocab'].unique()[0]

    return instrument_model, run_center


def preparations_for_run(run_path, sheet):
    """Given a run's path and sample sheet generates preparation files

    Parameters
    ----------
    run_path: str
        Path to the run folder
    sheet: sample_sheet.SampleSheet
        Sample sheet to convert

    Returns
    -------
    dict
        Filename to preparation file dictionary. Preparation files are
        represented as DataFrames.
    """

    _, run_id = os.path.split(os.path.normpath(run_path))
    run_date, instrument_code = parse_illumina_run_id(run_id)
    instrument_model, run_center = get_model_and_center(instrument_code)

    output = {}

    required_columns = {'well_description', 'sample_plate',
                        'sample_well', 'i7_index_id', 'index',
                        'i5_index_id', 'index2'}

    # TODO: How should we handle this stuff?
    not_present = required_columns - set(sheet.columns)
    if not_present:
        warnings.warn('These required columns were not found %s'
                      % ', '.join(not_present), UserWarning)
        for col in not_present:
            sheet[col] = 'MISSING_FROM_THE_SAMPLE_SHEET'

    for project, project_sheet in sheet.groupby('sample_project'):

        project_name = remove_qiita_id(project)
        qiita_id = project.replace(project_name + '_', '')

        # if the Qiita ID is not found then make for an easy find/replace
        if qiita_id == project:
            qiita_id = 'QIITA-ID'

        for lane, lane_sheet in project_sheet.groupby('lane'):

            # this is the portion of the loop that creates the prep
            data = []
            for sample_name, sample in lane_sheet.iterrows():
                run_prefix = get_run_prefix(run_path, project, sample_name,
                                            lane)

                # we don't care about the sample if there's no file
                if run_prefix is None:
                    continue

                row = {c: '' for c in PREP_COLUMNS}

                row["sample_name"] = sample.well_description
                row["experiment_design_description"] = "EXPERIMENT_DESC"
                row["library_construction_protocol"] = "LIBRARY_PROTOCOL"
                row["platform"] = "Illumina"
                row["run_center"] = "UCSDMI"
                row["run_date"] = run_date
                row["run_prefix"] = run_prefix
                row["sequencing_meth"] = "sequencing by synthesis"
                row["center_name"] = "CENTER_NAME"
                row["center_project_name"] = project_name

                row["instrument_model"] = instrument_model
                row["runid"] = run_id
                row["sample_plate"] = sample.sample_plate
                row["sample_well"] = sample.sample_well
                row["i7_index_id"] = sample['i7_index_id']
                row["index"] = sample['index']
                row["i5_index_id"] = sample['i5_index_id']
                row["index2"] = sample['index2']
                row["lane"] = lane
                row["sample_project"] = project
                row["well_description"] = '%s.%s.%s' % (sample.sample_plate,
                                                        sample.sample_name,
                                                        sample.sample_well)

                data.append(row)

            if not data:
                warnings.warn('Project %s and Lane %s have no data' %
                              (project, lane), UserWarning)

            # label the prep based on the run, project and lane
            name = run_id + '.' + project + '.' + lane
            output[name] = pd.DataFrame(columns=PREP_COLUMNS, data=data)

    return output
