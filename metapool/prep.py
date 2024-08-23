import re
import os
import gzip
import warnings
import pandas as pd

from glob import glob
from datetime import datetime
from string import ascii_letters, digits
from metapool.mp_strings import get_short_name_and_id
from metapool.plate import PlateReplication
from collections import Counter


REQUIRED_MF_COLUMNS = {'sample_name', 'barcode', 'primer', 'primer_plate',
                       'well_id_384', 'plating', 'extractionkit_lot',
                       'extraction_robot', 'tm1000_8_tool', 'primer_date',
                       'mastermix_lot', 'water_lot', 'processing_robot',
                       'tm300_8_tool', 'tm50_8_tool', 'sample_plate',
                       'project_name', 'orig_name', 'well_description',
                       'experiment_design_description', 'tm10_8_tool',
                       'library_construction_protocol', 'linker', 'platform',
                       'run_center', 'run_date', 'run_prefix', 'pcr_primers',
                       'sequencing_meth', 'target_gene', 'target_subfragment',
                       'center_name', 'center_project_name', 'well_id_96',
                       'instrument_model', 'runid'}


PREP_MF_COLUMNS = ['sample_name', 'barcode', 'center_name',
                   'center_project_name', 'experiment_design_description',
                   'instrument_model', 'lane', 'library_construction_protocol',
                   'platform', 'run_center', 'run_date', 'run_prefix', 'runid',
                   'sample_plate', 'sequencing_meth', 'linker', 'primer',
                   'primer_plate', 'well_id_384', 'plating',
                   'extractionkit_lot', 'extraction_robot', 'tm1000_8_tool',
                   'primer_date', 'mastermix_lot', 'water_lot',
                   'processing_robot', 'tm300_8_tool', 'tm50_8_tool',
                   'project_name', 'orig_name', 'well_description',
                   'pcr_primers', 'target_gene', 'tm10_8_tool',
                   'target_subfragment', 'well_id_96']

AMPLICON_PREP_COLUMN_RENAMER = {
    'Sample': 'sample_name',
    'Golay Barcode': 'barcode',
    '515FB Forward Primer (Parada)': 'primer',
    'Project Plate': 'project_plate',
    'Project Name': 'project_name',
    'Well': 'well',
    'Primer Plate #': 'primer_plate_number',
    'Plating': 'plating',
    'Extraction Kit Lot': 'extractionkit_lot',
    'Extraction Robot': 'extraction_robot',
    'TM1000 8 Tool': 'tm1000_8_tool',
    'Primer Date': 'primer_date',
    'MasterMix Lot': 'mastermix_lot',
    'Water Lot': 'water_lot',
    'Processing Robot': 'processing_robot',
    'sample sheet Sample_ID': 'well_description'
}

# put together by Gail, based on the instruments we know of
INSTRUMENT_LOOKUP = pd.DataFrame({
    'FS10001773': {'machine prefix': 'FS', 'Vocab': 'Illumina iSeq',
                   'Machine type': 'iSeq', 'run_center': 'KLM'},
    'A00953': {'machine prefix': 'A', 'Vocab': 'Illumina NovaSeq 6000',
               'Machine type': 'NovaSeq', 'run_center': 'IGM'},
    'A00169': {'machine prefix': 'A', 'Vocab': 'Illumina NovaSeq 6000',
               'Machine type': 'NovaSeq', 'run_center': 'LJI'},
    'M05314': {'machine prefix': 'M', 'Vocab': 'Illumina MiSeq',
               'Machine type': 'MiSeq', 'run_center': 'KLM'},
    'K00180': {'machine prefix': 'K', 'Vocab': 'Illumina HiSeq 4000',
               'Machine type': 'HiSeq', 'run_center': 'IGM'},
    'D00611': {'machine prefix': 'D', 'Vocab': 'Illumina HiSeq 2500',
               'Machine type': 'HiSeq/RR', 'run_center': 'IGM'},
    'LH00444': {'machine prefix': 'LH', 'Vocab': 'Illumina NovaSeq X',
                'Machine type': 'NovaSeq', 'run_center': 'IGM'},
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

    # Format should be YYMMDD_machinename_XXXX_FC
    # this regex has two groups, the first one is the date, and the second one
    # is the machine name + suffix. This URL shows some examples
    # tinyurl.com/rmy67kw
    matches6 = re.search(r'^(\d{6})_(\w*)', run_id)

    # iSeq uses YYYYMMDD and has a trailing -XXXX value, for example
    # 20220303_FS10001773_6_BRB11606-1914
    # So the regex is updated to allow 8 numbers for the date, and to handle
    # the trailing -XXXX piece
    matches8 = re.search(r'^(\d{8})_([\w-]*)', run_id)

    if matches6 is None and matches8 is None:
        raise ValueError('Unrecognized run identifier format "%s". The '
                         'expected format is either '
                         'YYMMDD_machinename_XXXX_FC or '
                         'YYYYMMDD_machinename_XXXX-XXXX' %
                         run_id)

    if matches6 is None:
        matches = matches8
        fmt = "%Y%m%d"
    else:
        matches = matches6
        fmt = "%y%m%d"

    if len(matches.groups()) != 2:
        raise ValueError('Unrecognized run identifier format "%s". The '
                         'expected format is either '
                         'YYMMDD_machinename_XXXX_FC or '
                         'YYYYMMDD_machinename_XXXX-XXXX' %
                         run_id)

    # convert illumina's format to qiita's format
    run_date = datetime.strptime(matches[1], fmt).strftime('%Y-%m-%d')

    return run_date, matches[2]


def is_nonempty_gz_file(name):
    """Taken from https://stackoverflow.com/a/37878550/379593"""
    with gzip.open(name, 'rb') as f:
        try:
            file_content = f.read(1)
            return len(file_content) > 0
        except Exception:
            return False


def remove_qiita_id(project_name):
    return get_short_name_and_id(project_name)[0]


def get_run_prefix(run_path, project, sample_id, lane):
    """For a sample find the run prefix

    Parameters
    ----------
    run_path: str
        Base path for the run
    project: str
        Name of the project
    sample_id: str
        Sample ID (was sample_name). Changed to reflect name used for files.
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

    qc = os.path.join(base, 'trimmed_sequences')
    hf = os.path.join(base, 'filtered_sequences')

    if _exists_and_has_files(qc) and _exists_and_has_files(hf):
        path = hf
    elif _exists_and_has_files(qc):
        path = qc
    elif _exists_and_has_files(hf):
        path = hf
    else:
        path = base

    search_me = '%s_S*_L*%s_R*.fastq.gz' % (sample_id, lane)

    results = glob(os.path.join(path, search_me))

    # at this stage there should only be two files forward and reverse
    if len(results) == 2:
        forward, reverse = sorted(results)
        if is_nonempty_gz_file(forward) and is_nonempty_gz_file(reverse):
            f, r = os.path.basename(forward), os.path.basename(reverse)
            if len(f) != len(r):
                raise ValueError("Forward and reverse sequences filenames "
                                 "don't match f:%s r:%s" % (f, r))

            # The first character that's different is the number in R1/R2. We
            # find this position this way because sometimes filenames are
            # written as _R1_.. or _R1.trimmed... and splitting on _R1 might
            # catch some substrings not part of R1/R2.
            for i in range(len(f)):
                if f[i] != r[i]:
                    i -= 2
                    break

            return f[:i]
        else:
            return None
    elif len(results) > 2:
        warnings.warn(('There are %d matches for sample "%s" in lane %s. Only'
                       ' two matches are allowed (forward and reverse): %s') %
                      (len(results),
                       sample_id,
                       lane,
                       ', '.join(sorted(results))))
    return None


def get_run_prefix_mf(run_path, project):
    search_path = os.path.join(run_path, project, 'amplicon',
                               '*_SMPL1_S*R?_*.fastq.gz')
    results = glob(search_path)

    # at this stage there should only be two files forward and reverse
    if len(results) == 2:
        forward, reverse = sorted(results)
        if is_nonempty_gz_file(forward) and is_nonempty_gz_file(reverse):
            f, r = os.path.basename(forward), os.path.basename(reverse)
            if len(f) != len(r):
                raise ValueError("Forward and reverse sequences filenames "
                                 "don't match f:%s r:%s" % (f, r))

            # The first character that's different is the number in R1/R2. We
            # find this position this way because sometimes filenames are
            # written as _R1_.. or _R1.trimmed... and splitting on _R1 might
            # catch some substrings not part of R1/R2.
            for i in range(len(f)):
                if f[i] != r[i]:
                    i -= 2
                    break

            return f[:i]
        else:
            return None
    elif len(results) > 2:
        warnings.warn("%d possible matches found for determining run-prefix. "
                      "Only two matches should be found (forward and "
                      "reverse): %s" % (len(results),
                                        ', '.join(sorted(results))))

    return None


def _file_list(path):
    return [f for f in os.listdir(path)
            if not os.path.isdir(os.path.join(path, f))]


def _exists_and_has_files(path):
    return os.path.exists(path) and len(_file_list(path))


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
        run_center = INSTRUMENT_LOOKUP.loc[instrument_model, 'run_center']
        instrument_model = INSTRUMENT_LOOKUP.loc[instrument_model, 'Vocab']
    else:
        instrument_prefix = get_machine_code(instrument_model)

        if instrument_prefix not in INSTRUMENT_LOOKUP['machine prefix']:
            ValueError('Unrecognized machine prefix %s' % instrument_prefix)

        instrument_model = INSTRUMENT_LOOKUP[
            INSTRUMENT_LOOKUP['machine prefix'] == instrument_prefix
            ]['Vocab'].unique()[0]

    return instrument_model, run_center


def agp_transform(frame, study_id):
    """If the prep belongs to the American Gut Project fill in some blanks

    Parameters
    ----------
    frame: pd.DataFrame
        The preparation file for a single project.
    study_id: str
        The Qiita study identifier for this preparations.

    Returns
    -------
    pd.DataFrame:
        If the study_id is "10317" then:
            - `center_name`, 'library_construction_protocol', and
              `experiment_design_description` columns are filled in with
              default values.
            - `sample_name` will be zero-filled no 9 digits.
        Otherwise no changes are made to `frame`.
    """
    if study_id == '10317':
        def zero_fill(name):
            if 'blank' not in name.lower() and name[0].isdigit():
                return name.zfill(9)
            return name

        frame['sample_name'] = frame['sample_name'].apply(zero_fill)
        frame['center_name'] = 'UCSDMI'
        frame['library_construction_protocol'] = 'Knight Lab KHP'
        frame['experiment_design_description'] = (
            'samples of skin, saliva and feces and other samples from the AGP')
    return frame


def _check_invalid_names(sample_names):
    # taken from qiita.qiita_db.metadata.util.get_invalid_sample_names
    valid = set(ascii_letters + digits + '.')

    def _has_invalid_chars(name):
        return bool(set(name) - valid)

    invalid = sample_names[sample_names.apply(_has_invalid_chars)]

    if len(invalid):
        warnings.warn('The following sample names have invalid '
                      'characters: %s' %
                      ', '.join(['"%s"' % i for i in invalid.values]))


def process_sample(sample, prep_columns, run_center, run_date, run_prefix,
                   project_name, instrument_model, run_id, lane):
    # initialize result
    result = {c: '' for c in prep_columns}

    # hard-coded columns
    result["platform"] = "Illumina"
    result["sequencing_meth"] = "sequencing by synthesis"
    result["center_name"] = "UCSD"

    # manually generated columns
    result["run_center"] = run_center
    result["run_date"] = run_date
    result["run_prefix"] = run_prefix
    result["center_project_name"] = project_name
    result["instrument_model"] = instrument_model
    result["runid"] = run_id

    # lane is extracted from sample-sheet but unlike the others is passed
    # to this function explicitly.
    result["lane"] = lane

    # handle multiple types of sample-sheets, where columns such
    # as 'syndna_pool_number' may or may not be present.
    additional_columns = ['syndna_pool_number', 'mass_syndna_input_ng',
                          'extracted_gdna_concentration_ng_ul',
                          'vol_extracted_elution_ul',
                          'total_rna_concentration_ng_ul',
                          "sample_name", "experiment_design_description",
                          "library_construction_protocol", "sample_plate",
                          "i7_index_id", "index", "i5_index_id", "index2",
                          "sample_project", 'well_id_384', 'sample_well'
                          ]

    for attribute in additional_columns:
        if attribute in sample:
            result[attribute] = sample[attribute]

    # sanity-checks
    if 'well_id_384' in sample and 'sample_well' in sample:
        # sample_well will be 'Sample_Well' to the user.
        raise ValueError("'well_id_384' and 'Sample_Well' are both defined"
                         "in sample.")

    if 'well_id_384' not in sample and 'sample_well' not in sample:
        raise ValueError("'well_id_384' and 'Sample_Well' are both undefined"
                         "in sample.")

    well_id = result['well_id_384'] if 'well_id_384' in sample else result[
        'sample_well']

    result["well_description"] = '%s.%s.%s' % (sample.sample_plate,
                                               sample.sample_name,
                                               well_id)

    # return unordered result. let caller re-order columns as they see fit.
    return result


def preparations_for_run(run_path, sheet, generated_prep_columns,
                         carried_prep_columns):
    """Given a run's path and sample sheet generates preparation files

    Parameters
    ----------
    run_path: str
        Path to the run folder
    sheet: dataFrame
        dataFrame() of SampleSheet contents
    generated_prep_columns: list
        List of required columns for output that are not expected in
        KLSampleSheet. It is expected that prep.py can derive these values
        appropriately.
    carried_prep_columns: list
        List of required columns for output that are expected in KLSampleSheet.
        Varies w/different versions of KLSampleSheet.

    Returns
    -------
    dict
        Dictionary keyed by run identifier, project name and lane. Values are
        preparations represented as DataFrames.
    """
    _, run_id = os.path.split(os.path.normpath(run_path))
    run_date, instrument_code = parse_illumina_run_id(run_id)
    instrument_model, run_center = get_model_and_center(instrument_code)

    output = {}

    # well_description is no longer a required column, since the sample-name
    #  is now taken from column sample_name, which is distinct from sample_id.
    # sample_id remains the version of sample_name acceptable to bcl2fastq/
    #  bcl-convert, while sample_name is now the original string. The
    #  well_description and description columns are no longer needed or
    #  required, but a warning will be generated if neither are present as
    #  they are customarily present.

    # if 'well_description' is defined instead as 'description', rename it.
    # well_description is a recommended column but is not required.
    if 'well_description' not in set(sheet.columns):
        warnings.warn("'well_description' is not present in sample-sheet. It "
                      "is not a required column but it is a recommended one.")
        if 'description' in sheet:
            warnings.warn("Using 'description' instead of 'well_description'"
                          " because that column isn't present", UserWarning)
            # copy and drop the original column
            sheet['well_description'] = sheet['description'].copy()
            sheet.drop('description', axis=1, inplace=True)

    not_present = set(carried_prep_columns) - set(sheet.columns)

    if not_present:
        raise ValueError("Required columns are missing: %s" %
                         ', '.join(not_present))

    all_columns = sorted(carried_prep_columns + generated_prep_columns)

    for project, project_sheet in sheet.groupby('sample_project'):
        project_name, qiita_id = get_short_name_and_id(project)

        # if the Qiita ID is not found then make for an easy find/replace
        if qiita_id == project:
            qiita_id = 'QIITA-ID'

        for lane, lane_sheet in project_sheet.groupby('lane'):
            # this is the portion of the loop that creates the prep
            data = []

            for well_id_col, sample in lane_sheet.iterrows():
                if isinstance(sample, pd.core.series.Series):
                    sample_id = well_id_col
                else:
                    sample_id = sample.sample_id

                run_prefix = get_run_prefix(run_path, project, sample_id, lane)

                # ignore the sample if there's no file
                if run_prefix is not None:
                    data.append(process_sample(sample, all_columns,
                                               run_center, run_date,
                                               run_prefix, project_name,
                                               instrument_model, run_id, lane))

            if not data:
                warnings.warn('Project %s and Lane %s have no data' %
                              (project, lane), UserWarning)

            # the American Gut Project is a special case. We'll likely continue
            # to grow this study with more and more runs. So we fill some of
            # the blanks if we can verify the study id corresponds to the AGP.
            # This was a request by Daniel McDonald and Gail
            prep = agp_transform(pd.DataFrame(columns=all_columns,
                                              data=data), qiita_id)

            _check_invalid_names(prep.sample_name)

            output[(run_id, project, lane)] = prep

    return output


def extract_run_date_from_run_id(run_id):
    # assume first segment of run_id will always be a valid date.
    tmp = run_id.split('_')[0]
    year = tmp[0:2]
    month = tmp[2:4]
    date = tmp[4:6]

    return ("20%s/%s/%s" % (year, month, date))


def preparations_for_run_mapping_file(run_path, mapping_file):
    """Given a run's path and mapping-file generates preparation files

        Parameters
        ----------
        run_path: str
            Path to the run folder
        mapping_file: pandas.DataFrame
            mapping-file to convert

        Returns
        -------
        dict
            Dictionary keyed by run identifier, project name and lane. Values
            are preparations represented as DataFrames.
        """

    # lowercase all columns
    mapping_file.columns = mapping_file.columns.str.lower()

    # add a faked value for 'lane' to preserve the original logic.
    # lane will always be '1' for amplicon runs.
    mapping_file['lane'] = pd.Series(
        ['1' for x in range(len(mapping_file.index))])

    # add a faked column for 'Sample_ID' to preserve the original logic.
    # count-related code will need to search run_directories based on
    # sample-id, not sample-name.
    mapping_file['Sample_ID'] = mapping_file.apply(
        # importing bcl_scrub_name led to a circular import
        lambda x: re.sub(r'[^0-9a-zA-Z\-\_]+', '_', x['sample_name']), axis=1)

    output = {}

    # lane is also technically a required column but since we generate it,
    # it's not included in REQUIRED_MF_COLUMNS.
    not_present = REQUIRED_MF_COLUMNS - set(mapping_file.columns)

    if not_present:
        raise ValueError("Required columns are missing: %s" %
                         ', '.join(not_present))

    for project, project_sheet in mapping_file.groupby('project_name'):
        _, qiita_id = get_short_name_and_id(project)

        # if the Qiita ID is not found then notify the user.
        if qiita_id == project:
            raise ValueError("Values in project_name must be appended with a"
                             " Qiita Study ID.")

        # note that run_prefix and run_id columns are required columns in
        # mapping-files. We expect these columns to be blank when seqpro is
        # run, however.
        run_prefix = get_run_prefix_mf(run_path, project)

        run_id = run_prefix.split('_SMPL1')[0]

        # return an Error if run_prefix could not be determined,
        # as it is vital for amplicon prep-info files. All projects will have
        # the same run_prefix.
        if run_prefix is None:
            raise ValueError("A run-prefix could not be determined.")

        for lane, lane_sheet in project_sheet.groupby('lane'):
            lane_sheet = lane_sheet.set_index('sample_name')

            # this is the portion of the loop that creates the prep
            data = []

            for sample_name, sample in lane_sheet.iterrows():
                row = {c: '' for c in PREP_MF_COLUMNS}

                row["sample_name"] = sample_name
                row["experiment_design_description"] = \
                    sample.experiment_design_description
                row["library_construction_protocol"] = \
                    sample.library_construction_protocol
                row["platform"] = sample.platform
                row["run_center"] = sample.run_center
                row["run_date"] = extract_run_date_from_run_id(run_id)
                row["run_prefix"] = run_prefix
                row["sequencing_meth"] = sample.sequencing_meth
                row["center_name"] = sample.center_name
                row["center_project_name"] = sample.center_project_name
                row["instrument_model"] = sample.instrument_model
                row["runid"] = run_id
                row["sample_plate"] = sample.sample_plate
                row["lane"] = lane
                row["barcode"] = sample.barcode
                row["linker"] = sample.linker
                row["primer"] = sample.primer
                row['primer_plate'] = sample.primer_plate
                if 'well_id_384' in sample:
                    row["well_id_384"] = sample.well_id_384
                elif 'sample_well' in sample:
                    row["sample_well"] = sample.sample_well
                row['well_id_96'] = sample.well_id_96
                row['plating'] = sample.plating
                row['extractionkit_lot'] = sample.extractionkit_lot
                row['extraction_robot'] = sample.extraction_robot
                row['tm1000_8_tool'] = sample.tm1000_8_tool
                row['primer_date'] = sample.primer_date
                row['mastermix_lot'] = sample.mastermix_lot
                row['water_lot'] = sample.water_lot
                row['processing_robot'] = sample.processing_robot
                row['tm300_8_tool'] = sample.tm300_8_tool
                row['tm50_8_tool'] = sample.tm50_8_tool
                row['tm10_8_tool'] = sample.tm10_8_tool
                row['project_name'] = sample.project_name
                row['orig_name'] = sample.orig_name
                row['well_description'] = sample.well_description
                row['pcr_primers'] = sample.pcr_primers
                row['target_gene'] = sample.target_gene
                row['target_subfragment'] = sample.target_subfragment
                data.append(row)

            if not data:
                warnings.warn('Project %s and Lane %s have no data' %
                              (project, lane), UserWarning)

            # the American Gut Project is a special case. We'll likely continue
            # to grow this study with more and more runs. So we fill some of
            # the blanks if we can verify the study id corresponds to the AGP.
            # This was a request by Daniel McDonald and Gail
            prep = agp_transform(pd.DataFrame(columns=PREP_MF_COLUMNS,
                                              data=data),
                                 qiita_id)

            _check_invalid_names(prep.sample_name)

            output[(run_id, project, lane)] = prep

    return output


def parse_prep(prep_path):
    """Parses a prep or pre-prep file into a DataFrame

    Parameters
    ----------
    prep_path: str or file-like
        Filepath to the preparation or pre-preparation file

    Returns
    -------
    pd.DataFrame
        Parsed file with the sample_name column set as the index.
    """
    prep = pd.read_csv(prep_path, sep='\t', dtype=str, keep_default_na=False,
                       na_values=[])

    prep.set_index('sample_name', verify_integrity=True, inplace=True)

    return prep


def generate_qiita_prep_file(platedf, seqtype):
    """Renames columns from prep sheet to common ones

    Parameters
    ----------
    platedf: pd.DataFrame
        dataframe that needs to be renamed and added to generate
        Qiita mapping file
    seqtype: str
        designates what type of amplicon sequencing

    Returns
    -------
    pd.DataFrame
        df formatted with all the Qiita prep info for the designated
        amplicon sequencing type
    """

    if seqtype == '16S':
        column_renamer = {
            'Sample': 'sample_name',
            'Golay Barcode': 'barcode',
            '515FB Forward Primer (Parada)': 'primer',
            'Project Name': 'project_name',
            'Project_Name': 'project_name',
            'Well': 'well_id_384',
            'Primer Plate #': 'primer_plate',
            'Plating': 'plating',
            'Extraction Kit Lot': 'extractionkit_lot',
            'Extraction Robot': 'extraction_robot',
            'TM1000 8 Tool': 'tm1000_8_tool',
            'TM300 8 Tool': 'tm300_8_tool',
            'TM50 8 Tool': 'tm50_8_tool',
            'TM10 8 Tool': 'tm10_8_tool',
            'Primer Date': 'primer_date',
            'MasterMix Lot': 'mastermix_lot',
            'Water Lot': 'water_lot',
            'Processing Robot': 'processing_robot',
            'Sample Plate': 'sample_plate',
            'Forward Primer Linker': 'linker',
            'Date': 'platemap_generation_date',
            'Project Abbreviation': 'project_abbreviation'
        }
    else:
        column_renamer = {
            'Sample': 'sample_name',
            'Golay Barcode': 'barcode',
            'Reverse complement of 3prime Illumina Adapter': 'primer',
            'Project Name': 'project_name',
            'Project_Name': 'project_name',
            'Well': 'well_id_384',
            'Primer Plate #': 'primer_plate',
            'Plating': 'plating',
            'Extraction Kit Lot': 'extractionkit_lot',
            'Extraction Robot': 'extraction_robot',
            'TM1000 8 Tool': 'tm1000_8_tool',
            'TM300 8 Tool': 'tm300_8_tool',
            'TM50 8 Tool': 'tm50_8_tool',
            'TM10 8 Tool': 'tm10_8_tool',
            'Primer Date': 'primer_date',
            'MasterMix Lot': 'mastermix_lot',
            'Water Lot': 'water_lot',
            'Processing Robot': 'processing_robot',
            'Sample Plate': 'sample_plate',
            'Reverse Primer Linker': 'linker',
            'Date': 'platemap_generation_date',
            'Project Abbreviation': 'project_abbreviation'
        }

    prep = platedf.copy()
    prep.rename(column_renamer, inplace=True, axis=1)

    primers_16S = 'FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT'
    primers_18S = 'FWD:GTACACACCGCCCGTC; REV:TGATCCTTCTGCAGGTTCACCTAC'
    primers_ITS = 'FWD:CTTGGTCATTTAGAGGAAGTAA; REV:GCTGCGTTCTTCATCGATGC'

    ptl_16S = 'Illumina EMP protocol 515fbc, 806r amplification of 16S rRNA V4'
    prtcl_18S = 'Illumina EMP 18S rRNA 1391f EukBr'
    prtcl_ITS = 'Illumina  EMP protocol amplification of ITS1fbc, ITS2r'

    prep['orig_name'] = prep['sample_name']
    prep['well_description'] = (prep['sample_plate'] + '.'
                                + prep['sample_name'] + '.'
                                + prep['well_id_384'])
    prep['center_name'] = 'UCSDMI'
    prep['run_center'] = 'UCSDMI'
    prep['platform'] = 'Illumina'
    prep['sequencing_meth'] = 'Sequencing by synthesis'

    if seqtype == '16S':
        prep['pcr_primers'] = primers_16S
        prep['target_subfragment'] = 'V4'
        prep['target_gene'] = '16S rRNA'
        prep['library_construction_protocol'] = ptl_16S
    elif seqtype == '18S':
        prep['pcr_primers'] = primers_18S
        prep['target_subfragment'] = 'V9'
        prep['target_gene'] = '18S rRNA'
        prep['library_construction_protocol'] = prtcl_18S
    elif seqtype == 'ITS':
        prep['pcr_primers'] = primers_ITS
        prep['target_subfragment'] = 'ITS_1_2'
        prep['target_gene'] = 'ITS'
        prep['library_construction_protocol'] = prtcl_ITS
    else:
        raise ValueError(f'Unrecognized value "{seqtype}" for seqtype')

    # Additional columns to add if not defined
    extra_cols = ['center_project_name', 'experiment_design_description',
                  'instrument_model', 'run_date', 'run_prefix', 'runid',
                  'tm10_8_tool', 'tm300_8_tool', 'tm50_8_tool']

    for c in extra_cols:
        if c not in prep.columns:
            prep[c] = ''

    if not prep.columns.is_unique:
        column_names = Counter(list(prep.columns))
        duplicates = [x for x in column_names if column_names[x] > 1]
        raise ValueError(
            "One or more columns in DataFrame resolve to the same"
            " column name(s): %s" % ", ".join(duplicates)
        )

    # the approved order of columns in the prep-file.
    column_order = ['sample_name', 'barcode', 'primer', 'primer_plate',
                    'well_id_384', 'plating', 'extractionkit_lot',
                    'extraction_robot', 'tm1000_8_tool', 'primer_date',
                    'mastermix_lot', 'water_lot', 'processing_robot',
                    'tm300_8_tool', 'tm50_8_tool', 'tm10_8_tool',
                    'sample_plate', 'project_name', 'orig_name',
                    'well_description', 'experiment_design_description',
                    'library_construction_protocol', 'linker', 'platform',
                    'run_center', 'run_date', 'run_prefix', 'pcr_primers',
                    'sequencing_meth', 'target_gene', 'target_subfragment',
                    'center_name', 'center_project_name', 'instrument_model',
                    'runid']

    # although some columns should be allowed to pass through, there are some
    # columns currently passed into metapool through the pre-prep file that
    # should definitely _not_ appear in the output. These include the
    # following:
    remove_these = {'Blank', 'Col', 'Compressed Plate Name', 'Plate Position',
                    'EMP Primer Plate Well', 'Forward Primer Pad', 'Name',
                    'Illumina 5prime Adapter', 'Original Name', 'Plate', 'Row',
                    'Primer For PCR', 'Project Plate', 'index',
                    'Plate elution volume', 'Plate map file', 'Time',
                    'RackID'}

    # only remove the columns from the above set that are actually present in
    # the prep dataframe to avoid possible errors.
    prep.drop(list(set(prep.columns) & remove_these), axis=1, inplace=True)

    # reorder the dataframe's columns according to the order in the list above.
    # Unrecognized columns in the input will become the right-most columns in
    # the output, sorted in alphabetical order.
    column_order += sorted(list(set(prep.columns) - set(column_order)))

    return prep[column_order]


def qiita_scrub_name(name):
    """Modifies a sample name to be Qiita compatible

    Parameters
    ----------
    name : str
        the sample name

    Returns
    -------
    str
        the sample name, formatted for qiita
    """
    return re.sub(r'[^0-9a-zA-Z\-\.]+', '.', name)


def pre_prep_needs_demuxing(pre_prep):
    """Returns True if pre-prep needs to be demultiplexed.

    Parameters
    ----------
    pre_prep: Dataframe representing pre-prep file.
        Object from where to extract the data.

    Returns
    -------
    bool
        True if pre-prep needs to be demultiplexed.
    """
    if 'contains_replicates' in pre_prep:
        contains_replicates = pre_prep.contains_replicates.apply(
            lambda x: x.lower() == 'true').unique()

        # By convention, all values in this column must either be True or
        # False.
        if len(contains_replicates) > 1:
            raise ValueError("all values in contains_replicates column must "
                             "either be True or False")

        # return either True or False, depending on the values found.
        return list(contains_replicates)[0]

    # legacy pre-prep does not handle replicates or no replicates were
    # found.
    return False


def demux_pre_prep(pre_prep):
    """Given a pre-prep file w/samples that are plate-replicates, generate new
       pre-prep files for each unique plate-replicate.

    Parameters
    ----------
    pre_prep: A Pandas dataframe representing a pre-prep file.
        Object from where to extract the data.

    Returns
    -------
    list of pre_preps
    """
    if not pre_prep_needs_demuxing(pre_prep):
        raise ValueError("pre_prep does not need to be demultiplexed")

    # use PlateReplication object to convert each sample's 384 well location
    # into a 96-well location + quadrant. Since replication is performed at
    # the plate-level, this will identify which replicants belong in which
    # new sample-sheet.
    plate = PlateReplication(None)

    pre_prep['quad'] = pre_prep.apply(lambda row:
                                      plate.get_96_well_location_and_quadrant(
                                          row.well_id_384)[0], axis=1)

    res = []

    for quad in sorted(pre_prep['quad'].unique()):
        # for each unique quadrant found, create a new dataframe that's a
        # subset containing only members of that quadrant. Delete the temporary
        # 'quad' column afterwards and reset the index to an integer value
        # starting at zero; the current-index will revert to a column named
        # 'sample_id'. Return the list of new dataframes.
        res.append(pre_prep[pre_prep['quad'] == quad].drop(['quad'], axis=1))

    return res
