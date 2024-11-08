from types import MappingProxyType
from metapool.mp_strings import get_qiita_id_from_project_name, \
    SAMPLE_NAME_KEY, SAMPLE_TYPE_KEY, PRIMARY_STUDY_KEY, \
    SECONDARY_STUDIES_KEY, PM_PROJECT_NAME_KEY, PM_PROJECT_PLATE_KEY, \
    PM_SAMPLE_KEY, QIITA_ID_KEY

_BLANK_ROOT = "BLANK"
_BLANK_SAMPLE_TYPE = "control blank"
_STUDY_ID_DELIMITER = ";"

SAMPLE_CONTEXT_COLS = MappingProxyType({
    SAMPLE_NAME_KEY: str,
    SAMPLE_TYPE_KEY: str,
    PRIMARY_STUDY_KEY: str,
    SECONDARY_STUDIES_KEY: str})


def is_blank(sample_name, sample_context=None):
    if sample_context is not None:
        sample_record_mask = \
            sample_context[SAMPLE_NAME_KEY] == sample_name

        # if the sample name IS in the context, check further
        if sample_record_mask.any():
            found_samples = sample_record_mask.value_counts()[True]
            if found_samples != 1:
                raise ValueError(
                    f"Expected exactly one record with sample name "
                    f"'{sample_name}', but found {found_samples} records")

            is_sample_blank = \
                (sample_context.loc[
                    sample_record_mask,
                    SAMPLE_TYPE_KEY]).iloc[0] == _BLANK_SAMPLE_TYPE
        else:
            # if the sample name is not in the context, it is not any kind of
            # control, let alone a blank
            is_sample_blank = False
    # if there is no context, we can only check by name
    else:
        is_sample_blank = _is_blank_by_name(sample_name)
    return is_sample_blank


def get_blank_root():
    return _BLANK_ROOT


# Generate the metadata dictionary SampleContext section, which looks for
# example like the below, and add it to the metadata dictionary
#
#    'SampleContext': [
#        {
#            'Sample_Name': 'BLANK.NPH.4.G11',
#            'PrimaryQiitaStudy': '12986',
#            'SecondaryQiitaStudies': '12000;10981',
#            'Sample_Type': 'control blank'
#        },
#        {
#            'Sample_Name': 'BLANK.NPH.8.G10',
#            'PrimaryQiitaStudy': '12986',
#            'SecondaryQiitaStudies': '',
#            'Sample_Type': 'control blank'
#        },
#        {
#            'Sample_Name': 'BLANK.NPH.18.G01',
#            'PrimaryQiitaStudy': '12986',
#            'SecondaryQiitaStudies': '10981',
#            'Sample_Type': 'control blank'
#        },
#    ]
def get_delimited_controls_details_from_compressed_plate(
        the_plate_df, blanks_mask=None):
    if blanks_mask is None:
        blanks_mask = the_plate_df[PM_SAMPLE_KEY].apply(is_blank)

    # get all the unique Project Plate values for the blanks
    names_of_plates_w_blanks = \
        the_plate_df.loc[blanks_mask, PM_PROJECT_PLATE_KEY].unique()

    # start building the SampleContext object
    blanks_context_df = the_plate_df.loc[
        blanks_mask,
        [PM_SAMPLE_KEY, PM_PROJECT_NAME_KEY, PM_PROJECT_PLATE_KEY]].copy()
    blanks_context_df.rename(columns={PM_SAMPLE_KEY: SAMPLE_NAME_KEY},
                             inplace=True)
    blanks_context_df[PRIMARY_STUDY_KEY] = \
        blanks_context_df[PM_PROJECT_NAME_KEY].apply(
            get_qiita_id_from_project_name)
    blanks_context_df[SAMPLE_TYPE_KEY] = _BLANK_SAMPLE_TYPE

    # for each plate in names_of_plates_w_blanks
    for name_of_curr_plate_w_blanks in names_of_plates_w_blanks:
        # note that this mask is for ALL of plate_df, not just blank rows
        curr_plate_mask = \
            the_plate_df[PM_PROJECT_PLATE_KEY] == name_of_curr_plate_w_blanks
        # get all the unique Project Name values for plate_df rows that have
        # Project Plate value == name_of_curr_plate_w_blanks
        projects_on_curr_plate = the_plate_df.loc[
            curr_plate_mask, PM_PROJECT_NAME_KEY].unique().tolist()
        # add projects_on_curr_plate to blanks_context_df for
        # this plate's blanks
        curr_plate_blanks_mask = \
            blanks_context_df[
                PM_PROJECT_PLATE_KEY] == name_of_curr_plate_w_blanks

        # for each row in blanks_context_df that has the current plate name,
        # remove that row's project name from projects_on_curr_plate list to
        # generate the list of secondary studies for that blank
        curr_plate_records_df = blanks_context_df[curr_plate_blanks_mask]
        curr_secondary_studies_strs = curr_plate_records_df.apply(
                _get_secondary_studies,
                all_projects_on_plate=projects_on_curr_plate,
                axis=1)
        blanks_context_df.loc[
            curr_plate_blanks_mask, SECONDARY_STUDIES_KEY] = \
            curr_secondary_studies_strs
    # next plate name of plate w blanks

    # clean up by removing the Project Plate and Project Name columns
    blanks_context_df.drop(columns=[PM_PROJECT_PLATE_KEY, PM_PROJECT_NAME_KEY],
                           inplace=True)

    # turn blanks_context_df into a list of dictionaries
    blanks_context_list = blanks_context_df.to_dict(orient="records")
    return blanks_context_list


def get_all_projects_in_context(sample_context):
    if sample_context is None or sample_context.empty:
        return []

    # get unique set of values in the PRIMARY_QIITA_STUDY_KEY col
    study_ids = sample_context[PRIMARY_STUDY_KEY].unique().tolist()

    # get unique set of values in the SECONDARY_QIITA_STUDY_KEY col
    secondaries = sample_context[SECONDARY_STUDIES_KEY].unique().tolist()
    # these are delimited strings, so we need to split them into lists,
    # append them, and then deduplicate
    for s in secondaries:
        study_ids.extend(_split_secondary_studies(s))
    study_ids = sorted(list(set(study_ids)))
    return study_ids


def get_controls_details_from_context(sample_context):
    if sample_context is None:
        return {}

    result = {}
    # convert the sample context dataframe to a list of dictionaries
    # and split the secondary studies strings into lists
    controls_dicts_list = sample_context.to_dict(orient="records")
    for curr_dict in controls_dicts_list:
        curr_sample_name = curr_dict[SAMPLE_NAME_KEY]
        curr_dict[SECONDARY_STUDIES_KEY] = sorted(_split_secondary_studies(
            curr_dict[SECONDARY_STUDIES_KEY]))
        result[curr_sample_name] = curr_dict
    return result


def denormalize_controls_details(controls_details):
    denormalized_controls_details = []
    if controls_details is None:
        return None

    def _denormalize_details_record(a_record, qiita_id):
        # remove the primary and secondary studies entries from the record,
        # replace with a single QIITA_ID_KEY entry; leave everything else alone
        denormalized_record = a_record.copy()
        denormalized_record.pop(SECONDARY_STUDIES_KEY)
        denormalized_record.pop(PRIMARY_STUDY_KEY)
        denormalized_record[QIITA_ID_KEY] = qiita_id
        return denormalized_record

    for curr_record in controls_details.values():
        for curr_secondary_study in curr_record[SECONDARY_STUDIES_KEY]:
            curr_denorm_secondary_record = _denormalize_details_record(
                curr_record, curr_secondary_study)
            denormalized_controls_details.append(curr_denorm_secondary_record)
        # next secondary study, if any

        curr_denorm_primary_record = _denormalize_details_record(
            curr_record, curr_record[PRIMARY_STUDY_KEY])
        denormalized_controls_details.append(curr_denorm_primary_record)
    # next controls_detail record

    denormalized_controls_details = sorted(
        denormalized_controls_details,
        key=lambda k: (k[SAMPLE_NAME_KEY], k[QIITA_ID_KEY]))
    return denormalized_controls_details


def make_manual_control_details(sample_name, primary_study,
                                secondary_studies=None, sample_type=None):

    if sample_type is None:
        sample_type = _BLANK_SAMPLE_TYPE

    if secondary_studies is None:
        secondary_studies = []

    details_dict = {
        SAMPLE_NAME_KEY: sample_name,
        SAMPLE_TYPE_KEY: sample_type,
        PRIMARY_STUDY_KEY: primary_study,
        SECONDARY_STUDIES_KEY: secondary_studies
    }
    return details_dict


def _is_blank_by_name(sample_name):
    return sample_name.startswith(_BLANK_ROOT)


def _split_secondary_studies(secondary_studies_str):
    secondaries = []
    if secondary_studies_str:
        secondaries = secondary_studies_str.split(_STUDY_ID_DELIMITER)
    return secondaries


def _get_secondary_studies(curr_row, all_projects_on_plate):
    row_project = curr_row[PM_PROJECT_NAME_KEY]

    # NB: below will error if row_project is not in
    # projects_on_curr_plate; this is as it should be because
    # that situation should never arise and if it does, we need to
    # stop and figure out why
    secondary_projects = all_projects_on_plate.copy()
    secondary_projects.remove(row_project)
    secondary_projects = sorted(secondary_projects)

    # now split each secondary_projects on _ and make a list of the
    # last element of each split
    secondary_qiita_ids = [get_qiita_id_from_project_name(x)
                           for x in secondary_projects]
    secondary_qiita_ids = sorted(secondary_qiita_ids)

    secondary_qiita_ids_str = _STUDY_ID_DELIMITER.join(secondary_qiita_ids)
    return secondary_qiita_ids_str
