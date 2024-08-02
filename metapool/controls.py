from types import MappingProxyType

BLANK_ROOT = "BLANK"
BLANK_SAMPLE_TYPE = "control blank"
STUDY_ID_DELIMITER = ";"
QIITA_SAMPLE_NAME_KEY = "sample_name"
SAMPLE_TYPE_KEY = "sample_type"
PRIMARY_STUDY_KEY = "primary_qiita_study"
SECONDARY_STUDIES_KEY = "secondary_qiita_studies"

SAMPLE_CONTEXT_COLS = MappingProxyType({
    QIITA_SAMPLE_NAME_KEY: str,
    SAMPLE_TYPE_KEY: str,
    PRIMARY_STUDY_KEY: str,
    SECONDARY_STUDIES_KEY: str})


def is_blank(sample_name, sample_context=None):
    if sample_context is not None:
        sample_record_mask = \
            sample_context[QIITA_SAMPLE_NAME_KEY] == sample_name

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
                    SAMPLE_TYPE_KEY]).iloc[0] == BLANK_SAMPLE_TYPE
        else:
            # if the sample name is not in the context, it is not any kind of
            # control, let alone a blank
            is_sample_blank = False
    # if there is no context, we can only check by name
    else:
        is_sample_blank = _is_blank_by_name(sample_name)
    return is_sample_blank


def get_all_projects_in_context(sample_context):
    if sample_context is None:
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
    controls_dicts_list = sample_context.to_dict(orient="records")
    for curr_dict in controls_dicts_list:
        curr_sample_name = curr_dict[QIITA_SAMPLE_NAME_KEY]
        curr_dict[SECONDARY_STUDIES_KEY] = _split_secondary_studies(
            curr_dict[SECONDARY_STUDIES_KEY])
        result[curr_sample_name] = curr_dict
    return result


def make_manual_control_details(sample_name, primary_study,
                                secondary_studies=None, sample_type=None):

    if sample_type is None:
        sample_type = BLANK_SAMPLE_TYPE

    if secondary_studies is None:
        secondary_studies = []

    details_dict = {
        QIITA_SAMPLE_NAME_KEY: sample_name,
        SAMPLE_TYPE_KEY: sample_type,
        PRIMARY_STUDY_KEY: primary_study,
        SECONDARY_STUDIES_KEY: secondary_studies
    }
    return details_dict


# def copy_controls_between_projects(controls_info):
#     for curr_control in controls_info:
#         curr_sample_name = curr_control[QIITA_SAMPLE_NAME_KEY]
#         curr_primary_qiita_study = curr_control[PRIMARY_STUDY_KEY]
#         curr_secondary_qiita_studies = curr_control[SECONDARY_STUDIES_KEY]
#         for curr_secondary_qiita_study in curr_secondary_qiita_studies:
#             self.copy_sequences(curr_sample_name, curr_primary_qiita_study,
#                                 curr_secondary_qiita_study)
#         # next secondary study
#     # next control


def _get_blanks_mask_by_name_format(a_df):
    return a_df[QIITA_SAMPLE_NAME_KEY].str.startswith(BLANK_ROOT)


def _is_blank_by_name(sample_name):
    return sample_name.startswith(BLANK_ROOT)


def _split_secondary_studies(secondary_studies_str):
    secondaries = []
    if secondary_studies_str:
        secondaries = secondary_studies_str.split(STUDY_ID_DELIMITER)
    return secondaries
