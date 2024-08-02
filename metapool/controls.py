from types import MappingProxyType

BLANK_ROOT = "BLANK"
BLANK_SAMPLE_TYPE = "control blank"
STUDY_ID_DELIMITER = ";"
SAMPLE_NAME_KEY = "sample_name"
SAMPLE_TYPE_KEY = "sample_type"
PRIMARY_STUDY_KEY = "primary_qiita_study"
SECONDARY_STUDIES_KEY = "secondary_qiita_studies"

SAMPLE_CONTEXT_COLS = MappingProxyType({
    SAMPLE_NAME_KEY: str,
    SAMPLE_TYPE_KEY: str,
    PRIMARY_STUDY_KEY: str,
    SECONDARY_STUDIES_KEY: str})


def is_blank(sample_name, sample_context=None):
    if sample_context is not None:
        sample_record_mask = sample_context[SAMPLE_NAME_KEY] == sample_name
        if not sample_record_mask.any():
            found_samples = 0
        else:
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


def _get_blanks_mask_by_name_format(a_df):
    return a_df[SAMPLE_NAME_KEY].str.startswith(BLANK_ROOT)


def _is_blank_by_name(sample_name):
    return sample_name.startswith(BLANK_ROOT)


def _split_secondary_studies(secondary_studies_str):
    secondaries = []
    if secondary_studies_str:
        secondaries = secondary_studies_str.split(STUDY_ID_DELIMITER)
    return secondaries
