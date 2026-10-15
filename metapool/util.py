import os
import pandas as pd
import warnings
from metapool.mp_strings import get_qiita_id_from_project_name, \
    SAMPLE_NAME_KEY, QIITA_ID_KEY, PM_PROJECT_NAME_KEY, PM_PROJECT_ABBREV_KEY

QIITA_STUDY_ID_KEY = 'qiita_study_id'
SET_SUFFIX = '_set'


def _missing_like(series, index):
    """An all-missing column to stand in for series in a frame lacking it.

    Text keeps its dtype (str under pandas 3), so a later concat with series
    stays text; anything else gets an object column of missing values, which
    holds any value.
    """
    if pd.api.types.is_string_dtype(series):
        return pd.Series(None, index=index, dtype=series.dtype)
    return pd.Series(None, index=index, dtype="object")


def join_dfs_from_files(input_fps, req_cols_to_extract,
                        opt_cols_to_extract=None, unique_cols=None,
                        dtype=str, sep="\t"):
    if opt_cols_to_extract is None:
        opt_cols_to_extract = []

    if unique_cols is None:
        unique_cols = req_cols_to_extract

    uniques_required = [x in req_cols_to_extract for x in unique_cols]
    if not all(uniques_required):
        raise ValueError("All unique_cols must be in req_cols_to_extract")

    growing_df = None
    for curr_fp in input_fps:
        if not os.path.isfile(curr_fp):
            raise ValueError(
                "Problem! %s is not a path to a valid file" % curr_fp)

        # load the file into a dataframe, then extract the columns of interest
        curr_df = pd.read_csv(curr_fp, sep=sep, dtype=dtype)
        missing_req_cols = [x for x in req_cols_to_extract if
                            x not in curr_df.columns]
        if len(missing_req_cols) > 0:
            raise ValueError(f"File '{curr_fp}' is missing required columns "
                             f"{missing_req_cols}")

        curr_opt_cols = \
            [x for x in opt_cols_to_extract if x in curr_df.columns]
        curr_cols_to_extract = req_cols_to_extract + curr_opt_cols
        curr_extracted_df = curr_df[curr_cols_to_extract]

        if growing_df is None:
            growing_df = curr_extracted_df
        else:
            not_in_old = \
                set(curr_extracted_df.columns) - set(growing_df.columns)
            not_in_new = \
                set(growing_df.columns) - set(curr_extracted_df.columns)
            for curr_missing in not_in_old:
                growing_df[curr_missing] = _missing_like(
                    curr_extracted_df[curr_missing], growing_df.index)
            for curr_missing in not_in_new:
                curr_extracted_df[curr_missing] = _missing_like(
                    growing_df[curr_missing], curr_extracted_df.index)
            growing_df = pd.concat([growing_df, curr_extracted_df],
                                   ignore_index=True)
        # end if this isn't the first file
    # next file path

    # if any value in the unique_col is duplicated, raise an error
    errs = []
    for curr_unique_col in unique_cols:
        if growing_df[curr_unique_col].duplicated().any():
            dupes = growing_df[curr_unique_col][
                growing_df[curr_unique_col].duplicated()].unique().tolist()
            errs.append(f"Duplicate '{curr_unique_col}' values found in "
                        f"files: {dupes}")
    # next unique_col
    if len(errs) > 0:
        raise ValueError('\n'.join(errs))

    return growing_df


# sample_accession_df should be combined across all studies and should
# have a 'sample_name' col.
# metadata_df should be combined across all studies and have
# 'sample_name' and 'qiita_study_id' cols.
# each entry in studies_info dict should have a key named 'Project Name'
# and one named 'Project Abbreviation', as shown in example below:
# studies_info = [
#     {
#     'Project Name': 'Celeste_Adaptation_12986', # PROJECTNAME_QIITAID
#     'Project Abbreviation': 'ADAPT', # PROJECTNAME
#     'sample_accession_fp': './test_data/Plate_Maps/sa_file_1.tsv',
#     'qiita_metadata_fp': './test_data/Plate_Maps/12986_20230314-090655.txt',
#     'experiment_design_description': 'isolate sequencing',
#     'HumanFiltering': 'False',
#     'Email': 'r@gmail.com'
#     },
#     <etc, etc>
# ]
def extend_sample_accession_df(sample_accession_df, studies_info, metadata_df):
    local_metadata_df = metadata_df.copy()
    local_metadata_df[QIITA_ID_KEY] = local_metadata_df[QIITA_STUDY_ID_KEY]

    # extract qiita study ids from the Project Name in studies_info entries
    local_studies_info = studies_info.copy()
    for curr_study in local_studies_info:
        curr_project_name = curr_study[PM_PROJECT_NAME_KEY]
        curr_qiita_id = get_qiita_id_from_project_name(curr_project_name)
        curr_study[QIITA_ID_KEY] = curr_qiita_id
    studies_df = pd.DataFrame(local_studies_info)

    # check for qiita_study_ids in studies_df that aren't in the metadata
    _check_for_missing_df_ids(studies_df, local_metadata_df, QIITA_ID_KEY,
                              'studies', 'metadata')
    # merge the metadata df with the studies_df on the qiita_study_id
    metadata_plus_df = pd.merge(local_metadata_df, studies_df, on=QIITA_ID_KEY)
    # pull the qiita study id off the sample name
    metadata_plus_df[SAMPLE_NAME_KEY] = metadata_plus_df.apply(
        lambda x: x[SAMPLE_NAME_KEY].replace(f"{x[QIITA_ID_KEY]}.", ""),
        axis=1
    )

    # now add the project name and abbreviation to the sample_accession_df
    # by merging on the qiita_id.sample_name
    extension_cols = \
        [SAMPLE_NAME_KEY, PM_PROJECT_NAME_KEY, PM_PROJECT_ABBREV_KEY]
    _check_for_missing_df_ids(
        sample_accession_df, metadata_plus_df, SAMPLE_NAME_KEY,
        'sample accession', 'metadata')
    sample_accession_plus_df = pd.merge(
        sample_accession_df, metadata_plus_df[extension_cols],
        on=SAMPLE_NAME_KEY)
    return sample_accession_plus_df


# add the project abbreviation to the compression layout so the user doesn't
# have to enter it in two places
def extend_compression_layout_info(compression_layout, studies_info):
    extended_compression_layout = compression_layout.copy()
    # for each dict in compression_layout
    for curr_plate in extended_compression_layout:
        # get its Project Name
        curr_project_name = curr_plate[PM_PROJECT_NAME_KEY]
        # find that Project Name in studies_info
        found_study = None
        for curr_study in studies_info:
            if curr_study[PM_PROJECT_NAME_KEY] == curr_project_name:
                found_study = curr_study
                break
        # next study

        # add the Project Abbreviation from the studies_info to the
        # compression layout dict
        if found_study is None:
            raise ValueError(f"{PM_PROJECT_NAME_KEY} '{curr_project_name}' "
                             f"in the compression layout is not in the "
                             f"studies info")
        curr_plate[PM_PROJECT_ABBREV_KEY] = found_study[PM_PROJECT_ABBREV_KEY]
    # next plate

    return extended_compression_layout


def _check_for_missing_df_ids(subset_df, superset_df, id_col,
                              sub_name, sup_name):
    # if there are any qiita_study_ids in the studies_df that aren't in the
    # metadata_df, raise an error
    found_mask = subset_df[id_col].isin(superset_df[id_col])
    if not found_mask.all():
        missing_vals = subset_df.loc[~found_mask, id_col].unique()
        raise ValueError(f'Some {id_col} values in the {sub_name} dataframe '
                         f'are not in the {sup_name} dataframe: '
                         f'{", ".join(missing_vals)}')


def get_set_fp(set_fbase, set_id, extension="txt"):
    return f"{set_fbase}{SET_SUFFIX}_{set_id}.{extension}"


def warn_if_fp_exists(fp):
    if os.path.isfile(fp):
        warnings.warn(f"Warning! This file exists already: {fp}.")


def drop_unnamed_nan_columns(df):
    """
    Drops columns whose names start with 'Unnamed:' and contain only NaN values
    """
    return df.drop(
        columns=[
            col for col in df.columns
            if col.startswith('Unnamed:') and df[col].isna().all()
        ]
    )


def convert_to_bool(x):
    """
    Converts a value to a boolean if true or false, otherwise returns None.

    Parameters
    ----------
    x: object
        The value to convert to a boolean.

    Returns
    -------
    bool or None
        For booleans or strings that a human would read as true or false
        (barring case, whitespace, etc.), returns True or False. For any other
        type or value, returns None.  Note that it does NOT convert
        integers or floats to booleans.

    """
    if type(x) is bool:
        # column type is already correct.
        return x
    elif type(x) is str:
        # strings should be converted to bool if possible.
        if x.strip().lower() == 'true':
            return True
        elif x.strip().lower() == 'false':
            return False

    return None


def set_unique_index(df, keys):
    """Return df indexed by keys; raise if that index has duplicates.

    Stands in for set_index(..., verify_integrity=True), which pandas 3
    deprecates.

    Parameters
    ----------
    df: pd.DataFrame
        The dataframe to index.
    keys: str or list of str
        The column(s) to use as the index.

    Returns
    -------
    pd.DataFrame
        A new dataframe indexed by keys.

    Raises
    ------
    ValueError
        If the index built from keys contains duplicate values.
    """
    result = df.set_index(keys)
    if not result.index.is_unique:
        dups = result.index[result.index.duplicated()].unique().tolist()
        raise ValueError(f"Index has duplicate keys: {dups}")
    return result
