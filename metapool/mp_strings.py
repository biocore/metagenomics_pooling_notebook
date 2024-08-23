# I am aware that "having a global constants class" is largely regarded as
# an antipattern in Python because it increases dependence between the
# modules using it and keeps pieces tightly coupled that should be loosely
# coupled.  However, the current codebase has some literal strings repeated
# *throughout* the code, and I need them all in one place while we get them
# replaced with constants so that we can (a) prevent re-definition of the same
# constant string in multiple places and (b) get a clear view of which of
# them are really cross-domain and which really do belong to a specific
# domain (where they can then be moved in the future).  Constants that I
# already *know* to be domain-specific are not included here.

# Note:
# A number of these string literals are also used in amplicon-related code
# (e.g., prep.py, scripts/seqpro.py, seqpro_mf.py, etc.).  However, I don't
# know if they are the same because they come from the same source, or if it is
# just coincidental.  For now, I have NOT changed the amplicon-related
# instances to use these constants.

import re

# These are the "standard"-format snake-case column names;
# if writing something NEW, use these
SAMPLE_NAME_KEY = "sample_name"
QIITA_ID_KEY = "qiita_id"
SAMPLE_PROJECT_KEY = "sample_project"
ORIG_NAME_KEY = "orig_name"
EXPT_DESIGN_DESC_KEY = 'experiment_design_description'
CONTAINS_REPLICATES_KEY = "contains_replicates"
PROJECT_SHORT_NAME_KEY = "short_project_name"
PROJECT_FULL_NAME_KEY = "full_project_name"
# NB: don't use this as the name for a column of sample identifiers.  It is
# the key for a list of sample details *objects* (e.g., dictionaries)
SAMPLES_DETAILS_KEY = "samples"

# Plate map (PM) column names
PM_SAMPLE_KEY = "Sample"
PM_PROJECT_PLATE_KEY = "Project Plate"
PM_PROJECT_NAME_KEY = "Project Name"
PM_PROJECT_ABBREV_KEY = 'Project Abbreviation'
PM_COMPRESSED_PLATE_NAME_KEY = "Compressed Plate Name"
PM_BLANK_KEY = "Blank"
PLATE_NAME_DELIMITER = "_"

# Currently used to describe sample context (only)
SAMPLE_TYPE_KEY = "sample_type"
PRIMARY_STUDY_KEY = "primary_qiita_study"
SECONDARY_STUDIES_KEY = "secondary_qiita_studies"


def parse_project_name(project_name):
    """
    Split fully-qualified project_name into a short project_name and a qiita-id
    :param project_name: A fully-qualified project name e.g: Feist_1161.
    :return: Dictionary of qiita_id, short_project_name, and full_project_name
    """
    if project_name is None or project_name == '':
        raise ValueError("project_name cannot be None or empty string")

    # project identifiers are digit groups at the end of the project name
    # preceded by an underscore CaporasoIllumina_550
    matches = re.search(r'^(.+)_(\d+)$', str(project_name))

    if matches is None:
        raise ValueError(f"'{project_name}' does not contain a Qiita-ID.")

    proj_info_dict = {
        QIITA_ID_KEY: matches[2],
        PROJECT_SHORT_NAME_KEY: matches[1],
        PROJECT_FULL_NAME_KEY: project_name
    }

    return proj_info_dict


def get_short_name_and_id(project_name):
    try:
        proj_info_dict = parse_project_name(project_name)
    except ValueError:
        return project_name, None

    return proj_info_dict[PROJECT_SHORT_NAME_KEY], proj_info_dict[QIITA_ID_KEY]


def get_qiita_id_from_project_name(project_name):
    return parse_project_name(project_name)[QIITA_ID_KEY]


def get_plate_num_from_plate_name(plate_name):
    return _split_plate_name(plate_name)[1]


def get_main_project_from_plate_name(plate_name):
    return _split_plate_name(plate_name)[0]


# Full names for uncompressed (i.e., 96-well) plates have the format
# <project_name>_Plate_<plate_number>, such as
# Celeste_Adaptation_12986_Plate_16.
def _split_plate_name(plate_name):
    # split on just the *last* delimiter.
    # Note that plates can have more than one project on them, and the
    # name embedded in the plate name is that of (only) the "main" project.
    if PLATE_NAME_DELIMITER not in plate_name:
        raise ValueError(f"Plate name '{plate_name}' is malformed.")
    main_project_info, plate_num = \
        plate_name.rsplit(PLATE_NAME_DELIMITER, maxsplit=1)

    # There may or may not be a "_Plate" in the name; remove if there
    main_project_info = main_project_info.replace(
        f"{PLATE_NAME_DELIMITER}Plate", "")
    return main_project_info, plate_num
