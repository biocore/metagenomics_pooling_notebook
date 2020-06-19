import re
from datetime import datetime


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
        Instrument model
    """

    # of the form
    # YYMMDD_machinename_XXXX_FC
    run_id = run_id.split('_')

    # convert illumina's format to qiita's format
    run_date = datetime.strptime(run_id[0], '%y%m%d').strftime('%Y-%m-%d')

    return run_date, run_id[1]


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
    pass


def remove_qiita_id(project_name):

    # project identifiers are digit groups at the end of the project name
    # preceded by an underscore CaporasoIllumina_550
    qiita_id_re = re.compile(r'(.+)_(\d+)$')

    # no matches
    matches = re.search(qiita_id_re, project_name)
    if matches is None:
        return project_name
    else:
        return matches[0]


def preparations_for_run(run_dir, sample_sheet):
    pass
