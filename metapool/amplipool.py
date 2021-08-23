import os
import pandas as pd

from metapool.plate import _decompress_well, _plate_position


def assign_emp_index(plate_df, metadata):
    """Assign an EMP index to wells based on their compressed position

    Parameters
    ----------
    plate_df: pd.DataFrame
        Object with a Well column, and other metadata variables (usually with
        384 rows).
    metadata: pd.DataFrame
        Object with all the plate metadata (usually with 1-4 rows).

    Returns
    -------
    pd.DataFrame
        A table resulting from joining the compressed plate, the plate
        metadata, and the EMP indices.
    """
    # dataframe of wells organized by plate with a barcode per well
    emp_indices = _load_emp_indices()

    # the Well column exists elsewhere already so we rename it to avoid
    # duplication in the merged table
    emp_indices.rename({'Well': 'EMP Primer Plate Well'}, axis=1, inplace=True)

    # helpful to perform joins
    plate_df['__decompressed_well__'] = plate_df.Well.apply(_decompress_well)
    plate_df['__plate_position__'] = plate_df.Well.apply(_plate_position)

    # merge the compressed plate with the metadata based on the plate position
    plate_df = plate_df.merge(metadata, how='left',
                              left_on='__plate_position__',
                              right_on='Plate Position')

    # merge tables based on the primer plate and the
    # 96-well-plate-ID
    plate_df = plate_df.merge(
        emp_indices, how='left',
        left_on=['Primer Plate #', '__decompressed_well__'],
        right_on=['Plate', 'EMP Primer Plate Well'])

    # remove all the helper columns
    plate_df.drop(['__plate_position__', '__decompressed_well__'], axis=1,
                  inplace=True)

    return plate_df


def _load_emp_indices():
    """Helper method to load EMP primer plates"""
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data',
                      'emp-16S-V4-515F-806R-parada-april.tsv')
    indices = pd.read_csv(fn, sep='\t', dtype=str, keep_default_na=False,
                          na_values=[])
    return indices
