import os
import pandas as pd

from metapool.plate import _decompress_well, _plate_position


def assign_emp_index(plate_df, metadata):
    """
    """
    plate_df = plate_df.copy()

    # dataframe of wells organized by plate with a barcode per well
    emp_indices = _load_emp_indices()
    emp_indices['__plate_well__'] = emp_indices['Plate'] + emp_indices['Well']

    # the Well column exists elsewhere already so we rename it to avoid
    # duplication in the merged table
    emp_indices['EMP Primer Plate Well'] = emp_indices['Well'].copy()
    emp_indices.drop(['Well'], axis=1, inplace=True)

    # helpful to perform joins
    plate_df['__decompressed_well__'] = plate_df.Well.apply(_decompress_well)
    plate_df['__plate_position__'] = plate_df.Well.apply(_plate_position)

    # merge the compressed plate with the metadata based on the plate position
    plate_df = plate_df.merge(metadata, how='left',
                              left_on='__plate_position__',
                              right_on='Plate Position')

    # create a new column concatenating the primer plate and the
    # 96-well-plate-ID. These are the two identifiers that we need to join the
    # EMP indices with the plated samples
    plate_df['__plate_well__'] = (plate_df['Primer Plate #'] +
                                  plate_df['__decompressed_well__'])
    plate_df = plate_df.merge(emp_indices, how='left', on='__plate_well__')

    # remove all the helper columns
    plate_df.drop(['__plate_position__', '__plate_well__',
                   '__decompressed_well__'], axis=1, inplace=True)

    return plate_df


def _load_emp_indices():
    """Helper method to load EMP primer plates"""
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data',
                      'emp-16S-V4-515F-806R-parada-april.tsv')
    indices = pd.read_csv(fn, sep='\t', dtype=str, keep_default_na=False,
                          na_values=[])
    return indices
