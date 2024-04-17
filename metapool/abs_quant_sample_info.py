import os
import pandas
import yaml

from typing import Optional

# Relative path to default config file
DEFAULT_CONFIG_FP = 'data/abs_quant_sample_info_calc.yml'

# TERMS for metadata columns:
# sample: raw sample, i.e., stool
# storage liquid: storage liquid used with sample, ie DNA/RNA shield
# storage tube: tube used to store sample and storage liquid;
#   storage tube also includes labels
# homogenate: sample + storage liquid, homogenously combined
# aliquot: refers to the result of aliquoting material into the matrix tube
#   for extraction

MASS_STORAGE_TUBE_ONLY_G_KEY = 'mass_storage_tube_only_g'
MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_G_KEY = \
    'mass_storage_tube_and_storage_liquid_before_sample_g'
MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_MG_KEY = \
    'mass_storage_tube_and_storage_liquid_before_sample_mg'
MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_G_KEY = \
    'mass_storage_tube_and_storage_liquid_after_sample_g'
MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_MG_KEY = \
    'mass_storage_tube_and_storage_liquid_after_sample_mg'
VOL_HOMOGENATE_ALIQUOT_INPUT_UL_KEY = 'vol_homogenate_aliquot_input_ul'
VOL_HOMOGENATE_ALIQUOT_INPUT_ML_KEY = 'vol_homogenate_aliquot_input_ml'
DENSITY_STORAGE_LIQUID_KG_L_KEY = 'density_storage_liquid_kg_l'
DENSITY_STORAGE_LIQUID_G_ML_KEY = 'density_storage_liquid_g_ml'
STORAGE_LIQUID_LOT_NUM_STR_KEY = 'storage_liquid_lot_number'
DENSITY_SAMPLE_G_ML_KEY = 'density_sample_g_ml'
DENSITY_SAMPLE_KG_L_KEY = 'density_sample_kg_l'
STORAGE_LIQUID_TYPE_KEY = 'storage_liquid_type'

DENSITY_STOOL_STANDARDIZED_G_ML_KEY = "density_stool_standardized_g_ml"

CALC_MASS_SAMPLE_IN_STORAGE_TUBE_MG_KEY = 'calc_mass_sample_in_storage_tube_mg'
CALC_MASS_SAMPLE_IN_STORAGE_TUBE_G_KEY = 'calc_mass_sample_in_storage_tube_g'
CALC_MASS_STORAGE_LIQUID_ONLY_G_KEY = 'calc_mass_storage_liquid_only_g'
CALC_VOL_STORAGE_LIQUID_ONLY_ML_KEY = 'calc_vol_storage_liquid_only_ml'
CALC_VOL_SAMPLE_IN_STORAGE_TUBE_ML_KEY = 'calc_vol_sample_in_storage_tube_ml'
CALC_VOL_HOMOGENATE_IN_STORAGE_TUBE_ML_KEY = \
    'calc_vol_homogenate_in_storage_tube_ml'
CALC_DENSITY_HOMOGENATE_G_ML_KEY = 'calc_density_homogenate_g_ml'
CALC_MASS_HOMOGENATE_ALIQUOT_INPUT_G_KEY = \
    'calc_mass_homogenate_aliquot_input_g'
CALC_MASS_SAMPLE_ALIQUOT_INPUT_G_KEY = 'calc_mass_sample_aliquot_input_g'
CALC_MASS_SAMPLE_ALIQUOT_INPUT_MG_KEY = 'calc_mass_sample_aliquot_input_mg'
CALC_MASS_STORAGE_LIQUID_ALIQUOT_INPUT_G = \
    'calc_mass_storage_liquid_aliquot_input_g'

INPUT_COLUMNS = [STORAGE_LIQUID_TYPE_KEY]

CONFIG_COLUMNS = [
    DENSITY_STORAGE_LIQUID_KG_L_KEY, DENSITY_STORAGE_LIQUID_G_ML_KEY,
    MASS_STORAGE_TUBE_ONLY_G_KEY,
    DENSITY_SAMPLE_G_ML_KEY, DENSITY_SAMPLE_KG_L_KEY
]

CORE_COLUMNS = [
    VOL_HOMOGENATE_ALIQUOT_INPUT_UL_KEY, VOL_HOMOGENATE_ALIQUOT_INPUT_ML_KEY,
    MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_G_KEY,
    MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_MG_KEY,
    MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_G_KEY,
    MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_MG_KEY,
    STORAGE_LIQUID_LOT_NUM_STR_KEY]

SLURRY_CALC_COLUMNS = INPUT_COLUMNS + CONFIG_COLUMNS + CORE_COLUMNS


def add_abs_quant_metadata(
        input_df: pandas.DataFrame, sample_density_g_ml_key: str,
        storage_liquid_type: str, config_fp: Optional[str] = None) \
        -> pandas.DataFrame:
    """
    Add metadata for absolute quantification calculations to the dataframe.

    Parameters
    ----------
    input_df: pandas.DataFrame
        Dataframe containing at least all the columns in CORE_COLUMNS.
    sample_density_g_ml_key: str
        Dictionary key for the density of the sample type for the samples in
        the input df (usually stool).
    storage_liquid_type: str
        Type of storage liquid used with the samples in input df.
    config_fp: str
        Filepath to the yaml configuration file.  If None, the default
        configuration file is used.

    Returns
    -------

    calc_df: pandas.DataFrame
        Extension of input calc_df with extensive additional columns for
        absolute quantification calculations.
    """

    _validate_input_df(input_df)
    calc_df = input_df.copy()

    # ensure that the lot "number" (really lot ID) is stored a string
    calc_df[STORAGE_LIQUID_LOT_NUM_STR_KEY] = calc_df[
        STORAGE_LIQUID_LOT_NUM_STR_KEY].astype(str)

    config_dict = _read_config_metadata(config_fp)
    calc_df = _add_config_metadata(
        calc_df, config_dict, sample_density_g_ml_key,
        storage_liquid_type)

    calc_df = _calc_abs_quant_metadata(calc_df)
    return calc_df


def _validate_input_df(input_df: pandas.DataFrame):
    """
    Validate that the input dataframe contains the required columns.

    Parameters
    ----------
    input_df: pandas.DataFrame
        Dataframe containing at least all the columns in CORE_COLUMNS.

    Raises
    -------
    ValueError
        If the input_df doesn't contain all the columns in CORE_COLUMNS.
    """

    input_cols = set(input_df.columns.tolist())
    core_cols = set(CORE_COLUMNS)
    if not core_cols.issubset(input_cols):
        raise ValueError(
            f"The input_df is missing the following required column(s): "
            f"{', '.join(core_cols - input_cols)}")


def _read_config_metadata(config_fp: str = None) -> dict:
    """
    Read the yaml configuration file into a dictionary.

    Parameters
    ----------
    config_fp: str
        Filepath to the yaml configuration file.  If None, the default
        configuration file is used.

    Returns
    -------
    config_dict: dict
        Dictionary containing the configuration metadata from the config file.

    """
    if config_fp is None:
        curr_dir = os.path.dirname(os.path.abspath(__file__))
        config_fp = os.path.join(curr_dir, DEFAULT_CONFIG_FP)

    with open(config_fp, "r") as f:
        config_dict = yaml.safe_load(f)
    return config_dict


def _add_config_metadata(
        calc_df: pandas.DataFrame, config_dict: dict,
        sample_density_g_ml_key: str, storage_liquid_type: str) \
        -> pandas.DataFrame:
    """
    Add metadata from the config file to the dataframe.

    Parameters
    ----------
    calc_df: pandas.DataFrame
        Dataframe of per-sample metadata containing at least
        a STORAGE_LIQUID_LOT_NUM_STR_KEY column.
    config_dict: dict
        Dictionary of metadata from the configuration.  The expected structure
        is:

        STORAGE_LIQUID_TYPE_KEY: {
            <storage_liquid_type>: {
                STORAGE_LIQUID_LOT_NUM_STR_KEY: {
                    <storage_liquid_lot_num>: {
                        DENSITY_STORAGE_LIQUID_G_ML_KEY: ##.##
                    },
            MASS_STORAGE_TUBE_ONLY_G_KEY: ##.##
            },
        <sample_density_g_ml_key>: ##.##
        }

    sample_density_g_ml_key: str
        Dictionary key for the density of the sample type for the samples in
        the input calc_df (usually stool).
    storage_liquid_type: str
        Type of storage liquid used with the samples in calc_df.

    Returns
    -------
    calc_df: pandas.DataFrame
        Extension of input calc_df with the following additional columns:
            STORAGE_LIQUID_TYPE_KEY,
            MASS_STORAGE_TUBE_ONLY_G_KEY
            STORAGE_LIQUID_LOT_NUM_STR_KEY
            DENSITY_STORAGE_LIQUID_G_ML_KEY
            DENSITY_STORAGE_LIQUID_G_ML_KEY
            DENSITY_SAMPLE_G_ML_KEY
            DENSITY_SAMPLE_KG_L_KEY
    """

    # yes, these steps could be chained together, but doing them separately
    # makes it easier to see where errors come from if a key is missing
    storage_liquids_dict = config_dict[STORAGE_LIQUID_TYPE_KEY]
    storage_liquid_dict = storage_liquids_dict[storage_liquid_type]

    calc_df[STORAGE_LIQUID_TYPE_KEY] = storage_liquid_type
    calc_df[MASS_STORAGE_TUBE_ONLY_G_KEY] = \
        storage_liquid_dict[MASS_STORAGE_TUBE_ONLY_G_KEY]

    lot_nums_dict = storage_liquid_dict.get(STORAGE_LIQUID_LOT_NUM_STR_KEY)

    calc_df[DENSITY_STORAGE_LIQUID_G_ML_KEY] = (
        calc_df[STORAGE_LIQUID_LOT_NUM_STR_KEY].apply(
            lambda x: lot_nums_dict.get(x).get(
                DENSITY_STORAGE_LIQUID_G_ML_KEY)))
    calc_df[DENSITY_STORAGE_LIQUID_KG_L_KEY] = \
        calc_df[DENSITY_STORAGE_LIQUID_G_ML_KEY]

    # different sample types have different densities, so have to specify
    # which one to use for this experiment
    calc_df[DENSITY_SAMPLE_G_ML_KEY] = config_dict[sample_density_g_ml_key]
    calc_df[DENSITY_SAMPLE_KG_L_KEY] = calc_df[DENSITY_SAMPLE_G_ML_KEY]

    return calc_df


def _calc_abs_quant_metadata(calc_df: pandas.DataFrame) -> pandas.DataFrame:
    """
    Calculate the mass of sample in each aliquot, with supporting values.

    Parameters
    ----------
    calc_df: pandas.DataFrame
        Dataframe containing at least the following columns:
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_MG_KEY,
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_MG_KEY,
            DENSITY_SAMPLE_G_ML_KEY,
            MASS_STORAGE_TUBE_ONLY_G_KEY,
            VOL_HOMOGENATE_ALIQUOT_INPUT_ML_KEY,
            DENSITY_STORAGE_LIQUID_G_ML_KEY

    Returns
    -------
    calc_df: pandas.DataFrame
        Extension of input dataframe with the following columns added:
            CALC_MASS_SAMPLE_IN_STORAGE_TUBE_MG_KEY,
            CALC_MASS_SAMPLE_IN_STORAGE_TUBE_G_KEY,
            CALC_MASS_STORAGE_LIQUID_ONLY_G_KEY,
            CALC_VOL_STORAGE_LIQUID_ONLY_ML_KEY,
            CALC_VOL_SAMPLE_IN_STORAGE_TUBE_ML_KEY,
            CALC_VOL_HOMOGENATE_IN_STORAGE_TUBE_ML_KEY,
            CALC_DENSITY_HOMOGENATE_G_ML_KEY,
            CALC_MASS_HOMOGENATE_ALIQUOT_INPUT_G_KEY,
            CALC_MASS_SAMPLE_ALIQUOT_INPUT_G_KEY,
            CALC_MASS_SAMPLE_ALIQUOT_INPUT_MG_KEY,
            CALC_MASS_STORAGE_LIQUID_ALIQUOT_INPUT_G

    """
    # Calculate mass of sample in storage tube and mass of the storage liquid
    calc_df[CALC_MASS_SAMPLE_IN_STORAGE_TUBE_MG_KEY] = \
        calc_df[MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_MG_KEY] - \
        calc_df[MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_MG_KEY]
    calc_df[CALC_MASS_SAMPLE_IN_STORAGE_TUBE_G_KEY] = \
        calc_df[CALC_MASS_SAMPLE_IN_STORAGE_TUBE_MG_KEY] / 1000
    calc_df[CALC_MASS_STORAGE_LIQUID_ONLY_G_KEY] = \
        calc_df[MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_G_KEY] - \
        calc_df[MASS_STORAGE_TUBE_ONLY_G_KEY]

    # Calculate volume of storage liquid only and
    # volume of sample in storage tube
    calc_df[CALC_VOL_STORAGE_LIQUID_ONLY_ML_KEY] = (
            calc_df[CALC_MASS_STORAGE_LIQUID_ONLY_G_KEY] /
            calc_df[DENSITY_STORAGE_LIQUID_G_ML_KEY])
    calc_df[CALC_VOL_SAMPLE_IN_STORAGE_TUBE_ML_KEY] = \
        ((
                 calc_df[
                     MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_G_KEY] -
                 calc_df[
                     MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_G_KEY]) /   # noqa: E501
         calc_df[DENSITY_SAMPLE_G_ML_KEY])

    # Calculate volume of homogenate in storage tube
    calc_df[CALC_VOL_HOMOGENATE_IN_STORAGE_TUBE_ML_KEY] = \
        calc_df[CALC_VOL_STORAGE_LIQUID_ONLY_ML_KEY] + \
        calc_df[CALC_VOL_SAMPLE_IN_STORAGE_TUBE_ML_KEY]

    # Calculate density of homogenate
    calc_df[CALC_DENSITY_HOMOGENATE_G_ML_KEY] = \
        ((
                 calc_df[
                     MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_G_KEY] -
                 calc_df[MASS_STORAGE_TUBE_ONLY_G_KEY]) /
         calc_df[CALC_VOL_HOMOGENATE_IN_STORAGE_TUBE_ML_KEY])

    # Calculate mass of the homogenate aliquotted into the matrix tube and
    # mass of the sample aliquotted into the matrix tube
    calc_df[CALC_MASS_HOMOGENATE_ALIQUOT_INPUT_G_KEY] = \
        calc_df[CALC_DENSITY_HOMOGENATE_G_ML_KEY] * \
        calc_df[VOL_HOMOGENATE_ALIQUOT_INPUT_ML_KEY]
    calc_df[CALC_MASS_SAMPLE_ALIQUOT_INPUT_G_KEY] = \
        ((
                 (calc_df[DENSITY_SAMPLE_G_ML_KEY] *
                  calc_df[DENSITY_STORAGE_LIQUID_G_ML_KEY] *
                  calc_df[VOL_HOMOGENATE_ALIQUOT_INPUT_ML_KEY]) -
                 (calc_df[DENSITY_SAMPLE_G_ML_KEY] *
                  calc_df[CALC_DENSITY_HOMOGENATE_G_ML_KEY] *
                  calc_df[VOL_HOMOGENATE_ALIQUOT_INPUT_ML_KEY])) /
         (calc_df[DENSITY_STORAGE_LIQUID_G_ML_KEY] -
          calc_df[DENSITY_SAMPLE_G_ML_KEY]))

    # Calculate mg of sample aliquotted into matrix tube and
    # mass of the storage liquid aliquotted into the matrix tube
    calc_df[CALC_MASS_SAMPLE_ALIQUOT_INPUT_MG_KEY] = \
        calc_df[CALC_MASS_SAMPLE_ALIQUOT_INPUT_G_KEY] * 1000
    calc_df[CALC_MASS_STORAGE_LIQUID_ALIQUOT_INPUT_G] = \
        calc_df[CALC_MASS_HOMOGENATE_ALIQUOT_INPUT_G_KEY] - \
        calc_df[CALC_MASS_SAMPLE_ALIQUOT_INPUT_G_KEY]

    return calc_df
