import os
import pandas
from pandas.testing import assert_frame_equal
from unittest import TestCase

from metapool.abs_quant_sample_info import \
    VOL_HOMOGENATE_ALIQUOT_INPUT_UL_KEY, VOL_HOMOGENATE_ALIQUOT_INPUT_ML_KEY, \
    MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_G_KEY, \
    MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_MG_KEY, \
    MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_G_KEY, \
    MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_MG_KEY, \
    MASS_STORAGE_TUBE_ONLY_G_KEY, DENSITY_STORAGE_LIQUID_KG_L_KEY, \
    DENSITY_STORAGE_LIQUID_G_ML_KEY, STORAGE_LIQUID_LOT_NUM_STR_KEY, \
    DENSITY_SAMPLE_G_ML_KEY, DENSITY_SAMPLE_KG_L_KEY, \
    STORAGE_LIQUID_TYPE_KEY, DENSITY_STOOL_STANDARDIZED_G_ML_KEY, \
    CALC_MASS_SAMPLE_IN_STORAGE_TUBE_MG_KEY, \
    CALC_MASS_SAMPLE_IN_STORAGE_TUBE_G_KEY, \
    CALC_MASS_STORAGE_LIQUID_ONLY_G_KEY, \
    CALC_VOL_STORAGE_LIQUID_ONLY_ML_KEY, \
    CALC_VOL_SAMPLE_IN_STORAGE_TUBE_ML_KEY, \
    CALC_VOL_HOMOGENATE_IN_STORAGE_TUBE_ML_KEY, \
    CALC_DENSITY_HOMOGENATE_G_ML_KEY, \
    CALC_MASS_HOMOGENATE_ALIQUOT_INPUT_G_KEY, \
    CALC_MASS_SAMPLE_ALIQUOT_INPUT_G_KEY, \
    CALC_MASS_SAMPLE_ALIQUOT_INPUT_MG_KEY, \
    CALC_MASS_STORAGE_LIQUID_ALIQUOT_INPUT_G, \
    _validate_input_df, _add_config_metadata, _calc_abs_quant_metadata, \
    add_abs_quant_metadata


class TestAbsQuantSampleInfo(TestCase):
    storage_liquid_type = "zymo_dna_rna_shield"

    def test_add_abs_quant_metadata(self):
        valid_core_input_dict = {
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_G_KEY:
                [16.9, 15.9],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_MG_KEY:
                [16900.0, 15900.0],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_MG_KEY:
                [20478, 19478],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_G_KEY:
                [20.478, 19.478],
            VOL_HOMOGENATE_ALIQUOT_INPUT_ML_KEY: [0.1, 0.1],
            VOL_HOMOGENATE_ALIQUOT_INPUT_UL_KEY: [100, 100],
            # NB: here lot numbers are being put in as numbers, which is
            # ok because the method under test converts them to strings
            STORAGE_LIQUID_LOT_NUM_STR_KEY: [123456789, 223456789],
            "some_other_column_we_can_ignore": ["blue", "green"]
        }

        adds_dict = {
            STORAGE_LIQUID_TYPE_KEY:
                [self.storage_liquid_type, self.storage_liquid_type],
            MASS_STORAGE_TUBE_ONLY_G_KEY: [7.18, 7.18],
            DENSITY_STORAGE_LIQUID_G_ML_KEY: [1.11, 1.01],
            DENSITY_STORAGE_LIQUID_KG_L_KEY: [1.11, 1.01],
            DENSITY_SAMPLE_G_ML_KEY: [1.06, 1.06],
            DENSITY_SAMPLE_KG_L_KEY: [1.06, 1.06],
            CALC_MASS_SAMPLE_IN_STORAGE_TUBE_MG_KEY: [3578.0, 3578.0],
            CALC_MASS_SAMPLE_IN_STORAGE_TUBE_G_KEY: [3.578, 3.578],
            CALC_MASS_STORAGE_LIQUID_ONLY_G_KEY: [9.72, 8.72],
            CALC_VOL_STORAGE_LIQUID_ONLY_ML_KEY: [8.756756757, 8.633663366],
            CALC_VOL_SAMPLE_IN_STORAGE_TUBE_ML_KEY: [3.375471698, 3.375471698],
            CALC_VOL_HOMOGENATE_IN_STORAGE_TUBE_ML_KEY:
                [12.13222845, 12.00913506],
            CALC_DENSITY_HOMOGENATE_G_ML_KEY: [1.096088822, 1.024053767],
            CALC_MASS_HOMOGENATE_ALIQUOT_INPUT_G_KEY:
                [0.109608882, 0.102405377],
            CALC_MASS_SAMPLE_ALIQUOT_INPUT_G_KEY: [0.029491697, 0.029793986],
            CALC_MASS_SAMPLE_ALIQUOT_INPUT_MG_KEY: [29.49169654, 29.79398583],
            CALC_MASS_STORAGE_LIQUID_ALIQUOT_INPUT_G:
                [0.080117186, 0.072611391]
        }

        expected_out_dict = valid_core_input_dict | adds_dict
        expected_out_dict[STORAGE_LIQUID_LOT_NUM_STR_KEY] = \
            [str(x) for x in expected_out_dict[STORAGE_LIQUID_LOT_NUM_STR_KEY]]
        expected_out_df = pandas.DataFrame.from_dict(expected_out_dict)

        curr_dir = os.path.dirname(os.path.abspath(__file__))
        test_config_fp = os.path.join(
            curr_dir, "data", "alt_abs_quant_sample_info_calc.yml")

        input_df = pandas.DataFrame.from_dict(valid_core_input_dict)
        output_df = add_abs_quant_metadata(
            input_df, DENSITY_STOOL_STANDARDIZED_G_ML_KEY,
            self.storage_liquid_type, test_config_fp)

        assert_frame_equal(expected_out_df, output_df)

    def test_add_abs_quant_metadata_w_error(self):
        invalid_core_input_dict = {
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_G_KEY:
                [20.478, 19.478],
            VOL_HOMOGENATE_ALIQUOT_INPUT_ML_KEY: [0.1, 0.1],
            VOL_HOMOGENATE_ALIQUOT_INPUT_UL_KEY: [100, 100],
            # NB: here lot numbers are being put in as numbers, which is
            # ok because the method under test converts them to strings
            STORAGE_LIQUID_LOT_NUM_STR_KEY: [123456789, 223456789],
            "some_other_column_we_can_ignore": ["blue", "green"]
        }

        curr_dir = os.path.dirname(os.path.abspath(__file__))
        test_config_fp = os.path.join(
            curr_dir, "data", "alt_abs_quant_sample_info_calc.yml")

        input_df = pandas.DataFrame.from_dict(invalid_core_input_dict)

        err_msg = "The input_df is missing the following required"
        with self.assertRaisesRegex(ValueError, err_msg):
            add_abs_quant_metadata(
                input_df, DENSITY_STOOL_STANDARDIZED_G_ML_KEY,
                self.storage_liquid_type, test_config_fp)

    def test__validate_input_df(self):
        valid_core_input_dict = {
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_G_KEY: [16.9],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_MG_KEY:
                [16900.0],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_MG_KEY: [20478],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_G_KEY: [20.478],
            VOL_HOMOGENATE_ALIQUOT_INPUT_ML_KEY: [0.1],
            VOL_HOMOGENATE_ALIQUOT_INPUT_UL_KEY: [100],
            STORAGE_LIQUID_LOT_NUM_STR_KEY: ["123456789"],
            "some_other_column_we_can_ignore": ["blue"]
        }

        input_df = pandas.DataFrame.from_dict(valid_core_input_dict)
        try:
            _validate_input_df(input_df)
        except ValueError:
            self.fail("Raised ValueError incorrectly")

    def test__validate_input_df_error(self):
        invalid_core_input_dict = {
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_MG_KEY: [20478],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_G_KEY: [20.478],
            VOL_HOMOGENATE_ALIQUOT_INPUT_UL_KEY: [100],
            "some_other_column_we_can_ignore": ["blue"]
        }

        input_df = pandas.DataFrame.from_dict(invalid_core_input_dict)
        err_msg = "The input_df is missing the following required"
        with self.assertRaisesRegex(ValueError, err_msg):
            _validate_input_df(input_df)

    def test__add_config_metadata(self):
        valid_core_input_dict = {
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_G_KEY:
                [16.9, 15.9],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_MG_KEY:
                [16900.0, 15900.0],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_MG_KEY:
                [20478, 19478],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_G_KEY:
                [20.478, 19.478],
            VOL_HOMOGENATE_ALIQUOT_INPUT_ML_KEY: [0.1, 0.1],
            VOL_HOMOGENATE_ALIQUOT_INPUT_UL_KEY: [100, 100],
            DENSITY_STORAGE_LIQUID_G_ML_KEY: [1.1027, 1.1027],
            DENSITY_STORAGE_LIQUID_KG_L_KEY: [1.1027, 1.1027],
            STORAGE_LIQUID_LOT_NUM_STR_KEY: ["123456789", "223456789"],
            "some_other_column_we_can_ignore": ["blue", "green"]
        }

        config_dict = {
            STORAGE_LIQUID_TYPE_KEY: {
                self.storage_liquid_type: {
                    STORAGE_LIQUID_LOT_NUM_STR_KEY: {
                        "123456789": {
                            DENSITY_STORAGE_LIQUID_G_ML_KEY: 1.11
                        },
                        "223456789": {
                            DENSITY_STORAGE_LIQUID_G_ML_KEY: 1.01
                        }
                    },
                    MASS_STORAGE_TUBE_ONLY_G_KEY: 7.18
                },
            },
            DENSITY_STOOL_STANDARDIZED_G_ML_KEY: 1.06
        }

        config_adds_dict = {
            STORAGE_LIQUID_TYPE_KEY:
                [self.storage_liquid_type, self.storage_liquid_type],
            MASS_STORAGE_TUBE_ONLY_G_KEY: [7.18, 7.18],
            DENSITY_STORAGE_LIQUID_G_ML_KEY: [1.11, 1.01],
            DENSITY_STORAGE_LIQUID_KG_L_KEY: [1.11, 1.01],
            DENSITY_SAMPLE_G_ML_KEY: [1.06, 1.06],
            DENSITY_SAMPLE_KG_L_KEY: [1.06, 1.06],
        }

        expected_out_dict = valid_core_input_dict | config_adds_dict
        expected_out_df = pandas.DataFrame.from_dict(expected_out_dict)

        input_df = pandas.DataFrame.from_dict(valid_core_input_dict)
        output_df = _add_config_metadata(
            input_df, config_dict, DENSITY_STOOL_STANDARDIZED_G_ML_KEY,
            self.storage_liquid_type)

        assert_frame_equal(expected_out_df, output_df)

    def test__calc_abs_quant_metadata(self):
        single_calc_input_dict = {
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_G_KEY:
                [17.0416],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_BEFORE_SAMPLE_MG_KEY:
                [17041.6],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_MG_KEY: [20478],
            MASS_STORAGE_TUBE_AND_STORAGE_LIQUID_AFTER_SAMPLE_G_KEY: [20.478],
            DENSITY_SAMPLE_G_ML_KEY: [1.06],
            MASS_STORAGE_TUBE_ONLY_G_KEY: [7],
            VOL_HOMOGENATE_ALIQUOT_INPUT_ML_KEY: [0.1],
            DENSITY_STORAGE_LIQUID_G_ML_KEY: [1.1027],
        }

        single_calc_output_adds_dict = {
            CALC_MASS_SAMPLE_IN_STORAGE_TUBE_MG_KEY: [3436.4],
            CALC_MASS_SAMPLE_IN_STORAGE_TUBE_G_KEY: [3.4364],
            CALC_MASS_STORAGE_LIQUID_ONLY_G_KEY: [10.0416],
            CALC_VOL_STORAGE_LIQUID_ONLY_ML_KEY: [9.106375261],
            CALC_VOL_SAMPLE_IN_STORAGE_TUBE_ML_KEY: [3.241886792],
            CALC_VOL_HOMOGENATE_IN_STORAGE_TUBE_ML_KEY: [12.34826205],
            CALC_DENSITY_HOMOGENATE_G_ML_KEY: [1.091489632],
            CALC_MASS_HOMOGENATE_ALIQUOT_INPUT_G_KEY: [0.109148963],
            CALC_MASS_SAMPLE_ALIQUOT_INPUT_G_KEY: [0.027829017],
            CALC_MASS_SAMPLE_ALIQUOT_INPUT_MG_KEY: [27.829017],
            CALC_MASS_STORAGE_LIQUID_ALIQUOT_INPUT_G: [0.081319946]
        }

        input_df = pandas.DataFrame.from_dict(single_calc_input_dict)
        output_df = _calc_abs_quant_metadata(input_df)

        expected_out_df = single_calc_input_dict | single_calc_output_adds_dict
        expected_out_df = pandas.DataFrame.from_dict(expected_out_df)

        assert_frame_equal(expected_out_df, output_df)
