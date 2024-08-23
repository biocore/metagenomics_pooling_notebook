import pandas as pd
import os

from unittest import TestCase
from metapool.mp_strings import SAMPLE_NAME_KEY, SAMPLE_TYPE_KEY, \
    PRIMARY_STUDY_KEY, SECONDARY_STUDIES_KEY
from metapool.controls import is_blank, get_all_projects_in_context, \
    get_controls_details_from_context, make_manual_control_details, \
    get_delimited_controls_details_from_compressed_plate


class ControlsTests(TestCase):
    def setUp(self):
        self.maxDiff = None
        self.path = os.path.dirname(__file__)

        self.compressed_plate_df = pd.DataFrame([
            {"Sample": "sample1", "Blank": True,
             "Project Name": "Study_1", "Project Plate": "Study_1_Plate_11"},
            {"Sample": "sample2", "Blank": False,
             "Project Name": "Study_1", "Project Plate": "Study_1_Plate_11"},
            {"Sample": "sample3", "Blank": False,
             "Project Name": "Study_4", "Project Plate": "Study_1_Plate_11"},
            {"Sample": "sample4", "Blank": False,
             "Project Name": "Study_5", "Project Plate": "Study_1_Plate_11"},
            {"Sample": "BLANK.2", "Blank": True,
             "Project Name": "Study_2", "Project Plate": "Study_2_Plate_21"},
            {"Sample": "sm1", "Blank": False,
             "Project Name": "Study_2", "Project Plate": "Study_2_Plate_21"},
            {"Sample": "sm2", "Blank": False,
             "Project Name": "Study_3", "Project Plate": "Study_2_Plate_21"},
            {"Sample": "sm3", "Blank": False,
             "Project Name": "Study_6", "Project Plate": "Study_2_Plate_21"},
            {"Sample": "BLANK.3", "Blank": True,
             "Project Name": "Study_3", "Project Plate": "Study_3_Plate_13"},
            {"Sample": "samp1", "Blank": False,
             "Project Name": "Study_3", "Project Plate": "Study_3_Plate_13"},
            {"Sample": "blank4", "Blank": True,
             "Project Name": "Study_10", "Project Plate": "Study_10_Plate_1"},
            {"Sample": "samples1", "Blank": False,
             "Project Name": "Study_10", "Project Plate": "Study_10_Plate_1"},
            {"Sample": "samples2", "Blank": False,
             "Project Name": "Study_11", "Project Plate": "Study_10_Plate_1"}
        ])

    def test_is_blank_w_context(self):
        sample_context = pd.DataFrame({
            SAMPLE_NAME_KEY: ["sample1", "BLANK.2", "BLANK.3", "blank4"],
            SAMPLE_TYPE_KEY: ["control blank", "positive control",
                              "control blank", "positive control"],
            PRIMARY_STUDY_KEY: ["1", "2", "3", "10"],
            SECONDARY_STUDIES_KEY: ["4;5", "6", "", "11"]
        })

        # this one would NOT be expected to be a blank if you looked at its
        # name, but IS because of its sample type
        self.assertTrue(is_blank("sample1", sample_context))
        # this one WOULD be expected to be a blank if you looked at its
        # name, but ISN'T because of its sample type
        self.assertFalse(is_blank("BLANK.2", sample_context))
        # this one is a blank either way you decide
        self.assertTrue(is_blank("BLANK.3", sample_context))
        # this one is NOT a blank either way you decide--just having the
        # word "blank" in the name doesn't make it a blank, has to be caps
        self.assertFalse(is_blank("blank4", sample_context))

    def test_is_blank_w_context_err(self):
        sample_context = pd.DataFrame({
            SAMPLE_NAME_KEY: ["sample1", "BLANK.2", "BLANK.3", "sample1"],
            SAMPLE_TYPE_KEY: ["control blank", "positive control",
                              "control blank", "positive control"],
            PRIMARY_STUDY_KEY: ["1", "2", "3", "10"],
            SECONDARY_STUDIES_KEY: ["4;5", "6", "", "11"]
        })

        # malformed context section
        with self.assertRaisesRegex(
                ValueError, "Expected exactly one record with sample name "
                            "'sample1', but found 2 records"):
            is_blank("sample1", sample_context)

    def test_is_blank_wo_context(self):
        names = ["sample1", "BLANK.2", "BLANK.3", "blank4"]
        exp_vals = [False, True, True, False]
        for i in range(len(names)):
            self.assertEqual(is_blank(names[i]), exp_vals[i])

    def test_get_all_projects_in_context(self):
        sample_context = pd.DataFrame({
            SAMPLE_NAME_KEY: ["sample1", "BLANK.2", "BLANK.3", "blank4"],
            SAMPLE_TYPE_KEY: ["control blank", "positive control",
                              "control blank", "positive control"],
            PRIMARY_STUDY_KEY: ["1", "2", "3", "10"],
            SECONDARY_STUDIES_KEY: ["4;5", "6;3", "", "11"]
        })

        self.assertEqual(get_all_projects_in_context(sample_context),
                         ["1",  "10", "11", "2", "3", "4", "5", "6"])

    def test_get_all_projects_in_context_none(self):
        self.assertEqual(get_all_projects_in_context(None), [])

    def test_get_controls_details_from_context(self):
        sample_context = pd.DataFrame({
            SAMPLE_NAME_KEY: ["sample1", "BLANK.2", "BLANK.3", "blank4"],
            SAMPLE_TYPE_KEY: ["control blank", "positive control",
                              "control blank", "positive control"],
            PRIMARY_STUDY_KEY: ["1", "2", "3", "10"],
            SECONDARY_STUDIES_KEY: ["4;5", "6;3", "", "11"]
        })

        exp_details = {
            "sample1": {
                SAMPLE_NAME_KEY: "sample1",
                SAMPLE_TYPE_KEY: "control blank",
                PRIMARY_STUDY_KEY: "1",
                SECONDARY_STUDIES_KEY: ["4", "5"]
            },
            "BLANK.2": {
                SAMPLE_NAME_KEY: "BLANK.2",
                SAMPLE_TYPE_KEY: "positive control",
                PRIMARY_STUDY_KEY: "2",
                SECONDARY_STUDIES_KEY: ["3", "6"]
            },
            "BLANK.3": {
                SAMPLE_NAME_KEY: "BLANK.3",
                SAMPLE_TYPE_KEY: "control blank",
                PRIMARY_STUDY_KEY: "3",
                SECONDARY_STUDIES_KEY: []
            },
            "blank4": {
                SAMPLE_NAME_KEY: "blank4",
                SAMPLE_TYPE_KEY: "positive control",
                PRIMARY_STUDY_KEY: "10",
                SECONDARY_STUDIES_KEY: ["11"]
            }
        }

        obs_details = get_controls_details_from_context(sample_context)
        self.assertEqual(exp_details, obs_details)

    def test_get_controls_details_from_context_none(self):
        self.assertEqual(get_controls_details_from_context(None), {})

    def test_make_manual_control_details(self):
        exp_details = {
                SAMPLE_NAME_KEY: "sample1",
                SAMPLE_TYPE_KEY: "control blank",
                PRIMARY_STUDY_KEY: "1",
                SECONDARY_STUDIES_KEY: ["4", "5"]
            }

        obs_details = make_manual_control_details(
            "sample1", "1", ["4", "5"], "control blank")
        self.assertEqual(obs_details, exp_details)

    def test_make_manual_control_details_w_defaults(self):
        exp_details = {
                SAMPLE_NAME_KEY: "sample1",
                SAMPLE_TYPE_KEY: "control blank",
                PRIMARY_STUDY_KEY: "1",
                SECONDARY_STUDIES_KEY: []
            }

        obs_details = make_manual_control_details("sample1", "1")
        self.assertEqual(obs_details, exp_details)

    def test_get_delimited_controls_details_from_compressed_plate_w_mask(self):
        exp_details = [
            {SAMPLE_NAME_KEY: "sample1", SAMPLE_TYPE_KEY: "control blank",
             PRIMARY_STUDY_KEY: "1", SECONDARY_STUDIES_KEY: "4;5"},
            {SAMPLE_NAME_KEY: "BLANK.2", SAMPLE_TYPE_KEY: "control blank",
             PRIMARY_STUDY_KEY: "2", SECONDARY_STUDIES_KEY: "3;6"},
            {SAMPLE_NAME_KEY: "BLANK.3", SAMPLE_TYPE_KEY: "control blank",
             PRIMARY_STUDY_KEY: "3", SECONDARY_STUDIES_KEY: ""},
            {SAMPLE_NAME_KEY: "blank4", SAMPLE_TYPE_KEY: "control blank",
             PRIMARY_STUDY_KEY: "10", SECONDARY_STUDIES_KEY: "11"}
        ]

        obs_details = get_delimited_controls_details_from_compressed_plate(
            self.compressed_plate_df, self.compressed_plate_df["Blank"])

        self.assertListEqual(exp_details, obs_details)

    def test_get_delimited_controls_details_from_compressed_plate_wo_mask(
            self):

        exp_details = [
            {SAMPLE_NAME_KEY: "BLANK.2", SAMPLE_TYPE_KEY: "control blank",
                PRIMARY_STUDY_KEY: "2", SECONDARY_STUDIES_KEY: "3;6"},
            {SAMPLE_NAME_KEY: "BLANK.3", SAMPLE_TYPE_KEY: "control blank",
                PRIMARY_STUDY_KEY: "3", SECONDARY_STUDIES_KEY: ""}
        ]

        obs_details = get_delimited_controls_details_from_compressed_plate(
            self.compressed_plate_df)

        self.assertListEqual(exp_details, obs_details)
