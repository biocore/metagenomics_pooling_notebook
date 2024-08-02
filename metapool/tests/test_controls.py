import pandas as pd
import os

from unittest import TestCase
from metapool.controls import SAMPLE_NAME_KEY, SAMPLE_TYPE_KEY, \
    PRIMARY_STUDY_KEY, SECONDARY_STUDIES_KEY, \
    is_blank, get_all_projects_in_context


class ControlsTests(TestCase):
    def setUp(self):
        self.maxDiff = None
        self.path = os.path.dirname(__file__)

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

        with self.assertRaisesRegex(
                ValueError, "Expected exactly one record with sample name "
                            "'sample1', but found 2 records"):
            is_blank("sample1", sample_context)

        with self.assertRaisesRegex(
                ValueError, "Expected exactly one record with sample name "
                            "'sample14', but found 0 records"):
            is_blank("sample14", sample_context)

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
