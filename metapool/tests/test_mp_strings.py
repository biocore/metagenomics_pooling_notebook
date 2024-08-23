from unittest import TestCase

from metapool.mp_strings import _split_plate_name, \
    get_main_project_from_plate_name, get_plate_num_from_plate_name, \
    parse_project_name, get_short_name_and_id, get_qiita_id_from_project_name


class TestMpStrings(TestCase):
    def test__split_plate_name_w_Plate(self):
        plate_name = "Celeste_Adaptation_12986_Plate_16"
        obs = _split_plate_name(plate_name)
        self.assertEqual(obs, ("Celeste_Adaptation_12986", "16"))

    def test__split_plate_name_wo_Plate(self):
        plate_name = "Celeste_Adaptation_12986_16"
        obs = _split_plate_name(plate_name)
        self.assertEqual(obs, ("Celeste_Adaptation_12986", "16"))

    def test__split_plate_name_malformed(self):
        plate_name = "CelesteAdaptation12986Plate16"
        err_msg = "Plate name 'CelesteAdaptation12986Plate16' is malformed."
        with self.assertRaisesRegex(ValueError, err_msg):
            _split_plate_name(plate_name)

    def test_get_main_project_from_plate_name_w_Plate(self):
        plate_name = "Celeste_Adaptation_12986_Plate_16"
        obs = get_main_project_from_plate_name(plate_name)
        self.assertEqual(obs, "Celeste_Adaptation_12986")

    def test_get_main_project_from_plate_name_wo_Plate(self):
        plate_name = "Celeste_Adaptation_12986_16"
        obs = get_main_project_from_plate_name(plate_name)
        self.assertEqual(obs, "Celeste_Adaptation_12986")

    def test_get_plate_num_from_plate_name_w_Plate(self):
        plate_name = "Celeste_Adaptation_12986_Plate_16"
        obs = get_plate_num_from_plate_name(plate_name)
        self.assertEqual(obs, "16")

    def test_get_plate_num_from_plate_name_wo_Plate(self):
        plate_name = "Celeste_Adaptation_12986_16"
        obs = get_plate_num_from_plate_name(plate_name)
        self.assertEqual(obs, "16")

    def test_parse_project_name(self):
        exp = {
            "qiita_id": "1161",
            "short_project_name": "A_Feist",
            "full_project_name": "A_Feist_1161"
        }

        obs = parse_project_name("A_Feist_1161")
        self.assertEqual(obs, exp)

    def test_parse_project_name_err_no_qiita_id(self):
        with self.assertRaisesRegex(
                ValueError, "'A_Feist' does not contain a Qiita-ID."):
            parse_project_name("A_Feist")

    def test_parse_project_name_err_missing(self):
        with self.assertRaisesRegex(
                ValueError, "project_name cannot be None or empty string"):
            parse_project_name("")

        with self.assertRaisesRegex(
                ValueError, "project_name cannot be None or empty string"):
            parse_project_name(None)

    def test_get_short_name_and_id(self):
        exp = ("A_Feist", "1161")
        obs = get_short_name_and_id("A_Feist_1161")
        self.assertEqual(obs, exp)

    def test_get_short_name_and_id_no_qiita_id(self):
        obs = get_short_name_and_id("A_Feist")
        self.assertEqual(obs, ("A_Feist", None))

    def test_get_qiita_id_from_project_name(self):
        obs = get_qiita_id_from_project_name("A_Feist_1161")
        self.assertEqual(obs, "1161")

    def test_get_qiita_id_from_project_name_err_no_qiita_id(self):
        with self.assertRaisesRegex(
                ValueError, "'A_Feist' does not contain a Qiita-ID."):
            get_qiita_id_from_project_name("A_Feist")

    def test_get_qiita_id_from_project_name_err_missing(self):
        with self.assertRaisesRegex(
                ValueError, "project_name cannot be None or empty string"):
            get_qiita_id_from_project_name("")

        with self.assertRaisesRegex(
                ValueError, "project_name cannot be None or empty string"):
            get_qiita_id_from_project_name(None)
