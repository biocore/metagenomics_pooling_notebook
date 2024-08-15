import pandas
from pandas.testing import assert_frame_equal
import os

from unittest import TestCase
from metapool.notebook_support import join_dfs_from_files


class NotebookSupportTests(TestCase):
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.sa_partial_1_fp = os.path.join(self.data_dir, 'sa_partial_1.csv')
        self.sa_partial_2_fp = os.path.join(self.data_dir, 'sa_partial_2.csv')
        self.sa_partial_2_w_qiita_id_fp = os.path.join(
            self.data_dir, 'sa_partial_2_w_qiita_id.csv')

        self.basic_join_df = pandas.DataFrame([
            {'sample_name': '41B.Month6.1', 'TubeCode': '0363132553'},
            {'sample_name': '41B.Month6.10', 'TubeCode': '0363132554'},
            {'sample_name': '41B.Month6.13', 'TubeCode': '0363132557'},
            {'sample_name': '20001CC.B.FW2', 'TubeCode': '0363134233'},
            {'sample_name': '20001CF.A.FW2', 'TubeCode': '0363134234'},
        ])

    def test_join_dfs_from_files_basic(self):
        # test that the function can join dataframes, no specifiied
        # optional columns, dtypes, unique columns, etc.
        # Happily ignores the qiita_id column, which is in the files and not
        # unique, because nobody asked for it.
        input_fps = [self.sa_partial_1_fp, self.sa_partial_2_w_qiita_id_fp]
        joined_df = join_dfs_from_files(
            input_fps, ["sample_name", "TubeCode"], sep=",")
        assert_frame_equal(joined_df, self.basic_join_df)

        # throw error if required column(s) are missing
        with self.assertRaisesRegex(
                ValueError, r"sa_partial_2.csv' is missing required "
                            r"columns \['qiita_id'\]"):
            join_dfs_from_files(
                [self.sa_partial_1_fp, self.sa_partial_2_fp],
                ["sample_name", "qiita_id"], sep=",")

    def test_join_dfs_from_files_deluxe(self):
        # specify everything--optional columns, dtypes, unique columns
        exp_df = pandas.DataFrame([
            {'sample_name': '41B.Month6.1', 'TubeCode': 363132553,
             'qiita_id': '12986'},
            {'sample_name': '41B.Month6.10', 'TubeCode': 363132554,
             'qiita_id': '12986'},
            {'sample_name': '41B.Month6.13', 'TubeCode': 363132557,
             'qiita_id': '14577'},
            {'sample_name': '20001CC.B.FW2', 'TubeCode': 363134233,
             'qiita_id': None},
            {'sample_name': '20001CF.A.FW2', 'TubeCode': 363134234,
             'qiita_id': None},
        ])
        joined_df = join_dfs_from_files(
            [self.sa_partial_1_fp, self.sa_partial_2_fp],
            ["sample_name", "TubeCode"],
            opt_cols_to_extract=["qiita_id"],
            unique_cols=["sample_name"],
            dtype={'sample_name': str, 'TubeCode': int, 'qiita_id': str},
            sep=",")
        assert_frame_equal(joined_df, exp_df)

    def test_join_dfs_from_files_opt_col(self):
        # test that the function can join dataframes, incl optional columns

        # first, show correct behavior when optional is in the first file
        # examined but not in subsequent one(s)
        exp_df = pandas.DataFrame([
            {'sample_name': '41B.Month6.1', 'TubeCode': '0363132553',
             'qiita_id': '12986'},
            {'sample_name': '41B.Month6.10', 'TubeCode': '0363132554',
             'qiita_id': '12986'},
            {'sample_name': '41B.Month6.13', 'TubeCode': '0363132557',
             'qiita_id': '14577'},
            {'sample_name': '20001CC.B.FW2', 'TubeCode': '0363134233',
             'qiita_id': None},
            {'sample_name': '20001CF.A.FW2', 'TubeCode': '0363134234',
             'qiita_id': None},
        ])
        joined_df = join_dfs_from_files(
            [self.sa_partial_1_fp, self.sa_partial_2_fp],
            ["sample_name", "TubeCode"],
            opt_cols_to_extract=["qiita_id"], sep=",")
        assert_frame_equal(joined_df, exp_df)

        # now show correct behavior when optional column is NOT in
        # the first file examined but IS in subsequent one(s)
        exp_df = pandas.DataFrame([
            {'sample_name': '20001CC.B.FW2', 'TubeCode': '0363134233',
             'qiita_id': None},
            {'sample_name': '20001CF.A.FW2', 'TubeCode': '0363134234',
             'qiita_id': None},
            {'sample_name': '41B.Month6.1', 'TubeCode': '0363132553',
             'qiita_id': '12986'},
            {'sample_name': '41B.Month6.10', 'TubeCode': '0363132554',
             'qiita_id': '12986'},
            {'sample_name': '41B.Month6.13', 'TubeCode': '0363132557',
             'qiita_id': '14577'}
        ])
        joined_df = join_dfs_from_files(
            [self.sa_partial_2_fp, self.sa_partial_1_fp],
            ["sample_name", "TubeCode"],
            opt_cols_to_extract=["qiita_id"], sep=",")
        assert_frame_equal(joined_df, exp_df)

        # if optional column is in none of the dataframes, it won't be in the
        # output dataframe either
        input_fps = [self.sa_partial_1_fp, self.sa_partial_2_fp]
        joined_df = join_dfs_from_files(
            input_fps, ["sample_name", "TubeCode"],
            opt_cols_to_extract=["tube_id"], sep=",")
        assert_frame_equal(joined_df, self.basic_join_df)

    def test_join_dfs_from_files_dtypes(self):
        # specify a dtype for only one column, let the other be inferred
        exp_df = pandas.DataFrame([
            {'sample_name': '41B.Month6.1', 'TubeCode': 363132553},
            {'sample_name': '41B.Month6.10', 'TubeCode': 363132554},
            {'sample_name': '41B.Month6.13', 'TubeCode': 363132557},
            {'sample_name': '20001CC.B.FW2', 'TubeCode': 363134233},
            {'sample_name': '20001CF.A.FW2', 'TubeCode': 363134234},
        ])

        input_fps = [self.sa_partial_1_fp, self.sa_partial_2_fp]
        joined_df = join_dfs_from_files(
            input_fps, ["sample_name", "TubeCode"],
            dtype={'sample_name': str}, sep=",")
        assert_frame_equal(joined_df, exp_df)

    def test_join_dfs_from_files_unique(self):
        # happy behavior if all the requested unique columns are present and,
        # in fact, unique, even if other columns aren't
        exp_df = pandas.DataFrame([
            {'sample_name': '41B.Month6.1', 'TubeCode': '0363132553',
             'qiita_id': '12986'},
            {'sample_name': '41B.Month6.10', 'TubeCode': '0363132554',
             'qiita_id': '12986'},
            {'sample_name': '41B.Month6.13', 'TubeCode': '0363132557',
             'qiita_id': '14577'},
            {'sample_name': '20001CC.B.FW2', 'TubeCode': '0363134233',
             'qiita_id': '14577'},
            {'sample_name': '20001CF.A.FW2', 'TubeCode': '0363134234',
             'qiita_id': '10317'},
        ])
        input_fps = [self.sa_partial_1_fp, self.sa_partial_2_w_qiita_id_fp]
        joined_df = join_dfs_from_files(
            input_fps, ["sample_name", "TubeCode", "qiita_id"],
            unique_cols=["sample_name", "TubeCode"], sep=",")
        assert_frame_equal(joined_df, exp_df)

        # throw error if specified unique column is not in the required cols
        with self.assertRaisesRegex(
                ValueError, "All unique_cols must be in req_cols_to_extract"):
            join_dfs_from_files(
                [self.sa_partial_1_fp, self.sa_partial_2_fp],
                ["sample_name"],
                unique_cols=["TubeCode"], sep=",")

        # throw error if there are duplicate values in the unique column(s)
        with self.assertRaisesRegex(
                ValueError, r"Duplicate 'qiita_id' values found in files: "
                            r"\['12986', '14577'\]"):
            join_dfs_from_files(
                [self.sa_partial_1_fp, self.sa_partial_2_w_qiita_id_fp],
                ["sample_name", "TubeCode", "qiita_id"],
                unique_cols=["TubeCode", "qiita_id"], sep=",")

    # def test__generate_sample_context(self):
    #     # note that test must be run from the root of the metapool module
    #     # for this path to resolve
    #     path = os.path.dirname(__file__)
    #     plate_df_path = os.path.join(path,
    #         '../../notebooks/test_output/QC/YYYY_MM_DD_Celeste_Adaptation_df_A.txt')
    #         # 'notebooks/test_output/QC/YYYY_MM_DD_Celeste_Adaptation_df_A.txt'
    #
    #     plate_df = pd.read_csv(plate_df_path, dtype="str", sep='\t')
    #
    #     generated_context = _generate_sample_context(plate_df)
    #     self.assertTrue(False)
