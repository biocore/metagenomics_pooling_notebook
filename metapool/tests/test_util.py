import pandas
from pandas.testing import assert_frame_equal
import os

from unittest import TestCase
from metapool.util import join_dfs_from_files, \
    extend_sample_accession_df, extend_compression_layout_info, \
    _check_for_missing_df_ids


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

        self.studies_info = [
            {
                'Project Name': 'Celeste_Adaptation_12986',
                'Project Abbreviation': 'ADAPT',
                'sample_accession_fp': './adapt_sa.tsv',
                'qiita_metadata_fp': './adapt_metadata.txt',
                'experiment_design_description': 'isolate sequencing',
                'HumanFiltering': 'False',
                'Email': 'r@gmail.com'
            },
            {
                'Project Name': 'CHILD_15510',
                'Project Abbreviation': 'CHILD',
                'sample_accession_fp': './child_sa.tsv',
                'qiita_metadata_fp': './child_metadata.txt',
                'experiment_design_description': 'whole genome sequencing',
                'HumanFiltering': 'True',
                'Email': 'l@ucsd.edu'
            },
            {
                'Project Name': 'Celeste_Marmoset_14577',
                'Project Abbreviation': 'MARMO',
                'sample_accession_fp': './marmo_sa.tsv',
                'qiita_metadata_fp': './marmo_metadata.txt',
                'experiment_design_description': 'whole genome sequencing',
                'HumanFiltering': 'False',
                'Email': 'c@ucsd.edu'
            }
        ]

        self.compress_layout = [
            {
                # top left plate
                'Plate Position': 1,  # as int
                'Plate map file': './16_plate_map.tsv',
                'Project Name': 'Celeste_Adaptation_12986',
                'Project Plate': 'Plate_16',  # Plate_#
                'Plate elution volume': 70
            },
            {
                # top right plate
                'Plate Position': 2,  # as int
                'Plate map file': './17_plate_map.tsv',
                'Project Name': 'Celeste_Adaptation_12986',
                'Project Plate': 'Plate_17',  # Plate_#
                'Plate elution volume': 70
            },
            {
                # bottom left plate
                'Plate Position': 3,  # as int
                'Plate map file': './18_plate_map.tsv',
                'Project Name': 'Celeste_Adaptation_12986',
                'Project Plate': 'Plate_18',  # Plate_#
                'Plate elution volume': 70
            },
            {
                # bottom right plate
                'Plate Position': 4,  # as int
                'Plate map file': './CHILD_1000_plate_map.tsv',
                'Project Name': 'CHILD_15510',
                'Project Plate': 'Plate_1000',  # Plate_#
                'Plate elution volume': 70
            },
        ]

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

    def test_extend_sample_accession_df(self):
        samp_acc_df = pandas.DataFrame([
            {'sample_name': '41B.Month6.1', 'TubeCode': '0363132553'},
            {'sample_name': '41B.Month6.10', 'TubeCode': '0363132554'},
            {'sample_name': '41B.Month6.13', 'TubeCode': '0363132557'},
            {'sample_name': '20001CC.B.FW2', 'TubeCode': '0363134233'},
            {'sample_name': '20001CF.A.FW2', 'TubeCode': '0363134234'}
        ])

        metadata_df = pandas.DataFrame([
            {'sample_name': '12986.41B.Month6.1', 'qiita_study_id': '12986'},
            {'sample_name': '12986.41B.Month6.10', 'qiita_study_id': '12986'},
            {'sample_name': '14577.41B.Month6.13', 'qiita_study_id': '14577'},
            {'sample_name': '15510.20001CC.B.FW2', 'qiita_study_id': '15510'},
            {'sample_name': '15510.20001CF.A.FW2', 'qiita_study_id': '15510'}
        ])

        exp_df = pandas.DataFrame([
            {'sample_name': '41B.Month6.1',
             'TubeCode': '0363132553',
             'Project Name': 'Celeste_Adaptation_12986',
             'Project Abbreviation': 'ADAPT'},
            {'sample_name': '41B.Month6.10',
             'TubeCode': '0363132554',
             'Project Name': 'Celeste_Adaptation_12986',
             'Project Abbreviation': 'ADAPT'},
            {'sample_name': '41B.Month6.13',
             'TubeCode': '0363132557',
             'Project Name': 'Celeste_Marmoset_14577',
             'Project Abbreviation': 'MARMO'},
            {'sample_name': '20001CC.B.FW2',
             'TubeCode': '0363134233',
             'Project Name': 'CHILD_15510',
             'Project Abbreviation': 'CHILD'},
            {'sample_name': '20001CF.A.FW2',
             'TubeCode': '0363134234',
             'Project Name': 'CHILD_15510',
             'Project Abbreviation': 'CHILD'}
            ])

        obs_df = extend_sample_accession_df(
            samp_acc_df, self.studies_info, metadata_df)

        assert_frame_equal(exp_df, obs_df)

    def test__check_for_missing_df_ids(self):
        samp_acc_df = pandas.DataFrame([
            {'sample_name': '41B.Month6.1', 'TubeCode': '0363132553'},
            {'sample_name': '41B.Month6.10', 'TubeCode': '0363132554'}
        ])

        metadata_df = pandas.DataFrame([
            {'sample_name': '41B.Month6.1', 'qiita_study_id': '12986'},
            {'sample_name': '41B.Month6.10', 'qiita_study_id': '12986'},
            {'sample_name': '41B.Month6.13', 'qiita_study_id': '14577'},
            {'sample_name': '20001CC.B.FW2', 'qiita_study_id': '15510'},
            {'sample_name': '20001CF.A.FW2', 'qiita_study_id': '15510'}
        ])

        # return nothing: all sample names in samp_acc_df are in metadata_df
        _check_for_missing_df_ids(samp_acc_df, metadata_df, 'sample_name',
                                  'sample_accession', 'metadata')

        # now remove a sample name from metadata_df; should raise error
        metadata_df = metadata_df.iloc[1:]
        err_msg = ("Some sample_name values in the sample_accession dataframe "
                   "are not in the metadata dataframe: 41B.Month6.1")
        with self.assertRaisesRegex(ValueError, err_msg):
            _check_for_missing_df_ids(samp_acc_df, metadata_df, 'sample_name',
                                      'sample_accession', 'metadata')

    def test_extend_compression_layout_info(self):
        exp_layout = [
            {
                # top left plate
                'Plate Position': 1,  # as int
                'Plate map file': './16_plate_map.tsv',
                'Project Name': 'Celeste_Adaptation_12986',
                'Project Abbreviation': 'ADAPT',
                'Project Plate': 'Plate_16',  # Plate_#
                'Plate elution volume': 70
            },
            {
                # top right plate
                'Plate Position': 2,  # as int
                'Plate map file': './17_plate_map.tsv',
                'Project Name': 'Celeste_Adaptation_12986',
                'Project Abbreviation': 'ADAPT',
                'Project Plate': 'Plate_17',  # Plate_#
                'Plate elution volume': 70
            },
            {
                # bottom left plate
                'Plate Position': 3,  # as int
                'Plate map file': './18_plate_map.tsv',
                'Project Name': 'Celeste_Adaptation_12986',
                'Project Abbreviation': 'ADAPT',
                'Project Plate': 'Plate_18',  # Plate_#
                'Plate elution volume': 70
            },
            {
                # bottom right plate
                'Plate Position': 4,  # as int
                'Plate map file': './CHILD_1000_plate_map.tsv',
                'Project Name': 'CHILD_15510',
                'Project Abbreviation': 'CHILD',
                'Project Plate': 'Plate_1000',  # Plate_#
                'Plate elution volume': 70
            },
        ]

        obs_layout = extend_compression_layout_info(
            self.compress_layout, self.studies_info)
        self.assertListEqual(exp_layout, obs_layout)

    def test_extend_compression_layout_info_err(self):
        incomplete_studies_info = self.studies_info.copy()
        incomplete_studies_info.pop(0)

        with self.assertRaisesRegex(
                ValueError, "Project Name 'Celeste_Adaptation_12986' in the "
                            "compression layout is not in the studies info"):
            extend_compression_layout_info(self.compress_layout,
                                           incomplete_studies_info)
