import os
import pandas as pd

from unittest import TestCase, main
from metapool.metapool import parse_sample_sheet
from metapool.prep import (preparations_for_run, remove_qiita_id,
                           get_run_prefix, is_nonempty_gz_file,
                           sample_sheet_to_dataframe, parse_illumina_run_id)


class Tests(TestCase):
    def setUp(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.good_run = os.path.join(data_dir, 'runs',
                                     '191103_D32611_0365_G00DHB5YXX')

        self.ss = os.path.join(self.good_run, 'sample-sheet.csv')

    def test_preparations_for_run(self):
        ss = sample_sheet_to_dataframe(parse_sample_sheet(self.ss))
        obs = preparations_for_run(self.good_run, ss)

        exp = pd.DataFrame()
        pd.testing.assert_frame_equal(obs['D32611_0365_G00DHB5YXX.Baz.1'],
                                      exp)

        exp = pd.DataFrame()
        pd.testing.assert_frame_equal(obs['D32611_0365_G00DHB5YXX.Baz.3'],
                                      exp)

        exp = pd.DataFrame()
        pd.testing.assert_frame_equal(
            obs['D32611_0365_G00DHB5YXX.FooBar_666.3'], exp)

    def test_remove_qiita_id(self):
        obs = remove_qiita_id('project_1')
        self.assertEqual(obs, 'project')

        obs = remove_qiita_id('project_00333333')
        self.assertEqual(obs, 'project')

        obs = remove_qiita_id('project')
        self.assertEqual(obs, 'project')

        obs = remove_qiita_id('project_')
        self.assertEqual(obs, 'project_')

    def test_get_run_prefix(self):
        # project 1
        obs = get_run_prefix(self.good_run, 'Baz', 'sample1', '1')
        self.assertEqual('sample1_S11_L001', obs)

        obs = get_run_prefix(self.good_run, 'Baz', 'sample1', '3')
        self.assertEqual('sample1_S11_L003', obs)

        obs = get_run_prefix(self.good_run, 'Baz', 'sample2', '1')
        self.assertEqual('sample2_S10_L001', obs)

        obs = get_run_prefix(self.good_run, 'Baz', 'sample2', '3')
        self.assertIsNone(obs)

        # project 2
        obs = get_run_prefix(self.good_run, 'FooBar_666', 'sample31', '3')
        self.assertEqual('sample31_S13_L003', obs)

        obs = get_run_prefix(self.good_run, 'FooBar_666', 'sample32', '3')
        self.assertEqual('sample32_S19_L003', obs)

        obs = get_run_prefix(self.good_run, 'FooBar_666', 'sample34', '3')
        self.assertEqual('sample34_S33_L003', obs)

    def test_is_non_empty_gz_file(self):
        non_empty = os.path.join(self.good_run, 'Baz', 'sample2_S10_L001_R1_00'
                                 '1.fastq.gz')
        self.assertTrue(is_nonempty_gz_file(non_empty))
        non_empty = os.path.join(self.good_run, 'Baz', 'sample2_S10_L001_R2_00'
                                 '1.fastq.gz')
        self.assertTrue(is_nonempty_gz_file(non_empty))

        empty = os.path.join(self.good_run, 'Baz/atropos_qc/sample2_S10_L003_R'
                             '1_001.fastq.gz')
        self.assertFalse(is_nonempty_gz_file(empty))
        empty = os.path.join(self.good_run, 'Baz/atropos_qc/sample2_S10_L003_R'
                             '2_001.fastq.gz')
        self.assertFalse(is_nonempty_gz_file(empty))

    def test_sample_sheet_to_dataframe(self):
        ss = parse_sample_sheet(self.ss)
        obs = sample_sheet_to_dataframe(ss)

        exp = pd.DataFrame()
        pd.testing.assert_frame_equal(obs, exp)

    def test_parse_illumina_run_id(self):
        date, rid = parse_illumina_run_id('161004_D00611_0365_AH2HJ5BCXY')
        self.assertEqual(date, '2016-10-04')
        self.assertEqual(rid, 'D00611_0365_AH2HJ5BCXY')

        date, rid = parse_illumina_run_id('160909_K00180_0244_BH7VNKBBXX')
        self.assertEqual(date, '2016-09-09')
        self.assertEqual(rid, 'K00180_0244_BH7VNKBBXX')


if __name__ == "__main__":
    main()
