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

        columns = ['lane', 'sample_name', 'sample_plate', 'sample_well',
                   'i7_index_id', 'index', 'i5_index_id', 'index2',
                   'sample_project', 'well_description']
        index = ['sample1', 'sample2', 'sample1', 'sample2', 'sample31',
                 'sample32', 'sample34', 'sample44']
        data = [['1', 'sample1', 'FooBar_666_p1', 'A1', 'iTru7_107_07',
                 'CCGACTAT', 'iTru5_01_A', 'ACCGACAA', 'FooBar_666',
                 'important-sample1'],
                ['1', 'sample2', 'FooBar_666_p1', 'A2', 'iTru7_107_08',
                 'CCGACTAT', 'iTru5_01_A', 'CTTCGCAA', 'FooBar_666',
                 'important-sample2'],
                ['3', 'sample1', 'FooBar_666_p1', 'A3', 'iTru7_107_09',
                 'GCCTTGTT', 'iTru5_01_A', 'AACACCAC', 'FooBar_666',
                 'important-sample1'],
                ['3', 'sample2', 'FooBar_666_p1', 'A4', 'iTru7_107_10',
                 'AACTTGCC', 'iTru5_01_A', 'CGTATCTC', 'FooBar_666',
                 'important-sample2'],
                ['3', 'sample31', 'Baz_p3', 'A5', 'iTru7_107_11', 'CAATGTGG',
                 'iTru5_01_A', 'GGTACGAA', 'Baz', 'important-sample31'],
                ['3', 'sample32', 'Baz_p3', 'B6', 'iTru7_107_12', 'AAGGCTGA',
                 'iTru5_01_A', 'CGATCGAT', 'Baz', 'important-sample32'],
                ['3', 'sample34', 'Baz_p3', 'B8', 'iTru7_107_13', 'TTACCGAG',
                 'iTru5_01_A', 'AAGACACC', 'Baz', 'important-sample34'],
                ['3', 'sample44', 'Baz_p3', 'B99', 'iTru7_107_14', 'GTCCTAAG',
                 'iTru5_01_A', 'CATCTGCT', 'Baz', 'important-sample44']]

        exp = pd.DataFrame(index=index, data=data, columns=columns)
        exp.index.name = 'sample_id'
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
