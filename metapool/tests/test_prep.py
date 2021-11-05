import os
import pandas as pd

from unittest import TestCase, main
from metapool.sample_sheet import KLSampleSheet, sample_sheet_to_dataframe
from metapool.prep import (preparations_for_run, remove_qiita_id,
                           get_run_prefix, is_nonempty_gz_file,
                           get_machine_code, get_model_and_center,
                           parse_illumina_run_id,
                           _check_invalid_names, agp_transform, parse_prep,
                           generate_qiita_prep_file, qiita_scrub_name)


class TestPrep(TestCase):
    def setUp(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.good_run = os.path.join(data_dir, 'runs',
                                     '191103_D32611_0365_G00DHB5YXX')
        self.good_run_new_version = os.path.join(
            data_dir, 'runs', '191104_D32611_0365_G00DHB5YXZ')
        self.OKish_run_new_version = os.path.join(
            data_dir, 'runs', '191104_D32611_0365_OK15HB5YXZ')

        self.ss = os.path.join(self.good_run, 'sample-sheet.csv')
        self.prep = os.path.join(data_dir, 'prep.tsv')

    def _check_run_191103_D32611_0365_G00DHB5YXX(self, obs):
        "Convenience method to check the output of a whole run"

        exp = {('191103_D32611_0365_G00DHB5YXX', 'Baz', '1'),
               ('191103_D32611_0365_G00DHB5YXX', 'Baz', '3'),
               ('191103_D32611_0365_G00DHB5YXX', 'FooBar_666', '3')}
        self.assertEqual(set(obs.keys()), exp)

        columns = ['sample_name', 'experiment_design_description',
                   'library_construction_protocol', 'platform', 'run_center',
                   'run_date', 'run_prefix', 'sequencing_meth', 'center_name',
                   'center_project_name', 'instrument_model', 'runid',
                   'sample_plate', 'sample_well', 'i7_index_id', 'index',
                   'i5_index_id', 'index2', 'lane', 'sample_project',
                   'well_description']

        data = [['importantsample1', 'Eqiiperiment',
                 'Knight Lab Kapa HP', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample1_S11_L003', 'sequencing by synthesis', 'CENTER_NAME',
                 'Baz', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'FooBar_666_p1', 'A3',
                 'iTru7_107_09', 'GCCTTGTT', 'iTru5_01_A', 'AACACCAC', '3',
                 'Baz', 'FooBar_666_p1.sample1.A3'],
                ['importantsample44', 'Eqiiperiment',
                 'Knight Lab Kapa HP', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample44_S14_L003', 'sequencing by synthesis', 'CENTER_NAME',
                 'Baz', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'Baz_p3', 'B99',
                 'iTru7_107_14', 'GTCCTAAG', 'iTru5_01_A', 'CATCTGCT', '3',
                 'Baz', 'Baz_p3.sample44.B99']]
        exp = pd.DataFrame(data=data, columns=columns)
        obs_df = obs[('191103_D32611_0365_G00DHB5YXX', 'Baz', '3')]

        # make sure the columns are in the same order before comparing
        obs_df = obs_df[exp.columns].copy()
        pd.testing.assert_frame_equal(obs_df, exp)

        data = [['importantsample1', 'Eqiiperiment',
                 'Knight Lab Kapa HP', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample1_S11_L001', 'sequencing by synthesis', 'CENTER_NAME',
                 'Baz', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'FooBar_666_p1', 'A1',
                 'iTru7_107_07', 'CCGACTAT', 'iTru5_01_A', 'ACCGACAA', '1',
                 'Baz', 'FooBar_666_p1.sample1.A1'],
                ['importantsample2', 'Eqiiperiment',
                 'Knight Lab Kapa HP', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample2_S10_L001', 'sequencing by synthesis', 'CENTER_NAME',
                 'Baz', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'FooBar_666_p1', 'A2',
                 'iTru7_107_08', 'CCGACTAT', 'iTru5_01_A', 'CTTCGCAA', '1',
                 'Baz', 'FooBar_666_p1.sample2.A2']]
        exp = pd.DataFrame(columns=columns, data=data)
        obs_df = obs[('191103_D32611_0365_G00DHB5YXX', 'Baz', '1')]

        # make sure the columns are in the same order before comparing
        obs_df = obs_df[exp.columns].copy()
        pd.testing.assert_frame_equal(obs_df, exp)

        data = [['importantsample31', 'SomethingWitty',
                 'Knight Lab Kapa HP', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample31_S13_L003', 'sequencing by synthesis',
                 'CENTER_NAME', 'FooBar', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'FooBar_666_p1', 'A5',
                 'iTru7_107_11', 'CAATGTGG', 'iTru5_01_A', 'GGTACGAA', '3',
                 'FooBar_666', 'FooBar_666_p1.sample31.A5'],
                ['importantsample32', 'SomethingWitty',
                 'Knight Lab Kapa HP', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample32_S19_L003', 'sequencing by synthesis', 'CENTER_NAME',
                 'FooBar', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'FooBar_666_p1', 'B6',
                 'iTru7_107_12', 'AAGGCTGA', 'iTru5_01_A', 'CGATCGAT', '3',
                 'FooBar_666', 'FooBar_666_p1.sample32.B6'],
                ['importantsample34', 'SomethingWitty',
                 'Knight Lab Kapa HP', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample34_S33_L003', 'sequencing by synthesis', 'CENTER_NAME',
                 'FooBar', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'FooBar_666_p1', 'B8',
                 'iTru7_107_13', 'TTACCGAG', 'iTru5_01_A', 'AAGACACC', '3',
                 'FooBar_666', 'FooBar_666_p1.sample34.B8']]
        exp = pd.DataFrame(columns=columns, data=data)
        obs_df = obs[('191103_D32611_0365_G00DHB5YXX', 'FooBar_666', '3')]

        # make sure the columns are in the same order before comparing
        obs_df = obs_df[exp.columns].copy()
        pd.testing.assert_frame_equal(obs_df, exp)

    def test_preparations_for_run(self):
        ss = sample_sheet_to_dataframe(KLSampleSheet(self.ss))
        obs = preparations_for_run(self.good_run, ss,
                                   pipeline='atropos-and-bowtie2')
        self._check_run_191103_D32611_0365_G00DHB5YXX(obs)

    def test_preparations_for_run_missing_columns(self):
        # Check that warnings are raised whenever we overwrite the
        # "well_description" column with the "description" column
        ss = sample_sheet_to_dataframe(KLSampleSheet(self.ss))
        ss['description'] = ss['well_description'].copy()
        ss.drop('well_description', axis=1, inplace=True)

        with self.assertWarns(UserWarning) as cm:
            obs = preparations_for_run(self.good_run, ss,
                                       pipeline='atropos-and-bowtie2')

            self.assertEqual(str(cm.warnings[0].message), 'These required '
                             'columns were not found well_description')
            self.assertEqual(str(cm.warnings[1].message), "Using 'description'"
                             " instead of 'well_description' because that "
                             "column isn't present")

        self._check_run_191103_D32611_0365_G00DHB5YXX(obs)

    def test_invalid_sample_names_show_warning(self):
        ss = sample_sheet_to_dataframe(KLSampleSheet(self.ss))

        ss['well_description'] = ss['well_description'].str.replace(
            'importantsample44', 'important-sample44')

        with self.assertWarns(UserWarning) as cm:
            _check_invalid_names(ss['well_description'])
            self.assertEqual(str(cm.warnings[0].message), 'The following '
                             'sample names have invalid characters: '
                             '"important-sample44"')

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
        obs = get_run_prefix(self.good_run, 'Baz', 'sample1', '1',
                             'atropos-and-bowtie2')
        self.assertEqual('sample1_S11_L001', obs)

        obs = get_run_prefix(self.good_run, 'Baz', 'sample1', '3',
                             'atropos-and-bowtie2')
        self.assertEqual('sample1_S11_L003', obs)

        obs = get_run_prefix(self.good_run, 'Baz', 'sample2', '1',
                             'atropos-and-bowtie2')
        self.assertEqual('sample2_S10_L001', obs)

        obs = get_run_prefix(self.good_run, 'Baz', 'sample2', '3',
                             'atropos-and-bowtie2')
        self.assertIsNone(obs)

        # project 2
        obs = get_run_prefix(self.good_run, 'FooBar_666', 'sample31', '3',
                             'atropos-and-bowtie2')
        self.assertEqual('sample31_S13_L003', obs)

        obs = get_run_prefix(self.good_run, 'FooBar_666', 'sample32', '3',
                             'atropos-and-bowtie2')
        self.assertEqual('sample32_S19_L003', obs)

        obs = get_run_prefix(self.good_run, 'FooBar_666', 'sample34', '3',
                             'atropos-and-bowtie2')
        self.assertEqual('sample34_S33_L003', obs)

    def test_get_run_prefix_fastp_minimap(self):
        obs = get_run_prefix(self.good_run_new_version, 'Baz', 'sample1', '1',
                             'fastp-and-minimap2')
        self.assertEqual('sample1_S11_L001', obs)

        obs = get_run_prefix(self.good_run_new_version, 'Baz', 'sample1', '3',
                             'fastp-and-minimap2')
        self.assertEqual('sample1_S11_L003', obs)

        obs = get_run_prefix(self.good_run_new_version, 'Baz', 'sample2', '1',
                             'fastp-and-minimap2')
        self.assertEqual('sample2_S10_L001', obs)

        obs = get_run_prefix(self.good_run_new_version, 'Baz', 'sample44', '3',
                             'fastp-and-minimap2')
        self.assertIsNone(obs)

        obs = get_run_prefix(self.good_run_new_version, 'Baz', 'sample2', '3',
                             'fastp-and-minimap2')
        self.assertIsNone(obs)

        # project 2
        obs = get_run_prefix(self.good_run_new_version, 'FooBar_666',
                             'sample31', '3', 'fastp-and-minimap2')
        self.assertEqual('sample31_S13_L003', obs)

        obs = get_run_prefix(self.good_run_new_version, 'FooBar_666',
                             'sample32', '3', 'fastp-and-minimap2')
        self.assertEqual('sample32_S19_L003', obs)

        obs = get_run_prefix(self.good_run_new_version, 'FooBar_666',
                             'sample34', '3', 'fastp-and-minimap2')
        self.assertIsNone(obs)

    def test_get_run_prefix_more_than_forward_and_reverse(self):
        message = (r'There are 3 matches for sample "sample31" in lane 3\. '
                   r'Only two matches are allowed \(forward and reverse\): '
                   '(.*)metapool/tests/data/runs/191104_D32611_0365_OK15HB5YXZ'
                   r'/FooBar_666/filtered_sequences/sample31_S13_L003_R1\.'
                   r'filtered\.fastq\.gz, '
                   '(.*)metapool/tests/data/runs/191104_D32611_0365_OK15HB5YXZ'
                   r'/FooBar_666/filtered_sequences/sample31_S13_L003_R2\.'
                   r'filtered\.fastq\.gz, '
                   '(.*)metapool/tests/data/runs/191104_D32611_0365_OK15HB5YXZ'
                   r'/FooBar_666/filtered_sequences/sample31_S14_L003_R1\.'
                   r'filtered\.fastq\.gz')
        # project 2
        with self.assertWarnsRegex(Warning, message):
            obs = get_run_prefix(self.OKish_run_new_version, 'FooBar_666',
                                 'sample31', '3', 'fastp-and-minimap2')
            self.assertIsNone(obs)

        obs = get_run_prefix(self.OKish_run_new_version, 'FooBar_666',
                             'sample32', '3', 'fastp-and-minimap2')
        self.assertEqual('sample32_S19_L003', obs)

        obs = get_run_prefix(self.OKish_run_new_version, 'FooBar_666',
                             'sample34', '3', 'fastp-and-minimap2')
        self.assertIsNone(obs)

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

    def test_parse_illumina_run_id(self):
        date, rid = parse_illumina_run_id('161004_D00611_0365_AH2HJ5BCXY')
        self.assertEqual(date, '2016-10-04')
        self.assertEqual(rid, 'D00611_0365_AH2HJ5BCXY')

        date, rid = parse_illumina_run_id('160909_K00180_0244_BH7VNKBBXX')
        self.assertEqual(date, '2016-09-09')
        self.assertEqual(rid, 'K00180_0244_BH7VNKBBXX')

    def test_machine_code(self):
        obs = get_machine_code('K00180')
        self.assertEqual(obs, 'K')

        obs = get_machine_code('D00611')
        self.assertEqual(obs, 'D')

        obs = get_machine_code('MN01225')
        self.assertEqual(obs, 'MN')

        with self.assertRaisesRegex(ValueError,
                                    'Cannot find a machine code. This '
                                    'instrument model is malformed 8675309. '
                                    'The machine code is a one or two '
                                    'character prefix.'):
            get_machine_code('8675309')

    def test_get_model_and_center(self):
        obs = get_model_and_center('D32611_0365_G00DHB5YXX')
        self.assertEqual(obs, ('Illumina HiSeq 2500', 'UCSDMI'))

        obs = get_model_and_center('A86753_0365_G00DHB5YXX')
        self.assertEqual(obs, ('Illumina NovaSeq 6000', 'UCSDMI'))

        obs = get_model_and_center('A00953_0032_AHWMGJDDXX')
        self.assertEqual(obs, ('Illumina NovaSeq 6000', 'IGM'))

        obs = get_model_and_center('A00169_8131_AHKXYNDHXX')
        self.assertEqual(obs, ('Illumina NovaSeq 6000', 'LJI'))

        obs = get_model_and_center('M05314_0255_000000000-J46T9')
        self.assertEqual(obs, ('Illumina MiSeq', 'KLM'))

        obs = get_model_and_center('K00180_0957_AHCYKKBBXY')
        self.assertEqual(obs, ('Illumina HiSeq 4000', 'IGM'))

        obs = get_model_and_center('D00611_0712_BH37W2BCX3_RKL0040')
        self.assertEqual(obs, ('Illumina HiSeq 2500', 'IGM'))

        obs = get_model_and_center('MN01225_0002_A000H2W3FY')
        self.assertEqual(obs, ('Illumina MiniSeq', 'CMI'))

    def test_agp_transform(self):
        columns = ['sample_name', 'experiment_design_description',
                   'library_construction_protocol', 'platform', 'run_center',
                   'run_date', 'run_prefix', 'sequencing_meth', 'center_name',
                   'center_project_name', 'instrument_model', 'runid',
                   'sample_plate', 'sample_well', 'i7_index_id', 'index',
                   'i5_index_id', 'index2', 'lane', 'sample_project',
                   'well_description']

        data = [['importantsample1', 'EXPERIMENT_DESC',
                 'LIBRARY_PROTOCOL', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample1_S11_L003', 'sequencing by synthesis', 'CENTER_NAME',
                 'Baz', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'FooBar_666_p1', 'A3',
                 'iTru7_107_09', 'GCCTTGTT', 'iTru5_01_A', 'AACACCAC', '3',
                 'Baz', 'FooBar_666_p1.sample1.A3'],
                ['importantsample44', 'EXPERIMENT_DESC',
                 'LIBRARY_PROTOCOL', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample44_S14_L003', 'sequencing by synthesis', 'CENTER_NAME',
                 'Baz', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'Baz_p3', 'B99',
                 'iTru7_107_14', 'GTCCTAAG', 'iTru5_01_A', 'CATCTGCT', '3',
                 'Baz', 'Baz_p3.sample44.B99']]
        obs = pd.DataFrame(data=data, columns=columns)
        exp = obs.copy()
        exp['center_name'] = 'UCSDMI'
        exp['library_construction_protocol'] = 'Knight Lab KHP'
        exp['experiment_design_description'] = (
            'samples of skin, saliva and feces and other samples from the AGP')

        pd.testing.assert_frame_equal(agp_transform(obs, '10317'), exp)

    def test_agp_transform_zfill(self):
        columns = ['sample_name', 'experiment_design_description',
                   'library_construction_protocol', 'platform', 'run_center',
                   'run_date', 'run_prefix', 'sequencing_meth', 'center_name',
                   'center_project_name', 'instrument_model', 'runid',
                   'sample_plate', 'sample_well', 'i7_index_id', 'index',
                   'i5_index_id', 'index2', 'lane', 'sample_project',
                   'well_description']

        # the first sample name should be padded with 2 zeros
        # the second number should be padded with 1 zero
        # the third should be ignored
        data = [['8675309', 'EXPERIMENT_DESC',
                 'LIBRARY_PROTOCOL', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample1_S11_L003', 'sequencing by synthesis', 'CENTER_NAME',
                 'Baz', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'FooBar_666_p1', 'A3',
                 'iTru7_107_09', 'GCCTTGTT', 'iTru5_01_A', 'AACACCAC', '3',
                 'Baz', 'FooBar_666_p1.sample1.A3'],
                ['867530.9', 'EXPERIMENT_DESC',
                 'LIBRARY_PROTOCOL', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample44_S14_L003', 'sequencing by synthesis', 'CENTER_NAME',
                 'Baz', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'Baz_p3', 'B99',
                 'iTru7_107_14', 'GTCCTAAG', 'iTru5_01_A', 'CATCTGCT', '3',
                 'Baz', 'Baz_p3.sample44.B99'],
                ['notanumber', 'EXPERIMENT_DESC',
                 'LIBRARY_PROTOCOL', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample44_S14_L003', 'sequencing by synthesis', 'CENTER_NAME',
                 'Baz', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'Baz_p3', 'B99',
                 'iTru7_107_14', 'GTCCTAAG', 'iTru5_01_A', 'CATCTGCT', '3',
                 'Baz', 'Baz_p3.sample44.B99']
                ]
        obs = pd.DataFrame(data=data, columns=columns)
        exp = obs.copy()
        exp['center_name'] = 'UCSDMI'
        exp['library_construction_protocol'] = 'Knight Lab KHP'
        exp['experiment_design_description'] = (
            'samples of skin, saliva and feces and other samples from the AGP')
        exp.loc[0, 'sample_name'] = '008675309'
        exp.loc[1, 'sample_name'] = '0867530.9'

        pd.testing.assert_frame_equal(agp_transform(obs, '10317'), exp)

    def test_agp_transform_no_agp(self):
        columns = ['sample_name', 'experiment_design_description',
                   'library_construction_protocol', 'platform', 'run_center',
                   'run_date', 'run_prefix', 'sequencing_meth', 'center_name',
                   'center_project_name', 'instrument_model', 'runid',
                   'sample_plate', 'sample_well', 'i7_index_id', 'index',
                   'i5_index_id', 'index2', 'lane', 'sample_project',
                   'well_description']

        data = [['importantsample1', 'EXPERIMENT_DESC',
                 'LIBRARY_PROTOCOL', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample1_S11_L003', 'sequencing by synthesis', 'CENTER_NAME',
                 'Baz', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'FooBar_666_p1', 'A3',
                 'iTru7_107_09', 'GCCTTGTT', 'iTru5_01_A', 'AACACCAC', '3',
                 'Baz', 'FooBar_666_p1.sample1.A3'],
                ['importantsample44', 'EXPERIMENT_DESC',
                 'LIBRARY_PROTOCOL', 'Illumina', 'UCSDMI', '2019-11-03',
                 'sample44_S14_L003', 'sequencing by synthesis', 'CENTER_NAME',
                 'Baz', 'Illumina HiSeq 2500',
                 '191103_D32611_0365_G00DHB5YXX', 'Baz_p3', 'B99',
                 'iTru7_107_14', 'GTCCTAAG', 'iTru5_01_A', 'CATCTGCT', '3',
                 'Baz', 'Baz_p3.sample44.B99']]
        obs = pd.DataFrame(data=data, columns=columns)
        exp = obs.copy()

        # there shouldn't be any changes
        pd.testing.assert_frame_equal(agp_transform(obs, '666'), exp)

    def test_parse_prep(self):
        columns = [
            'barcode', 'primer', 'project_name', 'well_id',
            'primer_plate', 'plating', 'extractionkit_lot',
            'extraction_robot', 'tm1000_8_tool', 'primer_date',
            'mastermix_lot', 'water_lot', 'processing_robot',
            'sample_plate', 'linker', 'orig_name', 'well_description',
            'pcr_primers', 'center_name', 'run_center', 'platform',
            'target_subfragment', 'target_gene', 'sequencing_meth',
            'library_construction_protocol']

        data = [
            ['AGCCTTCGTCGC', 'GTGYCAGCMGCCGCGGTAA', 'THDMI_10317', 'A1', '1',
             'SF', '166032128', 'Carmen_HOWE_KF3', '109379Z', '2021-08-17',
             '978215', 'RNBJ0628', 'Echo550', 'THDMI_UK_Plate_2', 'GT',
             'X00180471', 'THDMI_UK_Plate_2.X00180471.A1',
             'FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT', 'UCSDMI',
             'UCSDMI', 'Illumina', 'V4', '16S rRNA', 'Sequencing by synthesis',
             'Illumina EMP protocol 515fbc, 806r amplification of 16S rRNA V4'
             ],
            ['CGTATAAATGCG', 'GTGYCAGCMGCCGCGGTAA', 'THDMI_10317', 'C1', '1',
             'SF', '166032128', 'Carmen_HOWE_KF3', '109379Z', '2021-08-17',
             '978215', 'RNBJ0628', 'Echo550', 'THDMI_UK_Plate_2', 'GT',
             'X00180199', 'THDMI_UK_Plate_2.X00180199.C1',
             'FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT', 'UCSDMI',
             'UCSDMI', 'Illumina', 'V4', '16S rRNA',
             'Sequencing by synthesis',
             'Illumina EMP protocol 515fbc, 806r amplification of 16S rRNA V4']
        ]
        index = pd.Index(['X00180471', 'X00180199'], name='sample_name')
        exp = pd.DataFrame(data=data, columns=columns, index=index)

        obs = parse_prep(self.prep)

        pd.testing.assert_frame_equal(obs, exp)

    def test_generate_qiita_prep_file(self):
        columns1 = ['Sample', 'Row', 'Col', 'Blank', 'Project Plate',
                    'Project Name', 'Compressed Plate Name', 'Well',
                    'Plate Position', 'Primer Plate #', 'Sample Plate',
                    'Project_Name', 'Plating', 'Extraction Kit Lot',
                    'Extraction Robot', 'TM1000 8 Tool', 'Primer Date',
                    'MasterMix Lot', 'Water Lot', 'Processing Robot',
                    'Original Name', 'Plate', 'EMP Primer Plate Well', 'Name',
                    'Illumina 5prime Adapter', 'Golay Barcode',
                    'Forward Primer Pad', 'Forward Primer Linker',
                    '515FB Forward Primer (Parada)', 'Primer For PCR',
                    'sample sheet Sample_ID']
        columns2 = ['Sample', 'Row', 'Col', 'Blank', 'Project Plate',
                    'Project Name', 'Compressed Plate Name', 'Well',
                    'Plate Position', 'Primer Plate #', 'Sample Plate',
                    'Project_Name', 'Plating', 'Extraction Kit Lot',
                    'Extraction Robot', 'TM1000 8 Tool', 'Primer Date',
                    'MasterMix Lot', 'Water Lot', 'Processing Robot',
                    'Original Name', 'Plate', 'EMP Primer Plate Well', 'Name',
                    'Reverse complement of 3prime Illumina Adapter',
                    'Golay Barcode',
                    'Reverse Primer Pad', 'Reverse Primer Linker',
                    'Reverse primer (EukBr)', 'Primer For PCR',
                    'sample sheet Sample_ID']
        columns3 = ['Sample', 'Row', 'Col', 'Blank', 'Project Plate',
                    'Project Name', 'Compressed Plate Name', 'Well',
                    'Plate Position', 'Primer Plate #', 'Sample Plate',
                    'Project_Name', 'Plating', 'Extraction Kit Lot',
                    'Extraction Robot', 'TM1000 8 Tool', 'Primer Date',
                    'MasterMix Lot', 'Water Lot', 'Processing Robot',
                    'Original Name', 'Plate', 'EMP Primer Plate Well', 'Name',
                    'Reverse complement of 3prime Illumina Adapter',
                    'Golay Barcode',
                    'Reverse Primer Pad', 'Reverse Primer Linker',
                    'ITS2 Reverse Primer', 'Primer For PCR',
                    'sample sheet Sample_ID']
        data1 = [['X00180471', 'A', 1, False, 'THDMI_10317_PUK2',
                  'THDMI_10317', 'THDMI_10317_UK2-US6', 'A1',
                  '1', '1', 'THDMI_UK_Plate_2', 'THDMI UK',
                  'SF', '166032128', 'Carmen_HOWE_KF3',
                  '109379Z', '2021-08-17', '978215', 'RNBJ0628',
                  'Echo550', '', '1', 'A1', '515rcbc0',
                  'AATGATACGGCGACCACCGAGATCTACACGCT', 'AGCCTTCGTCGC',
                  'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
                  ('AATGATACGGCGACCACCGAGATCTACACGCTAGCCTTCGTCGCT'
                   'ATGGTAATTGTGTGYCAGCMGCCGCGGTAA'),
                  'X00180471'],
                 ['X00180199', 'C', 1, False, 'THDMI_10317_PUK2',
                  'THDMI_10317', 'THDMI_10317_UK2-US6', 'C1',
                  '1', '1', 'THDMI_UK_Plate_2', 'THDMI UK',
                  'SF', '166032128', 'Carmen_HOWE_KF3',
                  '109379Z', '2021-08-17', '978215', 'RNBJ0628',
                  'Echo550', '', '1', 'B1', '515rcbc12',
                  'AATGATACGGCGACCACCGAGATCTACACGCT',
                  'CGTATAAATGCG', 'TATGGTAATT', 'GT',
                  'GTGYCAGCMGCCGCGGTAA',
                  ('AATGATACGGCGACCACCGAGATCTACACGCTCGTATAAATGCG'
                   'TATGGTAATTGTGTGYCAGCMGCCGCGGTAA'),
                  'X00180199']]
        data2 = [['X00180471', 'A', 1, False, 'THDMI_10317_PUK2',
                  'THDMI_10317',
                  'THDMI_10317_UK2-US6', 'A1', '1', '1', 'THDMI_UK_Plate_2',
                  'THDMI UK', 'SF', '166032128', 'Carmen_HOWE_KF3', '109379Z',
                  '2021-08-17', '978215', 'RNBJ0628', 'Echo550', '', '1', 'A1',
                  'EukBr_Hiseq_0017', 'CAAGCAGAAGACGGCATACGAGAT',
                  'ACGAGACTGATT',
                  'AGTCAGTCAG', 'CA', 'TGATCCTTCTGCAGGTTCACCTAC',
                  ('CAAGCAGAAGACGGCATACGAGATACGAGACTGATTAGTCAGTCAGCAT'
                   'GATCCTTCTGCAGGTTCACCTAC'),
                  'X00180471'],
                 ['X00180199', 'C', 1, False, 'THDMI_10317_PUK2',
                  'THDMI_10317',
                  'THDMI_10317_UK2-US6', 'C1', '1', '1', 'THDMI_UK_Plate_2',
                  'THDMI UK', 'SF', '166032128', 'Carmen_HOWE_KF3', '109379Z',
                  '2021-08-17', '978215', 'RNBJ0628', 'Echo550', '', '1', 'B1',
                  'EukBr_Hiseq_0029', 'CAAGCAGAAGACGGCATACGAGAT',
                  'GAATACCAAGTC',
                  'AGTCAGTCAG', 'CA', 'TGATCCTTCTGCAGGTTCACCTAC',
                  ('CAAGCAGAAGACGGCATACGAGATGAATACCAAGTCAGTCAGTCAGCAT'
                   'GATCCTTCTGCAGGTTCACCTAC'),
                  'X00180199']]
        data3 = [['X00180471', 'A', 1, False, 'THDMI_10317_PUK2',
                  'THDMI_10317',
                  'THDMI_10317_UK2-US6', 'A1', '1', '1', 'THDMI_UK_Plate_2',
                  'THDMI UK', 'SF', '166032128', 'Carmen_HOWE_KF3', '109379Z',
                  '2021-08-17', '978215', 'RNBJ0628', 'Echo550', '', '1',
                  'A1',
                  'kabir_ITS2rcbc0', 'CAAGCAGAAGACGGCATACGAGAT',
                  'TCCCTTGTCTCC',
                  '', 'CG', 'GCTGCGTTCTTCATCGATGC',
                  ('CAAGCAGAAGACGGCATACGAGATTCCCTTGTCTCCC'
                   'GGCTGCGTTCTTCATCGATGC'),
                  'X00180471'],
                 ['X00180199', 'C', 1, False, 'THDMI_10317_PUK2',
                  'THDMI_10317',
                  'THDMI_10317_UK2-US6', 'C1', '1', '1', 'THDMI_UK_Plate_2',
                  'THDMI UK', 'SF', '166032128', 'Carmen_HOWE_KF3', '109379Z',
                  '2021-08-17', '978215', 'RNBJ0628', 'Echo550', '',
                  '1', 'B1',
                  'kabir_ITS2rcbc12', 'CAAGCAGAAGACGGCATACGAGAT',
                  'TGCATACACTGG',
                  '', 'CG', 'GCTGCGTTCTTCATCGATGC',
                  ('CAAGCAGAAGACGGCATACGAGATTGCATACACTGGC'
                   'GGCTGCGTTCTTCATCGATGC'),
                  'X00180199']]

        self.platedf1 = pd.DataFrame(columns=columns1, data=data1)
        self.platedf2 = pd.DataFrame(columns=columns2, data=data2)
        self.platedf3 = pd.DataFrame(columns=columns3, data=data3)

        obs1 = generate_qiita_prep_file(self.platedf1, '16S')
        obs2 = generate_qiita_prep_file(self.platedf2, '18S')
        obs3 = generate_qiita_prep_file(self.platedf3, 'ITS')

        common_col = ['sample_name', 'barcode', 'primer', 'project_name',
                      'well_id', 'primer_plate', 'plating',
                      'extractionkit_lot', 'extraction_robot',
                      'tm1000_8_tool', 'primer_date', 'mastermix_lot',
                      'water_lot', 'processing_robot', 'sample_plate',
                      'linker', 'orig_name', 'well_description',
                      'center_name', 'run_center', 'platform',
                      'sequencing_meth', 'pcr_primers',
                      'target_subfragment', 'target_gene',
                      'library_construction_protocol']

        exp1 = pd.DataFrame(columns=common_col,
                            data=[['X00180471', 'AGCCTTCGTCGC',
                                   'GTGYCAGCMGCCGCGGTAA', 'THDMI_10317',
                                   'A1', '1', 'SF', '166032128',
                                   'Carmen_HOWE_KF3', '109379Z',
                                   '2021-08-17', '978215', 'RNBJ0628',
                                   'Echo550', 'THDMI_UK_Plate_2', 'GT',
                                   'X00180471',
                                   'THDMI_UK_Plate_2.X00180471.A1',
                                   'UCSDMI', 'UCSDMI', 'Illumina',
                                   'Sequencing by synthesis',
                                   ('FWD:GTGYCAGCMGCCGCGGTAA; '
                                    'REV:GGACTACNVGGGTWTCTAAT'), 'V4',
                                   '16S rRNA',
                                   ('Illumina EMP protocol 515fbc, '
                                    '806r amplification of 16S rRNA V4')],
                                  ['X00180199', 'CGTATAAATGCG',
                                   'GTGYCAGCMGCCGCGGTAA', 'THDMI_10317',
                                   'C1', '1', 'SF', '166032128',
                                   'Carmen_HOWE_KF3', '109379Z',
                                   '2021-08-17', '978215', 'RNBJ0628',
                                   'Echo550', 'THDMI_UK_Plate_2', 'GT',
                                   'X00180199',
                                   'THDMI_UK_Plate_2.X00180199.C1',
                                   'UCSDMI', 'UCSDMI', 'Illumina',
                                   'Sequencing by synthesis',
                                   ('FWD:GTGYCAGCMGCCGCGGTAA; '
                                    'REV:GGACTACNVGGGTWTCTAAT'), 'V4',
                                   '16S rRNA',
                                   ('Illumina EMP protocol 515fbc, '
                                    '806r amplification of 16S rRNA V4')]])
        exp2 = pd.DataFrame(columns=common_col,
                            data=[['X00180471', 'ACGAGACTGATT',
                                   'CAAGCAGAAGACGGCATACGAGAT', 'THDMI_10317',
                                   'A1', '1', 'SF', '166032128',
                                   'Carmen_HOWE_KF3', '109379Z',
                                   '2021-08-17', '978215', 'RNBJ0628',
                                   'Echo550', 'THDMI_UK_Plate_2', 'CA',
                                   'X00180471',
                                   'THDMI_UK_Plate_2.X00180471.A1',
                                   'UCSDMI', 'UCSDMI', 'Illumina',
                                   'Sequencing by synthesis',
                                   ('FWD:GTACACACCGCCCGTC; '
                                    'REV:TGATCCTTCTGCAGGTTCACCTAC'), 'V9',
                                   '18S rRNA',
                                   'Illumina EMP 18S rRNA 1391f EukBr'],
                                  ['X00180199', 'GAATACCAAGTC',
                                   'CAAGCAGAAGACGGCATACGAGAT', 'THDMI_10317',
                                   'C1', '1', 'SF', '166032128',
                                   'Carmen_HOWE_KF3', '109379Z',
                                   '2021-08-17', '978215', 'RNBJ0628',
                                   'Echo550', 'THDMI_UK_Plate_2', 'CA',
                                   'X00180199',
                                   'THDMI_UK_Plate_2.X00180199.C1',
                                   'UCSDMI', 'UCSDMI', 'Illumina',
                                   'Sequencing by synthesis',
                                   ('FWD:GTACACACCGCCCGTC; '
                                    'REV:TGATCCTTCTGCAGGTTCACCTAC'), 'V9',
                                   '18S rRNA',
                                   'Illumina EMP 18S rRNA 1391f EukBr']])
        exp3 = pd.DataFrame(columns=common_col,
                            data=[['X00180471', 'TCCCTTGTCTCC',
                                   'CAAGCAGAAGACGGCATACGAGAT',
                                   'THDMI_10317',
                                   'A1', '1', 'SF', '166032128',
                                   'Carmen_HOWE_KF3', '109379Z',
                                   '2021-08-17', '978215', 'RNBJ0628',
                                   'Echo550', 'THDMI_UK_Plate_2', 'CG',
                                   'X00180471',
                                   'THDMI_UK_Plate_2.X00180471.A1',
                                   'UCSDMI', 'UCSDMI', 'Illumina',
                                   'Sequencing by synthesis',
                                   ('FWD:CTTGGTCATTTAGAGGAAGTAA; '
                                    'REV:GCTGCGTTCTTCATCGATGC'), 'ITS_1_2',
                                   'ITS',
                                   ('Illumina  EMP protocol amplification '
                                    'of ITS1fbc, ITS2r')],
                                  ['X00180199', 'TGCATACACTGG',
                                   'CAAGCAGAAGACGGCATACGAGAT',
                                   'THDMI_10317',
                                   'C1', '1', 'SF', '166032128',
                                   'Carmen_HOWE_KF3', '109379Z',
                                   '2021-08-17', '978215', 'RNBJ0628',
                                   'Echo550', 'THDMI_UK_Plate_2', 'CG',
                                   'X00180199',
                                   'THDMI_UK_Plate_2.X00180199.C1',
                                   'UCSDMI', 'UCSDMI', 'Illumina',
                                   'Sequencing by synthesis',
                                   ('FWD:CTTGGTCATTTAGAGGAAGTAA; '
                                    'REV:GCTGCGTTCTTCATCGATGC'), 'ITS_1_2',
                                   'ITS',
                                   ('Illumina  EMP protocol amplification '
                                    'of ITS1fbc, ITS2r')]])

        pd.testing.assert_frame_equal(obs1, exp1)
        pd.testing.assert_frame_equal(obs2, exp2)
        pd.testing.assert_frame_equal(obs3, exp3)

    def test_qiita_scrub_name(self):
        self.assertEqual(qiita_scrub_name('its my life'), 'its.my.life')
        self.assertEqual(qiita_scrub_name('7-11.10()'), '7-11.10.')
        self.assertEqual(qiita_scrub_name('{}/<>'), '.')
        self.assertEqual(qiita_scrub_name('th!s.has.happened'),
                         'th.s.has.happened')


if __name__ == "__main__":
    main()
