from os.path import join, dirname

import pandas
import pandas as pd

from unittest import TestCase, main
from metapool.sample_sheet import (MetagenomicSampleSheetv90,
                                   MetagenomicSampleSheetv100,
                                   sample_sheet_to_dataframe)
from metapool.prep import (preparations_for_run, remove_qiita_id,
                           get_run_prefix, is_nonempty_gz_file,
                           get_machine_code, get_model_and_center,
                           parse_illumina_run_id,
                           _check_invalid_names, agp_transform, parse_prep,
                           generate_qiita_prep_file, qiita_scrub_name,
                           preparations_for_run_mapping_file, demux_pre_prep,
                           pre_prep_needs_demuxing)


class TestPrep(TestCase):
    def setUp(self):
        data_dir = join(dirname(__file__), 'data')
        self.good_run = join(data_dir, 'runs', '191103_D32611_0365_G00DHB5YXX')
        self.good_run_new_version = join(
            data_dir, 'runs', '191104_D32611_0365_G00DHB5YXZ')
        self.OKish_run_new_version = join(
            data_dir, 'runs', '191104_D32611_0365_OK15HB5YXZ')

        self.amplicon_run = join(data_dir, 'runs',
                                 '230207_M05314_0346_000000000-KVMGL')

        self.ss = join(self.good_run, 'sample-sheet.csv')
        self.prep = join(data_dir, 'prep.tsv')
        self.legacy_prep = join(data_dir, "legacy_prep.tsv")
        self.mf = join(self.amplicon_run, 'sample_mapping_file.tsv')

    def _check_run_191103_D32611_0365_G00DHB5YXX(self, obs):
        "Convenience method to check the output of a whole run"

        exp = {('191103_D32611_0365_G00DHB5YXX', 'Baz_12345', '1'),
               ('191103_D32611_0365_G00DHB5YXX', 'Baz_12345', '3'),
               ('191103_D32611_0365_G00DHB5YXX', 'FooBar_666', '3')}
        self.assertEqual(set(obs.keys()), exp)

        columns = ['sample_name', 'experiment_design_description',
                   'library_construction_protocol', 'platform', 'run_center',
                   'run_date', 'run_prefix', 'sequencing_meth', 'center_name',
                   'center_project_name', 'instrument_model', 'runid',
                   'sample_plate', 'well_id_384', 'i7_index_id', 'index',
                   'i5_index_id', 'index2', 'lane', 'sample_project',
                   'well_description']

        data = [['sample.1', 'Eqiiperiment', 'Knight Lab Kapa HP',
                 'Illumina', 'UCSDMI', '2019-11-03', 'sample_1_S11_L003',
                 'sequencing by synthesis', 'UCSD', 'Baz',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'FooBar_666_p1', 'A3', 'iTru7_107_09', 'GCCTTGTT',
                 'iTru5_01_A', 'AACACCAC', '3', 'Baz_12345',
                 'FooBar_666_p1.sample.1.A3'],
                ['sample.44', 'Eqiiperiment', 'Knight Lab Kapa HP',
                 'Illumina', 'UCSDMI', '2019-11-03', 'sample_44_S14_L003',
                 'sequencing by synthesis', 'UCSD', 'Baz',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'Baz_12345_p3', 'B99', 'iTru7_107_14', 'GTCCTAAG',
                 'iTru5_01_A', 'CATCTGCT', '3', 'Baz_12345',
                 'Baz_12345_p3.sample.44.B99']]

        exp = pd.DataFrame(data=data, columns=columns)

        # obs is a dictionary of dataframes, where the keys are tuples of
        # strings, rather than ordinary strings. We'll use the dataframe
        # associated w/('191103_D32611_0365_G00DHB5YXX', 'Baz_12345', '3') for
        # our unittest.
        obs_df = obs[('191103_D32611_0365_G00DHB5YXX', 'Baz_12345', '3')]

        # make sure the columns are in the same order before comparing
        obs_df = obs_df[exp.columns].copy()

        pd.testing.assert_frame_equal(obs_df, exp)

        data = [['sample.1', 'Eqiiperiment', 'Knight Lab Kapa HP',
                 'Illumina', 'UCSDMI', '2019-11-03', 'sample_1_S11_L001',
                 'sequencing by synthesis', 'UCSD', 'Baz',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'FooBar_666_p1', 'A1', 'iTru7_107_07', 'CCGACTAT',
                 'iTru5_01_A', 'ACCGACAA', '1', 'Baz_12345',
                 'FooBar_666_p1.sample.1.A1'],
                ['sample.2', 'Eqiiperiment', 'Knight Lab Kapa HP',
                 'Illumina', 'UCSDMI', '2019-11-03', 'sample_2_S10_L001',
                 'sequencing by synthesis', 'UCSD', 'Baz',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'FooBar_666_p1', 'A2', 'iTru7_107_08', 'CCGACTAT',
                 'iTru5_01_A', 'CTTCGCAA', '1', 'Baz_12345',
                 'FooBar_666_p1.sample.2.A2']]
        exp = pd.DataFrame(columns=columns, data=data)
        obs_df = obs[('191103_D32611_0365_G00DHB5YXX', 'Baz_12345', '1')]

        # make sure the columns are in the same order before comparing
        obs_df = obs_df[exp.columns].copy()
        pd.testing.assert_frame_equal(obs_df, exp)

        data = [['sample.31', 'SomethingWitty', 'Knight Lab Kapa HP',
                 'Illumina', 'UCSDMI', '2019-11-03', 'sample_31_S13_L003',
                 'sequencing by synthesis', 'UCSD', 'FooBar',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'FooBar_666_p1', 'A5', 'iTru7_107_11', 'CAATGTGG',
                 'iTru5_01_A', 'GGTACGAA', '3', 'FooBar_666',
                 'FooBar_666_p1.sample.31.A5'],
                ['sample.32', 'SomethingWitty', 'Knight Lab Kapa HP',
                 'Illumina', 'UCSDMI', '2019-11-03', 'sample_32_S19_L003',
                 'sequencing by synthesis', 'UCSD', 'FooBar',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'FooBar_666_p1', 'B6', 'iTru7_107_12', 'AAGGCTGA',
                 'iTru5_01_A', 'CGATCGAT', '3', 'FooBar_666',
                 'FooBar_666_p1.sample.32.B6'],
                ['sample.34', 'SomethingWitty', 'Knight Lab Kapa HP',
                 'Illumina', 'UCSDMI', '2019-11-03', 'sample_34_S33_L003',
                 'sequencing by synthesis', 'UCSD', 'FooBar',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'FooBar_666_p1', 'B8', 'iTru7_107_13', 'TTACCGAG',
                 'iTru5_01_A', 'AAGACACC', '3', 'FooBar_666',
                 'FooBar_666_p1.sample.34.B8']]
        exp = pd.DataFrame(columns=columns, data=data)
        obs_df = obs[('191103_D32611_0365_G00DHB5YXX', 'FooBar_666', '3')]

        # make sure the columns are in the same order before comparing
        obs_df = obs_df[exp.columns].copy()
        pd.testing.assert_frame_equal(obs_df, exp)

    def _check_run_230207_M05314_0346_000000000_KVMGL(self, obs):
        # confirm correct keys are present for the output prep
        exp_keys = {('230207_M05314_0346_000000000-KVMGL',
                     'ABTX_20230208_ABTX_11052', '1')}

        self.assertEqual(set(obs.keys()), exp_keys)

        # confirm the observed prep-info output contains the expected
        # columns.
        obs_df = obs[('230207_M05314_0346_000000000-KVMGL',
                      'ABTX_20230208_ABTX_11052', '1')]

        exp_columns = ['sample_name', 'barcode', 'center_name', 'platform',
                       'center_project_name', 'experiment_design_description',
                       'instrument_model', 'lane', 'run_center', 'run_date',
                       'library_construction_protocol', 'run_prefix', 'runid',
                       'sample_plate', 'sequencing_meth', 'linker', 'primer',
                       'target_gene', 'pcr_primers', 'primer_plate',
                       'processing_robot', 'well_description',
                       'extraction_robot', 'tm300_8_tool', 'water_lot',
                       'extractionkit_lot', 'target_subfragment',
                       'well_id_384', 'project_name', 'tm1000_8_tool',
                       'orig_name', 'plating', 'primer_date', 'mastermix_lot',
                       'tm50_8_tool', 'well_id_96', 'tm10_8_tool']

        self.assertEqual(set(exp_columns), set(obs_df.columns))

        exp_data = [
            ["sample.1", "AGCCTTCGTCGC", "UCSDMI", "Illumina",
             "SOME_CENTER_PROJECT_NAME",
             "This is a description of the experiment design.",
             "SOME_INSTRUMENT_MODEL", "1", "UCSDMI", "2023/02/07",
             "Illumina EMP protocol 515fbc, 806r amplification of 16S rRNA V4",
             "230207_M05314_0346_000000000-KVMGL_SMPL1_S1_L001",
             "230207_M05314_0346_000000000-KVMGL",
             "ABTX_20230208_11052_Plate_238", "Sequencing by synthesis", "GT",
             "GTGTGYCAGCMGCCGCGGTAA", "16S rRNA",
             "FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT", 1,
             "Echo 550", "ABTX_20230208_11052_Plate_238_11.8.21.RK.FH_A1",
             float("nan"), float("nan"), 1317793, float("nan"), "V4", "A1",
             "ABTX_20230208_ABTX_11052", "108379Z", "sample.1", "HT", 122822,
             1331807, float("nan"), "A1", float("nan"), ],
            ["sample.2", "TCCATACCGGAA", "UCSDMI", "Illumina",
             "SOME_CENTER_PROJECT_NAME",
             "This is a description of the experiment design.",
             "SOME_INSTRUMENT_MODEL", "1", "UCSDMI", "2023/02/07",
             "Illumina EMP protocol 515fbc, 806r amplification of 16S rRNA V4",
             "230207_M05314_0346_000000000-KVMGL_SMPL1_S1_L001",
             "230207_M05314_0346_000000000-KVMGL",
             "ABTX_20230208_11052_Plate_238", "Sequencing by synthesis", "GT",
             "GTGTGYCAGCMGCCGCGGTAA", "16S rRNA",
             "FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT", 1,
             "Echo 550", "ABTX_20230208_11052_Plate_238_11.17.21.RK.FH_A2",
             float("nan"), float("nan"), 1317793, float("nan"), "V4", "A2",
             "ABTX_20230208_ABTX_11052", "108379Z", "sample.2", "HT", 122822,
             1331807, float("nan"), "A2", float("nan"), ],
        ]

        # confirm that the observed data in the prep-info output matches
        # what's expected.

        exp_df = pd.DataFrame(columns=exp_columns, data=exp_data)

        # ensure the column order for the observed dataframe is the same as
        # what's expected. (a canonical column ordering isn't required.)
        obs_df = obs_df[exp_df.columns].copy()

        pd.testing.assert_frame_equal(obs_df, exp_df)

    def test_preparations_for_run(self):
        sheet = MetagenomicSampleSheetv100(self.ss)

        # obs will be a dictionary of dataframes, with the keys being
        # a triplet of strings, rather than a single string.
        obs = preparations_for_run(self.good_run,
                                   sample_sheet_to_dataframe(sheet),
                                   sheet.GENERATED_PREP_COLUMNS,
                                   sheet.CARRIED_PREP_COLUMNS)
        self._check_run_191103_D32611_0365_G00DHB5YXX(obs)

    def test_preparations_for_run_missing_columns(self):
        # Check that warnings are raised whenever we overwrite the
        # "well_description" column with the "description" column
        sheet = MetagenomicSampleSheetv100(self.ss)
        ss = sample_sheet_to_dataframe(sheet)
        ss['description'] = ss['well_description'].copy()
        ss.drop('well_description', axis=1, inplace=True)

        with self.assertWarns(UserWarning) as cm:
            obs = preparations_for_run(self.good_run, ss,
                                       sheet.GENERATED_PREP_COLUMNS,
                                       sheet.CARRIED_PREP_COLUMNS)

            self.assertEqual(str(cm.warnings[0].message), "'well_description' "
                                                          "is not present in s"
                                                          "ample-sheet. It is "
                                                          "not a required colu"
                                                          "mn but it is a reco"
                                                          "mmended one.")
            self.assertEqual(str(cm.warnings[1].message), "Using 'description'"
                                                          " instead of 'well_d"
                                                          "escription' because"
                                                          " that column isn't "
                                                          "present")

        self._check_run_191103_D32611_0365_G00DHB5YXX(obs)

    def test_preparations_for_run_mf(self):
        # mf has a header w/mixed case. This will test whether mapping-file
        # headers are properly converted to all lower-case.
        mf = pandas.read_csv(self.mf, delimiter='\t')

        # obs will be a dictionary of dataframes, with the keys being
        # a triplet of strings, rather than a single string.
        obs = preparations_for_run_mapping_file(self.amplicon_run, mf)

        self._check_run_230207_M05314_0346_000000000_KVMGL(obs)

        # remove a column, simulating reading in a file with a column
        # missing.
        mf = mf.drop('project_name', axis=1)

        with self.assertRaisesRegex(ValueError,
                                    'Required columns are missing: '
                                    'project_name'):
            preparations_for_run_mapping_file(self.amplicon_run, mf)

    def test_invalid_sample_names_show_warning(self):
        ss = MetagenomicSampleSheetv90(self.ss)

        self.assertIsNotNone(ss)

        ss = sample_sheet_to_dataframe(ss)

        ss['well_description'] = ss['well_description'].str.replace(
            'importantsample44', 'important-sample44')

        with self.assertWarns(UserWarning) as cm:
            _check_invalid_names(ss['well_description'])
            self.assertEqual(str(cm.warnings[0].message), 'The following sampl'
                                                          'e names have invali'
                                                          'd characters: "impo'
                                                          'rtant-sample44"'
                             )

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
        obs = get_run_prefix(self.good_run, 'Baz_12345', 'sample_1', '1')
        self.assertEqual('sample_1_S11_L001', obs)

        obs = get_run_prefix(self.good_run, 'Baz_12345', 'sample_1', '3')
        self.assertEqual('sample_1_S11_L003', obs)

        obs = get_run_prefix(self.good_run, "Baz_12345", "sample_2", "1")
        self.assertEqual("sample_2_S10_L001", obs)

        # project 2
        obs = get_run_prefix(self.good_run, 'FooBar_666', 'sample_31', '3')
        self.assertEqual('sample_31_S13_L003', obs)

        obs = get_run_prefix(self.good_run, 'FooBar_666', 'sample_32', '3')
        self.assertEqual('sample_32_S19_L003', obs)

        obs = get_run_prefix(self.good_run, 'FooBar_666', 'sample_34', '3')
        self.assertEqual('sample_34_S33_L003', obs)

    def test_get_run_prefix_fastp_minimap(self):
        obs = get_run_prefix(self.good_run_new_version, 'Baz_12345', 'sample1',
                             '1')
        self.assertEqual('sample1_S11_L001', obs)

        obs = get_run_prefix(self.good_run_new_version, 'Baz_12345', 'sample1',
                             '3')
        self.assertEqual('sample1_S11_L003', obs)

        obs = get_run_prefix(self.good_run_new_version, 'Baz_12345', 'sample2',
                             '1')
        self.assertEqual('sample2_S10_L001', obs)

        obs = get_run_prefix(self.good_run_new_version, 'Baz_12345',
                             'sample44', '3')
        self.assertIsNone(obs)

        obs = get_run_prefix(self.good_run_new_version, 'Baz_12345', 'sample2',
                             '3')
        self.assertIsNone(obs)

        # project 2
        obs = get_run_prefix(self.good_run_new_version, 'FooBar_666',
                             'sample31', '3')
        self.assertEqual('sample31_S13_L003', obs)

        obs = get_run_prefix(self.good_run_new_version, 'FooBar_666',
                             'sample32', '3')
        self.assertEqual('sample32_S19_L003', obs)

        obs = get_run_prefix(self.good_run_new_version, 'FooBar_666',
                             'sample34', '3')
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
                                 'sample31', '3')
            self.assertIsNone(obs)

        obs = get_run_prefix(self.OKish_run_new_version, 'FooBar_666',
                             'sample32', '3')
        self.assertEqual('sample32_S19_L003', obs)

        obs = get_run_prefix(self.OKish_run_new_version, 'FooBar_666',
                             'sample34', '3')
        self.assertIsNone(obs)

    def test_is_non_empty_gz_file(self):
        non_empty = join(self.good_run, 'Baz_12345', 'sample_2_S10_L001_R1_0'
                                                     '01.fastq.gz')
        self.assertTrue(is_nonempty_gz_file(non_empty))
        non_empty = join(self.good_run, 'Baz_12345', 'sample_2_S10_L001_R2_0'
                                                     '01.fastq.gz')
        self.assertTrue(is_nonempty_gz_file(non_empty))
        empty = join(self.good_run, 'Baz_12345/atropos_qc/sample_2_S10_L003_'
                                    'R1_001.fastq.gz')
        self.assertFalse(is_nonempty_gz_file(empty))
        empty = join(self.good_run, 'Baz_12345/atropos_qc/sample_2_S10_L003_'
                                    'R2_001.fastq.gz')
        self.assertFalse(is_nonempty_gz_file(empty))

    def test_parse_illumina_run_id(self):
        date, rid = parse_illumina_run_id('161004_D00611_0365_AH2HJ5BCXY')
        self.assertEqual(date, '2016-10-04')
        self.assertEqual(rid, 'D00611_0365_AH2HJ5BCXY')

        date, rid = parse_illumina_run_id('160909_K00180_0244_BH7VNKBBXX')
        self.assertEqual(date, '2016-09-09')
        self.assertEqual(rid, 'K00180_0244_BH7VNKBBXX')

        date, rid = parse_illumina_run_id('20220303_FS10001773_6_BRB11606-1914')  # noqa
        self.assertEqual(date, '2022-03-03')
        self.assertEqual(rid, 'FS10001773_6_BRB11606-1914')

    def test_parse_illumina_run_id_malformed(self):
        bad_entries = ['0220303_FS10001773_6_BRB11606-1914',
                       'verybad', '123456-bad', '12345678-bad']
        for bad in bad_entries:
            with self.assertRaises(ValueError):
                parse_illumina_run_id(bad)

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
                   'sample_plate', 'well_id_384', 'i7_index_id', 'index',
                   'i5_index_id', 'index2', 'lane', 'sample_project',
                   'well_description']

        data = [['importantsample1', 'EXPERIMENT_DESC', 'LIBRARY_PROTOCOL',
                 'Illumina', 'UCSDMI', '2019-11-03', 'sample1_S11_L003',
                 'sequencing by synthesis', 'UCSD', 'Baz_12345',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'FooBar_666_p1', 'A3', 'iTru7_107_09', 'GCCTTGTT',
                 'iTru5_01_A', 'AACACCAC', '3', 'Baz_12345',
                 'FooBar_666_p1.sample1.A3'],
                ['importantsample44', 'EXPERIMENT_DESC', 'LIBRARY_PROTOCOL',
                 'Illumina', 'UCSDMI', '2019-11-03', 'sample44_S14_L003',
                 'sequencing by synthesis', 'UCSD', 'Baz_12345',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'Baz_12345_p3', 'B99', 'iTru7_107_14', 'GTCCTAAG',
                 'iTru5_01_A', 'CATCTGCT', '3', 'Baz_12345',
                 'Baz_12345_p3.sample44.B99']]
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
                   'sample_plate', 'well_id_384', 'i7_index_id', 'index',
                   'i5_index_id', 'index2', 'lane', 'sample_project',
                   'well_description']

        # the first sample name should be padded with 2 zeros
        # the second number should be padded with 1 zero
        # the third should be ignored
        data = [['8675309', 'EXPERIMENT_DESC', 'LIBRARY_PROTOCOL', 'Illumina',
                 'UCSDMI', '2019-11-03', 'sample1_S11_L003',
                 'sequencing by synthesis', 'UCSD', 'Baz_12345',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'FooBar_666_p1', 'A3', 'iTru7_107_09', 'GCCTTGTT',
                 'iTru5_01_A', 'AACACCAC', '3', 'Baz_12345',
                 'FooBar_666_p1.sample1.A3'],
                ['867530.9', 'EXPERIMENT_DESC', 'LIBRARY_PROTOCOL', 'Illumina',
                 'UCSDMI', '2019-11-03', 'sample44_S14_L003',
                 'sequencing by synthesis', 'UCSD', 'Baz_12345',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'Baz_12345_p3', 'B99', 'iTru7_107_14', 'GTCCTAAG',
                 'iTru5_01_A', 'CATCTGCT', '3', 'Baz_12345',
                 'Baz_12345_p3.sample44.B99'],
                ['notanumber', 'EXPERIMENT_DESC', 'LIBRARY_PROTOCOL',
                 'Illumina', 'UCSDMI', '2019-11-03', 'sample44_S14_L003',
                 'sequencing by synthesis', 'UCSD', 'Baz_12345',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'Baz_12345_p3', 'B99', 'iTru7_107_14', 'GTCCTAAG',
                 'iTru5_01_A', 'CATCTGCT', '3', 'Baz_12345',
                 'Baz_12345_p3.sample44.B99']]
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
                   'sample_plate', 'well_id_384', 'i7_index_id', 'index',
                   'i5_index_id', 'index2', 'lane', 'sample_project',
                   'well_description']

        data = [['importantsample1', 'EXPERIMENT_DESC', 'LIBRARY_PROTOCOL',
                 'Illumina', 'UCSDMI', '2019-11-03', 'sample1_S11_L003',
                 'sequencing by synthesis', 'UCSD', 'Baz_12345',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'FooBar_666_p1', 'A3', 'iTru7_107_09', 'GCCTTGTT',
                 'iTru5_01_A', 'AACACCAC', '3', 'Baz_12345',
                 'FooBar_666_p1.sample1.A3'],
                ['importantsample44', 'EXPERIMENT_DESC', 'LIBRARY_PROTOCOL',
                 'Illumina', 'UCSDMI', '2019-11-03', 'sample44_S14_L003',
                 'sequencing by synthesis', 'UCSD', 'Baz_12345',
                 'Illumina HiSeq 2500', '191103_D32611_0365_G00DHB5YXX',
                 'Baz_12345_p3', 'B99', 'iTru7_107_14', 'GTCCTAAG',
                 'iTru5_01_A', 'CATCTGCT', '3', 'Baz_12345',
                 'Baz_12345_p3.sample44.B99']]
        obs = pd.DataFrame(data=data, columns=columns)
        exp = obs.copy()

        # there shouldn't be any changes
        pd.testing.assert_frame_equal(agp_transform(obs, '666'), exp)

    def test_parse_prep(self):
        columns = [
            'barcode', 'primer', 'project_name', 'well_id_384',
            'primer_plate', 'plating', 'extractionkit_lot',
            'extraction_robot', 'tm1000_8_tool', 'primer_date',
            'mastermix_lot', 'water_lot', 'processing_robot',
            'sample_plate', 'linker', 'orig_name', 'well_description',
            'pcr_primers', 'center_name', 'run_center', 'platform',
            'target_subfragment', 'target_gene', 'sequencing_meth',
            'library_construction_protocol']

        additional_columns = [
            "vol_extracted_elution_ul",
            "platemap_generation_date",
            "project_abbreviation",
            "TubeCode",
            "Kathseq_RackID",
            "number_of_cells",
            "description",
        ]

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

        additional_data = ["70", "20230627", "ADAPT", "0363132553", "", "", ""]

        index = pd.Index(["X00180471", "X00180199"], name="sample_name")
        exp = pd.DataFrame(data=data, columns=columns, index=index)

        # prep-info files do not currently have versions, unlike sample-sheets.
        # Moreover, a prep-info file w/or w/out karathoseq columns will be
        # accepted by Qiita and are possible w/optional karathoseq columns in
        # MetagenomicSampleSheetv101(). Hence, we should test both versions.
        obs = parse_prep(self.legacy_prep)
        pd.testing.assert_frame_equal(obs, exp)

        # append karathoseq metadata and recreate the expected DataFrame.
        columns += additional_columns
        data[0] += additional_data
        data[1] += additional_data

        exp = pd.DataFrame(data=data, columns=columns, index=index)
        obs = parse_prep(self.prep)

        pd.testing.assert_frame_equal(obs, exp)

    def test_qiita_prep_file_duplicate_columns(self):
        columns1 = [
            "Sample",
            "Row",
            "Col",
            "Blank",
            "Project Plate",
            # will convert to project_name
            "Compressed Plate Name",
            "Well",
            "Plate Position",
            "Primer Plate #",
            "Sample Plate",
            # will convert to project_name
            "Project_Name",
            "Plating",
            "Extraction Kit Lot",
            "Extraction Robot",
            "TM1000 8 Tool",
            "Primer Date",
            "MasterMix Lot",
            "Water Lot",
            "Processing Robot",
            "Original Name",
            "Plate",
            "EMP Primer Plate Well",
            "Name",
            "Illumina 5prime Adapter",
            "Golay Barcode",
            "Forward Primer Pad",
            "Forward Primer Linker",
            "515FB Forward Primer (Parada)",
            "Primer For PCR",
            "sample sheet Sample_ID",
        ]

        data1 = [
            [
                "X00180471",
                "A",
                1,
                False,
                "THDMI_10317_PUK2",
                "THDMI_UK",
                "THDMI_10317_UK2-US6",
                "A1",
                "1",
                "1",
                "THDMI_UK_Plate_2",
                "THDMI_UK",
                "SF",
                "166032128",
                "Carmen_HOWE_KF3",
                "109379Z",
                "2021-08-17",
                "978215",
                "RNBJ0628",
                "Echo550",
                "",
                "1",
                "A1",
                "515rcbc0",
                "AATGATACGGCGACCACCGAGATCTACACGCT",
                "AGCCTTCGTCGC",
                "TATGGTAATT",
                "GT",
                "GTGYCAGCMGCCGCGGTAA",
                (
                    "AATGATACGGCGACCACCGAGATCTACACGCTAGCCTTCGTCGCT"
                    "ATGGTAATTGTGTGYCAGCMGCCGCGGTAA"
                ),
                "X00180471",
            ],
            [
                "X00180199",
                "C",
                1,
                False,
                "THDMI_10317_PUK2",
                "THDMI_UK",
                "THDMI_10317_UK2-US6",
                "C1",
                "1",
                "1",
                "THDMI_UK_Plate_2",
                "THDMI_UK",
                "SF",
                "166032128",
                "Carmen_HOWE_KF3",
                "109379Z",
                "2021-08-17",
                "978215",
                "RNBJ0628",
                "Echo550",
                "",
                "1",
                "B1",
                "515rcbc12",
                "AATGATACGGCGACCACCGAGATCTACACGCT",
                "CGTATAAATGCG",
                "TATGGTAATT",
                "GT",
                "GTGYCAGCMGCCGCGGTAA",
                (
                    "AATGATACGGCGACCACCGAGATCTACACGCTCGTATAAATGCG"
                    "TATGGTAATTGTGTGYCAGCMGCCGCGGTAA"
                ),
                "X00180199",
            ],
        ]

        with self.assertRaisesRegex(ValueError, ""):
            generate_qiita_prep_file(pd.DataFrame(columns=columns1,
                                                  data=data1), "16S")

    def test_generate_qiita_prep_file(self):
        columns1 = [
            "Sample",
            "Row",
            "Col",
            "Blank",
            "Project Plate",
            "Compressed Plate Name",
            "Well",
            "Plate Position",
            "Primer Plate #",
            "Sample Plate",
            "Project_Name",
            "Plating",
            "Extraction Kit Lot",
            "Extraction Robot",
            "TM1000 8 Tool",
            "Primer Date",
            "MasterMix Lot",
            "Water Lot",
            "Processing Robot",
            "Original Name",
            "Plate",
            "EMP Primer Plate Well",
            "Name",
            "Illumina 5prime Adapter",
            "Golay Barcode",
            "Forward Primer Pad",
            "Forward Primer Linker",
            "515FB Forward Primer (Parada)",
            "Primer For PCR",
            "sample sheet Sample_ID",
        ]
        columns2 = [
            "Sample",
            "Row",
            "Col",
            "Blank",
            "Project Plate",
            "Compressed Plate Name",
            "Well",
            "Plate Position",
            "Primer Plate #",
            "Sample Plate",
            "Project_Name",
            "Plating",
            "Extraction Kit Lot",
            "Extraction Robot",
            "TM1000 8 Tool",
            "Primer Date",
            "MasterMix Lot",
            "Water Lot",
            "Processing Robot",
            "Original Name",
            "Plate",
            "EMP Primer Plate Well",
            "Name",
            "Reverse complement of 3prime Illumina Adapter",
            "Golay Barcode",
            "Reverse Primer Pad",
            "Reverse Primer Linker",
            "Reverse primer (EukBr)",
            "Primer For PCR",
            "sample sheet Sample_ID",
        ]
        columns3 = [
            "Sample",
            "Row",
            "Col",
            "Blank",
            "Project Plate",
            "Compressed Plate Name",
            "Well",
            "Plate Position",
            "Primer Plate #",
            "Sample Plate",
            "Project_Name",
            "Plating",
            "Extraction Kit Lot",
            "Extraction Robot",
            "TM1000 8 Tool",
            "Primer Date",
            "MasterMix Lot",
            "Water Lot",
            "Processing Robot",
            "Original Name",
            "Plate",
            "EMP Primer Plate Well",
            "Name",
            "Reverse complement of 3prime Illumina Adapter",
            "Golay Barcode",
            "Reverse Primer Pad",
            "Reverse Primer Linker",
            "ITS2 Reverse Primer",
            "Primer For PCR",
            "sample sheet Sample_ID",
        ]

        data1 = [
            [
                "X00180471",
                "A",
                1,
                False,
                "THDMI_10317_PUK2",
                "THDMI_10317_UK2-US6",
                "A1",
                "1",
                "1",
                "THDMI_UK_Plate_2",
                "THDMI_UK",
                "SF",
                "166032128",
                "Carmen_HOWE_KF3",
                "109379Z",
                "2021-08-17",
                "978215",
                "RNBJ0628",
                "Echo550",
                "",
                "1",
                "A1",
                "515rcbc0",
                "AATGATACGGCGACCACCGAGATCTACACGCT",
                "AGCCTTCGTCGC",
                "TATGGTAATT",
                "GT",
                "GTGYCAGCMGCCGCGGTAA",
                (
                    "AATGATACGGCGACCACCGAGATCTACACGCTAGCCTTCGTCGCT"
                    "ATGGTAATTGTGTGYCAGCMGCCGCGGTAA"
                ),
                "X00180471",
            ],
            [
                "X00180199",
                "C",
                1,
                False,
                "THDMI_10317_PUK2",
                "THDMI_10317_UK2-US6",
                "C1",
                "1",
                "1",
                "THDMI_UK_Plate_2",
                "THDMI_UK",
                "SF",
                "166032128",
                "Carmen_HOWE_KF3",
                "109379Z",
                "2021-08-17",
                "978215",
                "RNBJ0628",
                "Echo550",
                "",
                "1",
                "B1",
                "515rcbc12",
                "AATGATACGGCGACCACCGAGATCTACACGCT",
                "CGTATAAATGCG",
                "TATGGTAATT",
                "GT",
                "GTGYCAGCMGCCGCGGTAA",
                (
                    "AATGATACGGCGACCACCGAGATCTACACGCTCGTATAAATGCG"
                    "TATGGTAATTGTGTGYCAGCMGCCGCGGTAA"
                ),
                "X00180199",
            ],
        ]

        data2 = [
            [
                "X00180471",
                "A",
                1,
                False,
                "THDMI_10317_PUK2",
                "THDMI_10317_UK2-US6",
                "A1",
                "1",
                "1",
                "THDMI_UK_Plate_2",
                "THDMI_UK",
                "SF",
                "166032128",
                "Carmen_HOWE_KF3",
                "109379Z",
                "2021-08-17",
                "978215",
                "RNBJ0628",
                "Echo550",
                "",
                "1",
                "A1",
                "EukBr_Hiseq_0017",
                "CAAGCAGAAGACGGCATACGAGAT",
                "ACGAGACTGATT",
                "AGTCAGTCAG",
                "CA",
                "TGATCCTTCTGCAGGTTCACCTAC",
                (
                    "CAAGCAGAAGACGGCATACGAGATACGAGACTGATTAGTCAGTCAGCAT"
                    "GATCCTTCTGCAGGTTCACCTAC"
                ),
                "X00180471",
            ],
            [
                "X00180199",
                "C",
                1,
                False,
                "THDMI_10317_PUK2",
                "THDMI_10317_UK2-US6",
                "C1",
                "1",
                "1",
                "THDMI_UK_Plate_2",
                "THDMI_UK",
                "SF",
                "166032128",
                "Carmen_HOWE_KF3",
                "109379Z",
                "2021-08-17",
                "978215",
                "RNBJ0628",
                "Echo550",
                "",
                "1",
                "B1",
                "EukBr_Hiseq_0029",
                "CAAGCAGAAGACGGCATACGAGAT",
                "GAATACCAAGTC",
                "AGTCAGTCAG",
                "CA",
                "TGATCCTTCTGCAGGTTCACCTAC",
                (
                    "CAAGCAGAAGACGGCATACGAGATGAATACCAAGTCAGTCAGTCAGCAT"
                    "GATCCTTCTGCAGGTTCACCTAC"
                ),
                "X00180199",
            ],
        ]

        data3 = [
            [
                "X00180471",
                "A",
                1,
                False,
                "THDMI_10317_PUK2",
                "THDMI_10317_UK2-US6",
                "A1",
                "1",
                "1",
                "THDMI_UK_Plate_2",
                "THDMI_UK",
                "SF",
                "166032128",
                "Carmen_HOWE_KF3",
                "109379Z",
                "2021-08-17",
                "978215",
                "RNBJ0628",
                "Echo550",
                "",
                "1",
                "A1",
                "kabir_ITS2rcbc0",
                "CAAGCAGAAGACGGCATACGAGAT",
                "TCCCTTGTCTCC",
                "",
                "CG",
                "GCTGCGTTCTTCATCGATGC",
                (
                    "CAAGCAGAAGACGGCATACGAGATTCCCTTGTC"
                    "TCCCGGCTGCGTTCTTCATCGATGC"
                ),
                "X00180471",
            ],
            [
                "X00180199",
                "C",
                1,
                False,
                "THDMI_10317_PUK2",
                "THDMI_10317_UK2-US6",
                "C1",
                "1",
                "1",
                "THDMI_UK_Plate_2",
                "THDMI_UK",
                "SF",
                "166032128",
                "Carmen_HOWE_KF3",
                "109379Z",
                "2021-08-17",
                "978215",
                "RNBJ0628",
                "Echo550",
                "",
                "1",
                "B1",
                "kabir_ITS2rcbc12",
                "CAAGCAGAAGACGGCATACGAGAT",
                "TGCATACACTGG",
                "",
                "CG",
                "GCTGCGTTCTTCATCGATGC",
                (
                    "CAAGCAGAAGACGGCATACGAGATTGCATACAC"
                    "TGGCGGCTGCGTTCTTCATCGATGC"
                ),
                "X00180199",
            ],
        ]

        self.platedf1 = pd.DataFrame(columns=columns1, data=data1)
        self.platedf2 = pd.DataFrame(columns=columns2, data=data2)
        self.platedf3 = pd.DataFrame(columns=columns3, data=data3)

        obs1 = generate_qiita_prep_file(self.platedf1, '16S')
        obs2 = generate_qiita_prep_file(self.platedf2, '18S')
        obs3 = generate_qiita_prep_file(self.platedf3, 'ITS')

        exp_columns1 = [
            "sample_name",
            "barcode",
            "primer",
            "primer_plate",
            "well_id_384",
            "plating",
            "extractionkit_lot",
            "extraction_robot",
            "tm1000_8_tool",
            "primer_date",
            "mastermix_lot",
            "water_lot",
            "processing_robot",
            "tm300_8_tool",
            "tm50_8_tool",
            "tm10_8_tool",
            "sample_plate",
            "project_name",
            "orig_name",
            "well_description",
            "experiment_design_description",
            "library_construction_protocol",
            "linker",
            "platform",
            "run_center",
            "run_date",
            "run_prefix",
            "pcr_primers",
            "sequencing_meth",
            "target_gene",
            "target_subfragment",
            "center_name",
            "center_project_name",
            "instrument_model",
            "runid",
            "sample sheet Sample_ID",
        ]

        exp_columns2 = [
            "sample_name",
            "barcode",
            "primer",
            "primer_plate",
            "well_id_384",
            "plating",
            "extractionkit_lot",
            "extraction_robot",
            "tm1000_8_tool",
            "primer_date",
            "mastermix_lot",
            "water_lot",
            "processing_robot",
            "tm300_8_tool",
            "tm50_8_tool",
            "tm10_8_tool",
            "sample_plate",
            "project_name",
            "orig_name",
            "well_description",
            "experiment_design_description",
            "library_construction_protocol",
            "linker",
            "platform",
            "run_center",
            "run_date",
            "run_prefix",
            "pcr_primers",
            "sequencing_meth",
            "target_gene",
            "target_subfragment",
            "center_name",
            "center_project_name",
            "instrument_model",
            "runid",
            "Reverse Primer Pad",
            "Reverse primer (EukBr)",
            "sample sheet Sample_ID",
        ]

        exp_columns3 = [
            "sample_name",
            "barcode",
            "primer",
            "primer_plate",
            "well_id_384",
            "plating",
            "extractionkit_lot",
            "extraction_robot",
            "tm1000_8_tool",
            "primer_date",
            "mastermix_lot",
            "water_lot",
            "processing_robot",
            "tm300_8_tool",
            "tm50_8_tool",
            "tm10_8_tool",
            "sample_plate",
            "project_name",
            "orig_name",
            "well_description",
            "experiment_design_description",
            "library_construction_protocol",
            "linker",
            "platform",
            "run_center",
            "run_date",
            "run_prefix",
            "pcr_primers",
            "sequencing_meth",
            "target_gene",
            "target_subfragment",
            "center_name",
            "center_project_name",
            "instrument_model",
            "runid",
            "ITS2 Reverse Primer",
            "Reverse Primer Pad",
            "sample sheet Sample_ID",
        ]

        exp1 = pd.DataFrame(
            columns=exp_columns1,
            data=[
                [
                    "X00180471",
                    "AGCCTTCGTCGC",
                    "GTGYCAGCMGCCGCGGTAA",
                    "1",
                    "A1",
                    "SF",
                    "166032128",
                    "Carmen_HOWE_KF3",
                    "109379Z",
                    "2021-08-17",
                    "978215",
                    "RNBJ0628",
                    "Echo550",
                    "",
                    "",
                    "",
                    "THDMI_UK_Plate_2",
                    "THDMI_UK",
                    "X00180471",
                    "THDMI_UK_Plate_2.X00180471.A1",
                    "",
                    (
                        "Illumina EMP protocol 515fbc, 806r "
                        "amplification of 16S rRNA V4"
                    ),
                    "GT",
                    "Illumina",
                    "UCSDMI",
                    "",
                    "",
                    ("FWD:GTGYCAGCMGCCGCGGTAA; " "REV:GGACTACNVGGGTWTCTAAT"),
                    "Sequencing by synthesis",
                    "16S rRNA",
                    "V4",
                    "UCSDMI",
                    "",
                    "",
                    "",
                    "X00180471",
                ],
                [
                    "X00180199",
                    "CGTATAAATGCG",
                    "GTGYCAGCMGCCGCGGTAA",
                    "1",
                    "C1",
                    "SF",
                    "166032128",
                    "Carmen_HOWE_KF3",
                    "109379Z",
                    "2021-08-17",
                    "978215",
                    "RNBJ0628",
                    "Echo550",
                    "",
                    "",
                    "",
                    "THDMI_UK_Plate_2",
                    "THDMI_UK",
                    "X00180199",
                    "THDMI_UK_Plate_2.X00180199.C1",
                    "",
                    (
                        "Illumina EMP protocol 515fbc, "
                        "806r amplification of 16S rRNA V4"
                    ),
                    "GT",
                    "Illumina",
                    "UCSDMI",
                    "",
                    "",
                    ("FWD:GTGYCAGCMGCCGCGGTAA; " "REV:GGACTACNVGGGTWTCTAAT"),
                    "Sequencing by synthesis",
                    "16S rRNA",
                    "V4",
                    "UCSDMI",
                    "",
                    "",
                    "",
                    "X00180199",
                ],
            ],
        )

        exp2 = pd.DataFrame(
            columns=exp_columns2,
            data=[
                [
                    "X00180471",
                    "ACGAGACTGATT",
                    "CAAGCAGAAGACGGCATACGAGAT",
                    "1",
                    "A1",
                    "SF",
                    "166032128",
                    "Carmen_HOWE_KF3",
                    "109379Z",
                    "2021-08-17",
                    "978215",
                    "RNBJ0628",
                    "Echo550",
                    "",
                    "",
                    "",
                    "THDMI_UK_Plate_2",
                    "THDMI_UK",
                    "X00180471",
                    "THDMI_UK_Plate_2.X00180471.A1",
                    "",
                    "Illumina EMP 18S rRNA 1391f EukBr",
                    "CA",
                    "Illumina",
                    "UCSDMI",
                    "",
                    "",
                    ("FWD:GTACACACCGCCCGTC; " "REV:TGATCCTTCTGCAGGTTCACCTAC"),
                    "Sequencing by synthesis",
                    "18S rRNA",
                    "V9",
                    "UCSDMI",
                    "",
                    "",
                    "",
                    "AGTCAGTCAG",
                    "TGATCCTTCTGCAGGTTCACCTAC",
                    "X00180471",
                ],
                [
                    "X00180199",
                    "GAATACCAAGTC",
                    "CAAGCAGAAGACGGCATACGAGAT",
                    "1",
                    "C1",
                    "SF",
                    "166032128",
                    "Carmen_HOWE_KF3",
                    "109379Z",
                    "2021-08-17",
                    "978215",
                    "RNBJ0628",
                    "Echo550",
                    "",
                    "",
                    "",
                    "THDMI_UK_Plate_2",
                    "THDMI_UK",
                    "X00180199",
                    "THDMI_UK_Plate_2.X00180199.C1",
                    "",
                    "Illumina EMP 18S rRNA 1391f EukBr",
                    "CA",
                    "Illumina",
                    "UCSDMI",
                    "",
                    "",
                    ("FWD:GTACACACCGCCCGTC; " "REV:TGATCCTTCTGCAGGTTCACCTAC"),
                    "Sequencing by synthesis",
                    "18S rRNA",
                    "V9",
                    "UCSDMI",
                    "",
                    "",
                    "",
                    "AGTCAGTCAG",
                    "TGATCCTTCTGCAGGTTCACCTAC",
                    "X00180199",
                ],
            ],
        )

        exp3 = pd.DataFrame(
            columns=exp_columns3,
            data=[
                [
                    "X00180471",
                    "TCCCTTGTCTCC",
                    "CAAGCAGAAGACGGCATACGAGAT",
                    "1",
                    "A1",
                    "SF",
                    "166032128",
                    "Carmen_HOWE_KF3",
                    "109379Z",
                    "2021-08-17",
                    "978215",
                    "RNBJ0628",
                    "Echo550",
                    "",
                    "",
                    "",
                    "THDMI_UK_Plate_2",
                    "THDMI_UK",
                    "X00180471",
                    "THDMI_UK_Plate_2.X00180471.A1",
                    "",
                    (
                        "Illumina  EMP protocol amplification of "
                        "ITS1fbc, ITS2r"
                    ),
                    "CG",
                    "Illumina",
                    "UCSDMI",
                    "",
                    "",
                    (
                        "FWD:CTTGGTCATTTAGAGGAAGTAA; "
                        "REV:GCTGCGTTCTTCATCGATGC"
                    ),
                    "Sequencing by synthesis",
                    "ITS",
                    "ITS_1_2",
                    "UCSDMI",
                    "",
                    "",
                    "",
                    "GCTGCGTTCTTCATCGATGC",
                    "",
                    "X00180471",
                ],
                [
                    "X00180199",
                    "TGCATACACTGG",
                    "CAAGCAGAAGACGGCATACGAGAT",
                    "1",
                    "C1",
                    "SF",
                    "166032128",
                    "Carmen_HOWE_KF3",
                    "109379Z",
                    "2021-08-17",
                    "978215",
                    "RNBJ0628",
                    "Echo550",
                    "",
                    "",
                    "",
                    "THDMI_UK_Plate_2",
                    "THDMI_UK",
                    "X00180199",
                    "THDMI_UK_Plate_2.X00180199.C1",
                    "",
                    (
                        "Illumina  EMP protocol amplification of "
                        "ITS1fbc, ITS2r"
                    ),
                    "CG",
                    "Illumina",
                    "UCSDMI",
                    "",
                    "",
                    (
                        "FWD:CTTGGTCATTTAGAGGAAGTAA; "
                        "REV:GCTGCGTTCTTCATCGATGC"
                    ),
                    "Sequencing by synthesis",
                    "ITS",
                    "ITS_1_2",
                    "UCSDMI",
                    "",
                    "",
                    "",
                    "GCTGCGTTCTTCATCGATGC",
                    "",
                    "X00180199",
                ],
            ],
        )

        pd.testing.assert_frame_equal(obs1, exp1)
        pd.testing.assert_frame_equal(obs2, exp2)
        pd.testing.assert_frame_equal(obs3, exp3)

    def test_qiita_scrub_name(self):
        self.assertEqual(qiita_scrub_name('its my life'), 'its.my.life')
        self.assertEqual(qiita_scrub_name('7-11.10()'), '7-11.10.')
        self.assertEqual(qiita_scrub_name('{}/<>'), '.')
        self.assertEqual(qiita_scrub_name('th!s.has.happened'),
                         'th.s.has.happened')


class TestPrePrepReplicates(TestCase):
    def setUp(self):
        self.data_dir = join('metapool', 'tests', 'data')
        self.prep_w_replicates_path = join(self.data_dir,
                                           'pre_prep_w_replicates.csv')
        self.prep_wo_replicates_path = join(self.data_dir,
                                            'pre_prep_wo_replicates.csv')
        self.prep_legacy_path = join(self.data_dir, 'prep.tsv')
        self.prep_w_mixed_replicates = join(self.data_dir,
                                            'bad_pre_prep.csv')
        self.replicate_output_paths = [join(self.data_dir,
                                            'replicate_output5.txt'),
                                       join(self.data_dir,
                                            'replicate_output6.txt'),
                                       join(self.data_dir,
                                            'replicate_output7.txt')]

    def test_pre_prep_needs_demuxing(self):
        # Confirm that a pre-prep file w/contains_replicates column defined
        # and containing replicate samples results in True being returned
        # from pre_prep_needs_demuxing().
        df = parse_prep(self.prep_w_replicates_path)
        self.assertTrue(pre_prep_needs_demuxing(df))

        # Confirm that a pre-prep file w/contains_replicates column defined
        # and containing no replicate samples results in False being returned
        # from pre_prep_needs_demuxing().
        df = parse_prep(self.prep_wo_replicates_path)
        self.assertFalse(pre_prep_needs_demuxing(df))

        # Confirm that a legacy pre-prep file w/out contains_replicates
        # column returns False from pre_prep_needs_demuxing().
        df = parse_prep(self.prep_legacy_path)
        self.assertFalse(pre_prep_needs_demuxing(df))

        # Confirm that a pre-prep file w/contains_replicates column defined
        # w/some samples expressing True while others are False returns an
        # Error.
        with self.assertRaisesRegex(ValueError, 'all values in contains_'
                                                'replicates column must either'
                                                ' be True or False'):
            sheet = parse_prep(self.prep_w_mixed_replicates)
            pre_prep_needs_demuxing(sheet)

    def test_demux_pre_prep(self):
        # parse_prep() extended to support pre-prep files as well.
        df = parse_prep(self.prep_w_replicates_path)
        results = demux_pre_prep(df)

        # assert that the proper number of dataframes were returned.
        self.assertEqual(len(results), 3)

        # assert that each pre-prep dataframe appears in the correct order and
        # matches known results.
        for replicate_output_path, obs in zip(self.replicate_output_paths,
                                              results):
            print(f"testing against {replicate_output_path}...")
            exp = parse_prep(replicate_output_path)
            self.assertTrue(obs.equals(exp))


if __name__ == "__main__":
    main()
