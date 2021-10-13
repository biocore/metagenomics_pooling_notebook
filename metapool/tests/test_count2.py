import os
import tempfile
import pandas as pd

from unittest import main, TestCase

from metapool import KLSampleSheet
from metapool.count import (_parse_fastp_counts, bcl2fastq_counts)


class TestCount2(TestCase):
    '''
    TestCount2 repeats specific tests from TestCount(), using a bcl-convert-
    generated Demultiplex_Stats.csv file for input instead of a bcl2fastq-
    generated Stats.json file.
    '''
    def setUp(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')

        # Use 200418_A00953_0082_AH5TWYDSXY, a sample dataset containing
        # bcl-convert-generated metadata.
        self.run_dir = os.path.join(data_dir, 'runs',
                                    '200418_A00953_0082_AH5TWYDSXY')
        self.ss = KLSampleSheet(os.path.join(self.run_dir, 'sample-sheet.csv'))

        self.stats = pd.DataFrame(RUN_STATS)
        # help make comparisons consistent
        self.stats.sort_index(inplace=True)

        self.stats.index.set_names(['Sample_ID', 'Lane'], inplace=True)

    def test_parse_fastp_counts(self):
        obs = _parse_fastp_counts(
            os.path.join(self.run_dir, 'Trojecp_666', 'json',
                         'sample3_S457_L003_R1_001.json'))

        self.assertEqual(obs, 4692)

    def test_bcl2fastq_no_stats_file(self):
        bad_dir = os.path.join(os.path.abspath(self.run_dir), 'Trojecp_666')
        with self.assertRaisesRegex(IOError, f"Cannot find Stats.json '"
                                             f"{bad_dir}/Stats/Stats.json' or "
                                             "Demultiplex_Stats.csv '"
                                             f"{bad_dir}/Reports/Demultiplex_"
                                             "Stats.csv' for this run"):
            bcl2fastq_counts(bad_dir, self.ss)

    def test_bcl2fastq_counts_malformed_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats = os.path.join(tmp, 'Reports')
            os.makedirs(stats)

            with open(os.path.join(stats, 'Demultiplex_Stats.csv'), 'w') as f:
                f.write('')

            with self.assertRaisesRegex(pd.errors.EmptyDataError,
                                        'No columns to parse from file'):
                bcl2fastq_counts(tmp, self.ss)

    def test_bcl2fastq_counts_malformed_lane(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats = os.path.join(tmp, 'Reports')
            os.makedirs(stats)

            with open(os.path.join(stats, 'Demultiplex_Stats.csv'), 'w') as f:
                f.write('SampleID,Index,# Reads,# Perfect Index Reads,# One '
                        'Mismatch Index Reads,# of >= Q30 Bases (PF),Mean Qu'
                        'ality Score (PF),well_description\n')
                f.write('sample1,AAAAAAAA-GGGGGGGG,10000,0,0,0,50.00,sample1'
                        '\n')

            with self.assertRaisesRegex(KeyError, r"\['Lane'\] not in index"):
                bcl2fastq_counts(tmp, self.ss)

    def test_bcl2fastq_counts(self):
        obs = bcl2fastq_counts(self.run_dir, self.ss)
        pd.testing.assert_frame_equal(obs.sort_index(),
                                      self.stats[['raw_reads']])


RUN_STATS = {
    'raw_reads': {('sample1', '1'): 10000, ('sample2', '1'): 100000,
                  ('sample1', '3'): 100000, ('sample2', '3'): 2300000,
                  ('sample3', '3'): 300000, ('sample4', '3'): 400000,
                  ('sample5', '3'): 567000}}


if __name__ == '__main__':
    main()
