import json
import os
import tempfile
import shutil
import pandas as pd

from sample_sheet import Sample
from unittest import main, TestCase

from metapool import KLSampleSheet
from metapool.count import (_extract_name_and_lane, _parse_samtools_counts,
                            _parse_fastp_counts, bcl2fastq_counts,
                            fastp_counts, minimap2_counts, run_counts,
                            _parsefier)


class TestCount(TestCase):
    def setUp(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.run_dir = os.path.join(data_dir, 'runs',
                                    '200318_A00953_0082_AH5TWYDSXY')
        self.ss = KLSampleSheet(os.path.join(self.run_dir, 'sample-sheet.csv'))

        self.stats = pd.DataFrame(RUN_STATS)
        # help make comparisons consistent
        self.stats.sort_index(inplace=True)

        self.stats.index.set_names(['Sample_ID', 'Lane'], inplace=True)

    def test_extract_name_and_lane(self):
        self.assertEqual(
            _extract_name_and_lane('33333_G2750L_S2031_L001_I1_001.fastq.gz'),
                                  ('33333_G2750L', '1'))
        self.assertEqual(
            _extract_name_and_lane('33333_G2750L_S2031_L001_R1_001.fastq.gz'),
                                  ('33333_G2750L', '1'))
        self.assertEqual(
            _extract_name_and_lane('33333_G2750L_S2031_L001_R2_001.fastq.gz'),
                                  ('33333_G2750L', '1'))
        self.assertEqual(
            _extract_name_and_lane('33333_G2751R_S2072_L009_R1_001.fastq.gz'),
                                  ('33333_G2751R', '9'))
        self.assertEqual(
            _extract_name_and_lane('33333_G2751R_S2072_L010_R1_001.fastq.gz'),
                                  ('33333_G2751R', '10'))

    def test_extract_name_and_lane_terrible_pattern(self):
        # this is likely to never happen but we label a sample with the same
        # scheme that Illumina would use to identify differnt cells, lanes, and
        # orientations
        self.assertEqual(
            _extract_name_and_lane('S2031_L001_R1_S2031_L001_I1_001.fastq.gz'),
                                  ('S2031_L001_R1', '1'))

    def test_parsefier_multiple_matches_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = os.path.join(tmp, 'funky-rerun-with-repeated-samples')
            shutil.copytree(self.run_dir, run)

            # sample 3 exists, but not with cell number S458, so this should
            # raise an error because if this happense something else went wrong
            fake = os.path.join(run, 'Trojecp_666', 'json',
                                'sample3_S458_L003_R1_001.json')
            with open(fake, 'w') as f:
                f.write(json.dumps({}))

            msg = ('Multiple matches found for the same samples in the same '
                   'lane, only one match is expected: sample3 in lane 3')
            with self.assertRaisesRegex(ValueError, msg):
                _parsefier(run, self.ss, 'json', '.json', 'halloween',
                           lambda x: 1)

    def test_parsefier_no_logs_warns(self):
        self.ss.add_sample(Sample({
            'Sample_ID': 'H20_Myers',
            'Lane': '1',
            'Sample_Name': 'H20_Myers',
            'index': 'ACTTTGTTGGAA',
            'index2': 'GGTTAATTGAGA',
            'Sample_Project': 'Trojecp_666'
        }))

        exp = pd.DataFrame(data=[[1], [1], [1], [1], [1], [1], [1]],
                           columns=['halloween'],
                           index=self.stats.index.copy())

        with self.assertWarnsRegex(UserWarning, 'No halloween log found for '
                                   'these samples: H20_Myers'):
            obs = _parsefier(self.run_dir, self.ss, 'json', '.json',
                             'halloween', lambda x: 1)

        pd.testing.assert_frame_equal(obs.sort_index(), exp)

    def test_parse_fastp_malformed(self):
        with tempfile.NamedTemporaryFile('w+') as tmp:
            tmp.write(json.dumps({}))
            tmp.seek(0)

            with self.assertRaisesRegex(ValueError, 'The fastp log for '
                                                    f'{tmp.name} is'
                                                    ' malformed'):
                _parse_fastp_counts(tmp.name)

            tmp.write(json.dumps({'summary': {}}))
            tmp.seek(0)

            with self.assertRaisesRegex(ValueError, 'The fastp log for '
                                                    f'{tmp.name} is'
                                                    ' malformed'):
                _parse_fastp_counts(tmp.name)

            tmp.write(json.dumps({'summary': {'after_filtering': {}}}))
            tmp.seek(0)

            with self.assertRaisesRegex(ValueError, 'The fastp log for '
                                                    f'{tmp.name} is'
                                                    ' malformed'):
                _parse_fastp_counts(tmp.name)

    def test_parse_fastp_counts(self):
        obs = _parse_fastp_counts(
            os.path.join(self.run_dir, 'Trojecp_666', 'json',
                         'sample3_S457_L003_R1_001.json'))

        self.assertEqual(obs, 4692)

    def test_parse_samtools_malformed(self):
        with tempfile.NamedTemporaryFile('w+') as tmp:
            tmp.write('[hey] we processed like 42 reads\n')
            tmp.seek(0)

            with self.assertRaisesRegex(ValueError, 'The samtools log for '
                                                    f'{tmp.name} is'
                                                    ' malformed'):
                _parse_samtools_counts(tmp.name)

    def test_parse_samtools_counts(self):
        obs = _parse_samtools_counts(
            os.path.join(self.run_dir, 'Trojecp_666', 'samtools',
                         'sample4_S369_L003_R1_001.log'))

        self.assertEqual(obs, 2777)

    def test_bcl2fastq_no_stats_file(self):
        bad_dir = os.path.join(os.path.abspath(self.run_dir), 'Trojecp_666')

        with self.assertRaisesRegex(IOError,
                                    rf"Cannot find stats file \({bad_dir}"
                                    r"/Stats/Stats.json\) for this run"):
            bcl2fastq_counts(bad_dir, self.ss)

    def test_bcl2fastq_counts_malformed_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats = os.path.join(tmp, 'Stats')
            os.makedirs(stats)

            with open(os.path.join(stats, 'Stats.json'), 'w') as f:
                f.write(json.dumps({}))

            with self.assertRaisesRegex(KeyError, 'bcl stats file is missing '
                                                  'ConversionResults '
                                                  'attribute'):
                bcl2fastq_counts(tmp, self.ss)

    def test_bcl2fastq_counts_malformed_lane(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats = os.path.join(tmp, 'Stats')
            os.makedirs(stats)

            with open(os.path.join(stats, 'Stats.json'), 'w') as f:
                f.write(json.dumps({'ConversionResults': [{}]}))

            with self.assertRaisesRegex(KeyError, 'bcl stats file is missing '
                                                  'DemuxResults '
                                                  'attribute'):
                bcl2fastq_counts(tmp, self.ss)

    def test_bcl2fastq_counts_malformed_lane_number(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats = os.path.join(tmp, 'Stats')
            os.makedirs(stats)

            with open(os.path.join(stats, 'Stats.json'), 'w') as f:
                f.write(json.dumps(
                    {'ConversionResults': [{'DemuxResults': {}}]}))

            with self.assertRaisesRegex(KeyError, 'bcl stats file is missing '
                                                  'LaneNumber attribute'):
                bcl2fastq_counts(tmp, self.ss)

    def test_bcl2fastq_counts(self):
        obs = bcl2fastq_counts(self.run_dir, self.ss)
        pd.testing.assert_frame_equal(obs.sort_index(),
                                      self.stats[['bcl_counts']])

    def test_fastp_counts(self):
        obs = fastp_counts(self.run_dir, self.ss)
        pd.testing.assert_frame_equal(obs.sort_index(),
                                      self.stats[['fastp_counts']])

    def test_minimap2_counts(self):
        obs = minimap2_counts(self.run_dir, self.ss)
        pd.testing.assert_frame_equal(obs.sort_index(),
                                      self.stats[['minimap2_counts']])

    def test_count_collector(self):
        obs = run_counts(self.run_dir, self.ss)
        pd.testing.assert_frame_equal(obs.sort_index(), self.stats)


RUN_STATS = {
    'bcl_counts': {('sample1', '1'): 10000, ('sample2', '1'): 100000,
                   ('sample1', '3'): 100000, ('sample2', '3'): 2300000,
                   ('sample3', '3'): 300000, ('sample4', '3'): 400000,
                   ('sample5', '3'): 567000},
    'fastp_counts': {('sample1', '1'): 10800, ('sample2', '1'): 61404,
                     ('sample1', '3'): 335996, ('sample2', '3'): 18374,
                     ('sample3', '3'): 4692, ('sample4', '3'): 960,
                     ('sample5', '3'): 30846196},
    'minimap2_counts': {('sample1', '1'): 111172.0, ('sample2', '1'): 277611.0,
                        ('sample1', '3'): 1168275.0, ('sample2', '3'): 1277.0,
                        ('sample3', '3'): 33162.0, ('sample4', '3'): 2777.0,
                        ('sample5', '3'): 4337654.0}
}


if __name__ == '__main__':
    main()
