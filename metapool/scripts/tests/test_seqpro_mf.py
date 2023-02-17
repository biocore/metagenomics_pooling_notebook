import os
import unittest
from click.testing import CliRunner
from metapool.scripts.seqpro_mf import format_preparation_files_mf
from shutil import copy, copytree, rmtree
from os.path import join


class SeqproAmpliconTests(unittest.TestCase):
    def setUp(self):
        # we need to get the test data directory in the parent directory
        # important to use abspath because we use CliRunner.isolated_filesystem
        tests_dir = os.path.abspath(os.path.dirname(__file__))
        tests_dir = os.path.dirname(os.path.dirname(tests_dir))
        self.data_dir = os.path.join(tests_dir, 'tests', 'data')

        self.fastp_run = os.path.join(self.data_dir, 'runs',
                                      '230207_M05314_0346_000000000-KVMGL')
        self.fastp_sheet = os.path.join(self.fastp_run,
                                        'sample_mapping_file.tsv')

        # before continuing, create a copy of
        # 230207_M05314_0346_000000000-KVMGL and replace Stats sub-dir with
        # Reports.
        self.temp_copy = self.fastp_run.replace('230207', '240207')
        copytree(self.fastp_run, self.temp_copy)
        rmtree(join(self.temp_copy, 'Stats'))
        os.makedirs(join(self.temp_copy, 'Reports'))
        copy(join(self.data_dir, 'Demultiplex_Stats.csv'),
             join(self.temp_copy, 'Reports', 'Demultiplex_Stats.csv'))

    def test_fastp_run(self):
        runner = CliRunner()

        with runner.isolated_filesystem():
            result = runner.invoke(format_preparation_files_mf,
                                   args=[self.temp_copy,
                                         join(self.temp_copy,
                                              'sample_mapping_file.tsv'),
                                         './'])

            self.assertEqual(result.output, '')
            self.assertEqual(result.exit_code, 0)

            # not exp_preps is just one file w/a long name, not two files.
            exp_preps = ['240207_M05314_0346_000000000-KVMGL.'
                         'ABTX_20230208_ABTX_11052.1.tsv']
            self.assertEqual(sorted(os.listdir('./')), exp_preps)

    def tearDown(self):
        rmtree(self.temp_copy)


if __name__ == '__main__':
    unittest.main()
