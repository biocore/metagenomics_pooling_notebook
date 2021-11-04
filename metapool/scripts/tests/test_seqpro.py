import os
import unittest
from click.testing import CliRunner
from metapool.scripts.seqpro import format_preparation_files
from shutil import copy, copytree, rmtree
from os.path import join


class SeqproTests(unittest.TestCase):
    def setUp(self):
        # we need to get the test data directory in the parent directory
        # important to use abspath because we use CliRunner.isolated_filesystem
        tests_dir = os.path.abspath(os.path.dirname(__file__))
        tests_dir = os.path.dirname(os.path.dirname(tests_dir))
        data_dir = os.path.join(tests_dir, 'tests', 'data')

        self.run = os.path.join(data_dir, 'runs',
                                '191103_D32611_0365_G00DHB5YXX')
        self.sheet = os.path.join(self.run, 'sample-sheet.csv')

        self.fastp_run = os.path.join(data_dir, 'runs',
                                      '200318_A00953_0082_AH5TWYDSXY')
        self.fastp_sheet = os.path.join(self.fastp_run, 'sample-sheet.csv')

    def test_atropos_run(self):
        runner = CliRunner()

        with runner.isolated_filesystem():
            result = runner.invoke(format_preparation_files,
                                   args=[self.run, self.sheet, './',
                                         '--pipeline', 'atropos-and-bowtie2'])

            self.assertEqual(result.output,
                             'Stats collection is not supported for pipeline '
                             'atropos-and-bowtie2\n')
            self.assertEqual(result.exit_code, 0)

            exp_preps = [
                '191103_D32611_0365_G00DHB5YXX.Baz.1.tsv',
                '191103_D32611_0365_G00DHB5YXX.Baz.3.tsv',
                '191103_D32611_0365_G00DHB5YXX.FooBar_666.3.tsv'
            ]

            self.assertEqual(sorted(os.listdir('./')), exp_preps)

            for prep, exp_lines in zip(exp_preps, [4, 4, 5]):
                with open(prep) as f:
                    self.assertEqual(len(f.read().split('\n')), exp_lines,
                                     'Assertion error in %s' % prep)

    def test_fastp_run(self):
        runner = CliRunner()

        with runner.isolated_filesystem():
            result = runner.invoke(format_preparation_files,
                                   args=[self.fastp_run, self.fastp_sheet,
                                         './', '--pipeline',
                                         'fastp-and-minimap2'])

            self.assertEqual(result.output, '')
            self.assertEqual(result.exit_code, 0)

            exp_preps = [
                '200318_A00953_0082_AH5TWYDSXY.Project_1111.1.tsv',
                '200318_A00953_0082_AH5TWYDSXY.Project_1111.3.tsv',
                '200318_A00953_0082_AH5TWYDSXY.Trojecp_666.3.tsv'
            ]

            self.assertEqual(sorted(os.listdir('./')), exp_preps)

            for prep, exp_lines in zip(exp_preps, [4, 4, 5]):
                with open(prep) as f:
                    self.assertEqual(len(f.read().split('\n')), exp_lines,
                                     'Assertion error in %s' % prep)


class SeqproBCLConvertTests(unittest.TestCase):
    def setUp(self):
        # we need to get the test data directory in the parent directory
        # important to use abspath because we use CliRunner.isolated_filesystem
        tests_dir = os.path.abspath(os.path.dirname(__file__))
        tests_dir = os.path.dirname(os.path.dirname(tests_dir))
        self.data_dir = os.path.join(tests_dir, 'tests', 'data')

        self.fastp_run = os.path.join(self.data_dir, 'runs',
                                      '200318_A00953_0082_AH5TWYDSXY')
        self.fastp_sheet = os.path.join(self.fastp_run, 'sample-sheet.csv')

        # before continuing, create a copy of 200318_A00953_0082_AH5TWYDSXY
        # and replace Stats sub-dir with Reports.
        self.temp_copy = self.fastp_run.replace('200318', '200418')
        copytree(self.fastp_run, self.temp_copy)
        rmtree(join(self.temp_copy, 'Stats'))
        os.makedirs(join(self.temp_copy, 'Reports'))
        copy(join(self.data_dir, 'Demultiplex_Stats.csv'),
             join(self.temp_copy, 'Reports', 'Demultiplex_Stats.csv'))

    def test_fastp_run(self):
        runner = CliRunner()

        with runner.isolated_filesystem():
            result = runner.invoke(format_preparation_files,
                                   args=[self.temp_copy, self.fastp_sheet,
                                         './', '--pipeline',
                                         'fastp-and-minimap2'])
            self.assertEqual(result.output, '')
            self.assertEqual(result.exit_code, 0)

            exp_preps = [
                '200418_A00953_0082_AH5TWYDSXY.Project_1111.1.tsv',
                '200418_A00953_0082_AH5TWYDSXY.Project_1111.3.tsv',
                '200418_A00953_0082_AH5TWYDSXY.Trojecp_666.3.tsv'
            ]

            self.assertEqual(sorted(os.listdir('./')), exp_preps)

            for prep, exp_lines in zip(exp_preps, [4, 4, 5]):
                with open(prep) as f:
                    self.assertEqual(len(f.read().split('\n')), exp_lines,
                                     'Assertion error in %s' % prep)

    def tearDown(self):
        rmtree(self.temp_copy)


if __name__ == '__main__':
    unittest.main()
