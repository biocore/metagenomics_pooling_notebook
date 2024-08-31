import os
import unittest

import pandas as pd
from click.testing import CliRunner
from metapool.scripts.seqpro_mf import format_preparation_files_mf
from shutil import copy, copytree, rmtree
from os.path import join, exists
import re
from subprocess import Popen, PIPE


class SeqproAmpliconTests(unittest.TestCase):
    def setUp(self):
        # we need to get the test data directory in the parent directory
        # important to use abspath because we use CliRunner.isolated_filesystem
        tests_dir = os.path.abspath(os.path.dirname(__file__))
        tests_dir = os.path.dirname(os.path.dirname(tests_dir))
        self.test_dir = os.path.join(tests_dir, 'tests')
        self.data_dir = os.path.join(self.test_dir, 'data')
        self.vf_test_dir = os.path.join(self.test_dir, 'VFTEST')

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

    def tearDown(self):
        rmtree(self.temp_copy)
        # this output-path isn't created for all tests. ignore error if it
        # does not exist.
        rmtree(self.vf_test_dir, ignore_errors=True)

    def test_run(self):
        runner = CliRunner()

        with runner.isolated_filesystem():
            result = runner.invoke(format_preparation_files_mf,
                                   args=[self.temp_copy,
                                         join(self.temp_copy,
                                              'sample_mapping_file.tsv'),
                                         './'])

            # assert seqpro_mf returned successfully
            self.assertEqual(result.exit_code, 0)

            obs_fp = ('./230207_M05314_0346_000000000-KVMGL.ABTX_20230208_'
                      'ABTX_11052.1.tsv')

            # assert prep-info-file output exists
            self.assertTrue(exists(obs_fp))

            obs_df = pd.read_csv(obs_fp, delimiter='\t')

            # assert sample_name does not contain any '_' characters
            names = list(obs_df['sample_name'])

            # generate a list of sample-names that contain characters other
            # than alphanumerics + '.'
            names = [x for x in names if not bool(re.match(r"^[\w\d.]*$", x))]

            # assert that all sample-names were of the proper form.
            self.assertEqual(names, [])

            # confirm correct header columns exist
            self.assertEqual({'sample_name', 'center_name',
                              'center_project_name',
                              'experiment_design_description',
                              'instrument_model', 'lane',
                              'library_construction_protocol', 'platform',
                              'run_center', 'run_date', 'run_prefix', 'runid',
                              'sample_plate', 'sequencing_meth', 'barcode',
                              'linker', 'primer', 'extraction_robot',
                              'extractionkit_lot', 'mastermix_lot',
                              'orig_name', 'pcr_primers', 'plating',
                              'primer_date', 'primer_plate',
                              'processing_robot', 'project_name',
                              'target_gene', 'target_subfragment',
                              'tm1000_8_tool', 'tm300_8_tool', 'tm50_8_tool',
                              'tm10_8_tool', 'well_id_96', 'water_lot',
                              'well_description', 'well_id_384'},
                             set(obs_df.columns))

    def test_verbose_flag(self):
        self.maxDiff = None
        sample_dir = (f"{self.data_dir}/runs/230207_M05314_0346_000000000"
                      f"-KVMGL")

        cmd = ['seqpro_mf', '--verbose',
               sample_dir,
               join(sample_dir, 'sample_mapping_file.tsv'),
               self.vf_test_dir]

        proc = Popen(' '.join(cmd), universal_newlines=True, shell=True,
                     stdout=PIPE, stderr=PIPE)

        stdout, stderr = proc.communicate()
        return_code = proc.returncode

        tmp = []

        # remove trailing whitespace before splitting each line into pairs.
        for line in stdout.strip().split('\n'):
            qiita_id, file_path = line.split('\t')
            # truncate full-path output to be file-system agnostic.
            file_path = re.sub('^.*metagenomics_pooling_notebook/',
                               'metagenomics_pooling_notebook/', file_path)
            tmp.append(f'{qiita_id}\t{file_path}')

        stdout = '\n'.join(tmp)

        self.assertEqual(('11052\tmetagenomics_pooling_notebook/metapool/tests'
                          '/VFTEST/230207_M05314_0346_000000000-KVMGL.ABTX_202'
                          '30208_ABTX_11052.1.tsv'), stdout)
        self.assertEqual('', stderr)
        self.assertEqual(0, return_code)


if __name__ == '__main__':
    unittest.main()
