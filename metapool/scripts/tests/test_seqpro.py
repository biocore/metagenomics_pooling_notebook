import os
import re
import unittest
from click.testing import CliRunner
from metapool.scripts.seqpro import format_preparation_files
from shutil import copy, copytree, rmtree
from os.path import join, exists
from subprocess import Popen, PIPE
import pandas as pd


class SeqproTests(unittest.TestCase):
    def setUp(self):
        # we need to get the test data directory in the parent directory
        # important to use abspath because we use CliRunner.isolated_filesystem
        tests_dir = os.path.abspath(os.path.dirname(__file__))
        tests_dir = os.path.dirname(os.path.dirname(tests_dir))
        self.test_dir = os.path.join(tests_dir, "tests")
        data_dir = os.path.join(self.test_dir, "data")
        self.vf_test_dir = os.path.join(tests_dir, "tests", "VFTEST")

        self.run = os.path.join(
            data_dir, "runs", "191103_D32611_0365_G00DHB5YXX"
        )
        self.sheet = os.path.join(self.run, "sample-sheet.csv")

        self.fastp_run = os.path.join(
            data_dir, "runs", "200318_A00953_0082_AH5TWYDSXY"
        )

        self.fastp_sheet = os.path.join(self.fastp_run, "sample-sheet.csv")
        self.v90_test_sheet = os.path.join(
            self.fastp_run, "mgv90_test_sheet.csv"
        )

    def tearDown(self):
        rmtree(self.vf_test_dir, ignore_errors=True)

    def test_fastp_run(self):
        runner = CliRunner()

        with runner.isolated_filesystem():
            result = runner.invoke(
                format_preparation_files,
                args=[
                    self.fastp_run,
                    self.fastp_sheet,
                    "./"
                ],
            )

            self.assertEqual(result.output, "")
            self.assertEqual(result.exit_code, 0)

            exp_preps = [
                "200318_A00953_0082_AH5TWYDSXY.Project_1111.1.tsv",
                "200318_A00953_0082_AH5TWYDSXY.Project_1111.3.tsv",
                "200318_A00953_0082_AH5TWYDSXY.Trojecp_666.3.tsv",
            ]

            # assert filenames are correct, and contents are correct,
            # including columns, column order, and empty values are not
            # present.
            exp = {
                "200318_A00953_0082_AH5TWYDSXY.Project_1111.1.tsv": {
                    0: {
                        "experiment_design_description": "Eqiiperiment",
                        "well_description": "FooBar_666_p1.sample1.A1",
                        "library_construction_protocol": "Knight Lab Kapa HP",
                        "platform": "Illumina",
                        "run_center": "IGM",
                        "run_date": "2020-03-18",
                        "run_prefix": "sample1_S333_L001",
                        "sequencing_meth": "sequencing by synthesis",
                        "center_name": "UCSD",
                        "center_project_name": "Project",
                        "instrument_model": "Illumina NovaSeq 6000",
                        "runid": "200318_A00953_0082_AH5TWYDSXY",
                        "lane": 1,
                        "sample_project": "Project",
                        "i5_index_id": "iTru5_01_A",
                        "index2": "ACCGACAA",
                        "sample_plate": "FooBar_666_p1",
                        "well_id_384": "A1",
                        "sample_name": "sample1",
                        "index": "CCGACTAT",
                        "i7_index_id": "iTru7_107_07",
                        "raw_reads_r1r2": 10000,
                        "quality_filtered_reads_r1r2": 16.0,
                        "total_biological_reads_r1r2": 10800.0,
                        "fraction_passing_quality_filter": 0.0016
                    },
                    1: {
                        "experiment_design_description": "Eqiiperiment",
                        "well_description": "FooBar_666_p1.sample2.A2",
                        "library_construction_protocol": "Knight Lab Kapa HP",
                        "platform": "Illumina",
                        "run_center": "IGM",
                        "run_date": "2020-03-18",
                        "run_prefix": "sample2_S404_L001",
                        "sequencing_meth": "sequencing by synthesis",
                        "center_name": "UCSD",
                        "center_project_name": "Project",
                        "instrument_model": "Illumina NovaSeq 6000",
                        "runid": "200318_A00953_0082_AH5TWYDSXY",
                        "lane": 1,
                        "sample_project": "Project",
                        "i5_index_id": "iTru5_01_A",
                        "index2": "CTTCGCAA",
                        "sample_plate": "FooBar_666_p1",
                        "well_id_384": "A2",
                        "sample_name": "sample2",
                        "index": "CCGACTAT",
                        "i7_index_id": "iTru7_107_08",
                        "raw_reads_r1r2": 100000,
                        "quality_filtered_reads_r1r2": 16.0,
                        "total_biological_reads_r1r2": 61404.0,
                        "fraction_passing_quality_filter": 0.00016
                    },
                },
                "200318_A00953_0082_AH5TWYDSXY.Project_1111.3.tsv": {
                    0: {
                        "experiment_design_description": "Eqiiperiment",
                        "well_description": "FooBar_666_p1.sample1.A3",
                        "library_construction_protocol": "Knight Lab Kapa HP",
                        "platform": "Illumina",
                        "run_center": "IGM",
                        "run_date": "2020-03-18",
                        "run_prefix": "sample1_S241_L003",
                        "sequencing_meth": "sequencing by synthesis",
                        "center_name": "UCSD",
                        "center_project_name": "Project",
                        "instrument_model": "Illumina NovaSeq 6000",
                        "runid": "200318_A00953_0082_AH5TWYDSXY",
                        "lane": 3,
                        "sample_project": "Project",
                        "i5_index_id": "iTru5_01_A",
                        "index2": "AACACCAC",
                        "sample_plate": "FooBar_666_p1",
                        "well_id_384": "A3",
                        "sample_name": "sample1",
                        "index": "GCCTTGTT",
                        "i7_index_id": "iTru7_107_09",
                        "raw_reads_r1r2": 100000,
                        "quality_filtered_reads_r1r2": 16.0,
                        "total_biological_reads_r1r2": 335996.0,
                        "fraction_passing_quality_filter": 0.00016
                    },
                    1: {
                        "experiment_design_description": "Eqiiperiment",
                        "well_description": "FooBar_666_p1.sample2.A4",
                        "library_construction_protocol": "Knight Lab Kapa HP",
                        "platform": "Illumina",
                        "run_center": "IGM",
                        "run_date": "2020-03-18",
                        "run_prefix": "sample2_S316_L003",
                        "sequencing_meth": "sequencing by synthesis",
                        "center_name": "UCSD",
                        "center_project_name": "Project",
                        "instrument_model": "Illumina NovaSeq 6000",
                        "runid": "200318_A00953_0082_AH5TWYDSXY",
                        "lane": 3,
                        "sample_project": "Project",
                        "i5_index_id": "iTru5_01_A",
                        "index2": "CGTATCTC",
                        "sample_plate": "FooBar_666_p1",
                        "well_id_384": "A4",
                        "sample_name": "sample2",
                        "index": "AACTTGCC",
                        "i7_index_id": "iTru7_107_10",
                        "raw_reads_r1r2": 2300000,
                        "quality_filtered_reads_r1r2": 16.0,
                        "total_biological_reads_r1r2": 18374.0,
                        "fraction_passing_quality_filter":
                            0.000006956521739130435
                    },
                },
                "200318_A00953_0082_AH5TWYDSXY.Trojecp_666.3.tsv": {
                    0: {
                        "experiment_design_description": "SomethingWitty",
                        "well_description": "FooBar_666_p1.sample3.A5",
                        "library_construction_protocol": "Knight Lab Kapa HP",
                        "platform": "Illumina",
                        "run_center": "IGM",
                        "run_date": "2020-03-18",
                        "run_prefix": "sample3_S457_L003",
                        "sequencing_meth": "sequencing by synthesis",
                        "center_name": "UCSD",
                        "center_project_name": "Trojecp",
                        "instrument_model": "Illumina NovaSeq 6000",
                        "runid": "200318_A00953_0082_AH5TWYDSXY",
                        "lane": 3,
                        "sample_project": "Trojecp",
                        "i5_index_id": "iTru5_01_A",
                        "index2": "GGTACGAA",
                        "sample_plate": "FooBar_666_p1",
                        "well_id_384": "A5",
                        "sample_name": "sample3",
                        "index": "CAATGTGG",
                        "i7_index_id": "iTru7_107_11",
                        "raw_reads_r1r2": 300000,
                        "quality_filtered_reads_r1r2": 16.0,
                        "total_biological_reads_r1r2": 4692.0,
                        "fraction_passing_quality_filter":
                        0.00005333333333333333
                    },
                    1: {
                        "experiment_design_description": "SomethingWitty",
                        "well_description": "FooBar_666_p1.sample4.B6",
                        "library_construction_protocol": "Knight Lab Kapa HP",
                        "platform": "Illumina",
                        "run_center": "IGM",
                        "run_date": "2020-03-18",
                        "run_prefix": "sample4_S369_L003",
                        "sequencing_meth": "sequencing by synthesis",
                        "center_name": "UCSD",
                        "center_project_name": "Trojecp",
                        "instrument_model": "Illumina NovaSeq 6000",
                        "runid": "200318_A00953_0082_AH5TWYDSXY",
                        "lane": 3,
                        "sample_project": "Trojecp",
                        "i5_index_id": "iTru5_01_A",
                        "index2": "CGATCGAT",
                        "sample_plate": "FooBar_666_p1",
                        "well_id_384": "B6",
                        "sample_name": "sample4",
                        "index": "AAGGCTGA",
                        "i7_index_id": "iTru7_107_12",
                        "raw_reads_r1r2": 400000,
                        "quality_filtered_reads_r1r2": 16.0,
                        "total_biological_reads_r1r2": 960.0,
                        "fraction_passing_quality_filter": 0.00004
                    },
                    2: {
                        "experiment_design_description": "SomethingWitty",
                        "well_description": "FooBar_666_p1.sample5.B8",
                        "library_construction_protocol": "Knight Lab Kapa HP",
                        "platform": "Illumina",
                        "run_center": "IGM",
                        "run_date": "2020-03-18",
                        "run_prefix": "sample5_S392_L003",
                        "sequencing_meth": "sequencing by synthesis",
                        "center_name": "UCSD",
                        "center_project_name": "Trojecp",
                        "instrument_model": "Illumina NovaSeq 6000",
                        "runid": "200318_A00953_0082_AH5TWYDSXY",
                        "lane": 3,
                        "sample_project": "Trojecp",
                        "i5_index_id": "iTru5_01_A",
                        "index2": "AAGACACC",
                        "sample_plate": "FooBar_666_p1",
                        "well_id_384": "B8",
                        "sample_name": "sample5",
                        "index": "TTACCGAG",
                        "i7_index_id": "iTru7_107_13",
                        "raw_reads_r1r2": 567000,
                        "quality_filtered_reads_r1r2": 16.0,
                        "total_biological_reads_r1r2": 30846196.0,
                        "fraction_passing_quality_filter":
                        0.000028218694885361552
                    },
                },
            }

            for prep in exp_preps:
                obs = pd.read_csv(prep, sep="\t").to_dict("index")
                self.assertDictEqual(obs, exp[prep])

    def test_verbose_flag(self):
        self.maxDiff = None
        sample_dir = "metapool/tests/data/runs/200318_A00953_0082_AH5TWYDSXY"

        cmd = [
            "seqpro",
            "--verbose",
            sample_dir,
            join(sample_dir, "sample-sheet.csv"),
            self.vf_test_dir,
        ]

        proc = Popen(
            " ".join(cmd),
            universal_newlines=True,
            shell=True,
            stdout=PIPE,
            stderr=PIPE,
        )

        stdout, stderr = proc.communicate()
        return_code = proc.returncode

        tmp = []

        # remove trailing whitespace before splitting each line into pairs.
        for line in stdout.strip().split("\n"):
            qiita_id, file_path = line.split("\t")
            # truncate full-path output to be file-system agnostic.
            file_path = re.sub(
                "^.*metagenomics_pooling_notebook/",
                "metagenomics_pooling_notebook/",
                file_path,
            )
            tmp.append(f"{qiita_id}\t{file_path}")

        stdout = "\n".join(tmp)

        self.assertEqual(
            (
                "1111\tmetagenomics_pooling_notebook/metapool/tests"
                "/VFTEST/200318_A00953_0082_AH5TWYDSXY.Project_1111"
                ".1.tsv\n1111\tmetagenomics_pooling_notebook/metapo"
                "ol/tests/VFTEST/200318_A00953_0082_AH5TWYDSXY.Proj"
                "ect_1111.3.tsv\n666\tmetagenomics_pooling_notebook"
                "/metapool/tests/VFTEST/200318_A00953_0082_AH5TWYDS"
                "XY.Trojecp_666.3.tsv"
            ),
            stdout,
        )
        self.assertEqual("", stderr)
        self.assertEqual(0, return_code)

    def test_legacy_run(self):
        runner = CliRunner()

        # confirm seqpro runs w/out error when using a legacy sample-sheet.
        # the sheet used is sample-sheet.csv converted to v90 metagenomic
        # spec.

        with runner.isolated_filesystem():
            result = runner.invoke(
                format_preparation_files,
                args=[
                    self.fastp_run,
                    self.v90_test_sheet,
                    "./"
                ],
            )

            # all we need to test is that seqpro didn't exit abnormally.
            # this implies that it successfully processed a sample-sheet w/
            # a mixed-case column ('Sample_Well').
            self.assertEqual(result.output, "")
            self.assertEqual(result.exit_code, 0)


class SeqproBCLConvertTests(unittest.TestCase):
    def setUp(self):
        # we need to get the test data directory in the parent directory
        # important to use abspath because we use CliRunner.isolated_filesystem
        tests_dir = os.path.abspath(os.path.dirname(__file__))
        tests_dir = os.path.dirname(os.path.dirname(tests_dir))
        self.data_dir = os.path.join(tests_dir, "tests", "data")

        self.fastp_run = os.path.join(
            self.data_dir, "runs", "200318_A00953_0082_AH5TWYDSXY"
        )
        self.fastp_sheet = os.path.join(self.fastp_run, "sample-sheet.csv")

        # before continuing, create a copy of 200318_A00953_0082_AH5TWYDSXY
        # and replace Stats sub-dir with Reports.
        self.temp_copy = self.fastp_run.replace("200318", "200418")
        copytree(self.fastp_run, self.temp_copy)
        rmtree(join(self.temp_copy, "Stats"))
        os.makedirs(join(self.temp_copy, "Reports"))
        copy(
            join(self.data_dir, "Demultiplex_Stats.csv"),
            join(self.temp_copy, "Reports", "Demultiplex_Stats.csv"),
        )

    def test_fastp_run(self):
        runner = CliRunner()

        with runner.isolated_filesystem():
            result = runner.invoke(
                format_preparation_files,
                args=[
                    self.temp_copy,
                    self.fastp_sheet,
                    "./"
                ],
            )
            self.assertEqual(result.output, "")
            self.assertEqual(result.exit_code, 0)

            exp_preps = [
                "200418_A00953_0082_AH5TWYDSXY.Project_1111.1.tsv",
                "200418_A00953_0082_AH5TWYDSXY.Project_1111.3.tsv",
                "200418_A00953_0082_AH5TWYDSXY.Trojecp_666.3.tsv",
            ]

            self.assertEqual(sorted(os.listdir("./")), exp_preps)

            for prep, exp_lines in zip(exp_preps, [4, 4, 5]):
                with open(prep) as f:
                    self.assertEqual(
                        len(f.read().split("\n")),
                        exp_lines,
                        "Assertion error in %s" % prep,
                    )

    def tearDown(self):
        rmtree(self.temp_copy)


class MetricsTest(unittest.TestCase):
    def setUp(self):
        tests_dir = os.path.abspath(os.path.dirname(__file__))
        tests_dir = os.path.dirname(os.path.dirname(tests_dir))
        self.test_dir = os.path.join(tests_dir, "tests")
        data_dir = os.path.join(self.test_dir, "data")
        self.run = os.path.join(
            data_dir, "runs", "240124_LH00444_0046_A223N22LT4"
        )
        self.sheet = os.path.join(data_dir, "good_metatv10_sheet.csv")
        self.delete_these = []

    def tearDown(self):
        for path in self.delete_these:
            if exists(path):
                rmtree(path)

    def test_metrics_generation(self):
        runner = CliRunner()

        out_dir = join(self.test_dir, "test_output")
        self.delete_these.append(out_dir)

        with runner.isolated_filesystem():
            result = runner.invoke(
                format_preparation_files, args=[self.run, self.sheet, out_dir]
            )

            # confirm seqpro exited successfully before continuing test.
            self.assertEqual(result.exit_code, 0)

        # confirm output prep-info file exists and is named as expected.
        obs = join(out_dir, "240124_LH00444_0046_A223N22LT4.NPH_15288.7.tsv")
        self.assertTrue(exists(obs))

        # read in observed result and determine whether it matches
        # expectations.
        obs = pd.read_csv(obs, sep="\t", dtype=str)

        exp = {
            "center_name": {0: "UCSD", 1: "UCSD"},
            "center_project_name": {0: "NPH", 1: "NPH"},
            "experiment_design_description": {
                0: "shotgun sequencing",
                1: "shotgun sequencing",
            },
            "i5_index_id": {0: "iTru5_01_A", 1: "iTru5_02_A"},
            "i7_index_id": {0: "iTru7_206_08", 1: "iTru7_206_09"},
            "index": {0: "AAGCGCAT", 1: "CACTGACA"},
            "index2": {0: "ACCGACAA", 1: "CTTCGCAA"},
            "instrument_model": {
                0: "Illumina NovaSeq X",
                1: "Illumina NovaSeq X",
            },
            "lane": {0: 7, 1: 7},
            "library_construction_protocol": {
                0: "Knight Lab Kapa HyperPlus",
                1: "Knight Lab Kapa HyperPlus",
            },
            "platform": {0: "Illumina", 1: "Illumina"},
            "run_center": {0: "IGM", 1: "IGM"},
            "run_date": {0: "2024-01-24", 1: "2024-01-24"},
            "run_prefix": {0: "359250488_S19_L007", 1: "359180426_S19_L007"},
            "runid": {
                0: "240124_LH00444_0046_A223N22LT4",
                1: "240124_LH00444_0046_A223N22LT4",
            },
            "sample_name": {0: 359250488, 1: 359180426},
            "sample_plate": {0: "NPH_15288_P1", 1: "NPH_15288_P1"},
            "sample_project": {0: "NPH", 1: "NPH"},
            "sequencing_meth": {
                0: "sequencing by synthesis",
                1: "sequencing by synthesis",
            },
            # values below should reflect columns in sample-sheet input.
            "total_rna_concentration_ng_ul": {0: 4.912, 1: 124.168},
            # values below should reflect columns in sample-sheet input.
            "vol_extracted_elution_ul": {0: 70, 1: 70},
            "well_description": {
                0: "NPH_15288_P1.359250488.A1",
                1: "NPH_15288_P1.359180426.C1",
            },
            "well_id_384": {0: "A1", 1: "C1"},
            # values below should reflect '# Reads' column found in
            # Demultiplex_Stats.csv * 2 for both forward and reverse reads.
            "raw_reads_r1r2": {0: 59709594, 1: 46883162},
            # values below should reflect values extracted from fastp-
            # generated json files.
            "total_biological_reads_r1r2": {0: 55785054.0, 1: 41738784.0},
            # values below should reflect # of sequences found in a fastq file
            # by seqtk.
            "quality_filtered_reads_r1r2": {0: 4.0, 1: 4.0},
            "fraction_passing_quality_filter": {
                0: 6.699090936709433e-08,
                1: 8.531847745252336e-08,
            }
        }

        exp = pd.DataFrame(exp, dtype=str)

        self.assertTrue(obs.equals(exp))


if __name__ == "__main__":
    unittest.main()
