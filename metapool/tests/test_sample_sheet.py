import unittest
import os
import warnings
import tempfile

import pandas as pd
import sample_sheet

from metapool.sample_sheet import (KLSampleSheet, validate_sample_sheet,
                                   sample_sheet_to_dataframe,
                                   _add_metadata_to_sheet, _add_data_to_sheet,
                                   _validate_sample_sheet_metadata,
                                   make_sample_sheet)


# The classes below share the same filepaths, so we use this dummy class
class BaseTests(unittest.TestCase):
    def setUp(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.ss = os.path.join(data_dir, 'runs',
                               '191103_D32611_0365_G00DHB5YXX',
                               'sample-sheet.csv')

        self.good_ss = os.path.join(data_dir, 'good-sample-sheet.csv')

        self.no_project_ss = os.path.join(data_dir,
                                          'no-project-name-sample-sheet.csv')

        # "valid" upfront but will have repeated values after scrubbing
        self.ok_ss = os.path.join(data_dir, 'ok-sample-sheet.csv')

        self.scrubbable_ss = os.path.join(data_dir,
                                          'scrubbable-sample-sheet.csv')

        self.bad_project_name_ss = os.path.join(
            data_dir, 'bad-project-name-sample-sheet.csv')

        # TODO: should we remove because we don't use comments anymore?
        self.no_comments_ss = os.path.join(data_dir,
                                           'no-comments-sample-sheet.csv')


class KLSampleSheetTests(BaseTests):
    def test_sample_sheet_roundtripping(self):
        # testing with all the sheets we have access to
        sheets = [self.ss, self.good_ss,
                  self.no_project_ss, self.ok_ss,
                  self.scrubbable_ss, self.bad_project_name_ss]
        sheets = {sheet: KLSampleSheet(sheet) for sheet in sheets}

        for filename, sheet in sheets.items():
            with tempfile.NamedTemporaryFile('w+') as tmp:
                sheet.write(tmp)
                tmp.seek(0)
                observed = tmp.read()

                with open(filename) as expected:
                    self.assertEqual(observed.split(), expected.read().split(),
                                     f'Problem found with {filename}')

    def test_parse(self):
        sheet = KLSampleSheet(self.ss)

        exp = {
            'IEMFileVersion': '4',
            'Investigator Name': 'Caballero',
            'Experiment Name': 'RKL0042',
            'Date': '2/26/20',
            'Workflow': 'GenerateFASTQ',
            'Application': 'FASTQ Only',
            'Assay': 'Metagenomics',
            'Description': '',
            'Chemistry': 'Default'
        }

        self.assertEqual(sheet.Header, exp)
        self.assertEqual(sheet.Reads, [150, 150])
        self.assertEqual(sheet.Settings, {'ReverseComplement': '0'})

        data = (
            '1,sample1,sample1,FooBar_666_p1,A1,iTru7_107_07,CCGACTAT,'
            'iTru5_01_A,ACCGACAA,Baz,importantsample1\n'
            '1,sample2,sample2,FooBar_666_p1,A2,iTru7_107_08,CCGACTAT,'
            'iTru5_01_A,CTTCGCAA,Baz,importantsample2\n'
            '3,sample1,sample1,FooBar_666_p1,A3,iTru7_107_09,GCCTTGTT,'
            'iTru5_01_A,AACACCAC,Baz,importantsample1\n'
            '3,sample2,sample2,FooBar_666_p1,A4,iTru7_107_10,AACTTGCC,'
            'iTru5_01_A,CGTATCTC,Baz,importantsample2\n'
            '3,sample31,sample31,FooBar_666_p1,A5,iTru7_107_11,CAATGTGG,'
            'iTru5_01_A,GGTACGAA,FooBar_666,importantsample31\n'
            '3,sample32,sample32,FooBar_666_p1,B6,iTru7_107_12,AAGGCTGA,'
            'iTru5_01_A,CGATCGAT,FooBar_666,importantsample32\n'
            '3,sample34,sample34,FooBar_666_p1,B8,iTru7_107_13,TTACCGAG,'
            'iTru5_01_A,AAGACACC,FooBar_666,importantsample34\n'
            '3,sample44,sample44,Baz_p3,B99,iTru7_107_14,GTCCTAAG,'
            'iTru5_01_A,CATCTGCT,Baz,importantsample44\n'
        )
        keys = ['Lane', 'Sample_ID', 'Sample_Name', 'Sample_Plate',
                'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                'Sample_Project', 'Well_description']

        for sample, line in zip(sheet.samples, data.split()):
            values = line.strip().split(',')
            exp = sample_sheet.Sample(dict(zip(keys, values)))

            self.assertEqual(sample, exp)

        # check for Bioinformatics
        exp = pd.DataFrame(
            columns=['Sample_Project', 'QiitaID', 'BarcodesAreRC',
                     'ForwardAdapter', 'ReverseAdapter', 'HumanFiltering'],
            data=[
                ['Baz', '100', 'False', 'AACC', 'GGTT', 'False'],
                ['FooBar_666', '666', 'False', 'AACC', 'GGTT', 'False']
            ]
        )
        pd.testing.assert_frame_equal(sheet.Bioinformatics, exp)

        # check for Contact
        exp = pd.DataFrame(
            columns=['Email', 'Sample_Project'],
            data=[
                ['test@lol.com', 'Baz'],
                ['tester@rofl.com', 'FooBar_666']
            ]
        )
        pd.testing.assert_frame_equal(sheet.Contact, exp)


class SampleSheetWorkflow(BaseTests):
    def test_make_sample_sheet(self):
        self.fail()
        make_sample_sheet()

    def test_add_data_to_sheet(self):
        self.fail()
        _add_data_to_sheet()

    def test_add_metadata_to_sheet(self):
        self.fail()
        _add_metadata_to_sheet()

    def test_validate_sample_sheet_metadata(self):
        self.fail()
        _validate_sample_sheet_metadata()


class ValidateSampleSheetTests(BaseTests):
    def test_validate_sample_sheet(self):
        sheet = KLSampleSheet(self.good_ss)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            _ = validate_sample_sheet(sheet)

            # check no warnings are raised
            self.assertTrue(len(w) == 0)

    def test_validate_sample_sheet_no_sample_project(self):
        sheet = KLSampleSheet(self.no_project_ss)

        with self.assertRaisesRegex(ValueError, 'The Sample_project column '):
            validate_sample_sheet(sheet)

    def test_validate_sample_sheet_repeated_sample_ids(self):
        sheet = KLSampleSheet(self.ok_ss)

        # the sample identifiers are only repeated after scrubbing
        with self.assertRaisesRegex(ValueError, 'After scrubbing samples for '
                                    'bcl2fastq compatibility there are '
                                    'repeated identifiers. The following names'
                                    ' in the Sample_ID column are listed '
                                    'multiple times:\nCDPH-SAL_Salmonella_'
                                    'Typhi_MDL_144: 2'):
            validate_sample_sheet(sheet)

    def test_validate_sample_sheet_scrubbed_names(self):
        sheet = KLSampleSheet(self.scrubbable_ss)

        message = ('The following sample names were scrubbed for bcl2fastq '
                   'compatibility:\nCDPH-SAL_Salmonella_Typhi_MDL.143, '
                   'CDPH-SAL_Salmonella_Typhi_MDL.144, CDPH-SAL_Salmonella_'
                   'Typhi_MDL.145, CDPH-SAL_Salmonella_Typhi_MDL.146, CDPH-'
                   'SAL_Salmonella_Typhi_MDL.147, CDPH-SAL_Salmonella_Typhi'
                   '_MDL.148, CDPH-SAL_Salmonella_Typhi_MDL.149, CDPH-SAL_S'
                   'almonella_Typhi_MDL.150, CDPH-SAL_Salmonella_Typhi_MDL.'
                   '151, CDPH-SAL_Salmonella_Typhi_MDL.152, CDPH-SAL_Salmon'
                   'ella_Typhi_MDL.153, CDPH-SAL_Salmonella_Typhi_MDL.154, '
                   'CDPH-SAL_Salmonella_Typhi_MDL.155, CDPH-SAL_Salmonella_'
                   'Typhi_MDL.156, CDPH-SAL_Salmonella_Typhi_MDL.157, CDPH-'
                   'SAL_Salmonella_Typhi_MDL.158, CDPH-SAL_Salmonella_Typhi'
                   '_MDL.159, CDPH-SAL_Salmonella_Typhi_MDL.160, CDPH-SAL_S'
                   'almonella_Typhi_MDL.161, CDPH-SAL_Salmonella_Typhi_MDL.'
                   '162, CDPH-SAL_Salmonella_Typhi_MDL.163, CDPH-SAL_Salmon'
                   'ella_Typhi_MDL.164, CDPH-SAL_Salmonella_Typhi_MDL.165, '
                   'CDPH-SAL_Salmonella_Typhi_MDL.166, CDPH-SAL_Salmonella_'
                   'Typhi_MDL.167, CDPH-SAL_Salmonella_Typhi_MDL.168, P21_E'
                   '.coli ELI344, P21_E.coli ELI345, P21_E.coli ELI347, P21'
                   '_E.coli ELI348, P21_E.coli ELI349, P21_E.coli ELI350, P'
                   '21_E.coli ELI351, P21_E.coli ELI352, P21_E.coli ELI353,'
                   ' P21_E.coli ELI354, P21_E.coli ELI355, P21_E.coli ELI35'
                   '7, P21_E.coli ELI358, P21_E.coli ELI359, P21_E.coli ELI'
                   '361, P21_E.coli ELI362, P21_E.coli ELI363, P21_E.coli '
                   'ELI364, P21_E.coli ELI365, P21_E.coli ELI366, P21_E.coli '
                   'ELI367, P21_E.coli ELI368, P21_E.coli ELI369')

        with self.assertWarnsRegex(UserWarning, message):
            validate_sample_sheet(sheet)

    def test_validate_sample_sheet_bad_project_names(self):
        sheet = KLSampleSheet(self.bad_project_name_ss)

        with self.assertWarnsRegex(UserWarning, 'The following project names '
                                   'in the Sample_project column are missing a'
                                   ' Qiita study identifier: Feist, Gerwick'):
            validate_sample_sheet(sheet)

    def test_sample_sheet_to_dataframe(self):
        ss = KLSampleSheet(self.ss)
        obs = sample_sheet_to_dataframe(ss)

        columns = ['lane', 'sample_name', 'sample_plate', 'sample_well',
                   'i7_index_id', 'index', 'i5_index_id', 'index2',
                   'sample_project', 'well_description']
        index = ['sample1', 'sample2', 'sample1', 'sample2', 'sample31',
                 'sample32', 'sample34', 'sample44']

        exp = pd.DataFrame(index=index, data=DF_DATA, columns=columns)
        exp.index.name = 'sample_id'
        pd.testing.assert_frame_equal(obs, exp)


DF_DATA = [
    ['1', 'sample1', 'FooBar_666_p1', 'A1', 'iTru7_107_07', 'CCGACTAT',
     'iTru5_01_A', 'ACCGACAA', 'Baz', 'importantsample1'],
    ['1', 'sample2', 'FooBar_666_p1', 'A2', 'iTru7_107_08', 'CCGACTAT',
     'iTru5_01_A', 'CTTCGCAA', 'Baz', 'importantsample2'],
    ['3', 'sample1', 'FooBar_666_p1', 'A3', 'iTru7_107_09', 'GCCTTGTT',
     'iTru5_01_A', 'AACACCAC', 'Baz', 'importantsample1'],
    ['3', 'sample2', 'FooBar_666_p1', 'A4', 'iTru7_107_10', 'AACTTGCC',
     'iTru5_01_A', 'CGTATCTC', 'Baz', 'importantsample2'],
    ['3', 'sample31', 'FooBar_666_p1', 'A5', 'iTru7_107_11', 'CAATGTGG',
     'iTru5_01_A', 'GGTACGAA', 'FooBar_666', 'importantsample31'],
    ['3', 'sample32', 'FooBar_666_p1', 'B6', 'iTru7_107_12', 'AAGGCTGA',
     'iTru5_01_A', 'CGATCGAT', 'FooBar_666', 'importantsample32'],
    ['3', 'sample34', 'FooBar_666_p1', 'B8', 'iTru7_107_13', 'TTACCGAG',
     'iTru5_01_A', 'AAGACACC', 'FooBar_666', 'importantsample34'],
    ['3', 'sample44', 'Baz_p3', 'B99', 'iTru7_107_14', 'GTCCTAAG',
     'iTru5_01_A', 'CATCTGCT', 'Baz', 'importantsample44']]


if __name__ == '__main__':
    unittest.main()
