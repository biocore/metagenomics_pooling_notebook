import sys
import unittest
import tempfile
from datetime import datetime
from os.path import join, dirname

import pandas as pd
import sample_sheet

from metapool.sample_sheet import (KLSampleSheet, AmpliconSampleSheet,
                                   MetagenomicSampleSheetv100,
                                   MetagenomicSampleSheetv90,
                                   MetatranscriptomicSampleSheet,
                                   AbsQuantSampleSheetv10,
                                   sample_sheet_to_dataframe,
                                   make_sample_sheet,
                                   demux_sample_sheet, sheet_needs_demuxing)
from metapool.plate import ErrorMessage, WarningMessage


# The classes below share the same filepaths, so we use this dummy class
class BaseTests(unittest.TestCase):
    def setUp(self):
        data_dir = join(dirname(__file__), 'data')
        self.ss = join(data_dir, 'runs', '191103_D32611_0365_G00DHB5YXX',
                       'sample-sheet.csv')

        self.alt_ss = join(data_dir,
                           'good-sample-sheet-with-alt-col-names.csv')

        self.good_ss = join(data_dir, 'good-sample-sheet.csv')
        self.with_comments = join(data_dir, 'good-sample-sheet-but-'
                                            'with-comments.csv')

        fp = 'good-sample-sheet-with-comments-and-new-lines.csv'
        self.with_comments_and_new_lines = join(data_dir, fp)

        self.with_new_lines = join(data_dir, 'good-sample-sheet-with-'
                                             'new-lines.csv')

        self.no_project_ss = join(data_dir, 'no-project-name-sample-sheet.csv')

        # "valid" upfront but will have repeated values after scrubbing
        self.ok_ss = join(data_dir, 'ok-sample-sheet.csv')

        self.scrubbable_ss = join(data_dir, 'scrubbable-sample-sheet.csv')

        self.bad_project_name_ss = join(data_dir,
                                        'bad-project-name-sample-sheet.csv')

        bfx = [
            {
             'Sample_Project': 'Koening_ITS_101',
             'QiitaID': '101',
             'BarcodesAreRC': 'False',
             'ForwardAdapter': 'GATACA',
             'ReverseAdapter': 'CATCAT',
             'HumanFiltering': 'False',
             'library_construction_protocol': 'Knight Lab Kapa HP',
             'experiment_design_description': 'Eqiiperiment'
            },
            {
             'Sample_Project': 'Yanomani_2008_10052',
             'QiitaID': '10052',
             'BarcodesAreRC': 'False',
             'ForwardAdapter': 'GATACA',
             'ReverseAdapter': 'CATCAT',
             'HumanFiltering': 'False',
             'library_construction_protocol': 'Knight Lab Kapa HP',
             'experiment_design_description': 'Eqiiperiment'
            }
        ]

        contact = [
            {
             'Sample_Project': 'Koening_ITS_101',
             'Email': 'yoshiki@compy.com,ilike@turtles.com'
            },
            {
             'Sample_Project': 'Yanomani_2008_10052',
             'Email': 'mgdb@gmail.com'
            }
        ]

        self.md_ampl = {
            'Bioinformatics': bfx,
            'Contact': contact,
            'Assay': 'TruSeq HT',
            'SheetType': 'standard_metag',
            'SheetVersion': '0'
        }

        self.md_metag = {
            'Bioinformatics': bfx,
            'Contact': contact,
            'Assay': 'Metagenomic',
            'SheetType': 'standard_metag',
            'SheetVersion': '100'
        }


class KLSampleSheetTests(BaseTests):
    def test_instantiation(self):
        # base class can no longer be instantiated
        with self.assertRaises(TypeError, msg="TypeError: only children of "
                                              "'KLSampleSheet' may be insta"
                                              "ntiated"):
            KLSampleSheet()

        # child class should instantiate successfully.
        sheet = MetagenomicSampleSheetv100(self.good_ss)
        self.assertIsNotNone(sheet)

    def test_sample_sheet_roundtripping(self):
        # testing with all the sheets we have access to
        sheets = [self.ss, self.good_ss,
                  self.no_project_ss, self.ok_ss,
                  self.scrubbable_ss, self.bad_project_name_ss,
                  self.with_comments, self.with_comments_and_new_lines,
                  self.with_new_lines]
        sheets = {sheet: MetagenomicSampleSheetv100(sheet) for sheet in sheets}

        for filename, sheet in sheets.items():
            # write each KLSampleSheet object out to disk and compare the text
            # against the original.
            with tempfile.NamedTemporaryFile('w+') as tmp:
                sheet.write(tmp)
                tmp.seek(0)
                observed = tmp.read()

                # the following sample-sheets are identical to self.good_ss,
                # except for comments and/or empty lines. For these files,
                # observed needs to be compared to self.good_ss, since
                # comments and empty lines are ignored by metapool.
                if filename in {self.with_comments,
                                self.with_new_lines,
                                self.with_comments_and_new_lines}:
                    expected = self.good_ss
                else:
                    expected = filename

                with open(expected) as f:
                    # if the assertion fails, metapool is not processing
                    # filename as intended.
                    print({filename})
                    self.assertEqual(observed.split(),
                                     f.read().split(),
                                     f'Problem found with {filename}')

    def test_empty_write(self):
        exp = [
            '[Header],',
            ',',
            '[Reads],',
            ',',
            '[Settings],',
            ',',
            '[Data],',
            ',',
            ',',
            '[Bioinformatics],',
            ',',
            '[Contact],',
            ',',
            '']

        empty = MetagenomicSampleSheetv100()
        with tempfile.NamedTemporaryFile('w+') as tmp:
            empty.write(tmp)
            tmp.seek(0)
            observed = tmp.read()

            self.assertEqual(observed.split('\n'), exp)

    def test_empty_read(self):
        empty = [
            '[Header],',
            ',',
            '[Reads],',
            ',',
            '[Settings],',
            ',',
            '[Data],',
            ',',
            ',',
            '[Bioinformatics],',
            ',',
            '[Contact],',
            ',']

        with tempfile.NamedTemporaryFile('w+') as tmp:
            for line in empty:
                tmp.write(line + '\n')

            sheet = MetagenomicSampleSheetv100(tmp.name)

            self.assertEqual(sheet.samples, [])
            self.assertEqual(sheet.Settings, {})
            self.assertEqual(sheet.Header, {})
            self.assertEqual(sheet.Reads, [])
            self.assertIsNone(sheet.Bioinformatics)
            self.assertIsNone(sheet.Contact)

    def test_parse(self):
        sheet = MetagenomicSampleSheetv90(self.ss)

        exp = {
            'IEMFileVersion': '4',
            'SheetType': 'standard_metag',
            'SheetVersion': '100',
            'Investigator Name': 'Caballero',
            'Experiment Name': 'RKL0042',
            'Date': '2/26/20',
            'Workflow': 'GenerateFASTQ',
            'Application': 'FASTQ Only',
            'Assay': 'Metagenomic',
            'Description': '',
            'Chemistry': 'Default'
        }

        self.assertEqual(sheet.Header, exp)
        self.assertEqual(sheet.Reads, [150, 150])
        self.assertEqual(sheet.Settings, {'ReverseComplement': '0'})

        data = (
            '1,sample_1,sample.1,FooBar_666_p1,A1,iTru7_107_07,CCGACTAT,'
            'iTru5_01_A,ACCGACAA,Baz_12345,pool1,importantsample1,'
            'KnightLabKapaHP,Eqiiperiment\n'
            '1,sample_2,sample.2,FooBar_666_p1,A2,iTru7_107_08,CCGACTAT,'
            'iTru5_01_A,CTTCGCAA,Baz_12345,pool1,importantsample2,'
            'KnightLabKapaHP,Eqiiperiment\n'
            '3,sample_1,sample.1,FooBar_666_p1,A3,iTru7_107_09,GCCTTGTT,'
            'iTru5_01_A,AACACCAC,Baz_12345,pool1,importantsample1,'
            'KnightLabKapaHP,Eqiiperiment\n'
            '3,sample_2,sample.2,FooBar_666_p1,A4,iTru7_107_10,AACTTGCC,'
            'iTru5_01_A,CGTATCTC,Baz_12345,pool1,importantsample2,'
            'KnightLabKapaHP,Eqiiperiment\n'
            '3,sample_31,sample.31,FooBar_666_p1,A5,iTru7_107_11,CAATGTGG,'
            'iTru5_01_A,GGTACGAA,FooBar_666,pool1,importantsample31,'
            'KnightLabKapaHP,SomethingWitty\n'
            '3,sample_32,sample.32,FooBar_666_p1,B6,iTru7_107_12,AAGGCTGA,'
            'iTru5_01_A,CGATCGAT,FooBar_666,pool1,importantsample32,'
            'KnightLabKapaHP,SomethingWitty\n'
            '3,sample_34,sample.34,FooBar_666_p1,B8,iTru7_107_13,TTACCGAG,'
            'iTru5_01_A,AAGACACC,FooBar_666,pool1,importantsample34,'
            'KnightLabKapaHP,SomethingWitty\n'
            '3,sample_44,sample.44,Baz_p3,B99,iTru7_107_14,GTCCTAAG,'
            'iTru5_01_A,CATCTGCT,Baz_12345,pool1,importantsample44,'
            'KnightLabKapaHP,Eqiiperiment\n'
        )
        keys = ['Lane', 'Sample_ID', 'Sample_Name', 'Sample_Plate',
                'well_id_384', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                'Sample_Project', 'syndna_pool_number', 'Well_description',
                'library_construction_protocol', 'experiment_design_protocol']

        for sample, line in zip(sheet.samples, data.split()):
            values = line.strip().split(',')
            exp = sample_sheet.Sample(dict(zip(keys, values)))
            self.assertEqual(sample, exp)

        # check for Bioinformatics
        exp = pd.DataFrame(
            columns=['Sample_Project', 'QiitaID', 'BarcodesAreRC',
                     'ForwardAdapter', 'ReverseAdapter', 'HumanFiltering',
                     'library_construction_protocol',
                     'experiment_design_description', 'contains_replicates'],
            data=[
                ['Baz_12345', '100', 'False', 'AACC', 'GGTT', 'False',
                 'Knight Lab Kapa HP', 'Eqiiperiment', 'False'],
                ['FooBar_666', '666', 'False', 'AACC', 'GGTT', 'False',
                 'Knight Lab Kapa HP', 'SomethingWitty', 'False']
            ]
        )

        pd.testing.assert_frame_equal(sheet.Bioinformatics, exp)

        # check for Contact
        exp = pd.DataFrame(
            columns=['Email', 'Sample_Project'],
            data=[
                ['test@lol.com', 'Baz_12345'],
                ['tester@rofl.com', 'FooBar_666']
            ]
        )
        pd.testing.assert_frame_equal(sheet.Contact, exp)

    def test_parse_with_comments(self):
        # the two sample sheets are identical except for the comments
        exp = MetagenomicSampleSheetv100(self.good_ss)
        with self.assertWarnsRegex(UserWarning, 'Comments at the beginning '):
            obs = MetagenomicSampleSheetv100(self.with_comments)

            self.assertEqual(obs.Header, exp.Header)
            self.assertEqual(obs.Settings, exp.Settings)
            self.assertEqual(obs.Reads, exp.Reads)

            for o_sample, e_sample in zip(obs.samples, exp.samples):
                self.assertEqual(o_sample, e_sample)

            pd.testing.assert_frame_equal(obs.Contact, exp.Contact)
            pd.testing.assert_frame_equal(obs.Bioinformatics,
                                          exp.Bioinformatics)

            self.assertEqual(len(obs), 783)

    def test_merge(self):
        base = MetagenomicSampleSheetv100()
        base.Reads = [151, 151]
        base.add_sample(sample_sheet.Sample({
            'Sample_ID': 'y',
            'index': 'GGTACA',
            'index2': 'GGCGCC',
            'Sample_Name': 'y.sample'
        }))
        base.Contact = pd.DataFrame(self.md_metag['Contact'])

        hugo = MetagenomicSampleSheetv100()
        hugo.Reads = [151, 151]
        hugo.add_sample(sample_sheet.Sample({
            'Sample_ID': 'a',
            'index': 'GATACA',
            'index2': 'GCCGCC',
            'Sample_Name': 'a.sample'
        }))
        hugo.Contact = pd.DataFrame(self.md_metag['Contact'])

        paco = MetagenomicSampleSheetv100()
        paco.Reads = [151, 151]
        paco.add_sample(sample_sheet.Sample({
            'Sample_ID': 'b',
            'index': 'GATAAA',
            'index2': 'GCCACC',
            'Sample_Name': 'b.sample'
        }))

        luis = MetagenomicSampleSheetv100()
        luis.Reads = [151, 151]
        luis.add_sample(sample_sheet.Sample({
                'Sample_ID': 'c',
                'index': 'GATATA',
                'index2': 'GCCTCC',
                'Sample_Name': 'c.sample'}))

        base.merge([hugo, paco, luis])

        self.assertEqual(base.Reads, [151, 151])

        exp_samples = [
            sample_sheet.Sample({
                'Sample_ID': 'y',
                'index': 'GGTACA',
                'index2': 'GGCGCC',
                'Sample_Name': 'y.sample'}
            ),
            sample_sheet.Sample({
                'Sample_ID': 'a',
                'index': 'GATACA',
                'index2': 'GCCGCC',
                'Sample_Name': 'a.sample'}
            ),
            sample_sheet.Sample({
                'Sample_ID': 'b',
                'index': 'GATAAA',
                'index2': 'GCCACC',
                'Sample_Name': 'b.sample'}
            ),
            sample_sheet.Sample({
                'Sample_ID': 'c',
                'index': 'GATATA',
                'index2': 'GCCTCC',
                'Sample_Name': 'c.sample'}
            ),
        ]

        for obs, exp in zip(base.samples, exp_samples):
            self.assertEqual(obs, exp)

        # checks the items haven't been repeated
        contact = self.md_metag['Contact']
        pd.testing.assert_frame_equal(base.Contact, pd.DataFrame(contact))

    def test_merge_bioinformatics(self):
        base = MetagenomicSampleSheetv100()
        base.Reads = [151, 151]
        base.add_sample(sample_sheet.Sample({
            'Sample_ID': 'y',
            'index': 'GGTACA',
            'index2': 'GGCGCC',
            'Sample_Name': 'y.sample'
        }))
        base.Bioinformatics = pd.DataFrame(self.md_metag['Bioinformatics'])

        hugo = MetagenomicSampleSheetv100()
        hugo.Reads = [151, 151]
        hugo.add_sample(sample_sheet.Sample({
            'Sample_ID': 'a',
            'index': 'GATACA',
            'index2': 'GCCGCC',
            'Sample_Name': 'a.sample'
        }))
        hugo.Bioinformatics = pd.DataFrame(self.md_metag['Bioinformatics'])

        paco = MetagenomicSampleSheetv100()
        paco.Reads = [151, 151]
        paco.add_sample(sample_sheet.Sample({
            'Sample_ID': 'b',
            'index': 'GATAAA',
            'index2': 'GCCACC',
            'Sample_Name': 'b.sample'
        }))
        paco.Bioinformatics = pd.DataFrame(self.md_metag['Bioinformatics'])
        paco.Bioinformatics['Sample_Project'] = (
                'paco_' + paco.Bioinformatics['Sample_Project'])

        base.merge([hugo, paco])

        self.assertEqual(base.Reads, [151, 151])

        exp_samples = [
            sample_sheet.Sample({
                'Sample_ID': 'y',
                'index': 'GGTACA',
                'index2': 'GGCGCC',
                'Sample_Name': 'y.sample'}
            ),
            sample_sheet.Sample({
                'Sample_ID': 'a',
                'index': 'GATACA',
                'index2': 'GCCGCC',
                'Sample_Name': 'a.sample'}
            ),
            sample_sheet.Sample({
                'Sample_ID': 'b',
                'index': 'GATAAA',
                'index2': 'GCCACC',
                'Sample_Name': 'b.sample'}
            ),
        ]

        for obs, exp in zip(base.samples, exp_samples):
            self.assertEqual(obs, exp)

        self.assertIsNone(base.Contact)

        # check for Bioinformatics
        exp = pd.DataFrame(
            columns=['Sample_Project', 'QiitaID', 'BarcodesAreRC',
                     'ForwardAdapter', 'ReverseAdapter', 'HumanFiltering',
                     'library_construction_protocol',
                     'experiment_design_description'],
            data=[
                ['Koening_ITS_101', '101', 'False', 'GATACA', 'CATCAT',
                 'False', 'Knight Lab Kapa HP', 'Eqiiperiment'],
                ['Yanomani_2008_10052', '10052', 'False', 'GATACA', 'CATCAT',
                 'False', 'Knight Lab Kapa HP', 'Eqiiperiment'],
                ['paco_Koening_ITS_101', '101', 'False', 'GATACA', 'CATCAT',
                 'False', 'Knight Lab Kapa HP', 'Eqiiperiment'],
                ['paco_Yanomani_2008_10052', '10052', 'False', 'GATACA',
                 'CATCAT', 'False', 'Knight Lab Kapa HP', 'Eqiiperiment']
            ]
        )

        # checks the items haven't been repeated
        pd.testing.assert_frame_equal(base.Bioinformatics, exp)

    def test_merge_error(self):
        base = MetagenomicSampleSheetv100()
        base.Reads = [151, 151]
        base.Settings = {'ReverseComplement': 0,
                         'SomethingElse': '100'}

        hugo = MetagenomicSampleSheetv100()
        hugo.add_sample(sample_sheet.Sample({
            'Sample_ID': 'a',
            'index': 'GATACA',
            'index2': 'GCCGCC',
            'Sample_Name': 'a.sample'
        }))

        with self.assertRaisesRegex(ValueError, 'The Settings section is '
                                    'different for sample sheet 1'):
            base.merge([hugo])

    def test_merge_different_dates(self):
        base = MetagenomicSampleSheetv100()
        base.Header['Date'] = '08-01-1989'
        base.Settings = {'ReverseComplement': 0}

        hugo = MetagenomicSampleSheetv100()
        hugo.Header['Date'] = '04-26-2021'
        hugo.Settings = {'ReverseComplement': 0}

        hugo.add_sample(sample_sheet.Sample({
            'Sample_ID': 'a',
            'index': 'GATACA',
            'index2': 'GCCGCC',
            'Sample_Name': 'a.sample'
        }))

        base.merge([hugo])

        # keeps base's date
        self.assertEqual(dict(base.Header), {'Date': '08-01-1989'})

        # there should only be one sample
        self.assertEqual(len(base.samples), 1)
        self.assertEqual(base.samples[0],
                         sample_sheet.Sample({'Sample_ID': 'a',
                                              'index': 'GATACA',
                                              'index2': 'GCCGCC',
                                              'Sample_Name': 'a.sample'}))

    def test_validate(self):
        sheet = AmpliconSampleSheet()
        obs = sheet._validate_sample_sheet_metadata(self.md_ampl)
        self.assertEqual(obs, [])

    def test_more_attributes(self):
        sheet = AmpliconSampleSheet()
        self.md_ampl['Ride'] = 'the lightning'

        obs = sheet._validate_sample_sheet_metadata(self.md_ampl)
        exp = [ErrorMessage('These metadata keys are not supported: Ride')]
        self.assertEqual(obs, exp)

    def test_validate_missing_assay(self):
        sheet = AmpliconSampleSheet()
        self.md_ampl['Assay'] = 'NewAssayType'

        obs = sheet._validate_sample_sheet_metadata(self.md_ampl)
        exp = [ErrorMessage('NewAssayType is not a supported Assay')]
        self.assertEqual(obs, exp)

    def test_validate_missing_bioinformatics_data(self):
        sheet = AmpliconSampleSheet()
        del self.md_ampl['Bioinformatics']

        obs = sheet._validate_sample_sheet_metadata(self.md_ampl)
        exp = [ErrorMessage('Bioinformatics is a required attribute')]
        self.assertEqual(obs, exp)

    def test_validate_missing_column_in_bioinformatics(self):
        sheet = AmpliconSampleSheet()
        del self.md_ampl['Bioinformatics'][0]['Sample_Project']
        exp = [ErrorMessage('In the Bioinformatics section Project #1 does not'
                            ' have exactly these keys BarcodesAreRC, '
                            'ForwardAdapter, HumanFiltering, QiitaID, '
                            'ReverseAdapter, Sample_Project, '
                            'experiment_design_description, '
                            'library_construction_protocol')]
        obs = sheet._validate_sample_sheet_metadata(self.md_ampl)
        self.assertEqual(str(obs[0]), str(exp[0]))

    def test_alt_sample_sheet(self):
        # testing with all the sheets we have access to
        obs = MetagenomicSampleSheetv90(self.alt_ss).all_sample_keys

        exp = ['Lane',
               'Sample_ID',
               'Sample_Name',
               'Sample_Plate',
               'well_id_384',
               'I7_Index_ID',
               'index',
               'I5_Index_ID',
               'index2',
               'Sample_Project',
               'syndna_pool_number',
               'Well_description']

        self.assertEqual(set(obs), set(exp))


class SampleSheetWorkflow(BaseTests):
    def setUp(self):
        super().setUp()

        self.sheet = AmpliconSampleSheet()
        self.sheet.Header['IEM4FileVersion'] = 4
        self.sheet.Header['Investigator Name'] = 'Knight'
        self.sheet.Header['Experiment Name'] = 'RKO_experiment'
        self.sheet.Header['Date'] = '2021-08-17'
        self.sheet.Header['Workflow'] = 'GenerateFASTQ'
        self.sheet.Header['Application'] = 'FASTQ Only'
        self.sheet.Header['Assay'] = 'TruSeq HT'
        self.sheet.Header['Description'] = ''
        self.sheet.Header['Chemistry'] = 'Default'
        self.sheet.Reads = [151, 151]
        self.sheet.Settings['ReverseComplement'] = 0

        self.sheet.Bioinformatics = pd.DataFrame(
            columns=['Sample_Project', 'QiitaID', 'BarcodesAreRC',
                     'ForwardAdapter', 'ReverseAdapter', 'HumanFiltering'],
            data=[
                ['THDMI_10317', '10317', 'False', 'AACC', 'GGTT', 'False']
            ]
        )

        # check for Contact
        self.sheet.Contact = pd.DataFrame(
            columns=['Email', 'Sample_Project'],
            data=[
                ['daniel@tmi.com', 'THDMI_10317'],
            ]
        )

        data = [
            ['X00180471',
             'X00180471', 'A', 1, False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'A1', '1', '1', 'SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_2', 'THDMI UK', '', '1', 'A1',
             '515rcbc0', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'AGCCTTCGTCGC',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             'AATGATACGGCGACCACCGAGATCTACACGCTAGCCTTCGTCGCTATGGTAATTGTGTGYCAG'
             'CMGCCGCGGTAA', 'pool1'],
            ['X00180199',
             'X00180199', 'C', 1, False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'C1', '1', '1', 'SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_2', 'THDMI UK', '', '1', 'B1',
             '515rcbc12', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'CGTATAAATGCG',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             'AATGATACGGCGACCACCGAGATCTACACGCTCGTATAAATGCGTATGGTAATTGTGTGYCAG'
             'CMGCCGCGGTAA', 'pool1'],
            ['X00179789',
             'X00179789', 'E', 1, False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'E1', '1', '1', 'SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_2', 'THDMI UK', '', '1', 'C1',
             '515rcbc24', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'TGACTAATGGCC',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             'AATGATACGGCGACCACCGAGATCTACACGCTTGACTAATGGCCTATGGTAATTGTGTGYCAG'
             'CMGCCGCGGTAA', 'pool1'],
        ]

        self.table = pd.DataFrame(
            columns=['sample sheet Sample_ID',
                     'Sample', 'Row', 'Col', 'Blank', 'Project Plate',
                     'Project Name', 'Compressed Plate Name', 'Well',
                     'Plate Position', 'Primer Plate #', 'Plating',
                     'Extraction Kit Lot', 'Extraction Robot', 'TM1000 8 Tool',
                     'Primer Date', 'MasterMix Lot', 'Water Lot',
                     'Processing Robot', 'Sample Plate', 'Project_Name',
                     'Original Name', 'Plate', 'EMP Primer Plate Well', 'Name',
                     "Illumina 5' Adapter", 'Golay Barcode',
                     'Forward Primer Pad', 'Forward Primer Linker',
                     '515FB Forward Primer (Parada)', 'Primer For PCR',
                     'syndna_pool_number'],
            data=data
        )

    def test_validate_sample_sheet_metadata_empty(self):
        sheet = AmpliconSampleSheet()
        messages = sheet._validate_sample_sheet_metadata({})

        exp = [
            ErrorMessage('Assay is a required attribute'),
            ErrorMessage('Bioinformatics is a required attribute'),
            ErrorMessage('Contact is a required attribute'),
        ]

        self.assertEqual(messages, exp)

    def test_validate_sample_sheet_metadata_not_supported(self):
        sheet = AmpliconSampleSheet()
        self.md_ampl['Rush'] = 'XYZ'
        messages = sheet._validate_sample_sheet_metadata(self.md_ampl)

        exp = [
                ErrorMessage('These metadata keys are not supported: Rush'),
        ]

        self.assertEqual(messages, exp)

    def test_validate_sample_sheet_metadata_good(self):
        # self.md_ampl is patterned after legacy amplicon sample-sheet.
        sheet = AmpliconSampleSheet()
        messages = sheet._validate_sample_sheet_metadata(self.md_ampl)
        self.assertEqual(messages, [])

        # test _validate_sample_sheet_metadata() against a
        # MetagenomicSampleSheetv100 object which defines an extra column
        # (contains_replicates) in the Bioinformatics section. Since
        # self.metadata does not contain this extra column, ErrorMessage()s
        # should be returned saying as much.
        sheet = MetagenomicSampleSheetv100()
        messages = sheet._validate_sample_sheet_metadata(self.md_metag)

        exp_msgs = ['In the Bioinformatics section Project #1 does not have '
                    'exactly these keys BarcodesAreRC, ForwardAdapter, Human'
                    'Filtering, QiitaID, ReverseAdapter, Sample_Project, '
                    'contains_replicates, experiment_design_description, '
                    'library_construction_protocol',
                    'In the Bioinformatics section Project #2 does not have '
                    'exactly these keys BarcodesAreRC, ForwardAdapter, Human'
                    'Filtering, QiitaID, ReverseAdapter, Sample_Project, '
                    'contains_replicates, experiment_design_description, '
                    'library_construction_protocol']

        self.assertEqual(messages[0].message, exp_msgs[0])
        self.assertEqual(messages[1].message, exp_msgs[1])

    def test_validate_sample_sheet_metadata_bad_assay_types(self):
        sheet = AmpliconSampleSheet()

        invalid_types = ['SomeType', 'Metagenomics', 'Metatranscriptomics']

        for invalid_type in invalid_types:
            self.md_ampl['Assay'] = invalid_type
            messages = sheet._validate_sample_sheet_metadata(self.md_ampl)
            exp = f'ErrorMessage: {invalid_type} is not a supported Assay'
            self.assertEqual(str(messages[0]), exp)

    def test_make_sample_sheet(self):
        exp_bfx = pd.DataFrame(self.md_ampl['Bioinformatics'])
        exp_contact = pd.DataFrame(self.md_ampl['Contact'])

        # for amplicon we expect the following three columns to not be there
        message = (r'The column (I5_Index_ID|index2|Well_description) '
                   r'in the sample sheet is empty')
        with self.assertWarnsRegex(UserWarning, message):
            obs = make_sample_sheet(self.md_ampl, self.table, 'HiSeq4000',
                                    [5, 7], strict=False)

        self.assertIsInstance(obs, AmpliconSampleSheet)

        self.assertEqual(obs.Reads, [151, 151])
        self.assertEqual(obs.Settings, {'ReverseComplement': '0'})

        pd.testing.assert_frame_equal(obs.Bioinformatics, exp_bfx)
        pd.testing.assert_frame_equal(obs.Contact, exp_contact)

        header = {
            'IEMFileVersion': '4',
            'SheetType': 'standard_metag',
            'SheetVersion': '0',
            'Date': datetime.today().strftime('%Y-%m-%d'),
            'Workflow': 'GenerateFASTQ',
            'Application': 'FASTQ Only',
            'Assay': 'TruSeq HT',
            'Description': '',
            'Chemistry': 'Default',
        }

        self.assertEqual(obs.Header, header)
        self.assertEqual(len(obs.samples), 6)

        data = (
            [5, 'X00180471', 'X00180471', 'THDMI_10317_PUK2', 'A1', '515rcbc0',
             'AGCCTTCGTCGC', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00180471.A1'],
            [5, 'X00180199', 'X00180199', 'THDMI_10317_PUK2', 'C1',
             '515rcbc12', 'CGTATAAATGCG', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00180199.C1'],
            [5, 'X00179789', 'X00179789', 'THDMI_10317_PUK2', 'E1',
             '515rcbc24', 'TGACTAATGGCC', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00179789.E1'],
            [7, 'X00180471', 'X00180471', 'THDMI_10317_PUK2', 'A1', '515rcbc0',
             'AGCCTTCGTCGC', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00180471.A1'],
            [7, 'X00180199', 'X00180199', 'THDMI_10317_PUK2', 'C1',
             '515rcbc12', 'CGTATAAATGCG', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00180199.C1'],
            [7, 'X00179789', 'X00179789', 'THDMI_10317_PUK2', 'E1',
             '515rcbc24', 'TGACTAATGGCC', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00179789.E1'],
        )
        keys = ['Lane', 'Sample_ID', 'Sample_Name', 'Sample_Plate',
                'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                'Sample_Project', 'Well_description']

        for sample, row in zip(obs.samples, data):
            exp = sample_sheet.Sample(dict(zip(keys, row)))
            self.assertEqual(dict(sample), dict(exp))

    def test_column_alternatives(self):
        # confirm standard 'Well_description' column name behaved as intended.
        table2 = self.table.copy(deep=True)
        table2['Well_description'] = ['Row A', 'Row B', 'Row C']

        # allow 'Well_description' column to pass through to obs.
        obs = make_sample_sheet(self.md_ampl,
                                table2,
                                'HiSeq4000',
                                [5, 7],
                                strict=False)

        self.assertIsNotNone(obs, msg="make_sample_sheet() failed")
        self.assertIsInstance(obs, AmpliconSampleSheet)

        data = (
            [5, 'X00180471', 'X00180471', 'THDMI_10317_PUK2', 'A1', '515rcbc0',
             'AGCCTTCGTCGC', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00180471.A1'],
            [5, 'X00180199', 'X00180199', 'THDMI_10317_PUK2', 'C1',
             '515rcbc12', 'CGTATAAATGCG', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00180199.C1'],
            [5, 'X00179789', 'X00179789', 'THDMI_10317_PUK2', 'E1',
             '515rcbc24', 'TGACTAATGGCC', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00179789.E1'],
            [7, 'X00180471', 'X00180471', 'THDMI_10317_PUK2', 'A1', '515rcbc0',
             'AGCCTTCGTCGC', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00180471.A1'],
            [7, 'X00180199', 'X00180199', 'THDMI_10317_PUK2', 'C1',
             '515rcbc12', 'CGTATAAATGCG', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00180199.C1'],
            [7, 'X00179789', 'X00179789', 'THDMI_10317_PUK2', 'E1',
             '515rcbc24', 'TGACTAATGGCC', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00179789.E1'],
        )
        keys = ['Lane', 'Sample_ID', 'Sample_Name', 'Sample_Plate',
                'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                'Sample_Project', 'Well_description']

        for sample, row in zip(obs.samples, data):
            exp = sample_sheet.Sample(dict(zip(keys, row)))
            self.assertEqual(dict(sample), dict(exp))

        # Try making sample-sheet w/an alternate column name and confirm that
        # the results continue to be as expected.
        table2.rename({'Well_description': 'well_description'},
                      axis=1, inplace=True)

        obs = make_sample_sheet(self.md_ampl,
                                table2,
                                'HiSeq4000',
                                [5, 7],
                                strict=False)

        for sample, row in zip(obs.samples, data):
            exp = sample_sheet.Sample(dict(zip(keys, row)))
            self.assertEqual(dict(sample), dict(exp))

        # Try w/another alternate column name
        table2.rename({'well_description': 'description'},
                      axis=1, inplace=True)

        obs = make_sample_sheet(self.md_ampl,
                                table2,
                                'HiSeq4000',
                                [5, 7],
                                strict=False)

        for sample, row in zip(obs.samples, data):
            exp = sample_sheet.Sample(dict(zip(keys, row)))
            self.assertEqual(dict(sample), dict(exp))

    def test_remap_table_amplicon(self):
        columns = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well',
                   'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                   'Sample_Project', 'Well_description']

        data = [
            ['X00180471', 'X00180471', 'THDMI_10317_PUK2', 'A1', '515rcbc0',
             'AGCCTTCGTCGC', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00180471.A1'],
            ['X00180199', 'X00180199', 'THDMI_10317_PUK2', 'C1', '515rcbc12',
             'CGTATAAATGCG', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00180199.C1'],
            ['X00179789', 'X00179789', 'THDMI_10317_PUK2', 'E1', '515rcbc24',
             'TGACTAATGGCC', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00179789.E1'],
        ]

        exp = pd.DataFrame(columns=columns, data=data)

        # for amplicon we expect the following three columns to not be there.
        message = (r'The column (I5_Index_ID|index2) '
                   r'in the sample sheet is empty')
        with self.assertWarnsRegex(UserWarning, message):
            # because obs is generated from self.table (a pre-prep df), we
            # expect 'Well_description' to be empty since it is created and
            # populated before _remap_table() is called.
            sheet = AmpliconSampleSheet()
            obs = sheet._remap_table(self.table, strict=False)
            self.assertEqual(len(obs), 3)
            pd.testing.assert_frame_equal(obs, exp, check_like=True)

    def test_remap_table_metagenomics(self):
        data = [
            ['33-A1', 'A', 1, True, 'A1', 0, 0, 'AACGCACACTCGTCTT',
             'iTru5_19_A', 'AACGCACA', 'A1', 'iTru5_plate', 'iTru7_109_01',
             'CTCGTCTT', 'A22', 'iTru7_plate', '33-A1', 'pool1',
             'The_plate.33-A1.A1'],
            ['820072905-2', 'C', 1, False, 'C1', 1, 1, 'ATGCCTAGCGAACTGT',
             'iTru5_19_B', 'ATGCCTAG', 'B1', 'iTru5_plate', 'iTru7_109_02',
             'CGAACTGT', 'B22', 'iTru7_plate', '820072905-2',
             'pool1', 'The_plate.820072905-2.C1'],
            ['820029517-3', 'E', 1, False, 'E1', 2, 2, 'CATACGGACATTCGGT',
             'iTru5_19_C', 'CATACGGA', 'C1', 'iTru5_plate', 'iTru7_109_03',
             'CATTCGGT', 'C22', 'iTru7_plate', '820029517-3',
             'pool1', 'The_plate.820029517-3.E1']
        ]
        columns = ['Sample', 'Row', 'Col', 'Blank', 'Well', 'index',
                   'index combo', 'index combo seq', 'i5 name', 'i5 sequence',
                   'i5 well', 'i5 plate', 'i7 name', 'i7 sequence', 'i7 well',
                   'i7 plate', 'sample sheet Sample_ID', 'syndna_pool_number',
                   'Well_description']
        self.table = pd.DataFrame(data=data, columns=columns)
        self.table['Project Name'] = 'Tst_project_1234'
        self.table['Project Plate'] = 'The_plate'

        columns = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                   'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                   'Sample_Project', 'Well_description']
        data = [
            ['33-A1', '33-A1', 'The_plate', 'A1', 'iTru7_109_01',
             'CTCGTCTT', 'iTru5_19_A', 'AACGCACA', 'Tst_project_1234',
             'The_plate.33-A1.A1'],
            ['820072905-2', '820072905-2', 'The_plate', 'C1', 'iTru7_109_02',
             'CGAACTGT', 'iTru5_19_B', 'ATGCCTAG', 'Tst_project_1234',
             'The_plate.820072905-2.C1'],
            ['820029517-3', '820029517-3', 'The_plate', 'E1', 'iTru7_109_03',
             'CATTCGGT', 'iTru5_19_C', 'CATACGGA', 'Tst_project_1234',
             'The_plate.820029517-3.E1'],
        ]

        exp = pd.DataFrame(columns=columns, data=data)

        sheet = MetagenomicSampleSheetv100()

        obs = sheet._remap_table(self.table, strict=False)

        self.assertEqual(len(obs), 3)

        pd.testing.assert_frame_equal(obs, exp, check_like=True)

    def test_remap_table_metatranscriptomics(self):
        # note that Well_description is now included because it's expected
        # to be inserted by the function that calls _remap_table().
        data = [
            ['33-A1', 'A', 1, True, 'A1', 0, 0, 'AACGCACACTCGTCTT',
             'iTru5_19_A', 'AACGCACA', 'A1', 'iTru5_plate', 'iTru7_109_01',
             'CTCGTCTT', 'A22', 'iTru7_plate', '33-A1', 'The_plate.33-A1.A1'],
            ['820072905-2', 'C', 1, False, 'C1', 1, 1, 'ATGCCTAGCGAACTGT',
             'iTru5_19_B', 'ATGCCTAG', 'B1', 'iTru5_plate', 'iTru7_109_02',
             'CGAACTGT', 'B22', 'iTru7_plate', '820072905-2',
             'The_plate.820072905-2.C1'],
            ['820029517-3', 'E', 1, False, 'E1', 2, 2, 'CATACGGACATTCGGT',
             'iTru5_19_C', 'CATACGGA', 'C1', 'iTru5_plate', 'iTru7_109_03',
             'CATTCGGT', 'C22', 'iTru7_plate', '820029517-3',
             'The_plate.820029517-3.E1']
        ]
        columns = ['Sample', 'Row', 'Col', 'Blank', 'Well', 'index',
                   'index combo', 'index combo seq', 'i5 name', 'i5 sequence',
                   'i5 well', 'i5 plate', 'i7 name', 'i7 sequence', 'i7 well',
                   'i7 plate', 'sample sheet Sample_ID',
                   'Well_description']
        self.table = pd.DataFrame(data=data, columns=columns)
        self.table['Project Name'] = 'Tst_project_1234'
        self.table['Project Plate'] = 'The_plate'

        columns = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                   'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                   'Sample_Project', 'Well_description']
        data = [
            ['33-A1', '33-A1', 'The_plate', 'A1', 'iTru7_109_01', 'CTCGTCTT',
             'iTru5_19_A', 'AACGCACA', 'Tst_project_1234',
             'The_plate.33-A1.A1'],
            ['820072905-2', '820072905-2', 'The_plate', 'C1', 'iTru7_109_02',
             'CGAACTGT', 'iTru5_19_B', 'ATGCCTAG', 'Tst_project_1234',
             'The_plate.820072905-2.C1'],
            ['820029517-3', '820029517-3', 'The_plate', 'E1', 'iTru7_109_03',
             'CATTCGGT', 'iTru5_19_C', 'CATACGGA', 'Tst_project_1234',
             'The_plate.820029517-3.E1'],
        ]

        exp = pd.DataFrame(columns=columns, data=data)

        sheet = MetatranscriptomicSampleSheet()

        obs = sheet._remap_table(self.table, strict=False)
        obs = obs[['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                   'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                   'Sample_Project', 'Well_description']]

        self.assertEqual(len(obs), 3)
        pd.testing.assert_frame_equal(obs, exp, check_like=True)

    def test_add_data_to_sheet(self):
        # for amplicon we expect the following three columns to not be there
        message = (r'The column (I5_Index_ID|index2|Well_description) '
                   r'in the sample sheet is empty')

        with self.assertWarnsRegex(UserWarning, message):
            self.sheet._add_data_to_sheet(self.table, 'HiSeq4000', [1],
                                          'TruSeq HT', strict=False)

        self.assertEqual(len(self.sheet), 3)

        data = (
            [1, 'X00180471', 'X00180471', 'THDMI_10317_PUK2', 'A1', '515rcbc0',
             'AGCCTTCGTCGC', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00180471.A1'],
            [1, 'X00180199', 'X00180199', 'THDMI_10317_PUK2', 'C1',
             '515rcbc12', 'CGTATAAATGCG', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00180199.C1'],
            [1, 'X00179789', 'X00179789', 'THDMI_10317_PUK2', 'E1',
             '515rcbc24', 'TGACTAATGGCC', '', '', 'THDMI_10317',
             'THDMI_10317_PUK2.X00179789.E1'],
        )
        keys = ['Lane', 'Sample_ID', 'Sample_Name', 'Sample_Plate',
                'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                'Sample_Project', 'Well_description']

        for sample, row in zip(self.sheet.samples, data):
            exp = sample_sheet.Sample(dict(zip(keys, row)))
            self.assertEqual(dict(sample), dict(exp))

    def test_add_metadata_to_sheet_all_defaults_amplicon(self):
        sheet = AmpliconSampleSheet()

        self.md_ampl['Assay'] = 'TruSeq HT'
        exp_bfx = pd.DataFrame(self.md_ampl['Bioinformatics'])
        exp_contact = pd.DataFrame(self.md_ampl['Contact'])

        obs = sheet._add_metadata_to_sheet(self.md_ampl, 'HiSeq4000')

        self.assertEqual(obs.Reads, [151, 151])

        settings = {
            'ReverseComplement': '0',
        }
        self.assertEqual(obs.Settings, settings)

        pd.testing.assert_frame_equal(obs.Bioinformatics, exp_bfx)
        pd.testing.assert_frame_equal(obs.Contact, exp_contact)

        header = {
            'IEMFileVersion': '4',
            'SheetType': 'standard_metag',
            'SheetVersion': '0',
            'Date': datetime.today().strftime('%Y-%m-%d'),
            'Workflow': 'GenerateFASTQ',
            'Application': 'FASTQ Only',
            'Assay': 'TruSeq HT',
            'Description': '',
            'Chemistry': 'Default',
        }

        self.assertEqual(obs.Header, header)
        self.assertEqual(len(obs.samples), 0)

    def test_add_metadata_to_sheet_most_defaults(self):
        sheet = MetagenomicSampleSheetv100()
        exp_bfx = pd.DataFrame(self.md_metag['Bioinformatics'])
        exp_contact = pd.DataFrame(self.md_metag['Contact'])

        obs = sheet._add_metadata_to_sheet(self.md_metag,
                                           'HiSeq4000')

        self.assertEqual(obs.Reads, [151, 151])

        settings = {
            'ReverseComplement': '0',
            'MaskShortReads': '1',
            'OverrideCycles': 'Y151;I8N2;I8N2;Y151'
        }
        self.assertEqual(obs.Settings, settings)

        pd.testing.assert_frame_equal(obs.Bioinformatics, exp_bfx)
        pd.testing.assert_frame_equal(obs.Contact, exp_contact)

        header = {
            'IEMFileVersion': '4',
            'SheetType': 'standard_metag',
            'SheetVersion': '100',
            'Investigator Name': 'Knight',
            'Experiment Name': 'RKL_experiment',
            'Date': datetime.today().strftime('%Y-%m-%d'),
            'Workflow': 'GenerateFASTQ',
            'Application': 'FASTQ Only',
            'Assay': 'Metagenomic',
            'Description': '',
            'Chemistry': 'Default',
        }

        self.assertEqual(obs.Header, header)
        self.assertEqual(len(obs.samples), 0)

    def test_add_metadata_to_sheet_some_defaults(self):
        sheet = MetagenomicSampleSheetv100()

        # add a sample to make sure we can keep data around
        sheet.add_sample(sample_sheet.Sample({
            'Sample_ID': 'thy_sample',
            'Sample_Name': 'the_name_is_sample',
            'index': 'CCGACTAT',
            'index2': 'ACCGACCA',
        }))

        exp_bfx = pd.DataFrame(self.md_metag['Bioinformatics'])
        exp_contact = pd.DataFrame(self.md_metag['Contact'])
        self.md_metag['Date'] = '1970-01-01'

        obs = sheet._add_metadata_to_sheet(self.md_metag, 'HiSeq4000')

        self.assertEqual(obs.Reads, [151, 151])
        self.assertDictEqual(dict(obs.Settings),
                             {'ReverseComplement': '0',
                              'MaskShortReads': '1',
                              'OverrideCycles': 'Y151;I8N2;I8N2;Y151'})

        pd.testing.assert_frame_equal(obs.Bioinformatics, exp_bfx)
        pd.testing.assert_frame_equal(obs.Contact, exp_contact)

        header = {
            'IEMFileVersion': '4',
            'SheetType': 'standard_metag',
            'SheetVersion': '100',
            'Date': '1970-01-01',
            'Workflow': 'GenerateFASTQ',
            'Application': 'FASTQ Only',
            'Assay': 'Metagenomic',
            'Description': '',
            'Chemistry': 'Default',
            'Investigator Name': 'Knight',
            'Experiment Name': 'RKL_experiment'
        }

        self.assertDictEqual(dict(obs.Header), header)
        self.assertEqual(len(obs.samples), 1)

    def test_remove_options_for_iseq(self):
        sheet = MetagenomicSampleSheetv100()
        self.md_metag['Assay'] = 'Metagenomic'
        obs = sheet._add_metadata_to_sheet(self.md_metag, 'iSeq')

        settings = {
            'ReverseComplement': '0'
        }

        self.assertEqual(obs.Settings, settings)


class ValidateSampleSheetTests(BaseTests):
    def assertStdOutEqual(self, expected):
        # testing stdout: https://stackoverflow.com/a/12683001
        observed = sys.stdout.getvalue().strip()
        self.assertEqual(observed, expected)

    def test_validate_and_scrub_sample_sheet(self):
        sheet = MetagenomicSampleSheetv100(self.good_ss)
        # no errors
        self.assertTrue(sheet.validate_and_scrub_sample_sheet())

    def test_quiet_validate_and_scrub_sample_sheet(self):
        sheet = MetagenomicSampleSheetv100(self.good_ss)
        msgs = sheet.quiet_validate_and_scrub_sample_sheet()
        # no errors
        self.assertStdOutEqual('')
        self.assertEqual(msgs, [])

    def test_validate_and_scrub_sample_sheet_no_sample_project(self):
        sheet = MetagenomicSampleSheetv100(self.no_project_ss)
        self.assertFalse(sheet.validate_and_scrub_sample_sheet())

        self.assertStdOutEqual('ErrorMessage: The Sample_Project column in the'
                               ' Data section is missing')

    def test_quiet_validate_and_scrub_sample_sheet_no_sample_project(self):
        sheet = MetagenomicSampleSheetv100(self.no_project_ss)
        msgs = sheet.quiet_validate_and_scrub_sample_sheet()

        self.assertStdOutEqual('')
        self.assertEqual(msgs, [ErrorMessage('The Sample_Project column in '
                                             'the Data section is missing')])

    def test_validate_and_scrub_sample_sheet_missing_bioinformatics(self):
        sheet = MetagenomicSampleSheetv100(self.good_ss)
        sheet.Bioinformatics = None
        self.assertFalse(sheet.validate_and_scrub_sample_sheet())

        self.assertStdOutEqual('ErrorMessage: The Bioinformatics section '
                               'cannot be empty')

    def test_quiet_validate_scrub_sample_sheet_missing_bioinformatics(self):
        sheet = MetagenomicSampleSheetv100(self.good_ss)
        sheet.Bioinformatics = None
        msgs = sheet.quiet_validate_and_scrub_sample_sheet()

        self.assertStdOutEqual('')
        self.assertEqual(msgs, [ErrorMessage('The Bioinformatics section '
                                             'cannot be empty')])

    def test_validate_and_scrub_sample_sheet_missing_contact(self):
        sheet = MetagenomicSampleSheetv100(self.good_ss)
        sheet.Contact = None
        self.assertFalse(sheet.validate_and_scrub_sample_sheet())

        self.assertStdOutEqual('ErrorMessage: The Contact section '
                               'cannot be empty')

    def test_validate_and_scrub_sample_sheet_scrubbed_names(self):
        sheet = MetagenomicSampleSheetv100(self.scrubbable_ss)

        message = ('WarningMessage: '
                   'The following sample names were scrubbed for bcl2fastq '
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

        self.assertTrue(sheet.validate_and_scrub_sample_sheet())
        self.assertStdOutEqual(message)

    def test_quiet_validate_and_scrub_sample_sheet_scrubbed_names(self):
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
        message = WarningMessage(message)

        sheet = MetagenomicSampleSheetv100(self.scrubbable_ss)
        msgs = sheet.quiet_validate_and_scrub_sample_sheet()
        self.assertStdOutEqual('')
        self.assertEqual(msgs, [message])

    def test_validate_and_scrub_sample_sheet_scrubbed_project_names(self):
        sheet = MetagenomicSampleSheetv100(self.good_ss)

        remapper = {
            'NYU_BMS_Melanoma_13059': "NYU's Tisch Art Microbiome 13059",
            'Feist_11661': "The x.x microbiome project 1337"
        }

        for sample in sheet:
            sample['Sample_Project'] = remapper.get(sample['Sample_Project'],
                                                    sample['Sample_Project'])

        sheet.Contact.Sample_Project.replace(remapper, inplace=True)
        sheet.Bioinformatics.Sample_Project.replace(remapper, inplace=True)

        sheet.validate_and_scrub_sample_sheet()

        message = (
            'WarningMessage: The following project names were scrubbed for '
            'bcl2fastq compatibility. If the same invalid characters are also '
            'found in the Bioinformatics and Contacts sections those will be '
            'automatically scrubbed too:\n'
            "NYU's Tisch Art Microbiome 13059, The x.x microbiome project 1337"
        )
        self.assertStdOutEqual(message)

        scrubbed = {
            'NYU_s_Tisch_Art_Microbiome_13059',
            'The_x_x_microbiome_project_1337',
            'Gerwick_6123'
        }

        for sample in sheet:
            self.assertTrue(sample['Sample_Project'] in scrubbed,
                            sample['Sample_Project'])

        for project in sheet.Bioinformatics.Sample_Project:
            self.assertTrue(project in scrubbed)

        for project in sheet.Contact.Sample_Project:
            self.assertTrue(project in scrubbed)

    def test_validate_and_scrub_sample_sheet_bad_project_names(self):
        sheet = MetagenomicSampleSheetv100(self.bad_project_name_ss)

        message = ('ErrorMessage: The following project names in the '
                   'Sample_Project column are missing a Qiita study '
                   'identifier: Feist, Gerwick')

        self.assertFalse(sheet.validate_and_scrub_sample_sheet())
        self.assertStdOutEqual(message)

    def test_validate_and_scrub_sample_sheet_project_missing_lane(self):
        sheet = MetagenomicSampleSheetv100(self.good_ss)

        # set the lane value as empty for one of the two projects
        for sample in sheet.samples:
            if sample.Sample_Project == 'Feist_11661':
                sample.Lane = ' '

        self.assertFalse(sheet.validate_and_scrub_sample_sheet())
        message = ('ErrorMessage: The following projects are missing a Lane '
                   'value: Feist_11661')
        self.assertStdOutEqual(message)

    def test_sample_sheet_to_dataframe(self):
        ss = MetagenomicSampleSheetv100(self.ss)
        obs = sample_sheet_to_dataframe(ss)

        columns = ['lane', 'sample_name', 'sample_plate', 'well_id_384',
                   'i7_index_id', 'index', 'i5_index_id', 'index2',
                   'sample_project', 'well_description',
                   'library_construction_protocol',
                   'experiment_design_description']
        index = ['sample_1', 'sample_2', 'sample_1', 'sample_2', 'sample_31',
                 'sample_32', 'sample_34', 'sample_44']

        exp = pd.DataFrame(index=index, data=DF_DATA, columns=columns)
        exp.index.name = 'sample_id'
        pd.testing.assert_frame_equal(obs, exp)


class ProfileTests(BaseTests):
    def test_profile(self):
        sheet = AbsQuantSampleSheetv10()

        # confirm that AbsQuantSampleSheetv10() contains the right values
        # for sheet-type and sheet-version, not the default values inherited
        # from its parent.
        self.assertEqual(sheet._HEADER['SheetType'], 'abs_quant_metag')
        self.assertEqual(sheet._HEADER['SheetVersion'], '10')
        self.assertIn('mass_syndna_input_ng', sheet.data_columns)


class DemuxReplicatesTests(BaseTests):
    def setUp(self):
        self.data_dir = join('metapool', 'tests', 'data')
        self.sheet_w_replicates_path = join(self.data_dir,
                                            'good_sheet_w_replicates.csv')

        # bad_sheet_w_replicates.csv contains two projects, one of which
        # doesn't contain replicates. By convention, all projects in the sheet
        # must either contain replicates or not contain replicates.
        self.bad_sht_w_replicates_path = join(self.data_dir,
                                              'bad_sheet_w_replicates.csv')

        self.sheet_wo_replicates_path = join(self.data_dir,
                                             'sheet_wo_replicates.csv')

        self.legacy_sheet_path = join(self.data_dir, 'good-sample-sheet.csv')

        self.replicate_output_paths = [join(self.data_dir,
                                            'replicate_output1.csv'),
                                       join(self.data_dir,
                                            'replicate_output2.csv'),
                                       join(self.data_dir,
                                            'replicate_output3.csv')]

    def test_sheet_needs_demuxing(self):
        # confirm legacy sample-sheets w/out contains_replicates column will
        # return False, instead of raising an Error. For processing purposes,
        # it's only critical to know whether the sheet needs demuxing or not.
        sheet = MetagenomicSampleSheetv90(self.legacy_sheet_path)
        self.assertFalse(sheet_needs_demuxing(sheet))

        # confirm bad sample-sheet raises a ValueError for containing projects
        # that contain replicates and those that don't.
        with self.assertRaisesRegex(ValueError, "all projects in "
                                                "Bioinformatics section must "
                                                "either contain replicates or "
                                                "not."):
            sheet = MetagenomicSampleSheetv100(self.bad_sht_w_replicates_path)
            sheet_needs_demuxing(sheet)

        # test a valid sample-sheet with replicates.
        sheet = MetagenomicSampleSheetv100(self.sheet_w_replicates_path)
        self.assertTrue(sheet_needs_demuxing(sheet))

        # test a valid sample-sheet w/out replicates.
        sheet = MetagenomicSampleSheetv100(self.sheet_wo_replicates_path)
        self.assertFalse(sheet_needs_demuxing(sheet))

    def test_demux_sample_sheet(self):
        # we don't want to demux legacy sample-sheets. sheet_needs_demuxing()
        # should be used to determine if demux_sample_sheet() should be
        # called.
        with self.assertRaisesRegex(ValueError, "sample-sheet does not contain"
                                                " replicates"):
            sheet = MetagenomicSampleSheetv90(self.legacy_sheet_path)
            demux_sample_sheet(sheet)

        # by convention, all replication is done at the plate level, and all
        # projects in a sample-sheet will either contain replicates, or all of
        # them will not. Hence, a sample-sheet with both True and False in
        # the contains_replicates column in the [Bioinformatics] section should
        # raise an error.
        with self.assertRaisesRegex(ValueError, "all projects in Bioinfor"
                                                "matics section must either "
                                                "contain replicates or not."):
            sheet = MetagenomicSampleSheetv100(self.bad_sht_w_replicates_path)
            demux_sample_sheet(sheet)

        # as mentioned above, sheet_needs_demuxing() should be used to
        # determine if demux_sample_sheet() should be called. If a sample-
        # sheet is passed to demux_sample_sheet() and all projects are False,
        # an Error should be raised to alert the user of an unexpected
        # condition, rather than silently allow as a degenerative case.
        with self.assertRaisesRegex(ValueError, "all projects in Bioinfor"
                                                "matics section do not contain"
                                                " replicates"):
            sheet = MetagenomicSampleSheetv100(self.sheet_wo_replicates_path)
            demux_sample_sheet(sheet)

        # this test will need to compare the four completed sample-sheets
        # made using self.sheet_w_replicates_path against an expected result.

        # test sample-sheet w/both projects w/replicates and not.
        sheet = MetagenomicSampleSheetv100(self.sheet_w_replicates_path)
        results = demux_sample_sheet(sheet)

        # assert that the proper number of KLSampleSheets were returned.
        self.assertEqual(len(results), len(self.replicate_output_paths))

        # assert that each sample-sheet appears in the correct order and
        # matches known results.
        for replicate_output_path in self.replicate_output_paths:
            exp = MetagenomicSampleSheetv100(replicate_output_path)
            obs = results.pop(0)
            self.assertEqual(obs.Header, exp.Header)
            self.assertEqual(obs.Reads, exp.Reads)
            self.assertEqual(obs.Settings, exp.Settings)
            self.assertTrue(obs.Bioinformatics.equals(exp.Bioinformatics))
            self.assertTrue(obs.Contact.equals(exp.Contact))
            for o_sample, e_sample in zip(obs.samples, exp.samples):
                self.assertEqual(o_sample, e_sample)


DF_DATA = [
    ['1', 'sample.1', 'FooBar_666_p1', 'A1', 'iTru7_107_07', 'CCGACTAT',
     'iTru5_01_A', 'ACCGACAA', 'Baz_12345', 'importantsample1',
     'Knight Lab Kapa HP', 'Eqiiperiment'],
    ['1', 'sample.2', 'FooBar_666_p1', 'A2', 'iTru7_107_08', 'CCGACTAT',
     'iTru5_01_A', 'CTTCGCAA', 'Baz_12345', 'importantsample2',
     'Knight Lab Kapa HP', 'Eqiiperiment'],
    ['3', 'sample.1', 'FooBar_666_p1', 'A3', 'iTru7_107_09', 'GCCTTGTT',
     'iTru5_01_A', 'AACACCAC', 'Baz_12345', 'importantsample1',
     'Knight Lab Kapa HP', 'Eqiiperiment'],
    ['3', 'sample.2', 'FooBar_666_p1', 'A4', 'iTru7_107_10', 'AACTTGCC',
     'iTru5_01_A', 'CGTATCTC', 'Baz_12345', 'importantsample2',
     'Knight Lab Kapa HP', 'Eqiiperiment'],
    ['3', 'sample.31', 'FooBar_666_p1', 'A5', 'iTru7_107_11', 'CAATGTGG',
     'iTru5_01_A', 'GGTACGAA', 'FooBar_666', 'importantsample31',
     'Knight Lab Kapa HP', 'SomethingWitty'],
    ['3', 'sample.32', 'FooBar_666_p1', 'B6', 'iTru7_107_12', 'AAGGCTGA',
     'iTru5_01_A', 'CGATCGAT', 'FooBar_666', 'importantsample32',
     'Knight Lab Kapa HP', 'SomethingWitty'],
    ['3', 'sample.34', 'FooBar_666_p1', 'B8', 'iTru7_107_13', 'TTACCGAG',
     'iTru5_01_A', 'AAGACACC', 'FooBar_666', 'importantsample34',
     'Knight Lab Kapa HP', 'SomethingWitty'],
    ['3', 'sample.44', 'Baz_12345_p3', 'B99', 'iTru7_107_14', 'GTCCTAAG',
     'iTru5_01_A', 'CATCTGCT', 'Baz_12345', 'importantsample44',
     'Knight Lab Kapa HP', 'Eqiiperiment']]


if __name__ == '__main__':
    assert not hasattr(sys.stdout, "getvalue")
    unittest.main(module=__name__, buffer=True, exit=False)
