import sys
import unittest
import tempfile
from datetime import datetime
from os.path import join, dirname
from io import StringIO
from contextlib import redirect_stdout

import pandas as pd
from pandas.testing import assert_frame_equal
import sample_sheet

from metapool.mp_strings import (
    QIITA_ID_KEY, PROJECT_SHORT_NAME_KEY, PROJECT_FULL_NAME_KEY,
    CONTAINS_REPLICATES_KEY, SAMPLES_DETAILS_KEY, SAMPLE_PROJECT_KEY,
    ORIG_NAME_KEY, SAMPLE_NAME_KEY, SAMPLE_TYPE_KEY, PRIMARY_STUDY_KEY,
    SECONDARY_STUDIES_KEY, DESTINATION_WELL_384_KEY, SS_SAMPLE_ID_KEY)
from metapool.metapool import TUBECODE_KEY
from metapool.sample_sheet import (KLSampleSheet, AmpliconSampleSheet,
                                   MetagenomicSampleSheetv102,
                                   MetagenomicSampleSheetv101,
                                   MetagenomicSampleSheetv100,
                                   MetagenomicSampleSheetv90,
                                   MetatranscriptomicSampleSheetv0,
                                   MetatranscriptomicSampleSheetv10,
                                   AbsQuantSampleSheetv10,
                                   AbsQuantSampleSheetv11,
                                   TellseqMetagSampleSheetv10,
                                   TellseqAbsquantMetagSampleSheetv10,
                                   PacBioMetagSampleSheetv10,
                                   PacBioAbsquantSampleSheetv10,
                                   PacBioMetagSampleSheetv11,
                                   PacBioAbsquantSampleSheetv11,
                                   PacBioAbsquantSampleSheetv12,
                                   sample_sheet_to_dataframe,
                                   make_sample_sheet, load_sample_sheet,
                                   demux_sample_sheet, sheet_needs_demuxing,
                                   make_sections_dict,
                                   _ASSAY_KEY, _SHEET_VERSION_KEY,
                                   _SHEET_TYPE_KEY, _BIOINFORMATICS_KEY,
                                   _CONTACT_KEY, _SAMPLE_CONTEXT_KEY)
from metapool.plate import WarningMessage
from metapool.metapool import generate_override_cycles_value


# Default KLSampleSheet objects don't have a `contains_replicates`
# key in the Bioinformatics section, but metag v100 and later sheets do, so
# sometimes we need to expand the base dummy info
def _add_contains_replicates(source_bfx_list):
    out_bfx_list = []
    for x in source_bfx_list:
        curr_dict = x.copy()
        curr_dict['contains_replicates'] = False
        out_bfx_list.append(curr_dict)
    return out_bfx_list


# The classes below share the same filepaths, so we use this dummy class
class BaseTests(unittest.TestCase):
    def setUp(self):
        data_dir = join(dirname(__file__), 'data')
        self.data_dir = data_dir

        self.alt_ss = join(data_dir,
                           'good-sample-sheet-with-alt-col-names.csv')

        self.good_ss = join(data_dir, 'good-sample-sheet.csv')
        self.good_metag_ss_w_context = \
            join(data_dir, "good-sample-sheet_w_sample_context.csv")
        self.with_comments = join(data_dir, 'good-sample-sheet-but-'
                                            'with-comments.csv')

        self.good_w_bools = join(data_dir, 'good-sheet-w-odd-bools.csv')

        fp = 'good-sample-sheet-with-comments-and-new-lines.csv'
        self.with_comments_and_new_lines = join(data_dir, fp)

        self.with_new_lines = join(data_dir, 'good-sample-sheet-with-'
                                             'new-lines.csv')

        self.no_project_ss = join(data_dir, 'no-project-name-sample-sheet.csv')

        self.scrubbable_ss = join(data_dir, 'scrubbable-sample-sheet.csv')

        self.bad_project_name_ss = join(data_dir,
                                        'bad-project-name-sample-sheet.csv')

        self.good_run_info = join(data_dir, "runinfo_files/RunInfo1.xml")

        bfx = [
            {
                'Sample_Project': 'Koening_ITS_101',
                'QiitaID': '101',
                'BarcodesAreRC': False,
                'ForwardAdapter': 'GATACA',
                'ReverseAdapter': 'CATCAT',
                'HumanFiltering': False,
                'library_construction_protocol': 'Knight Lab Kapa HP',
                'experiment_design_description': 'Eqiiperiment'
            },
            {
                'Sample_Project': 'Yanomani_2008_10052',
                'QiitaID': '10052',
                'BarcodesAreRC': False,
                'ForwardAdapter': 'GATACA',
                'ReverseAdapter': 'CATCAT',
                'HumanFiltering': False,
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
            'Investigator Name': 'a PI',
            'Experiment Name': 'an experiment name',
            'Bioinformatics': bfx,
            'Contact': contact,
            'Assay': 'TruSeq HT',
            'SheetType': 'dummy_amp',
            'SheetVersion': '0'
        }

        self.md_metag = {
            'Bioinformatics': bfx,
            'Contact': contact,
            'Assay': 'Metagenomic',
            'SheetType': 'standard_metag',
            'SheetVersion': '100'
        }

    def _help_test_csv_files_exact_text_match(self, file_1, file_2):
        with open(file_1, 'r', encoding='utf-8') as f1, \
                open(file_2, 'r', encoding='utf-8') as f2:
            text1 = f1.read()
            text2 = f2.read()
        self.assertMultiLineEqual(text1, text2)


class KLSampleSheetTests(BaseTests):
    def test_instantiation(self):
        # base class can no longer be instantiated
        with self.assertRaises(TypeError, msg="TypeError: only children of "
                                              "'KLSampleSheet' may be insta"
                                              "ntiated"):
            KLSampleSheet()

        # child class should instantiate successfully.
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)
        self.assertIsNotNone(sheet)

    def test__get_section_exists(self):
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)
        section = sheet._get_section(_BIOINFORMATICS_KEY)
        self.assertEqual(pd.DataFrame, type(section))
        self.assertEqual((3, 8), section.shape)

    def test__get_section_missing(self):
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)
        delattr(sheet, _BIOINFORMATICS_KEY)
        section = sheet._get_section(_BIOINFORMATICS_KEY)
        self.assertIsNone(section)

    def test__get_section_missing_err(self):
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)
        delattr(sheet, _BIOINFORMATICS_KEY)
        err_msg = "Section 'Bioinformatics' does not exist in sample sheet."
        with self.assertRaisesRegex(ValueError, err_msg):
            _ = sheet._get_section(_BIOINFORMATICS_KEY, err_if_missing=True)

    def test__get_data_section_to_df(self):
        # use good_w_bools sheet because its data section is mercifully short
        exp = pd.DataFrame([
            {'Lane': '1', 'Sample_ID': 'CDPH-SAL_Salmonella_Typhi_MDL-143',
             'Sample_Name': 'CDPH-SAL_Salmonella_Typhi_MDL-143',
             'Sample_Plate': 'ProjectF_11661_P40', 'well_id_384': 'A1',
             'I7_Index_ID': 'iTru7_107_07', 'index': 'CCGACTAT',
             'I5_Index_ID': 'iTru5_01_A', 'index2': 'ACCGACAA',
             'Sample_Project': 'ProjectF_11661',
             'Well_description': 'Desc_for_CDPH-SAL_Salmonella Typhi_MDL-143'},
            {'Lane': '1', 'Sample_ID': '3A', 'Sample_Name': '3A',
             'Sample_Plate': 'ProjectG_tubes', 'well_id_384': 'I23',
             'I7_Index_ID': 'iTru7_201_03', 'index': 'GATAGGCT',
             'I5_Index_ID': 'iTru5_09_H', 'index2': 'AGAAGGAC',
             'Sample_Project': 'ProjectG_6123',
             'Well_description': 'Desc_for_3A'},
            {'Lane': '1', 'Sample_ID': 'LP127890A01',
             'Sample_Name': 'LP127890A01',
             'Sample_Plate': 'ProjectN_13059_P1', 'well_id_384': 'I3',
             'I7_Index_ID': 'iTru7_108_09', 'index': 'TCTCTAGG',
             'I5_Index_ID': 'iTru5_01_B', 'index2': 'AGTGGCAA',
             'Sample_Project': 'ProjectN_13059',
             'Well_description': 'Desc_for_LP127890A01'}])
        sheet = load_sample_sheet(self.good_w_bools, defer_validate=False)
        obs = sheet._get_data_section_to_df()
        assert_frame_equal(exp, obs)

    def test_sample_sheet_roundtripping(self):
        # testing with various good sheets we have access to
        sheets = [self.good_ss,
                  self.scrubbable_ss,
                  self.with_comments, self.with_comments_and_new_lines,
                  self.with_new_lines]

        for filename in sheets:
            sheet = load_sample_sheet(filename, defer_validate=False)

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
                    self.assertEqual(observed.split(),
                                     f.read().split(),
                                     f'Problem found with {filename}')

    def test_write_w_lane(self):
        test_fp = self.good_metag_ss_w_context.replace(
            ".csv", "_lane_overwritten.csv")
        sheet = MetagenomicSampleSheetv101(test_fp, defer_validate=False)

        with tempfile.NamedTemporaryFile('w+') as tmp:
            sheet.write(tmp, lane=3)
            tmp.seek(0)
            observed = tmp.read()

        with open(test_fp) as f:
            expected = f.read()

            # Normalize line endings to '\n' for both files
        expected1 = expected.replace('\r\n', '\n').replace('\r', '\n')
        observed1 = observed.replace('\r\n', '\n').replace('\r', '\n')
        self.assertEqual(expected1, observed1)

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

            sheet = MetagenomicSampleSheetv100(tmp.name, defer_validate=True)

            self.assertEqual(sheet.samples, [])
            self.assertEqual(sheet.Settings, {})
            self.assertEqual(sheet.Header, {})
            self.assertEqual(sheet.Reads, [])
            self.assertIsNone(sheet.Bioinformatics)
            self.assertIsNone(sheet.Contact)

    def test_parse(self):
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)

        exp = {
            'IEMFileVersion': '4',
            'SheetType': 'standard_metag',
            'SheetVersion': '90',
            'Investigator Name': 'Knight',
            'Experiment Name': 'RKL0042',
            'Date': '2/26/20',
            'Workflow': 'GenerateFASTQ',
            'Application': 'FASTQ Only',
            'Assay': 'Metagenomic',
            'Description': '',
            'Chemistry': 'Default'}

        self.assertEqual(sheet.Header, exp)
        self.assertEqual(sheet.Reads, [150, 150])

        exp = {'ReverseComplement': '0', 'MaskShortReads': '1',
               'OverrideCycles': 'Y151;I8N2;I8N2;Y151'}

        self.assertEqual(sheet.Settings, exp)

        # a selection of samples from good-sample-sheet.csv.
        # one sample from each project.

        exp = [
            {
                "Lane": "1",
                "Sample_ID": "3A",
                "Sample_Name": "3A",
                "Sample_Plate": "ProjectG_tubes",
                "Sample_Well": "I23",
                "I7_Index_ID": "iTru7_201_03",
                "index": "GATAGGCT",
                "I5_Index_ID": "iTru5_09_H",
                "index2": "AGAAGGAC",
                "Sample_Project": "ProjectG_6123",
                "Well_description": "Desc_for_3A"
            }, {
                "Lane": "1",
                "Sample_ID": "JM-MEC__Staphylococcus_aureusstrain_BERTI-"
                "B0387",
                "Sample_Name": "JM-MEC__Staphylococcus_aureusstrain_BERTI-"
                "B0387",
                "Sample_Plate": "ProjectF_11661_P43",
                "Sample_Well": "B10",
                "I7_Index_ID": "iTru7_102_03",
                "index": "CGTTGCAA",
                "I5_Index_ID": "iTru5_121_C",
                "index2": "CTATGCCT",
                "Sample_Project": "ProjectF_11661",
                "Well_description": "Desc_for_JM-MEC__Staphylococcus "
                "aureusstrain BERTI-B0387"
            }, {
                "Lane": "1",
                "Sample_ID": "EP981129A02",
                "Sample_Name": "EP981129A02",
                "Sample_Plate": "ProjectN_13059_P4",
                "Sample_Well": "H4",
                "I7_Index_ID": "iTru7_115_08",
                "index": "TGGTACAG",
                "I5_Index_ID": "iTru5_124_A",
                "index2": "GTACGATC",
                "Sample_Project": "ProjectN_13059",
                "Well_description": "Desc_for_EP981129A02"
            }
        ]

        obs = {}

        for sample in sheet.samples:
            # sample-names can be duplicated across projects, hence the need
            # to store each sample in a truly unique value.
            key = f"{sample.Sample_Name}_{sample.Sample_Project}"
            obs[key] = sample.to_json()

        for sample in exp:
            key = f'{sample["Sample_Name"]}_{sample["Sample_Project"]}'
            self.assertDictEqual(obs[key], sample)

        # check for Bioinformatics
        exp = pd.DataFrame(
            columns=['Sample_Project', 'QiitaID', 'BarcodesAreRC',
                     'ForwardAdapter', 'ReverseAdapter', 'HumanFiltering',
                     'library_construction_protocol',
                     'experiment_design_description'],
            data=[
                ['ProjectN_13059', '13059', False, 'AACC', 'GGTT',
                 False, 'Knight Lab Kapa HP', 'Eqiiperiment'],
                ['ProjectF_11661', '11661', False, 'AACC', 'GGTT', False,
                 'Knight Lab Kapa HP', 'Eqiiperiment'],
                ['ProjectG_6123', '6123', False, 'AACC', 'GGTT', False,
                 'Knight Lab Kapa HP', 'Eqiiperiment']
            ]
        )

        pd.testing.assert_frame_equal(sheet.Bioinformatics, exp)

        # check for Contact
        exp = pd.DataFrame(
            columns=['Email', 'Sample_Project'],
            data=[
                ['test@lol.com', 'ProjectF_11661']
            ]
        )

        pd.testing.assert_frame_equal(sheet.Contact, exp)

    def test_parse_with_comments(self):
        # the two sample sheets are identical except for the comments
        exp = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)
        with self.assertWarnsRegex(UserWarning, 'Comments at the beginning '):
            obs = MetagenomicSampleSheetv90(
                self.with_comments, defer_validate=False)

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
                ['Koening_ITS_101', '101', False, 'GATACA', 'CATCAT',
                 False, 'Knight Lab Kapa HP', 'Eqiiperiment'],
                ['Yanomani_2008_10052', '10052', False, 'GATACA', 'CATCAT',
                 False, 'Knight Lab Kapa HP', 'Eqiiperiment'],
                ['paco_Koening_ITS_101', '101', False, 'GATACA', 'CATCAT',
                 False, 'Knight Lab Kapa HP', 'Eqiiperiment'],
                ['paco_Yanomani_2008_10052', '10052', False, 'GATACA',
                 'CATCAT', False, 'Knight Lab Kapa HP', 'Eqiiperiment']
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
        hugo.Reads = [151, 151]
        hugo.add_sample(sample_sheet.Sample({
            'Sample_ID': 'a',
            'index': 'GATACA',
            'index2': 'GCCGCC',
            'Sample_Name': 'a.sample'
        }))

        err = "The Settings section is different for sample sheet 1"
        with self.assertRaisesRegex(ValueError, err):
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
        obs = sheet._validate_metadata_dict(self.md_ampl)
        self.assertEqual(obs, [])

    def test_more_attributes(self):
        sheet = AmpliconSampleSheet()
        self.md_ampl['Ride'] = 'the lightning'

        obs = sheet._validate_metadata_dict(self.md_ampl)
        exp = ['ErrorMessage: These metadata keys are not supported: Ride']
        msg_strs = [str(msg) for msg in obs]
        self.assertEqual(msg_strs, exp)

    def test_validate_missing_assay(self):
        sheet = AmpliconSampleSheet()
        self.md_ampl['Assay'] = 'NewAssayType'

        obs = sheet._validate_metadata_dict(self.md_ampl)
        exp = ["ErrorMessage: 'NewAssayType' is not a supported Assay"]
        msg_strs = [str(msg) for msg in obs]
        self.assertEqual(msg_strs, exp)

    def test_validate_missing_bioinformatics_data(self):
        sheet = AmpliconSampleSheet()
        del self.md_ampl['Bioinformatics']

        obs = sheet._validate_metadata_dict(self.md_ampl)
        exp = ['ErrorMessage: Bioinformatics is a required attribute']
        msg_strs = [str(msg) for msg in obs]
        self.assertEqual(msg_strs, exp)

    def test_validate_missing_column_in_bioinformatics(self):
        sheet = AmpliconSampleSheet()
        del self.md_ampl['Bioinformatics'][0]['Sample_Project']
        exp = ['ErrorMessage: In the Bioinformatics section Project #1 does '
               'not have exactly these keys BarcodesAreRC, '
               'ForwardAdapter, HumanFiltering, QiitaID, '
               'ReverseAdapter, Sample_Project, '
               'experiment_design_description, '
               'library_construction_protocol']
        obs = sheet._validate_metadata_dict(self.md_ampl)
        self.assertEqual(str(obs[0]), exp[0])

    def test_set_override_cycles(self):
        sheet = load_sample_sheet(self.good_ss, defer_validate=False)

        # assert that the original value of the sheet is as expected.
        self.assertEqual("Y151;I8N2;I8N2;Y151",
                         sheet.Settings['OverrideCycles'])

        # generate a known value that is different from above using a known
        # sample-sheet. Assume that adapters are of length 8.
        new_value = generate_override_cycles_value(self.good_run_info, 8)

        # assert that the new value is as expected.
        self.assertEqual("Y151;I8N4;Y151", new_value)

        # use set_override_cycles() to change the value and assert that it
        # is now different.
        sheet.set_override_cycles(new_value)
        self.assertEqual("Y151;I8N4;Y151", sheet.Settings['OverrideCycles'])

    def test_sample_is_a_blank_wo_context(self):
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)
        # NB: the sample names and sample ids in the test spreadsheet are the
        # same.  FWIW, the intention is to use the sample names here.
        self.assertFalse(sheet.sample_is_a_blank(
            'CDPH-SAL_Salmonella_Typhi_MDL-143'))
        self.assertTrue(sheet.sample_is_a_blank('BLANK_40_12G'))

        # sending in a sample name that doesn't exist in the Data section
        # raises an error
        with self.assertRaisesRegex(
                ValueError,
                "Sample 'blank_40_12g' not found in the Data section"):
            sheet.sample_is_a_blank('blank_40_12g')

    def test_sample_is_a_blank_w_context(self):
        sheet = MetagenomicSampleSheetv101(
            self.good_metag_ss_w_context, defer_validate=False)
        # NB: the sample names and sample ids in the test spreadsheet are the
        # same.  FWIW, the intention is to use the sample names here.
        self.assertFalse(sheet.sample_is_a_blank(
            'CDPH-SAL_Salmonella_Typhi_MDL-143'))
        self.assertTrue(sheet.sample_is_a_blank('BLANK_40_12G'))

        # sending in a sample name that doesn't exist in the Data section
        # raises an error
        with self.assertRaisesRegex(
                ValueError,
                "Sample 'blank_40_12g' not found in the Data section"):
            sheet.sample_is_a_blank('blank_40_12g')

    def test_get_controls_details_w_context(self):
        exp_blank_names = [
            'BLANK_40_12G', 'BLANK_40_12H', 'BLANK_41_12G', 'BLANK_41_12H',
            'BLANK_42_12G', 'BLANK_42_12H', 'BLANK_43_12G', 'BLANK_43_12H',
            'BLANK1_1A', 'BLANK1_1B', 'BLANK1_1C', 'BLANK1_1D', 'BLANK1_1E',
            'BLANK1_1F', 'BLANK1_1G', 'BLANK1_1H', 'BLANK2_2A', 'BLANK2_2B',
            'BLANK2_2C', 'BLANK2_2D', 'BLANK2_2E', 'BLANK2_2F', 'BLANK2_2G',
            'BLANK2_2H', 'BLANK3_3A', 'BLANK3_3B', 'BLANK3_3C', 'BLANK3_3D',
            'BLANK3_3E', 'BLANK3_3F', 'BLANK3_3G', 'BLANK3_3H', 'BLANK4_4A',
            'BLANK4_4B', 'BLANK4_4C', 'BLANK4_4D', 'BLANK4_4E', 'BLANK4_4F',
            'BLANK4_4G', 'BLANK4_4H']
        sheet = MetagenomicSampleSheetv101(
            self.good_metag_ss_w_context, defer_validate=False)
        obs_details = sheet.get_controls_details()
        self.assertEqual(len(exp_blank_names), len(obs_details))
        for curr_blank_name in exp_blank_names:
            self.assertTrue(curr_blank_name in obs_details)
            curr_details = obs_details[curr_blank_name]
            self.assertEqual(len(curr_details), 4)
            self.assertEqual(curr_details[SAMPLE_NAME_KEY],
                             curr_blank_name)
            self.assertEqual(curr_details[SAMPLE_TYPE_KEY], 'control blank')

            if curr_blank_name.startswith("BLANK_4"):
                expected_qiita_id = "11661"
            else:
                expected_qiita_id = "13059"
            self.assertEqual(curr_details[PRIMARY_STUDY_KEY],
                             expected_qiita_id)

            expected_secondary_studies = []
            if curr_blank_name in ["BLANK_41_12G", "BLANK_41_12H"]:
                expected_secondary_studies = ["10317", "11223"]
            self.assertEqual(curr_details[SECONDARY_STUDIES_KEY],
                             expected_secondary_studies)
        # next expected blank name

    def test_get_controls_details_wo_context(self):
        exp_blank_names = [
            'BLANK_40_12G', 'BLANK_40_12H', 'BLANK_41_12G', 'BLANK_41_12H',
            'BLANK_42_12G', 'BLANK_42_12H', 'BLANK_43_12G', 'BLANK_43_12H',
            'BLANK1_1A', 'BLANK1_1B', 'BLANK1_1C', 'BLANK1_1D', 'BLANK1_1E',
            'BLANK1_1F', 'BLANK1_1G', 'BLANK1_1H', 'BLANK2_2A', 'BLANK2_2B',
            'BLANK2_2C', 'BLANK2_2D', 'BLANK2_2E', 'BLANK2_2F', 'BLANK2_2G',
            'BLANK2_2H', 'BLANK3_3A', 'BLANK3_3B', 'BLANK3_3C', 'BLANK3_3D',
            'BLANK3_3E', 'BLANK3_3F', 'BLANK3_3G', 'BLANK3_3H', 'BLANK4_4A',
            'BLANK4_4B', 'BLANK4_4C', 'BLANK4_4D', 'BLANK4_4E', 'BLANK4_4F',
            'BLANK4_4G', 'BLANK4_4H']
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)
        obs_details = sheet.get_controls_details()
        self.assertEqual(len(exp_blank_names), len(obs_details))
        for curr_blank_name in exp_blank_names:
            self.assertTrue(curr_blank_name in obs_details)
            curr_details = obs_details[curr_blank_name]
            self.assertEqual(len(curr_details), 4)
            self.assertEqual(curr_details[SAMPLE_NAME_KEY],
                             curr_blank_name)
            self.assertEqual(curr_details[SAMPLE_TYPE_KEY], 'control blank')

            if curr_blank_name == "BLANK_41_12G":
                expected_qiita_id = "6123"
            elif curr_blank_name.startswith("BLANK_4"):
                expected_qiita_id = "11661"
            else:
                expected_qiita_id = "13059"
            self.assertEqual(curr_details[PRIMARY_STUDY_KEY],
                             expected_qiita_id)
            self.assertEqual(curr_details[SECONDARY_STUDIES_KEY], [])
        # next expected blank name

    def test_get_denormalized_controls_list(self):
        exp_details_list = [
            ('BLANK_40_12G', '11661', 'ProjectF_11661'),
            ('BLANK_40_12H', '11661', 'ProjectF_11661'),
            ('BLANK_41_12G', '11661', 'ProjectF_11661'),
            ('BLANK_41_12G', '10317', 'TMI_10317'),
            ('BLANK_41_12G', '11223', 'Other_11223'),
            ('BLANK_41_12H', '11661', 'ProjectF_11661'),
            ('BLANK_41_12H', '10317', 'TMI_10317'),
            ('BLANK_41_12H', '11223', 'Other_11223'),
            ('BLANK_42_12G', '11661', 'ProjectF_11661'),
            ('BLANK_42_12H', '11661', 'ProjectF_11661'),
            ('BLANK_43_12G', '11661', 'ProjectF_11661'),
            ('BLANK_43_12H', '11661', 'ProjectF_11661')]
        nyu_blanks = [
            'BLANK1_1A', 'BLANK1_1B', 'BLANK1_1C', 'BLANK1_1D', 'BLANK1_1E',
            'BLANK1_1F', 'BLANK1_1G', 'BLANK1_1H', 'BLANK2_2A', 'BLANK2_2B',
            'BLANK2_2C', 'BLANK2_2D', 'BLANK2_2E', 'BLANK2_2F', 'BLANK2_2G',
            'BLANK2_2H', 'BLANK3_3A', 'BLANK3_3B', 'BLANK3_3C', 'BLANK3_3D',
            'BLANK3_3E', 'BLANK3_3F', 'BLANK3_3G', 'BLANK3_3H', 'BLANK4_4A',
            'BLANK4_4B', 'BLANK4_4C', 'BLANK4_4D', 'BLANK4_4E', 'BLANK4_4F',
            'BLANK4_4G', 'BLANK4_4H']
        for curr_blank_name in nyu_blanks:
            exp_details_list.append(
                (curr_blank_name, '13059', 'ProjectN_13059'))
        exp_details_list = sorted(exp_details_list, key=lambda k: (k[0], k[1]))

        sheet = MetagenomicSampleSheetv101(
            self.good_metag_ss_w_context, defer_validate=False)
        obs_details_list = sheet.get_denormalized_controls_list()
        self.assertEqual(len(exp_details_list), len(obs_details_list))
        for i in range(len(exp_details_list)):
            curr_exp_details = exp_details_list[i]
            curr_obs_details = obs_details_list[i]

            self.assertEqual(len(curr_obs_details), 4)
            self.assertEqual(
                curr_obs_details[SAMPLE_NAME_KEY], curr_exp_details[0])
            self.assertEqual(
                curr_obs_details[SAMPLE_TYPE_KEY], 'control blank')
            self.assertEqual(
                curr_obs_details[QIITA_ID_KEY], curr_exp_details[1])
            self.assertEqual(
                curr_obs_details[PROJECT_FULL_NAME_KEY], curr_exp_details[2])
        # next expected blank

    def test_get_projects_details(self):
        exp_details = {
            'ProjectN_13059': {
                QIITA_ID_KEY: '13059',
                PROJECT_SHORT_NAME_KEY: 'ProjectN',
                PROJECT_FULL_NAME_KEY: 'ProjectN_13059',
                CONTAINS_REPLICATES_KEY: False,
                SAMPLES_DETAILS_KEY: {
                    'LP127890A01': {
                        SAMPLE_NAME_KEY: 'LP127890A01',
                        SAMPLE_PROJECT_KEY: 'ProjectN_13059',
                        SS_SAMPLE_ID_KEY: 'LP127890A01'
                    }
                },
            },
            'ProjectF_11661': {
                QIITA_ID_KEY: '11661',
                PROJECT_SHORT_NAME_KEY: 'ProjectF',
                PROJECT_FULL_NAME_KEY: 'ProjectF_11661',
                CONTAINS_REPLICATES_KEY: False,
                SAMPLES_DETAILS_KEY: {
                    'CDPH-SAL_Salmonella_Typhi_MDL-143': {
                        SAMPLE_NAME_KEY: 'CDPH-SAL_Salmonella_Typhi_MDL-143',
                        SAMPLE_PROJECT_KEY: 'ProjectF_11661',
                        SS_SAMPLE_ID_KEY: 'CDPH-SAL_Salmonella_Typhi_MDL-143'
                    }
                },
            },
            'ProjectG_6123': {
                QIITA_ID_KEY: '6123',
                PROJECT_SHORT_NAME_KEY: 'ProjectG',
                PROJECT_FULL_NAME_KEY: 'ProjectG_6123',
                CONTAINS_REPLICATES_KEY: False,
                SAMPLES_DETAILS_KEY: {
                    '3A': {
                        SAMPLE_NAME_KEY: '3A',
                        SAMPLE_PROJECT_KEY: 'ProjectG_6123',
                        SS_SAMPLE_ID_KEY: '3A'
                    }
                }
            }
        }

        sheet = MetagenomicSampleSheetv100(
            self.good_w_bools, defer_validate=False)
        obs_details = sheet.get_projects_details()
        self.assertEqual(exp_details, obs_details)

    def test_get_projects_details_w_orig_name(self):
        good_replicates_ss_fp = join(
            self.data_dir, 'good_standard_metagv100_w_replicates.csv')
        sheet = MetagenomicSampleSheetv100(
            good_replicates_ss_fp, defer_validate=False)

        # not going to check the whole thing, just checking one sample to
        # make sure the original name is being added correctly
        obs_details = sheet.get_projects_details()
        self.assertTrue("ProjectF_11661" in obs_details)
        an_obs_samples = obs_details["ProjectF_11661"][SAMPLES_DETAILS_KEY]
        self.assertTrue("BLANK.43.12G.A1" in an_obs_samples)
        an_obs_sample = an_obs_samples["BLANK.43.12G.A1"]
        self.assertTrue(ORIG_NAME_KEY in an_obs_sample)
        self.assertEqual(an_obs_sample[ORIG_NAME_KEY], "BLANK.43.12G")


class SampleSheetWorkflow(BaseTests):
    def setUp(self):
        super().setUp()
        self.maxDiff = None

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
            columns=[SS_SAMPLE_ID_KEY,
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

    # Default KLSampleSheet objects don't have a `contains_replicates`
    # key in the Bioinformatics section, but metag v100 sheets do, and
    # several tests here use v100, so they need to have a `contains_replicates`
    # entry added to the base metadata in order to pass
    def _add_replicates_to_md_metag(self):
        input_dict = {}
        for key, value in self.md_metag.items():
            mod_value = value
            if key == 'Bioinformatics':
                mod_value = _add_contains_replicates(value)
            input_dict[key] = mod_value
        return input_dict

    def test_validate_sample_sheet_metadata_empty(self):
        sheet = AmpliconSampleSheet()
        messages = sheet._validate_metadata_dict({})

        exp = [
            'ErrorMessage: Assay is a required attribute',
            'ErrorMessage: Bioinformatics is a required attribute',
            'ErrorMessage: Contact is a required attribute',
        ]

        msg_strs = [str(msg) for msg in messages]
        self.assertEqual(msg_strs, exp)

    def test_validate_sample_sheet_metadata_not_supported(self):
        sheet = AmpliconSampleSheet()
        self.md_ampl['Rush'] = 'XYZ'
        messages = sheet._validate_metadata_dict(self.md_ampl)

        exp = ['ErrorMessage: These metadata keys are not supported: Rush']
        msg_strs = [str(msg) for msg in messages]
        self.assertEqual(msg_strs, exp)

    def test_validate_sample_sheet_metadata_good(self):
        # self.md_ampl is patterned after legacy amplicon sample-sheet.
        sheet = AmpliconSampleSheet()
        messages = sheet._validate_metadata_dict(self.md_ampl)
        self.assertEqual(messages, [])

        # test _validate_metadata_dict() against a
        # MetagenomicSampleSheetv100 object which defines an extra column
        # (contains_replicates) in the Bioinformatics section. Since
        # self.metadata does not contain this extra column, ErrorMessage()s
        # should be returned saying as much.
        sheet = MetagenomicSampleSheetv100()
        messages = sheet._validate_metadata_dict(self.md_metag)

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
            messages = sheet._validate_metadata_dict(self.md_ampl)
            exp = f"ErrorMessage: '{invalid_type}' is not a supported Assay"
            self.assertEqual(str(messages[0]), exp)

    def test_make_sample_sheet(self):
        exp_bfx = pd.DataFrame(self.md_ampl['Bioinformatics'])
        exp_bfx['BarcodesAreRC'] = exp_bfx['BarcodesAreRC'].astype('bool')
        exp_bfx['HumanFiltering'] = exp_bfx['HumanFiltering'].astype('bool')

        exp_contact = pd.DataFrame(self.md_ampl['Contact'])

        # for amplicon we expect the following three columns to not be there
        message = (r'The column (I5_Index_ID|index2|Well_description) '
                   r'in the sample sheet is empty')

        message2_defer_true_prefix = r"ErrorMessage"
        message2_defer_false_prefix = r"Sample sheet instantiation failed"
        message2 = (r": The following "
                    r"projects need to be in the Data and Bioinformatics "
                    r"sections: Koening_ITS_101, THDMI_10317, "
                    r"Yanomani_2008_10052")

        with self.assertWarnsRegex(UserWarning, message):
            table2 = self.table.copy(deep=True)

            # assert that make_sample_sheet() raises an Error when the
            # projects are improperly defined; if defer_validate is True,
            # the error message starts with "ErrorMessage"
            msg2_defer_true = message2_defer_true_prefix + message2
            with self.assertRaisesRegex(ValueError, msg2_defer_true):
                make_sample_sheet(self.md_ampl, table2, 'HiSeq4000', [5, 7],
                                  defer_validate=True)

            # alternately, if defer_validate is False, the error message
            # starts with "Sample sheet instantiation failed".
            msg2_defer_false = message2_defer_false_prefix + message2
            with self.assertRaisesRegex(ValueError, msg2_defer_false):
                make_sample_sheet(self.md_ampl, table2, 'HiSeq4000', [5, 7],
                                  defer_validate=False)

            # second, correct the errors in the [Data] section.
            table2['Project Name'] = ['Koening_ITS_101', 'Yanomani_2008_10052',
                                      'Yanomani_2008_10052']

            obs = make_sample_sheet(self.md_ampl, table2, 'HiSeq4000',
                                    [5, 7])

        self.assertIsInstance(obs, AmpliconSampleSheet)

        self.assertEqual(obs.Reads, [151, 151])
        self.assertEqual(obs.Settings, {'ReverseComplement': '0'})

        pd.testing.assert_frame_equal(obs.Bioinformatics, exp_bfx)
        pd.testing.assert_frame_equal(obs.Contact, exp_contact)

        header = {
            'IEMFileVersion': '4',
            'SheetType': 'dummy_amp',
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
             'AGCCTTCGTCGC', '', '', 'Koening_ITS_101',
             'THDMI_10317_PUK2.X00180471.A1'],
            [5, 'X00180199', 'X00180199', 'THDMI_10317_PUK2', 'C1',
             '515rcbc12', 'CGTATAAATGCG', '', '', 'Yanomani_2008_10052',
             'THDMI_10317_PUK2.X00180199.C1'],
            [5, 'X00179789', 'X00179789', 'THDMI_10317_PUK2', 'E1',
             '515rcbc24', 'TGACTAATGGCC', '', '', 'Yanomani_2008_10052',
             'THDMI_10317_PUK2.X00179789.E1'],
            [7, 'X00180471', 'X00180471', 'THDMI_10317_PUK2', 'A1', '515rcbc0',
             'AGCCTTCGTCGC', '', '', 'Koening_ITS_101',
             'THDMI_10317_PUK2.X00180471.A1'],
            [7, 'X00180199', 'X00180199', 'THDMI_10317_PUK2', 'C1',
             '515rcbc12', 'CGTATAAATGCG', '', '', 'Yanomani_2008_10052',
             'THDMI_10317_PUK2.X00180199.C1'],
            [7, 'X00179789', 'X00179789', 'THDMI_10317_PUK2', 'E1',
             '515rcbc24', 'TGACTAATGGCC', '', '', 'Yanomani_2008_10052',
             'THDMI_10317_PUK2.X00179789.E1'],
        )
        keys = ['Lane', 'Sample_ID', 'Sample_Name', 'Sample_Plate',
                'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                'Sample_Project', 'Well_description']

        for sample, row in zip(obs.samples, data):
            exp = sample_sheet.Sample(dict(zip(keys, row)))
            self.assertEqual(dict(sample), dict(exp))

    def _help_make_test_column_alternatives_table(self):
        table2 = self.table.copy(deep=True)

        table2['Well_description'] = ['Row A', 'Row B', 'Row C']

        table2['Project Name'] = ['Koening_ITS_101', 'Yanomani_2008_10052',
                                  'Yanomani_2008_10052']
        return table2

    def _help_test_column_alternatives(self, obs):
        data = (
            [5, 'X00180471', 'X00180471', 'THDMI_10317_PUK2', 'A1', '515rcbc0',
             'AGCCTTCGTCGC', '', '', 'Koening_ITS_101',
             'THDMI_10317_PUK2.X00180471.A1'],
            [5, 'X00180199', 'X00180199', 'THDMI_10317_PUK2', 'C1',
             '515rcbc12', 'CGTATAAATGCG', '', '', 'Yanomani_2008_10052',
             'THDMI_10317_PUK2.X00180199.C1'],
            [5, 'X00179789', 'X00179789', 'THDMI_10317_PUK2', 'E1',
             '515rcbc24', 'TGACTAATGGCC', '', '', 'Yanomani_2008_10052',
             'THDMI_10317_PUK2.X00179789.E1'],
            [7, 'X00180471', 'X00180471', 'THDMI_10317_PUK2', 'A1', '515rcbc0',
             'AGCCTTCGTCGC', '', '', 'Koening_ITS_101',
             'THDMI_10317_PUK2.X00180471.A1'],
            [7, 'X00180199', 'X00180199', 'THDMI_10317_PUK2', 'C1',
             '515rcbc12', 'CGTATAAATGCG', '', '', 'Yanomani_2008_10052',
             'THDMI_10317_PUK2.X00180199.C1'],
            [7, 'X00179789', 'X00179789', 'THDMI_10317_PUK2', 'E1',
             '515rcbc24', 'TGACTAATGGCC', '', '', 'Yanomani_2008_10052',
             'THDMI_10317_PUK2.X00179789.E1'],
        )
        keys = ['Lane', 'Sample_ID', 'Sample_Name', 'Sample_Plate',
                'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                'Sample_Project', 'Well_description']

        for sample, row in zip(obs.samples, data):
            exp = sample_sheet.Sample(dict(zip(keys, row)))
            self.assertEqual(dict(sample), dict(exp))

    def test_column_alternatives_default(self):
        # confirm standard 'Well_description' column name behaved as intended.
        table2 = self._help_make_test_column_alternatives_table()

        # 'Well_description' column in the input does not cause an error--
        # but it IS silently overwritten!
        obs = make_sample_sheet(self.md_ampl,
                                table2,
                                'HiSeq4000',
                                [5, 7])

        self.assertIsNotNone(obs, msg="make_sample_sheet() failed")
        self.assertIsInstance(obs, AmpliconSampleSheet)

        self._help_test_column_alternatives(obs)

    def test_column_alternatives_err_duplicates(self):
        # Try making sample-sheet w/an alternate column name and confirm that
        # the results continue to be as expected.
        table2 = self._help_make_test_column_alternatives_table()

        table2.rename({'Well_description': 'well_description'},
                      axis=1, inplace=True)
        err_msg = "The remapped sample sheet column names contain duplicates: "
        with self.assertRaisesRegex(ValueError, err_msg):
            _ = make_sample_sheet(self.md_ampl,
                                  table2,
                                  'HiSeq4000',
                                  [5, 7])

    def test_column_alternatives(self):
        # Try making sample-sheet w/an alternate column name and confirm that
        # the results continue to be as expected.
        table2 = self._help_make_test_column_alternatives_table()

        table2.rename({'index2': 'i5 sequence'},
                      axis=1, inplace=True)

        obs = make_sample_sheet(self.md_ampl,
                                table2,
                                'HiSeq4000',
                                [5, 7])

        self._help_test_column_alternatives(obs)

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

            # functionality that handles empty I5_Index_ID and index2 columns,
            # as well as generates Well_description column was migrated up
            # to _remap_table()'s caller, _add_data_to_sheet(). Hence, call
            # this method to ensure that the observed table remains as
            # expected.
            obs = sheet._add_data_to_sheet(self.table, 'HiSeq4000', [1],
                                           'TruSeq HT')
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
                   'i7 plate', SS_SAMPLE_ID_KEY, 'syndna_pool_number',
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

        obs = sheet._remap_table(self.table)

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
                   'i7 plate', SS_SAMPLE_ID_KEY,
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

        sheet = MetatranscriptomicSampleSheetv0()

        obs = sheet._remap_table(self.table)
        obs = obs[['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                   'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                   'Sample_Project', 'Well_description']]

        self.assertEqual(len(obs), 3)
        pd.testing.assert_frame_equal(obs, exp, check_like=True)

    def test_remap_table_metatranscriptomicsv10(self):
        # note that Well_description is now included because it's expected
        # to be inserted by the function that calls _remap_table().
        data = [
            ['33-A1', 'A', 1, True, 'A1', 0, 0, 'AACGCACACTCGTCTT',
             'iTru5_19_A', 'AACGCACA', 'A1', 'iTru5_plate', 'iTru7_109_01',
             'CTCGTCTT', 'A22', 'iTru7_plate', '33-A1', 'The_plate.33-A1.A1',
             '1.2', '1.1'],
            ['820072905-2', 'C', 1, False, 'C1', 1, 1, 'ATGCCTAGCGAACTGT',
             'iTru5_19_B', 'ATGCCTAG', 'B1', 'iTru5_plate', 'iTru7_109_02',
             'CGAACTGT', 'B22', 'iTru7_plate', '820072905-2',
             'The_plate.820072905-2.C1', '1.4', '1.3'],
            ['820029517-3', 'E', 1, False, 'E1', 2, 2, 'CATACGGACATTCGGT',
             'iTru5_19_C', 'CATACGGA', 'C1', 'iTru5_plate', 'iTru7_109_03',
             'CATTCGGT', 'C22', 'iTru7_plate', '820029517-3',
             'The_plate.820029517-3.E1', '1.6', '1.5']
        ]
        columns = ['Sample', 'Row', 'Col', 'Blank', 'Well', 'index',
                   'index combo', 'index combo seq', 'i5 name', 'i5 sequence',
                   'i5 well', 'i5 plate', 'i7 name', 'i7 sequence', 'i7 well',
                   'i7 plate', SS_SAMPLE_ID_KEY, 'Well_description',
                   'vol_extracted_elution_ul', 'total_rna_concentration_ng_ul']
        self.table = pd.DataFrame(data=data, columns=columns)
        self.table['Project Name'] = 'Tst_project_1234'
        self.table['Project Plate'] = 'The_plate'

        columns = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                   'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                   'Sample_Project', 'total_rna_concentration_ng_ul',
                   'vol_extracted_elution_ul', 'Well_description']
        data = [
            ['33-A1', '33-A1', 'The_plate', 'A1', 'iTru7_109_01', 'CTCGTCTT',
             'iTru5_19_A', 'AACGCACA', 'Tst_project_1234', '1.1', '1.2',
             'The_plate.33-A1.A1'],
            ['820072905-2', '820072905-2', 'The_plate', 'C1', 'iTru7_109_02',
             'CGAACTGT', 'iTru5_19_B', 'ATGCCTAG', 'Tst_project_1234', '1.3',
             '1.4', 'The_plate.820072905-2.C1'],
            ['820029517-3', '820029517-3', 'The_plate', 'E1', 'iTru7_109_03',
             'CATTCGGT', 'iTru5_19_C', 'CATACGGA', 'Tst_project_1234', '1.5',
             '1.6', 'The_plate.820029517-3.E1'],
        ]

        exp = pd.DataFrame(columns=columns, data=data)

        sheet = MetatranscriptomicSampleSheetv10()

        obs = sheet._remap_table(self.table)
        obs = obs[['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                   'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                   'Sample_Project', 'total_rna_concentration_ng_ul',
                   'vol_extracted_elution_ul', 'Well_description']]

        self.assertEqual(len(obs), 3)
        pd.testing.assert_frame_equal(obs, exp, check_like=True)

    def test_add_data_to_sheet(self):
        # for amplicon we expect the following three columns to not be there
        message = (r'The column (I5_Index_ID|index2|Well_description) '
                   r'in the sample sheet is empty')

        with self.assertWarnsRegex(UserWarning, message):
            self.sheet._add_data_to_sheet(self.table, 'HiSeq4000', [1],
                                          'TruSeq HT')

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
            'SheetType': 'dummy_amp',
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
        md_metag_w_contains_reps = self._add_replicates_to_md_metag()

        exp_bfx_list = _add_contains_replicates(
            self.md_metag['Bioinformatics'])
        exp_bfx = pd.DataFrame(exp_bfx_list)
        exp_contact = pd.DataFrame(self.md_metag['Contact'])

        obs = sheet._add_metadata_to_sheet(md_metag_w_contains_reps,
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
        md_metag_w_contains_reps = self._add_replicates_to_md_metag()
        md_metag_w_contains_reps['Date'] = '1970-01-01'

        exp_bfx = pd.DataFrame(_add_contains_replicates(
            self.md_metag['Bioinformatics']))
        exp_contact = pd.DataFrame(self.md_metag['Contact'])

        obs = sheet._add_metadata_to_sheet(
            md_metag_w_contains_reps, 'HiSeq4000')

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
        md_metag_w_contains_reps = self._add_replicates_to_md_metag()
        obs = sheet._add_metadata_to_sheet(md_metag_w_contains_reps, 'iSeq')

        settings = {
            'ReverseComplement': '0'
        }

        self.assertEqual(obs.Settings, settings)

    def test_make_sections_dict(self):
        compressed_plate_df = pd.DataFrame([
            {"Sample": "sample1", "Blank": True, "contains_replicates": False,
             "Project Name": "Study_1", "Project Plate": "Study_1_Plate_11"},
            {"Sample": "sample2", "Blank": False, "contains_replicates": False,
             "Project Name": "Study_1", "Project Plate": "Study_1_Plate_11"},
            {"Sample": "sample3", "Blank": False, "contains_replicates": False,
             "Project Name": "Study_4", "Project Plate": "Study_1_Plate_11"},
            {"Sample": "sample4", "Blank": False, "contains_replicates": False,
             "Project Name": "Study_5", "Project Plate": "Study_1_Plate_11"},
            {"Sample": "BLANK.2", "Blank": True, "contains_replicates": False,
             "Project Name": "Study_2", "Project Plate": "Study_2_Plate_21"},
            {"Sample": "sm1", "Blank": False, "contains_replicates": False,
             "Project Name": "Study_2", "Project Plate": "Study_2_Plate_21"},
            {"Sample": "sm2", "Blank": False, "contains_replicates": False,
             "Project Name": "Study_3", "Project Plate": "Study_2_Plate_21"},
            {"Sample": "sm3", "Blank": False, "contains_replicates": False,
             "Project Name": "Study_6", "Project Plate": "Study_2_Plate_21"},
            {"Sample": "BLANK.3", "Blank": True, "contains_replicates": False,
             "Project Name": "Study_3", "Project Plate": "Study_3_Plate_13"},
            {"Sample": "samp1", "Blank": False, "contains_replicates": False,
             "Project Name": "Study_3", "Project Plate": "Study_3_Plate_13"},
            {"Sample": "blank4", "Blank": True, "contains_replicates": False,
             "Project Name": "Study_10", "Project Plate": "Study_10_Plate_1"},
            {"Sample": "samples1", "Blank": False,
             "contains_replicates": False,
             "Project Name": "Study_10", "Project Plate": "Study_10_Plate_1"},
            {"Sample": "samples2", "Blank": False,
             "contains_replicates": False,
             "Project Name": "Study_11", "Project Plate": "Study_10_Plate_1"}
        ])

        studies_info = [
            {
                'Project Name': 'Study_1',
                'Project Abbreviation': 'ADAPT',
                'sample_accession_fp': './adapt_sa.tsv',
                'qiita_metadata_fp': './adapt_metadata.txt',
                'experiment_design_description': 'isolate sequencing',
                'HumanFiltering': 'False',
                'Email': 'r@gmail.com'
            },
            {
                'Project Name': 'Study_2',
                'Project Abbreviation': 'CHILD',
                'sample_accession_fp': './child_sa.tsv',
                'qiita_metadata_fp': './child_metadata.txt',
                'experiment_design_description': 'whole genome sequencing',
                'HumanFiltering': 'True',
                'Email': 'l@ucsd.edu'
            },
            {
                'Project Name': 'Study_3',
                'Project Abbreviation': 'MARMO',
                'sample_accession_fp': './marmo_sa.tsv',
                'qiita_metadata_fp': './marmo_metadata.txt',
                'experiment_design_description': 'whole genome sequencing',
                'HumanFiltering': 'False',
                'Email': 'c@ucsd.edu'
            },
            {
                'Project Name': 'Study_4',
                'Project Abbreviation': 'MEME',
                'sample_accession_fp': './meme_sa.tsv',
                'qiita_metadata_fp': './meme_metadata.txt',
                'experiment_design_description': 'whole genome sequencing',
                'HumanFiltering': 'False',
                'Email': 'b@ucsd.edu'
            }
        ]

        bioinfo_base = {
            'ForwardAdapter': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
            'ReverseAdapter': 'GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT',
            'library_construction_protocol': 'Knight Lab Kapa HyperPlus',
            'BarcodesAreRC': 'True'
        }

        exp = {
            'Experiment Name': 'RKL001',
            'SheetType': 'standard_metag',
            'SheetVersion': '101',
            'Assay': 'Metagenomic',
            'Bioinformatics': [
                {'ForwardAdapter': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                 'ReverseAdapter': 'GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT',
                 'library_construction_protocol': 'Knight Lab Kapa HyperPlus',
                 'BarcodesAreRC': 'True', 'Sample_Project': 'Study_1',
                 'QiitaID': '1', 'HumanFiltering': 'False',
                 'experiment_design_description': 'isolate sequencing',
                 'contains_replicates': False},
                {'ForwardAdapter': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                 'ReverseAdapter': 'GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT',
                 'library_construction_protocol': 'Knight Lab Kapa HyperPlus',
                 'BarcodesAreRC': 'True', 'Sample_Project': 'Study_2',
                 'QiitaID': '2', 'HumanFiltering': 'True',
                 'experiment_design_description': 'whole genome sequencing',
                 'contains_replicates': False},
                {'ForwardAdapter': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                 'ReverseAdapter': 'GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT',
                 'library_construction_protocol': 'Knight Lab Kapa HyperPlus',
                 'BarcodesAreRC': 'True', 'Sample_Project': 'Study_3',
                 'QiitaID': '3', 'HumanFiltering': 'False',
                 'experiment_design_description': 'whole genome sequencing',
                 'contains_replicates': False},
                {'ForwardAdapter': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                 'ReverseAdapter': 'GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT',
                 'library_construction_protocol': 'Knight Lab Kapa HyperPlus',
                 'BarcodesAreRC': 'True', 'Sample_Project': 'Study_4',
                 'QiitaID': '4', 'HumanFiltering': 'False',
                 'experiment_design_description': 'whole genome sequencing',
                 'contains_replicates': False}
            ],
            'Contact': [
                {'Sample_Project': 'Study_1', 'Email': 'r@gmail.com'},
                {'Sample_Project': 'Study_2', 'Email': 'l@ucsd.edu'},
                {'Sample_Project': 'Study_3', 'Email': 'c@ucsd.edu'},
                {'Sample_Project': 'Study_4', 'Email': 'b@ucsd.edu'}
            ],

            # when there are multiple secondary studies, they are delimited
            # by a ";" (no spaces).  When there are NO secondary studies,
            # the value of the `secondary_qiita_studies` key is an empty string
            # (NOT a None).
            'SampleContext': [
                {'sample_name': 'sample1', 'primary_qiita_study': '1',
                 'sample_type': 'control blank',
                 'secondary_qiita_studies': '4;5'},
                {'sample_name': 'BLANK.2', 'primary_qiita_study': '2',
                 'sample_type': 'control blank',
                 'secondary_qiita_studies': '3;6'},
                {'sample_name': 'BLANK.3', 'primary_qiita_study': '3',
                 'sample_type': 'control blank',
                 'secondary_qiita_studies': ''},
                {'sample_name': 'blank4', 'primary_qiita_study': '10',
                 'sample_type': 'control blank',
                 'secondary_qiita_studies': '11'}
            ]
        }

        obs = make_sections_dict(
            compressed_plate_df, studies_info, "RKL001",
            'standard_metag', '101', bioinfo_base)
        self.assertDictEqual(exp, obs)


class ValidateSampleSheetTests(BaseTests):
    maxDiff = None

    # def test_init_defer_validate_warning(self):
    #     # Assure DeprecationWarning for sheet init w/o defer_validate param
    #
    #     with self.assertWarns(DeprecationWarning):
    #         _ = MetagenomicSampleSheetv90(self.good_ss)
    #
    # def test_load_sample_sheet_defer_validate_warning(self):
    #     # Assure DeprecationWarning for load_sample_sheet w/o defer_validate
    #
    #     with self.assertWarns(DeprecationWarning):
    #         _ = load_sample_sheet(self.good_ss)

    def test_init_default_validate_err(self):
        # Using default validation setting raises error
        err_msg = ("Sample sheet instantiation failed: The "
                   "Sample_Project column in the Data section is missing")
        with self.assertRaisesRegex(ValueError, err_msg):
            _ = MetagenomicSampleSheetv100(self.no_project_ss)

    def test_init_explicit_validate_err(self):
        # Explicitly requesting validation raises error
        err_msg = ("Sample sheet instantiation failed: The "
                   "Sample_Project column in the Data section is missing")
        with self.assertRaisesRegex(ValueError, err_msg):
            _ = MetagenomicSampleSheetv100(self.no_project_ss,
                                           defer_validate=False)

    def test_init_explicit_no_validate(self):
        # Explicitly deferring validation allows load without error
        sheet = MetagenomicSampleSheetv100(self.no_project_ss,
                                           defer_validate=True)
        self.assertEqual(MetagenomicSampleSheetv100, type(sheet))

    def test_validate_and_scrub_sample_sheet(self):
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)
        # no errors
        self.assertTrue(sheet.validate_and_scrub_sample_sheet())

    def test_quiet_validate_and_scrub_sample_sheet(self):
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)

        buffer = StringIO()
        with redirect_stdout(buffer):
            msgs = sheet.quiet_validate_and_scrub_sample_sheet()

        # no errors
        self.assertEqual(buffer.getvalue().strip(), '')
        self.assertEqual(msgs, [])

    def test_quiet_validate_and_scrub_sample_sheet_w_context(self):
        sheet = MetagenomicSampleSheetv101(
            self.good_metag_ss_w_context, defer_validate=False)

        buffer = StringIO()
        with redirect_stdout(buffer):
            msgs = sheet.quiet_validate_and_scrub_sample_sheet()

        # no errors
        self.assertEqual(buffer.getvalue().strip(), '')
        self.assertEqual(msgs, [])

    def test_validate_and_scrub_sample_sheet_no_sample_project(self):
        sheet = MetagenomicSampleSheetv100(
            self.no_project_ss, defer_validate=True)

        buffer = StringIO()
        with redirect_stdout(buffer):
            self.assertFalse(sheet.validate_and_scrub_sample_sheet())

        self.assertEqual(buffer.getvalue().strip(),
                         'ErrorMessage: The Sample_Project column in the'
                         ' Data section is missing')

    def test_quiet_validate_and_scrub_sample_sheet_no_sample_project(self):
        sheet = MetagenomicSampleSheetv100(
            self.no_project_ss, defer_validate=True)

        buffer = StringIO()
        with redirect_stdout(buffer):
            msgs = sheet.quiet_validate_and_scrub_sample_sheet()

        self.assertEqual(buffer.getvalue().strip(), '')
        msg_strs = [str(msg) for msg in msgs]
        self.assertEqual(msg_strs, ['ErrorMessage: The Sample_Project column'
                                    ' in the Data section is missing'])

    def test_validate_and_scrub_sample_sheet_missing_bioinformatics(self):
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=True)
        sheet.Bioinformatics = None

        buffer = StringIO()
        with redirect_stdout(buffer):
            self.assertFalse(sheet.validate_and_scrub_sample_sheet())

        self.assertEqual(
            buffer.getvalue().strip(),
            'ErrorMessage: The Bioinformatics section cannot be missing')

    def test_quiet_validate_scrub_sample_sheet_missing_bioinformatics(self):
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=True)
        sheet.Bioinformatics = None

        buffer = StringIO()
        with redirect_stdout(buffer):
            msgs = sheet.quiet_validate_and_scrub_sample_sheet()

        self.assertEqual(buffer.getvalue().strip(), '')
        msg_strs = [str(msg) for msg in msgs]
        self.assertEqual(
            msg_strs,
            ['ErrorMessage: The Bioinformatics section cannot be missing'])

    def test_validate_and_scrub_sample_sheet_missing_contact(self):
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=True)
        sheet.Contact = None

        buffer = StringIO()
        with redirect_stdout(buffer):
            self.assertFalse(sheet.validate_and_scrub_sample_sheet())

        self.assertEqual(buffer.getvalue().strip(),
                         'ErrorMessage: The Contact section cannot be missing')

    def test_validate_and_scrub_sample_sheet_scrubbed_names(self):
        sheet = MetagenomicSampleSheetv90(
            self.scrubbable_ss, defer_validate=True)

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

        buffer = StringIO()
        with redirect_stdout(buffer):
            self.assertTrue(sheet.validate_and_scrub_sample_sheet())

        self.assertEqual(buffer.getvalue().strip(), message)

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

        sheet = MetagenomicSampleSheetv90(
            self.scrubbable_ss, defer_validate=False)

        buffer = StringIO()
        with redirect_stdout(buffer):
            msgs = sheet.quiet_validate_and_scrub_sample_sheet()

        self.assertEqual(buffer.getvalue().strip(), '')
        self.assertEqual(msgs, [message])

    def test_validate_and_scrub_sample_sheet_scrubbed_project_names(self):
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=True)

        remapper = {
            'ProjectN_13059': "NYU's Tisch Art Microbiome 13059",
            'ProjectF_11661': "The x.x microbiome project 1337"
        }

        for sample in sheet:
            sample['Sample_Project'] = remapper.get(sample['Sample_Project'],
                                                    sample['Sample_Project'])

        # new pandas won't let you set value inplace on a slice
        project_remapper = {'Sample_Project': remapper}
        sheet.Contact.replace(project_remapper, inplace=True)
        sheet.Bioinformatics.replace(project_remapper, inplace=True)

        message = (
            'WarningMessage: The following project names were scrubbed for '
            'bcl2fastq compatibility. If the same invalid characters are also '
            'found in the Bioinformatics and Contact sections, those will be '
            'automatically scrubbed too:\n'
            "NYU's Tisch Art Microbiome 13059, The x.x microbiome project 1337"
        )

        buffer = StringIO()
        with redirect_stdout(buffer):
            sheet.validate_and_scrub_sample_sheet()

        self.assertEqual(buffer.getvalue().strip(), message)

        scrubbed = {
            'NYU_s_Tisch_Art_Microbiome_13059',
            'The_x_x_microbiome_project_1337',
            'ProjectG_6123'
        }

        for sample in sheet:
            self.assertTrue(sample['Sample_Project'] in scrubbed,
                            sample['Sample_Project'])

        for project in sheet.Bioinformatics.Sample_Project:
            self.assertTrue(project in scrubbed)

        for project in sheet.Contact.Sample_Project:
            self.assertTrue(project in scrubbed)

    def test_validate_and_scrub_sample_sheet_bad_project_names(self):
        sheet = MetagenomicSampleSheetv100(
            self.bad_project_name_ss, defer_validate=True)

        message = ('ErrorMessage: The following project names in the '
                   'Sample_Project column are missing a Qiita study '
                   'identifier: ProjectF, ProjectG')

        buffer = StringIO()
        with redirect_stdout(buffer):
            self.assertFalse(sheet.validate_and_scrub_sample_sheet())

        self.assertEqual(buffer.getvalue().strip(), message)

    def test_validate_and_scrub_sample_sheet_project_missing_lane(self):
        sheet = MetagenomicSampleSheetv90(self.good_ss, defer_validate=True)

        # set the lane value as empty for one of the two projects
        for sample in sheet.samples:
            if sample.Sample_Project == 'ProjectF_11661':
                sample.Lane = ' '

        message = ('ErrorMessage: The following projects are missing a Lane '
                   'value: ProjectF_11661')

        buffer = StringIO()
        with redirect_stdout(buffer):
            self.assertFalse(sheet.validate_and_scrub_sample_sheet())

        self.assertEqual(buffer.getvalue().strip(), message)

    def test_validate_and_scrub_sample_sheet_missing_project_names(self):
        sheet = MetagenomicSampleSheetv101(
            self.good_metag_ss_w_context, defer_validate=True)
        # pick a random entry in the sample context section and set its
        # qiita study id to something that doesn't exist in the other metadata
        a_blank_mask = sheet.SampleContext[SAMPLE_NAME_KEY] == "BLANK1_1A"
        sheet.SampleContext.loc[a_blank_mask, PRIMARY_STUDY_KEY] = "123456"

        message = ("ErrorMessage: The following projects were only found in "
                   "the SampleContext section: 123456. Projects need to be "
                   "listed in the Data and Bioinformatics section in order to "
                   "be included in the SampleContext section.")

        buffer = StringIO()
        with redirect_stdout(buffer):
            self.assertFalse(sheet.validate_and_scrub_sample_sheet())

        self.assertEqual(buffer.getvalue().strip(), message)

    def test_sample_sheet_to_dataframe(self):
        ss = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)
        obs = sample_sheet_to_dataframe(ss)

        # first, confirm that the function returns a DataFrame
        self.assertTrue(isinstance(obs, pd.DataFrame))

        # confirm that function returns the [Data] section of the sample-sheet
        # as a DataFrame.
        exp_columns = {'lane', 'sample_name', 'sample_plate', 'sample_well',
                       'i7_index_id', 'index', 'i5_index_id', 'index2',
                       'sample_project', 'well_description',
                       'library_construction_protocol',
                       'experiment_design_description'}
        self.assertEqual(set(obs.columns), exp_columns)

        # since good-sample-sheet.csv contains many samples, just check for
        # the below three.
        exp_sample_names = ["3A", "EP981129A02",
                            "JM-MEC__Staphylococcus_aureusstrain_BERTI-B0387"]
        obs_names = obs['sample_name'].to_list()
        for exp_name in exp_sample_names:
            self.assertIn(exp_name, obs_names)

    def test_sample_sheet_to_dataframe_no_lcase(self):
        ss = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)
        obs = sample_sheet_to_dataframe(ss, lcase_cols=False)

        # first, confirm that the function returns a DataFrame
        self.assertTrue(isinstance(obs, pd.DataFrame))

        # confirm that function returns the [Data] section of the sample-sheet
        # as a DataFrame.
        exp_columns = {'Lane', 'Sample_Name', 'Sample_Plate', 'Sample_Well',
                       'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                       'Sample_Project', 'Well_description',
                       'library_construction_protocol',
                       'experiment_design_description'}
        self.assertEqual(set(obs.columns), exp_columns)

        # since good-sample-sheet.csv contains many samples, just check for
        # the below three.
        exp_sample_names = ["3A", "EP981129A02",
                            "JM-MEC__Staphylococcus_aureusstrain_BERTI-B0387"]
        obs_names = obs['Sample_Name'].to_list()
        for exp_name in exp_sample_names:
            self.assertIn(exp_name, obs_names)

    def test_sample_sheet_to_dataframe_no_lcase_no_protocols(self):
        ss = MetagenomicSampleSheetv90(self.good_ss, defer_validate=False)
        obs = sample_sheet_to_dataframe(
            ss, lcase_cols=False, add_protocol_info=False)

        # first, confirm that the function returns a DataFrame
        self.assertTrue(isinstance(obs, pd.DataFrame))

        # confirm that function returns the [Data] section of the sample-sheet
        # as a DataFrame.
        exp_columns = {'Lane', 'Sample_Name', 'Sample_Plate', 'Sample_Well',
                       'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                       'Sample_Project', 'Well_description'}
        self.assertEqual(set(obs.columns), exp_columns)

        # since good-sample-sheet.csv contains many samples, just check for
        # the below three.
        exp_sample_names = ["3A", "EP981129A02",
                            "JM-MEC__Staphylococcus_aureusstrain_BERTI-B0387"]
        obs_names = obs['Sample_Name'].to_list()
        for exp_name in exp_sample_names:
            self.assertIn(exp_name, obs_names)

    def test_boolean_column_handling(self):
        sheet = MetagenomicSampleSheetv100(
            self.good_w_bools, defer_validate=False)

        # self.good_w_bools contains a [Bioinformatics] section w/multiple
        # projects and values for BarcodesAreRC and HumanFiltering columns that
        # are strings of mixed-case. Demonstrate that no matter the case, the
        # values are properly converted to bools, False for the former and
        # True for the latter.

        obs = set(sheet.Bioinformatics['BarcodesAreRC'])
        self.assertEqual(obs, {False})

        obs = set(sheet.Bioinformatics['HumanFiltering'])
        self.assertEqual(obs, {True})


class ProfileTests(BaseTests):
    def test_profile(self):
        sheet = AbsQuantSampleSheetv10()

        # confirm that AbsQuantSampleSheetv10() contains the right values
        # for sheet-type and sheet-version, not the default values inherited
        # from its parent.
        self.assertEqual(sheet._HEADER['SheetType'], 'abs_quant_metag')
        self.assertEqual(sheet._HEADER['SheetVersion'], '10')
        self.assertIn('mass_syndna_input_ng', sheet._data_columns)
        self.assertNotIn('SampleContext', sheet._ordered_section_keys)

    def test_profile_absquant_11(self):
        sheet = AbsQuantSampleSheetv11()

        # confirm that AbsQuantSampleSheetv10() contains the right values
        # for sheet-type and sheet-version, not the default values inherited
        # from its parent.
        self.assertEqual(sheet._HEADER['SheetType'], 'abs_quant_metag')
        self.assertEqual(sheet._HEADER['SheetVersion'], '11')
        self.assertIn('mass_syndna_input_ng', sheet._data_columns)
        self.assertIn('SampleContext', sheet._ordered_section_keys)


class DemuxReplicatesTests(BaseTests):
    def setUp(self):
        self.data_dir = join(dirname(__file__), 'data')
        self.sheet_w_replicates_path = \
            join(self.data_dir, 'good_standard_metagv100_w_replicates.csv')

        # bad_sheet_w_replicates.csv contains two projects, one of which
        # doesn't contain replicates. By convention, all projects in the sheet
        # must either contain replicates or not contain replicates.
        self.bad_sht_w_replicates_path = join(self.data_dir,
                                              'bad_sheet_w_replicates.csv')

        self.sheet_wo_replicates_path = join(self.data_dir,
                                             'sheet_wo_replicates.csv')

        self.legacy_sheet_path = \
            join(self.data_dir, 'good_standard_metagv90.csv')

        self.replicate_output_paths = [join(self.data_dir,
                                            'replicate_output1.csv'),
                                       join(self.data_dir,
                                            'replicate_output2.csv'),
                                       join(self.data_dir,
                                            'replicate_output3.csv')]

    def _help_test_demux_sample_sheet(
            self, sheet_class, input_path, output_paths):

        # read in a sample sheet containing replicates
        sheet = sheet_class(input_path, defer_validate=False)

        # demux and write out the demuxed sample-sheets
        results = demux_sample_sheet(sheet)

        # assert that the proper number of demuxed sample sheets were returned
        num_expected = len(output_paths)
        self.assertEqual(len(results), num_expected)

        # compare the completed sample-sheets against an expected results.
        # among other things confirm that 'orig_name' is not in the
        # output replicate csvs, indicating it has become the 'Sample_Name'
        # column for that replicate sample-sheet.
        self.maxDiff = None
        for curr_path_num in range(num_expected):
            with tempfile.NamedTemporaryFile(mode='w+') as tmp:
                results[curr_path_num].write(tmp)
                tmp.flush()
                self._help_test_csv_files_exact_text_match(
                    output_paths[curr_path_num], tmp.name)

    def test_sheet_needs_demuxing(self):
        # confirm legacy sample-sheets w/out contains_replicates column will
        # return False, instead of raising an Error. For processing purposes,
        # it's only critical to know whether the sheet needs demuxing or not.
        sheet = MetagenomicSampleSheetv90(
            self.legacy_sheet_path, defer_validate=False)
        self.assertFalse(sheet_needs_demuxing(sheet))

        # confirm bad sample-sheet raises a ValueError for containing projects
        # that contain replicates and those that don't.
        with self.assertRaisesRegex(ValueError, "All projects in "
                                                "Bioinformatics section must "
                                                "either contain replicates or "
                                                "not."):
            sheet = MetagenomicSampleSheetv100(
                self.bad_sht_w_replicates_path, defer_validate=False)
            sheet_needs_demuxing(sheet)

        # test a valid sample-sheet with replicates.
        sheet = MetagenomicSampleSheetv100(
            self.sheet_w_replicates_path, defer_validate=False)
        self.assertTrue(sheet_needs_demuxing(sheet))

        # test a valid sample-sheet w/out replicates.
        sheet = MetagenomicSampleSheetv100(
            self.sheet_wo_replicates_path, defer_validate=False)
        self.assertFalse(sheet_needs_demuxing(sheet))

    def test_sheet_needs_demuxing_no_replicates_support(self):
        # test valid sheet of a type that doesn't even support replicates.
        metat_fp = join(
            self.data_dir,
            MetatranscriptomicSampleSheetv10CreationTests.sample_sheet_name)

        sheet = MetatranscriptomicSampleSheetv10(
            metat_fp, defer_validate=False)
        self.assertFalse(sheet_needs_demuxing(sheet))

    def test_demux_sample_sheet_err_no_contains_replicates(self):
        # we don't want to demux legacy sample-sheets. sheet_needs_demuxing()
        # should be used to determine if demux_sample_sheet() should be
        # called.
        sheet = MetagenomicSampleSheetv90(
            self.legacy_sheet_path, defer_validate=False)
        err_msg = ("sample sheet does not have a 'contains_replicates' "
                   "column in the 'Bioinformatics' section.")
        with self.assertRaisesRegex(ValueError, err_msg):
            demux_sample_sheet(sheet)

    def test_demux_sample_sheet_err_contains_replicates_inconsistent(self):
        # by convention, all replication is done at the plate level, and all
        # projects in a sample-sheet will either contain replicates, or all of
        # them will not. Hence, a sample-sheet with both True and False in
        # the contains_replicates column in the [Bioinformatics] section should
        # raise an error. Note that this error would be caught at the
        # point the sheet was created if the defer_validate flag was False!
        sheet = MetagenomicSampleSheetv100(
            self.bad_sht_w_replicates_path, defer_validate=True)
        err = ("All projects in Bioinformatics section must either contain"
               " replicates or not.")
        with self.assertRaisesRegex(ValueError, err):
            demux_sample_sheet(sheet)

    def test_demux_sample_sheet_err_contains_replicates_false(self):
        # as mentioned above, sheet_needs_demuxing() should be used to
        # determine if demux_sample_sheet() should be called. If a sample
        # sheet is passed to demux_sample_sheet() and all projects are False,
        # an Error should be raised to alert the user of an unexpected
        # condition, rather than silently allow as a degenerative case.
        sheet = MetagenomicSampleSheetv100(
            self.sheet_wo_replicates_path, defer_validate=False)
        err = "No projects in Bioinformatics section contain replicates"
        with self.assertRaisesRegex(ValueError, err):
            demux_sample_sheet(sheet)

    def test_demux_sample_sheet(self):
        self._help_test_demux_sample_sheet(
            MetagenomicSampleSheetv100, self.sheet_w_replicates_path,
            self.replicate_output_paths)

    def test_demux_sample_sheet_w_context(self):
        # this test will need to compare the four completed sample-sheets
        # made using self.sheet_w_replicates_path against an expected result.
        demux_sheet_w_context_path = join(
            self.data_dir, "good_standard_metagv101_w_replicates.csv")
        context_output_paths = [x.replace('.csv', '_w_context.csv') for
                                x in self.replicate_output_paths]
        self._help_test_demux_sample_sheet(
            MetagenomicSampleSheetv101, demux_sheet_w_context_path,
            context_output_paths)


class AdditionalSampleSheetCreationTests(BaseTests):
    def _fill_test_metagenomic_sheet(self, sheet, bfx=None, data=None):
        sheet.Header['IEMFileVersion'] = 4
        sheet.Header['SheetType'] = 'standard_metag'
        sheet.Header['SheetVersion'] = '100'
        sheet.Header['Investigator Name'] = 'Knight'
        sheet.Header['Experiment Name'] = 'RKO_experiment'
        sheet.Header['Date'] = '2021-08-17'
        sheet.Header['Workflow'] = 'GenerateFASTQ'
        sheet.Header['Application'] = 'FASTQ Only'
        sheet.Header['Assay'] = 'Metagenomic'
        sheet.Header['Description'] = ''
        sheet.Header['Chemistry'] = 'Default'
        sheet.Reads = [151, 151]
        sheet.Settings['ReverseComplement'] = 0

        if not bfx:
            bfx = [
                ['Project1_99999', '99999', 'False', 'AACC', 'GGTT', 'False',
                 'False', 'protocol_1', 'a designed experiment']
            ]

        sheet.Bioinformatics = pd.DataFrame(
            columns=['Sample_Project', 'QiitaID', 'BarcodesAreRC',
                     'ForwardAdapter', 'ReverseAdapter', 'HumanFiltering',
                     'contains_replicates', 'library_construction_protocol',
                     'experiment_design_description'], data=bfx)

        sheet.Contact = pd.DataFrame(columns=['Email', 'Sample_Project'],
                                     data=[['c2cowart@ucsd.edu',
                                            'Project1_99999'],])

        header = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                  'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                  'Sample_Project', 'Well_description']

        if not data:
            data = [
                ['sample_1', 'sample.1', 'sample_plate_1', 'A1',
                 'iTru7_107_07', 'CCGACTAT', 'iTru5_01_A', 'ACCGACAA',
                 'Project1_99999', 'desc'],
                ['sample_2', 'sample.2', 'sample_plate_1', 'A2',
                 'iTru7_107_07', 'CCGACTAC', 'iTru5_01_A', 'ACCGACAT',
                 'Project1_99999', 'desc'],
                ['sample_3', 'sample.3', 'sample_plate_1', 'A3',
                 'iTru7_107_07', 'CCGACTAG', 'iTru5_01_A', 'ACCGACAG',
                 'Project1_99999', 'desc'],
            ]

        for row in data:
            # Add each row as a Sample() object. Each Sample() object takes
            # a dict as its initializer.
            sheet.add_sample(sample_sheet.Sample(dict(zip(header, row))))

        return sheet

    def test_metatranscriptomic_sheet_creation(self):
        # create a Metatranscriptomic-type sample-sheet from scratch and
        # manually populate the required fields.
        sheet = MetatranscriptomicSampleSheetv0()
        sheet.Header['IEMFileVersion'] = 4
        sheet.Header['SheetType'] = 'standard_metag'
        sheet.Header['SheetVersion'] = '0'
        sheet.Header['Investigator Name'] = 'Knight'
        sheet.Header['Experiment Name'] = 'RKO_experiment'
        sheet.Header['Date'] = '2021-08-17'
        sheet.Header['Workflow'] = 'GenerateFASTQ'
        sheet.Header['Application'] = 'FASTQ Only'
        sheet.Header['Assay'] = 'Metatranscriptomic'
        sheet.Header['Description'] = ''
        sheet.Header['Chemistry'] = 'Default'
        sheet.Reads = [151, 151]
        sheet.Settings['ReverseComplement'] = 0

        data = [
            ['Project1_99999', '99999', 'False', 'AACC', 'GGTT', 'False',
             'False', 'protocol_1', 'a designed experiment']
        ]

        sheet.Bioinformatics = pd.DataFrame(
            columns=['Sample_Project', 'QiitaID', 'BarcodesAreRC',
                     'ForwardAdapter', 'ReverseAdapter', 'HumanFiltering',
                     'contains_replicates', 'library_construction_protocol',
                     'experiment_design_description'], data=data)

        sheet.Contact = pd.DataFrame(columns=['Email', 'Sample_Project'],
                                     data=[['c2cowart@ucsd.edu',
                                            'Project1_99999'],])

        header = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                  'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                  'Sample_Project', 'Well_description']

        data = [
            ['sample_1', 'sample.1', 'sample_plate_1', 'A1', 'iTru7_107_07',
             'CCGACTAT', 'iTru5_01_A', 'ACCGACAA', 'Project1_99999', 'desc'],
            ['sample_2', 'sample.2', 'sample_plate_1', 'A2', 'iTru7_107_07',
             'CCGACTAC', 'iTru5_01_A', 'ACCGACAT', 'Project1_99999', 'desc'],
            ['sample_3', 'sample.3', 'sample_plate_1', 'A3', 'iTru7_107_07',
             'CCGACTAG', 'iTru5_01_A', 'ACCGACAG', 'Project1_99999', 'desc'],
        ]

        for row in data:
            # Add each row as a Sample() object. Each Sample() object takes
            # a dict as its initializer.
            sheet.add_sample(sample_sheet.Sample(dict(zip(header, row))))

        # Once sheet has been manually populated, validate it.
        self.assertTrue(sheet.validate_and_scrub_sample_sheet())

    def test_metatranscriptomic_sheet_creationv10(self):
        # create a Metatranscriptomic-type sample-sheet from scratch and
        # manually populate the required fields.
        sheet = MetatranscriptomicSampleSheetv10()
        sheet.Header['IEMFileVersion'] = 4
        sheet.Header['SheetType'] = 'standard_metat'
        sheet.Header['SheetVersion'] = '10'
        sheet.Header['Investigator Name'] = 'Knight'
        sheet.Header['Experiment Name'] = 'RKO_experiment'
        sheet.Header['Date'] = '2021-08-17'
        sheet.Header['Workflow'] = 'GenerateFASTQ'
        sheet.Header['Application'] = 'FASTQ Only'
        sheet.Header['Assay'] = 'Metatranscriptomic'
        sheet.Header['Description'] = ''
        sheet.Header['Chemistry'] = 'Default'
        sheet.Reads = [151, 151]
        sheet.Settings['ReverseComplement'] = 0

        data = [
            ['Project1_99999', '99999', 'False', 'AACC', 'GGTT', 'False',
             'False', 'protocol_1', 'a designed experiment']
        ]

        sheet.Bioinformatics = pd.DataFrame(
            columns=['Sample_Project', 'QiitaID', 'BarcodesAreRC',
                     'ForwardAdapter', 'ReverseAdapter', 'HumanFiltering',
                     'contains_replicates', 'library_construction_protocol',
                     'experiment_design_description'], data=data)

        sheet.Contact = pd.DataFrame(columns=['Email', 'Sample_Project'],
                                     data=[['c2cowart@ucsd.edu',
                                            'Project1_99999'],])

        header = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
                  'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                  'Sample_Project', 'total_rna_concentration_ng_ul',
                  'vol_extracted_elution_ul', 'Well_description']

        data = [
            ['sample_1', 'sample.1', 'sample_plate_1', 'A1', 'iTru7_107_07',
             'CCGACTAT', 'iTru5_01_A', 'ACCGACAA', 'Project1_99999', '1.0',
             '1.1', 'desc'],
            ['sample_2', 'sample.2', 'sample_plate_1', 'A2', 'iTru7_107_07',
             'CCGACTAC', 'iTru5_01_A', 'ACCGACAT', 'Project1_99999', '1.0',
             '1.1', 'desc'],
            ['sample_3', 'sample.3', 'sample_plate_1', 'A3', 'iTru7_107_07',
             'CCGACTAG', 'iTru5_01_A', 'ACCGACAG', 'Project1_99999', '1.0',
             '1.1', 'desc'],
        ]

        for row in data:
            # Add each row as a Sample() object. Each Sample() object takes
            # a dict as its initializer.
            sheet.add_sample(sample_sheet.Sample(dict(zip(header, row))))

        # Once sheet has been manually populated, validate it.
        self.assertTrue(sheet.validate_and_scrub_sample_sheet())

    def test_metatranscriptomic_sheet_load(self):
        metat_fp = join(
            self.data_dir,
            MetatranscriptomicSampleSheetv10CreationTests.sample_sheet_name)

        # confirm manual loading is w/out error.
        sheet = MetatranscriptomicSampleSheetv10(
            metat_fp, defer_validate=False)
        self.assertTrue(sheet.validate_and_scrub_sample_sheet())

        # Metat v10 should NOT include contains_replicates
        self.assertFalse(
            CONTAINS_REPLICATES_KEY in sheet.Bioinformatics.columns)

        # confirm load_sample_sheet() returns the correct child class of
        # KLSampleSheet.
        sheet = load_sample_sheet(metat_fp, defer_validate=False)
        self.assertIsInstance(sheet, MetatranscriptomicSampleSheetv10)

    def test_metagenomic_sheet_creation(self):
        # create a Metagenomic-type sample-sheet from scratch and manually
        # populate the required fields.
        sheet = MetagenomicSampleSheetv100()
        sheet = self._fill_test_metagenomic_sheet(sheet)
        sheet.Header['SheetVersion'] = '100'

        # Once sheet has been manually populated, validate it.
        self.assertTrue(sheet.validate_and_scrub_sample_sheet())

        # Insert a few errors into the sample-sheet to ensure it fails
        # validation.
        del (sheet.Header['Workflow'])
        sheet.Header['Assay'] = 'NotMetagenomic'

        obs = sheet.quiet_validate_and_scrub_sample_sheet()

        # convert ErrorMessages and WarningMessages into text strings for
        # testing.
        obs = set([str(msg) for msg in obs])

        exp = {"ErrorMessage: 'Workflow' is not declared in Header section",
               "ErrorMessage: 'Assay' value is not 'Metagenomic'"}

        self.assertEqual(obs, exp)

    def test_metagenomic_sheet_w_context_creation(self):
        # create a Metagenomic-type sample sheet that contains a
        # SampleContext section from scratch and manually
        # populate the required fields.
        data = [
            ['sample_1', 'sample.1', 'sample_plate_1', 'A1', 'iTru7_107_07',
             'CCGACTAT', 'iTru5_01_A', 'ACCGACAA', 'Project1_99999', 'desc'],
            ['sample_2', 'sample.2', 'sample_plate_1', 'A2', 'iTru7_107_07',
             'CCGACTAC', 'iTru5_01_A', 'ACCGACAT', 'Project2_14577', 'desc'],
            ['sample_3', 'sample.3', 'sample_plate_1', 'A3', 'iTru7_107_07',
             'CCGACTAG', 'iTru5_01_A', 'ACCGACAG', 'Project3_10317', 'desc'],
            ['BLANK_CHILD_1000_G4', 'BLANK.CHILD.1000.G4',
             'sample_plate_1', 'G4', 'iTru7_107_07',
             'CCGACTAA', 'iTru5_01_A', 'ACCGACAG',
             'Project1_99999', 'desc'],
            ['sample_4', 'sample.4', 'sample_plate_2', 'A3', 'iTru7_107_07',
             'CCGACTCT', 'iTru5_01_A', 'ACCGACAG', 'Project4_11223', 'desc'],
            ['BLANK_OTHER_10_H4', 'BLANK.OTHER.10.H4',
             'sample_plate_2', 'H4', 'iTru7_107_07',
             'CCGACTCC', 'iTru5_01_A', 'ACCGACAG',
             'Project4_11223', 'desc'],
        ]

        bfx = [
            ['Project1_99999', '99999', 'False', 'AACC', 'GGTT', 'False',
             'False', 'protocol_1', 'a designed experiment'],
            ['Project2_14577', '14577', 'False', 'AACC', 'GGTT', 'False',
             'False', 'protocol_1', 'a designed experiment'],
            ['Project3_10317', '10317', 'False', 'AACC', 'GGTT', 'False',
             'False', 'protocol_1', 'a designed experiment'],
            ['Project4_11223', '11223', 'False', 'AACC', 'GGTT', 'False',
             'False', 'protocol_1', 'a designed experiment']
        ]

        sample_context = [['BLANK.CHILD.1000.G4', 'control blank',
                           '99999', '14577;10317'],
                          ['BLANK.OTHER.10.H4', 'control blank',
                           '14577', ""],
                          ]

        sheet = MetagenomicSampleSheetv101()
        sheet = self._fill_test_metagenomic_sheet(sheet, bfx=bfx, data=data)
        sheet.Header['SheetVersion'] = '101'
        sheet.SampleContext = pd.DataFrame(
            columns=[SAMPLE_NAME_KEY, SAMPLE_TYPE_KEY,
                     PRIMARY_STUDY_KEY, SECONDARY_STUDIES_KEY],
            data=sample_context)

        # Once sheet has been manually populated, validate it.
        self.assertTrue(sheet.validate_and_scrub_sample_sheet())

        # Insert a few errors into the sample-sheet to ensure it fails
        # validation.
        del (sheet.Header['Workflow'])
        sheet.Header['Assay'] = 'NotMetagenomic'

        obs = sheet.quiet_validate_and_scrub_sample_sheet()

        # convert ErrorMessages and WarningMessages into text strings for
        # testing.
        obs = set([str(msg) for msg in obs])

        exp = {"ErrorMessage: 'Workflow' is not declared in Header section",
               "ErrorMessage: 'Assay' value is not 'Metagenomic'"}

        self.assertEqual(obs, exp)

    def test_metagenomic_sheet_w_context_load(self):
        # confirm manual loading is w/out error.
        sheet = MetagenomicSampleSheetv101(
            self.good_metag_ss_w_context, defer_validate=False)
        self.assertTrue(sheet.validate_and_scrub_sample_sheet())

        # confirm load_sample_sheet() returns the correct child class of
        # KLSampleSheet.
        sheet = load_sample_sheet(
            self.good_metag_ss_w_context, defer_validate=False)
        self.assertIsInstance(sheet, MetagenomicSampleSheetv101)


class SampleSheetLoadMakeAndLoadTests(BaseTests):
    sheet_class = KLSampleSheet
    sample_sheet_name = ""

    _INPUT_COLS = [
        SS_SAMPLE_ID_KEY, 'Sample', 'Row', 'Col', 'Blank',
        'Well', 'Project Plate', 'i7 name', 'i7 sequence',
        'i5 name', 'i5 sequence', 'Project Name']

    _INPUT_DATA = [
        ['sample_1', 'sample.1', '1', '1', 'False',
         'A1', 'sample_plate_1', 'iTru7_107_07', 'CCGACTAT',
         'iTru5_01_A', 'ACCGACAA', 'MyProject_99999'],
        ['sample_2', 'sample.2', '2', '1', 'False',
         'A2', 'sample_plate_1', 'iTru7_107_07', 'CCGACTAC',
         'iTru5_01_A', 'ACCGACAT', 'MyProject_99999'],
        ['sample_3', 'sample.3', '3', '1', 'False',
         'A3', 'sample_plate_1', 'iTru7_107_07', 'CCGACTAG',
         'iTru5_01_A', 'ACCGACAG', 'MyProject_99999'],
    ]

    _OUTPUT_COLS = [
        'Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
        'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
        'Sample_Project', 'Well_description']

    _BIOINFORMATICS = [
        {
            'Sample_Project': 'MyProject_99999',
            'QiitaID': '99999',
            'BarcodesAreRC': 'False',
            'ForwardAdapter': 'AACC',
            'ReverseAdapter': 'GGTT',
            'HumanFiltering': 'False',
            'library_construction_protocol': 'some protocol',
            'experiment_design_description': 'some description',
            'contains_replicates': 'False'
        }
    ]

    _CONTACTS = [
        {
            'Sample_Project': 'MyProject_99999',
            'Email': 'foo@bar.org'
        }
    ]

    _SAMPLE_CONTEXT = None

    data_dir = join(dirname(__file__), 'data')

    @staticmethod
    def _make_metadata(a_test_self):
        metadata = {
            _BIOINFORMATICS_KEY:
                [x.copy() for x in a_test_self._BIOINFORMATICS],
            _CONTACT_KEY: [x.copy() for x in a_test_self._CONTACTS],
            _ASSAY_KEY: a_test_self.sheet_class._HEADER[_ASSAY_KEY],
            _SHEET_TYPE_KEY: a_test_self.sheet_class._HEADER[_SHEET_TYPE_KEY],
            _SHEET_VERSION_KEY:
                a_test_self.sheet_class._HEADER[_SHEET_VERSION_KEY]}

        if a_test_self._SAMPLE_CONTEXT is not None:
            metadata[_SAMPLE_CONTEXT_KEY] = \
                [x.copy() for x in a_test_self._SAMPLE_CONTEXT]
        # endif _SAMPLE_CONTEXT

        return metadata

    @property
    def sample_sheet_fp(self):
        return join(self.data_dir, self.sample_sheet_name)

    def _help_test_instantiate_sample_sheet_from_path(
            self, sheet_class, sample_sheet_fp=None, output_cols=None):
        if sample_sheet_fp is None:
            sample_sheet_fp = self.sample_sheet_fp
        if output_cols is None:
            output_cols = self._OUTPUT_COLS

        sheet = sheet_class(sample_sheet_fp, defer_validate=False)

        obs = sheet._get_expected_data_columns()

        self.assertEqual(obs, tuple(output_cols))

        self.assertTrue(sheet.validate_and_scrub_sample_sheet())

    def _help_test_make_sample_sheet(
            self, sheet_class, input_cols=None, input_data=None,
            output_cols=None, sequencer=None):
        if input_cols is None:
            input_cols = self._INPUT_COLS
        if input_data is None:
            input_data = self._INPUT_DATA
        if output_cols is None:
            output_cols = self._OUTPUT_COLS
        if sequencer is None:
            sequencer = 'iSeq'

        table = pd.DataFrame(columns=input_cols, data=input_data)

        metadata = self._make_metadata(self)
        self._help_test_make_sample_sheet_from_metadata(
            sheet_class, metadata, table, output_cols, sequencer)

    def _help_test_make_sample_sheet_from_metadata(
            self, sheet_class, metadata, table, output_cols, sequencer):

        sheet = make_sample_sheet(
            metadata, table, sequencer, [1])

        self.assertIsNotNone(sheet)
        self.assertIsInstance(sheet, sheet_class)
        obs_columns = set(sheet.samples[0].to_json().keys())
        exp_columns = set(output_cols) | {'Lane'}
        self.assertEqual(exp_columns, obs_columns)

        self.assertTrue(sheet.validate_and_scrub_sample_sheet())

    def _help_test_load_sample_sheet(
            self, sheet_class, sample_sheet_fp=None, output_cols=None):
        if sample_sheet_fp is None:
            sample_sheet_fp = self.sample_sheet_fp
        if output_cols is None:
            output_cols = self._OUTPUT_COLS

        sheet1 = load_sample_sheet(sample_sheet_fp, defer_validate=False)
        self.assertEqual(type(sheet1), sheet_class)

        obs = sheet1._get_expected_data_columns()
        self.assertEqual(obs, tuple(output_cols))

        self.assertTrue(sheet1.validate_and_scrub_sample_sheet())

    def _help_test_roundtrip_sample_sheet(
            self, sheet_class, sample_sheet_fp=None):
        if sample_sheet_fp is None:
            sample_sheet_fp = self.sample_sheet_fp
        sheet1 = load_sample_sheet(sample_sheet_fp, defer_validate=False)
        self.assertEqual(type(sheet1), sheet_class)

        self.maxDiff = None
        # write the sample-sheet to a temporary file
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp:
            sheet1.write(tmp)
            tmp.flush()

            # confirm that the written sample-sheet matches the original
            self._help_test_csv_files_exact_text_match(
                sample_sheet_fp, tmp.name)


class MetagenomicSampleSheetv90CreationTests(SampleSheetLoadMakeAndLoadTests):
    sheet_class = MetagenomicSampleSheetv90
    sample_sheet_name = "good_standard_metagv90.csv"

    _OUTPUT_COLS = [
        'Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well',
        'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
        'Sample_Project', 'Well_description']

    @staticmethod
    def _make_metadata(a_test_self):
        metadata = SampleSheetLoadMakeAndLoadTests._make_metadata(a_test_self)
        # remove "contains_replicates" key from metadata[_BIOINFORMATICS_KEY]
        for bfx in metadata[_BIOINFORMATICS_KEY]:
            if CONTAINS_REPLICATES_KEY in bfx:
                del bfx[CONTAINS_REPLICATES_KEY]
        return metadata

    def test_MetagenomicSampleSheetv90_instantiate_from_path(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_MetagenomicSampleSheetv90_make_sample_sheet(self):
        self._help_test_make_sample_sheet(self.sheet_class)

    def test_MetagenomicSampleSheetv90_load_sample_sheet(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_MetagenomicSampleSheetv90_roundtrip(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)


class MetagenomicSampleSheetv100CreationTests(SampleSheetLoadMakeAndLoadTests):
    sheet_class = MetagenomicSampleSheetv100
    sample_sheet_name = "good_standard_metagv100_wo_replicates.csv"
    replicates_sheet_name = "good_standard_metagv100_w_replicates.csv"

    _REP_INPUT_COLS = \
        SampleSheetLoadMakeAndLoadTests._INPUT_COLS.copy() + \
        [ORIG_NAME_KEY, "Library Well"]

    _REP_OUTPUT_COLS = \
        SampleSheetLoadMakeAndLoadTests._OUTPUT_COLS.copy() + \
        [ORIG_NAME_KEY, DESTINATION_WELL_384_KEY]

    def setUp(self):
        self.data_dir = join(dirname(__file__), 'data')
        self.good_w_reps_sheet_fp = join(
            self.data_dir, self.replicates_sheet_name)

    def test_MetagenomicSampleSheetv100_instantiate_from_path_wo_reps(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_MetagenomicSampleSheetv100_make_sample_sheet_wo_reps(self):
        self._help_test_make_sample_sheet(self.sheet_class)

    def test_MetagenomicSampleSheetv100_load_sample_sheet_wo_reps(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_MetagenomicSampleSheetv100_roundtrip_wo_reps(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)

    def test_MetagenomicSampleSheetv100_instantiate_from_path_w_reps(self):
        self._help_test_instantiate_sample_sheet_from_path(
            self.sheet_class, self.good_w_reps_sheet_fp, self._REP_OUTPUT_COLS)

    def test_MetagenomicSampleSheetv100_make_sample_sheet_w_reps(self):
        metadata = self._make_metadata(self)
        bioinfo = [
            {'Sample_Project': 'ProjectF_11661', 'QiitaID': '11661',
             'BarcodesAreRC': False, 'ForwardAdapter': 'AACC',
             'ReverseAdapter': 'GGTT', 'HumanFiltering': False,
             'library_construction_protocol': 'Nextera',
             'experiment_design_description': 'Equipment',
             'contains_replicates': True},
            {'Sample_Project': 'ProjectN_13059', 'QiitaID': '13059',
             'BarcodesAreRC': False, 'ForwardAdapter': 'AACC',
             'ReverseAdapter': 'GGTT', 'HumanFiltering': False,
             'library_construction_protocol': 'Knight Lab Kapa HP',
             'experiment_design_description': 'Equipment',
             'contains_replicates': True}]
        contact = [
            {'Email': 'person@domain.edu', 'Sample_Project': 'ProjectF_11661'},
            {'Email': 'another_person@domain.edu',
             'Sample_Project': 'ProjectN_13059'}]
        metadata[_BIOINFORMATICS_KEY] = bioinfo
        metadata[_CONTACT_KEY] = contact

        data = {
            'orig_name': [
                'BLANK.43.12G', 'BLANK.43.12H', 'RMA.KHP.rpoS.Mage.Q97D',
                'RMA.KHP.rpoS.Mage.Q97L', 'RMA.KHP.rpoS.Mage.Q97N',
                'RMA.KHP.rpoS.Mage.Q97E', 'JBI.KHP.HGL.021', 'JBI.KHP.HGL.022',
                'JBI.KHP.HGL.023', 'JBI.KHP.HGL.024', 'AP581451B02',
                'EP256645B01', 'EP112567B02', 'EP337425B01', 'LP127890A01',
                'EP159692B04', 'EP987683A01', 'AP959450A03', 'SP464350A04',
                'EP121011B01', 'BLANK.43.12G', 'BLANK.43.12H',
                'RMA.KHP.rpoS.Mage.Q97D', 'RMA.KHP.rpoS.Mage.Q97L',
                'RMA.KHP.rpoS.Mage.Q97N', 'RMA.KHP.rpoS.Mage.Q97E',
                'JBI.KHP.HGL.021', 'JBI.KHP.HGL.022', 'JBI.KHP.HGL.023',
                'JBI.KHP.HGL.024', 'AP581451B02', 'EP256645B01', 'EP112567B02',
                'EP337425B01', 'LP127890A01', 'EP159692B04', 'EP987683A01',
                'AP959450A03', 'SP464350A04', 'EP121011B01', 'BLANK.43.12G',
                'BLANK.43.12H', 'RMA.KHP.rpoS.Mage.Q97D',
                'RMA.KHP.rpoS.Mage.Q97L', 'RMA.KHP.rpoS.Mage.Q97N',
                'RMA.KHP.rpoS.Mage.Q97E', 'JBI.KHP.HGL.021', 'JBI.KHP.HGL.022',
                'JBI.KHP.HGL.023', 'JBI.KHP.HGL.024', 'AP581451B02',
                'EP256645B01', 'EP112567B02', 'EP337425B01', 'LP127890A01',
                'EP159692B04', 'EP987683A01', 'AP959450A03', 'SP464350A04',
                'EP121011B01'],
            SS_SAMPLE_ID_KEY: [
                'BLANK_43_12G_A1', 'BLANK_43_12H_A3',
                'RMA_KHP_rpoS_Mage_Q97D_A5', 'RMA_KHP_rpoS_Mage_Q97L_A7',
                'RMA_KHP_rpoS_Mage_Q97N_A9', 'RMA_KHP_rpoS_Mage_Q97E_A11',
                'JBI_KHP_HGL_021_A13', 'JBI_KHP_HGL_022_A15',
                'JBI_KHP_HGL_023_A17', 'JBI_KHP_HGL_024_A19',
                'AP581451B02_A21', 'EP256645B01_A23', 'EP112567B02_C1',
                'EP337425B01_C3', 'LP127890A01_C5', 'EP159692B04_C7',
                'EP987683A01_C9', 'AP959450A03_C11', 'SP464350A04_C13',
                'EP121011B01_C15', 'BLANK_43_12G_A2', 'BLANK_43_12H_A4',
                'RMA_KHP_rpoS_Mage_Q97D_A6', 'RMA_KHP_rpoS_Mage_Q97L_A8',
                'RMA_KHP_rpoS_Mage_Q97N_A10', 'RMA_KHP_rpoS_Mage_Q97E_A12',
                'JBI_KHP_HGL_021_A14', 'JBI_KHP_HGL_022_A16',
                'JBI_KHP_HGL_023_A18', 'JBI_KHP_HGL_024_A20',
                'AP581451B02_A22', 'EP256645B01_A24', 'EP112567B02_C2',
                'EP337425B01_C4', 'LP127890A01_C6', 'EP159692B04_C8',
                'EP987683A01_C10', 'AP959450A03_C12', 'SP464350A04_C14',
                'EP121011B01_C16', 'BLANK_43_12G_B2', 'BLANK_43_12H_B4',
                'RMA_KHP_rpoS_Mage_Q97D_B6', 'RMA_KHP_rpoS_Mage_Q97L_B8',
                'RMA_KHP_rpoS_Mage_Q97N_B10', 'RMA_KHP_rpoS_Mage_Q97E_B12',
                'JBI_KHP_HGL_021_B14', 'JBI_KHP_HGL_022_B16',
                'JBI_KHP_HGL_023_B18', 'JBI_KHP_HGL_024_B20',
                'AP581451B02_B22', 'EP256645B01_B24', 'EP112567B02_D2',
                'EP337425B01_D4', 'LP127890A01_D6', 'EP159692B04_D8',
                'EP987683A01_D10', 'AP959450A03_D12', 'SP464350A04_D14',
                'EP121011B01_D16'],
            'Sample': [
                'BLANK.43.12G.A1', 'BLANK.43.12H.A3',
                'RMA.KHP.rpoS.Mage.Q97D.A5', 'RMA.KHP.rpoS.Mage.Q97L.A7',
                'RMA.KHP.rpoS.Mage.Q97N.A9', 'RMA.KHP.rpoS.Mage.Q97E.A11',
                'JBI.KHP.HGL.021.A13', 'JBI.KHP.HGL.022.A15',
                'JBI.KHP.HGL.023.A17', 'JBI.KHP.HGL.024.A19',
                'AP581451B02.A21', 'EP256645B01.A23', 'EP112567B02.C1',
                'EP337425B01.C3', 'LP127890A01.C5', 'EP159692B04.C7',
                'EP987683A01.C9', 'AP959450A03.C11', 'SP464350A04.C13',
                'EP121011B01.C15', 'BLANK.43.12G.A2', 'BLANK.43.12H.A4',
                'RMA.KHP.rpoS.Mage.Q97D.A6', 'RMA.KHP.rpoS.Mage.Q97L.A8',
                'RMA.KHP.rpoS.Mage.Q97N.A10', 'RMA.KHP.rpoS.Mage.Q97E.A12',
                'JBI.KHP.HGL.021.A14', 'JBI.KHP.HGL.022.A16',
                'JBI.KHP.HGL.023.A18', 'JBI.KHP.HGL.024.A20',
                'AP581451B02.A22', 'EP256645B01.A24', 'EP112567B02.C2',
                'EP337425B01.C4', 'LP127890A01.C6', 'EP159692B04.C8',
                'EP987683A01.C10', 'AP959450A03.C12', 'SP464350A04.C14',
                'EP121011B01.C16', 'BLANK.43.12G.B2', 'BLANK.43.12H.B4',
                'RMA.KHP.rpoS.Mage.Q97D.B6', 'RMA.KHP.rpoS.Mage.Q97L.B8',
                'RMA.KHP.rpoS.Mage.Q97N.B10', 'RMA.KHP.rpoS.Mage.Q97E.B12',
                'JBI.KHP.HGL.021.B14', 'JBI.KHP.HGL.022.B16',
                'JBI.KHP.HGL.023.B18', 'JBI.KHP.HGL.024.B20',
                'AP581451B02.B22', 'EP256645B01.B24', 'EP112567B02.D2',
                'EP337425B01.D4', 'LP127890A01.D6', 'EP159692B04.D8',
                'EP987683A01.D10', 'AP959450A03.D12', 'SP464350A04.D14',
                'EP121011B01.D16'],
            'extra_carried_col': [
                'A1', 'A3', 'A5', 'A7', 'A9', 'A11', 'A13', 'A15', 'A17',
                'A19', 'A21', 'A23', 'C1', 'C3', 'C5', 'C7', 'C9', 'C11',
                'C13', 'C15', 'A1', 'A3', 'A5', 'A7', 'A9', 'A11',
                'A13', 'A15', 'A17', 'A19', 'A21', 'A23', 'C1', 'C3', 'C5',
                'C7', 'C9', 'C11', 'C13', 'C15', 'A1', 'A3', 'A5', 'A7',
                'A9', 'A11', 'A13', 'A15', 'A17', 'A19', 'A21', 'A23', 'C1',
                'C3', 'C5', 'C7', 'C9', 'C11', 'C13', 'C15'],
            'Library Well': [
                'A1', 'A3', 'A5', 'A7', 'A9', 'A11', 'A13', 'A15', 'A17',
                'A19', 'A21', 'A23', 'C1', 'C3', 'C5', 'C7', 'C9', 'C11',
                'C13', 'C15', 'A2', 'A4', 'A6', 'A8', 'A10', 'A12',
                'A14', 'A16', 'A18', 'A20', 'A22', 'A24', 'C2', 'C4', 'C6',
                'C8', 'C10', 'C12', 'C14', 'C16', 'B2', 'B4', 'B6', 'B8',
                'B10', 'B12', 'B14', 'B16', 'B18', 'B20', 'B22', 'B24',
                'D2', 'D4', 'D6', 'D8', 'D10', 'D12', 'D14', 'D16'],
            'Project Plate': [
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectF_11661_P43', 'ProjectF_11661_P43',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectN_13059_P1', 'ProjectN_13059_P1',
                'ProjectN_13059_P1', 'ProjectN_13059_P1'],
            'Well': [
                'A1', 'A3', 'A5', 'A7', 'A9', 'A11', 'A13', 'A15', 'A17',
                'A19', 'A21', 'A23', 'C1', 'C3', 'C5', 'C7', 'C9', 'C11',
                'C13', 'C15', 'A2', 'A4', 'A6', 'A8', 'A10', 'A12', 'A14',
                'A16', 'A18', 'A20', 'A22', 'A24', 'C2', 'C4', 'C6', 'C8',
                'C10', 'C12', 'C14', 'C16', 'B2', 'B4', 'B6', 'B8', 'B10',
                'B12', 'B14', 'B16', 'B18', 'B20', 'B22', 'B24', 'D2', 'D4',
                'D6', 'D8', 'D10', 'D12', 'D14', 'D16'],
            'i7 name': [
                'iTru7_114_08', 'iTru7_114_09', 'iTru7_114_10',
                'iTru7_114_11', 'iTru7_114_12', 'iTru7_201_01',
                'iTru7_201_02', 'iTru7_201_03', 'iTru7_201_04',
                'iTru7_201_05', 'iTru7_108_05', 'iTru7_108_06',
                'iTru7_108_07', 'iTru7_108_08', 'iTru7_108_09',
                'iTru7_108_10', 'iTru7_108_11', 'iTru7_108_12',
                'iTru7_109_01', 'iTru7_109_04', 'iTru7_114_08',
                'iTru7_114_09', 'iTru7_114_10', 'iTru7_114_11',
                'iTru7_114_12', 'iTru7_201_01', 'iTru7_201_02',
                'iTru7_201_03', 'iTru7_201_04', 'iTru7_201_05',
                'iTru7_108_05', 'iTru7_108_06', 'iTru7_108_07',
                'iTru7_108_08', 'iTru7_108_09', 'iTru7_108_10',
                'iTru7_108_11', 'iTru7_108_12', 'iTru7_109_01',
                'iTru7_109_04', 'iTru7_114_08', 'iTru7_114_09',
                'iTru7_114_10', 'iTru7_114_11', 'iTru7_114_12',
                'iTru7_201_01', 'iTru7_201_02', 'iTru7_201_03',
                'iTru7_201_04', 'iTru7_201_05', 'iTru7_108_05',
                'iTru7_108_06', 'iTru7_108_07', 'iTru7_108_08',
                'iTru7_108_09', 'iTru7_108_10', 'iTru7_108_11',
                'iTru7_108_12', 'iTru7_109_01', 'iTru7_109_04'],
            'i7 sequence': [
                'CCGACTAT', 'ACCGACAA', 'CCGACTAT', 'CTTCGCAA', 'GCCTTGTT',
                'AACACCAC', 'AACTTGCC', 'CGTATCTC',
                'CAATGTGG', 'GGTACGAA', 'TCTGAGAG', 'ACCGCATA', 'GAAGTACC',
                'CAGGTATC', 'TCTCTAGG', 'AAGCACTG',
                'CCAAGCAA', 'TGTTCGAG', 'CTCGTCTT', 'TCGGTTAC', 'TCTGAGAG',
                'CAATAGCC', 'ACCGCATA', 'CATTCGTC',
                'GAAGTACC', 'AGTGGCAA', 'CAGGAATC', 'GTGGAATG', 'TCTCAAGG',
                'TGAGATGT', 'TCTGAAAG', 'ACAGCATA',
                'GAAATACC', 'CAGATATC', 'TCTATAGG', 'AAGTTATG', 'ACAAGCAA',
                'AGTTCGAG', 'ATCGTCTT', 'ACGGTTAC',
                'AATTCGGT', 'TCGGACTT', 'TCGGTAAC', 'CATGTGTG', 'AAGTCGAG',
                'TGCCTCAA', 'TATCGGTC', 'ATCTGACC',
                'TATTCGCC', 'CACAGACT', 'GCTGAGAG', 'GCCGCATA', 'GGAGTACC',
                'CGGGTATC', 'TGTCTAGG', 'AGGCACTG',
                'GCAAGCAA', 'GGTTCGAG', 'GTCGTCTT', 'GCGGTTAC'],
            'i5 name': [
                'iTru5_01_A', 'iTru5_02_A', 'iTru5_03_A', 'iTru5_04_A',
                'iTru5_05_A', 'iTru5_06_A', 'iTru5_07_A', 'iTru5_08_A',
                'iTru5_09_A', 'iTru5_10_A', 'iTru5_09_A', 'iTru5_10_A',
                'iTru5_11_A', 'iTru5_12_A', 'iTru5_01_B', 'iTru5_02_B',
                'iTru5_03_B', 'iTru5_04_B', 'iTru5_05_B', 'iTru5_08_B',
                'iTru5_01_A', 'iTru5_02_A', 'iTru5_03_A', 'iTru5_04_A',
                'iTru5_05_A', 'iTru5_06_A', 'iTru5_07_A', 'iTru5_08_A',
                'iTru5_09_A', 'iTru5_10_A', 'iTru5_09_A', 'iTru5_10_A',
                'iTru5_11_A', 'iTru5_12_A', 'iTru5_01_B', 'iTru5_02_B',
                'iTru5_03_B', 'iTru5_04_B', 'iTru5_05_B', 'iTru5_08_B',
                'iTru5_01_A', 'iTru5_02_A', 'iTru5_03_A', 'iTru5_04_A',
                'iTru5_05_A', 'iTru5_06_A', 'iTru5_07_A', 'iTru5_08_A',
                'iTru5_09_A', 'iTru5_10_A', 'iTru5_09_A', 'iTru5_10_A',
                'iTru5_11_A', 'iTru5_12_A', 'iTru5_01_B', 'iTru5_02_B',
                'iTru5_03_B', 'iTru5_04_B', 'iTru5_05_B', 'iTru5_08_B'],
            'i5 sequence': [
                'AAGGCTGA', 'CGATCGAT', 'TTACCGAG', 'AAGACACC', 'GTCCTAAG',
                'CATCTGCT', 'GAAGGTTC', 'CTCTCAGA',
                'GAAGAGGT', 'TCGTCTGA', 'CTCTCAGA', 'TCGTCTGA', 'CAATAGCC',
                'CATTCGTC', 'AGTGGCAA', 'GTGGTATG',
                'TGAGCTGT', 'CGTCAAGA', 'AAGCATCG', 'ACCTCTTC', 'AAGCACTG',
                'CGTCAAGA', 'CCAAGCAA', 'AAGCATCG',
                'TGTTCGAG', 'TACTCCAG', 'CTCGTCTT', 'GATACCTG', 'CGAACTGT',
                'ACCTCTTC', 'ATCTCAGA', 'TGGTCTGA',
                'CATTAGCC', 'CAGTCGTC', 'AGTCGCAA', 'GTGGAATG', 'TGAGCTGT',
                'CGTCCAGA', 'AAGAATCG', 'ACCACTTC',
                'GTATTAGC', 'CACTGAAG', 'AGTCGCTT', 'CACAGGAA', 'TGGCACTA',
                'CCATGAAC', 'GGTTGTCA', 'GCCAATAC',
                'AACCTCCT', 'AGCTACCA', 'CTCTCAGA', 'TCGTCTGA', 'CAATAGCC',
                'CATTCGTC', 'AGTGGCAA', 'GTGGTATG',
                'TGAGCTGT', 'CGTCAAGA', 'AAGCATCG', 'ACCTCTTC'],
            'Project Name': [
                'ProjectF_11661', 'ProjectF_11661', 'ProjectF_11661',
                'ProjectF_11661', 'ProjectF_11661',
                'ProjectF_11661', 'ProjectF_11661', 'ProjectF_11661',
                'ProjectF_11661', 'ProjectF_11661',
                'ProjectN_13059', 'ProjectN_13059', 'ProjectN_13059',
                'ProjectN_13059', 'ProjectN_13059',
                'ProjectN_13059', 'ProjectN_13059', 'ProjectN_13059',
                'ProjectN_13059', 'ProjectN_13059',
                'ProjectF_11661', 'ProjectF_11661', 'ProjectF_11661',
                'ProjectF_11661', 'ProjectF_11661',
                'ProjectF_11661', 'ProjectF_11661', 'ProjectF_11661',
                'ProjectF_11661', 'ProjectF_11661',
                'ProjectN_13059', 'ProjectN_13059', 'ProjectN_13059',
                'ProjectN_13059', 'ProjectN_13059',
                'ProjectN_13059', 'ProjectN_13059', 'ProjectN_13059',
                'ProjectN_13059', 'ProjectN_13059',
                'ProjectF_11661', 'ProjectF_11661', 'ProjectF_11661',
                'ProjectF_11661', 'ProjectF_11661',
                'ProjectF_11661', 'ProjectF_11661', 'ProjectF_11661',
                'ProjectF_11661', 'ProjectF_11661',
                'ProjectN_13059', 'ProjectN_13059', 'ProjectN_13059',
                'ProjectN_13059', 'ProjectN_13059',
                'ProjectN_13059', 'ProjectN_13059', 'ProjectN_13059',
                'ProjectN_13059', 'ProjectN_13059']}
        table = pd.DataFrame(data)

        self._help_test_make_sample_sheet_from_metadata(
            self.sheet_class, metadata, table, self._REP_OUTPUT_COLS, 'iSeq')

    def test_MetagenomicSampleSheetv100_load_sample_sheet_w_reps(self):
        self._help_test_load_sample_sheet(
            self.sheet_class, self.good_w_reps_sheet_fp, self._REP_OUTPUT_COLS)

    def test_MetagenomicSampleSheetv100_roundtrip_w_reps(self):
        self._help_test_roundtrip_sample_sheet(
            self.sheet_class, self.good_w_reps_sheet_fp)


class MetagenomicSampleSheetv101CreationTests(SampleSheetLoadMakeAndLoadTests):
    sheet_class = MetagenomicSampleSheetv101
    sample_sheet_name = "good-sample-sheet_w_sample_context.csv"

    _SAMPLE_CONTEXT = [
        {'sample_name': 'sample.3',
         'sample_type': 'control blank',
         'primary_qiita_study': '99999',
         'secondary_qiita_studies': ''
         }
    ]

    def test_MetagenomicSampleSheetv101_instantiate_from_path(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_MetagenomicSampleSheetv101_make_sample_sheet(self):
        self._help_test_make_sample_sheet(self.sheet_class)

    def test_MetagenomicSampleSheetv101_load_sample_sheet(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_MetagenomicSampleSheetv101_roundtrip(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)


class PacBioMetagSampleSheetv10CreationTests(SampleSheetLoadMakeAndLoadTests):
    sheet_class = PacBioMetagSampleSheetv10
    sample_sheet_name = "good_pacbio_metagv10.csv"

    _INPUT_COLS = [
        SS_SAMPLE_ID_KEY, 'Sample', 'Row', 'Col', 'Blank',
        'Well', 'barcode_id', 'Project Plate', 'Project Name']

    _INPUT_DATA = [
        ['sample_1', 'sample.1', '1', '1', 'False',
         'A1', 'bc3011', 'sample_plate_1', 'MyProject_99999'],
        ['sample_2', 'sample.2', '2', '1', 'False',
         'A2', 'bc0112', 'sample_plate_1', 'MyProject_99999'],
        ['sample_3', 'sample.3', '3', '1', 'False',
         'A3', 'bc9992', 'sample_plate_1', 'MyProject_99999'],
    ]

    _OUTPUT_COLS = [
        'Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well',
        'barcode_id', 'Sample_Project', 'Well_description']

    _BIOINFORMATICS = [
        {
            'Sample_Project': 'MyProject_99999',
            'QiitaID': '99999',
            'HumanFiltering': 'False',
            'library_construction_protocol': 'some protocol',
            'experiment_design_description': 'some description',
            'contains_replicates': 'False'
        }
    ]

    _SAMPLE_CONTEXT = MetagenomicSampleSheetv101CreationTests._SAMPLE_CONTEXT

    def test_PacBioMetagSampleSheetv10_instantiate_from_path(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_PacBioMetagSampleSheetv10_make_sample_sheet(self):
        self._help_test_make_sample_sheet(self.sheet_class, sequencer="Revio")

    def test_PacBioMetagSampleSheetv10_load_sample_sheet(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_PacBioMetagSampleSheetv10_roundtrip(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)


class PacBioAbsquantSampleSheetv10CreationTests(
        SampleSheetLoadMakeAndLoadTests):
    sheet_class = PacBioAbsquantSampleSheetv10
    sample_sheet_name = "good_pacbio_absquantv10.csv"

    _INPUT_COLS = PacBioMetagSampleSheetv10CreationTests._INPUT_COLS.copy() + \
        ['mass_syndna_input_ng', 'extracted_gdna_concentration_ng_ul',
         'vol_extracted_elution_ul', 'syndna_pool_number']

    _INPUT_DATA = [
        ['sample_1', 'sample.1', '1', '1', 'False',
         'A1', 'bc3011', 'sample_plate_1', 'MyProject_99999',
         '0.2', '1.0', '1.1', '1'],
        ['sample_2', 'sample.2', '2', '1', 'False',
         'A2', 'bc0112', 'sample_plate_1', 'MyProject_99999',
         '0.22', '1.0', '1.1', '1'],
        ['sample_3', 'sample.3', '3', '1', 'False',
         'A3', 'bc9992', 'sample_plate_1', 'MyProject_99999',
         '0.25', '1.0', '1.1', '1'],
    ]

    _OUTPUT_COLS = \
        PacBioMetagSampleSheetv10CreationTests._OUTPUT_COLS.copy() + \
        ['mass_syndna_input_ng', 'extracted_gdna_concentration_ng_ul',
         'vol_extracted_elution_ul', 'syndna_pool_number']

    _BIOINFORMATICS = PacBioMetagSampleSheetv10CreationTests._BIOINFORMATICS

    _SAMPLE_CONTEXT = MetagenomicSampleSheetv101CreationTests._SAMPLE_CONTEXT

    def test_PacBioAbsquantSampleSheetv10_instantiate_from_path(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_PacBioAbsquantSampleSheetv10_make_sample_sheet(self):
        self._help_test_make_sample_sheet(self.sheet_class, sequencer="Revio")

    def test_PacBioAbsquantSampleSheetv10_load_sample_sheet(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_PacBioAbsquantSampleSheetv10_roundtrip(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)


class PacBioMetagSampleSheetv11CreationTests(SampleSheetLoadMakeAndLoadTests):
    sheet_class = PacBioMetagSampleSheetv11
    sample_sheet_name = "good_pacbio_metagv11.csv"

    _INPUT_COLS = [
        SS_SAMPLE_ID_KEY, 'Sample', 'Row', 'Col', 'Blank',
        'Well', 'barcode_id', 'twist_adaptor_id',
        'Project Plate', 'Project Name']

    _INPUT_DATA = [
        ['sample_1', 'sample.1', '1', '1', 'False',
         'A1', 'bc3011', '16_UDI_1_A01_F--16_UDI_1_A01_R',
         'sample_plate_1', 'MyProject_99999'],
        ['sample_2', 'sample.2', '2', '1', 'False',
         'A2', 'bc0112', '16_UDI_2_B01_F--16_UDI_2_B01_R',
         'sample_plate_1', 'MyProject_99999'],
        ['sample_3', 'sample.3', '3', '1', 'False',
         'A3', 'bc9992', '16_UDI_5_E01_F--16_UDI_5_E01_R',
         'sample_plate_1', 'MyProject_99999'],
    ]

    _OUTPUT_COLS = [
        'Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well',
        'barcode_id', 'twist_adaptor_id',
        'Sample_Project', 'Well_description']

    _BIOINFORMATICS = [
        {
            'Sample_Project': 'MyProject_99999',
            'QiitaID': '99999',
            'HumanFiltering': 'False',
            'library_construction_protocol': 'some protocol',
            'experiment_design_description': 'some description',
            'contains_replicates': 'False'
        }
    ]

    _SAMPLE_CONTEXT = MetagenomicSampleSheetv101CreationTests._SAMPLE_CONTEXT

    def test_PacBioMetagSampleSheetv11_instantiate_from_path(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_PacBioMetagSampleSheetv11_make_sample_sheet(self):
        self._help_test_make_sample_sheet(self.sheet_class, sequencer="Revio")

    def test_PacBioMetagSampleSheetv11_load_sample_sheet(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_PacBioMetagSampleSheetv11_roundtrip(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)


class PacBioAbsquantSampleSheetv11CreationTests(
        SampleSheetLoadMakeAndLoadTests):
    sheet_class = PacBioAbsquantSampleSheetv11
    sample_sheet_name = "good_pacbio_absquantv11.csv"

    _INPUT_COLS = PacBioMetagSampleSheetv11CreationTests._INPUT_COLS.copy() + \
        ['mass_syndna_input_ng', 'extracted_gdna_concentration_ng_ul',
         'vol_extracted_elution_ul', 'syndna_pool_number', 'syndna_is_twisted']

    _INPUT_DATA = [
        ['sample_1', 'sample.1', '1', '1', 'False',
         'A1', 'bc3011', '16_UDI_1_A01_F--16_UDI_1_A01_R',
         'sample_plate_1', 'MyProject_99999',
         '0.2', '1.0', '1.1', '1', 'False'],
        ['sample_2', 'sample.2', '2', '1', 'False',
         'A2', 'bc0112', '16_UDI_2_B01_F--16_UDI_2_B01_R',
         'sample_plate_1', 'MyProject_99999',
         '0.22', '1.0', '1.1', '1', 'False'],
        ['sample_3', 'sample.3', '3', '1', 'False',
         'A3', 'bc9992', '16_UDI_5_E01_F--16_UDI_5_E01_R',
         'sample_plate_1', 'MyProject_99999',
         '0.25', '1.0', '1.1', '1', 'False'],
    ]

    _OUTPUT_COLS = \
        PacBioMetagSampleSheetv11CreationTests._OUTPUT_COLS.copy() + [
            'mass_syndna_input_ng', 'extracted_gdna_concentration_ng_ul',
            'vol_extracted_elution_ul', 'syndna_pool_number',
            'syndna_is_twisted']

    _BIOINFORMATICS = PacBioMetagSampleSheetv11CreationTests._BIOINFORMATICS

    _SAMPLE_CONTEXT = MetagenomicSampleSheetv101CreationTests._SAMPLE_CONTEXT

    def test_PacBioAbsquantSampleSheetv11_instantiate_from_path(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_PacBioAbsquantSampleSheetv11_make_sample_sheet(self):
        self._help_test_make_sample_sheet(self.sheet_class, sequencer="Revio")

    def test_PacBioAbsquantSampleSheetv11_load_sample_sheet(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_PacBioAbsquantSampleSheetv11_roundtrip(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)


class PacBioAbsquantSampleSheetv12CreationTests(
        SampleSheetLoadMakeAndLoadTests):
    sheet_class = PacBioAbsquantSampleSheetv12
    sample_sheet_name = "good_pacbio_absquantv12.csv"

    _INPUT_COLS = PacBioMetagSampleSheetv11CreationTests._INPUT_COLS.copy() + \
        ['mass_syndna_input_ng', 'extracted_gdna_concentration_ng_ul',
         'vol_extracted_elution_ul', 'syndna_pool_number',
         'sequenced_sample_gdna_mass_ng', 'calc_mass_sample_aliquot_input_g',
         'sample_volume_ul', 'sample_surface_area_cm2', 'syndna_is_twisted']

    _INPUT_DATA = [
        ['sample_1', 'sample.1', '1', '1', 'False',
         'A1', 'bc3011', '16_UDI_1_A01_F--16_UDI_1_A01_R',
         'sample_plate_1', 'MyProject_99999',
         '0.2', '1.0', '1.1', '1', '0.15', '0.02', '2.5', '', 'False'],
        ['sample_2', 'sample.2', '2', '1', 'False',
         'A2', 'bc0112', '16_UDI_2_B01_F--16_UDI_2_B01_R',
         'sample_plate_1', 'MyProject_99999',
         '0.22', '1.0', '1.1', '1', '0.18', '0.012', '', '1.2', 'False'],
        ['sample_3', 'sample.3', '3', '1', 'False',
         'A3', 'bc9992', '16_UDI_5_E01_F--16_UDI_5_E01_R',
         'sample_plate_1', 'MyProject_99999',
         '0.25', '1.0', '1.1', '1', '0.21', '0.105', '2.5', '1.2', 'False'],
    ]

    _OUTPUT_COLS = \
        PacBioMetagSampleSheetv11CreationTests._OUTPUT_COLS.copy() + [
            'mass_syndna_input_ng', 'extracted_gdna_concentration_ng_ul',
            'vol_extracted_elution_ul', 'syndna_pool_number',
            'sequenced_sample_gdna_mass_ng',
            'calc_mass_sample_aliquot_input_g',
            'sample_volume_ul', 'sample_surface_area_cm2', 'syndna_is_twisted']

    _BIOINFORMATICS = PacBioMetagSampleSheetv11CreationTests._BIOINFORMATICS

    _SAMPLE_CONTEXT = MetagenomicSampleSheetv101CreationTests._SAMPLE_CONTEXT

    def test_PacBioAbsquantSampleSheetv12_instantiate_from_path(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_PacBioAbsquantSampleSheetv12_make_sample_sheet(self):
        self._help_test_make_sample_sheet(self.sheet_class, sequencer="Revio")

    def test_PacBioAbsquantSampleSheetv12_load_sample_sheet(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_PacBioAbsquantSampleSheetv12_roundtrip(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)


class AbsQuantSampleSheetv10CreationTests(SampleSheetLoadMakeAndLoadTests):
    sheet_class = AbsQuantSampleSheetv10
    sample_sheet_name = "good_abs_quant_metagv10.csv"

    _INPUT_COLS = [
        SS_SAMPLE_ID_KEY, 'Sample', 'Row', 'Col', 'Blank',
        'Well', 'Project Plate', 'i7 name', 'i7 sequence',
        'i5 name', 'i5 sequence', 'Project Name',
        'mass_syndna_input_ng', 'extracted_gdna_concentration_ng_ul',
        'vol_extracted_elution_ul', 'syndna_pool_number'
    ]

    _INPUT_DATA = [
        ['sample_1', 'sample.1', '1', '1', 'sample_plate_1', 'False',
         'A1', 'iTru7_107_07', 'CCGACTAT',
         'iTru5_01_A', 'ACCGACAA', 'MyProject_99999',
         '0.2', '1.0', '1.1', '1'],
        ['sample_2', 'sample.2', '2', '1', 'sample_plate_1', 'False',
         'A2', 'iTru7_107_07', 'CCGACTAC',
         'iTru5_01_A', 'ACCGACAT', 'MyProject_99999',
         '0.22', '1.0', '1.1', '1'],
        ['sample_3', 'sample.3', '3', '1', 'sample_plate_1', 'False',
         'A3', 'iTru7_107_07', 'CCGACTAG',
         'iTru5_01_A', 'ACCGACAG', 'MyProject_99999',
         '0.25', '1.0', '1.1', '1'],
    ]

    _OUTPUT_COLS = [
        'Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
        'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
        'Sample_Project', 'Well_description',
        'mass_syndna_input_ng', 'extracted_gdna_concentration_ng_ul',
        'vol_extracted_elution_ul', 'syndna_pool_number'
    ]

    def test_AbsQuantSampleSheetv10_instantiate_from_path(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_AbsQuantSampleSheetv10_make_sample_sheet(self):
        self._help_test_make_sample_sheet(self.sheet_class)

    def test_AbsQuantSampleSheetv10_load_sample_sheet(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_AbsQuantSampleSheetv10_roundtrip(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)


class AbsQuantSampleSheetv11CreationTests(SampleSheetLoadMakeAndLoadTests):
    sheet_class = AbsQuantSampleSheetv11
    sample_sheet_name = "good_abs_quant_metagv11.csv"

    _INPUT_COLS = AbsQuantSampleSheetv10CreationTests._INPUT_COLS

    _INPUT_DATA = AbsQuantSampleSheetv10CreationTests._INPUT_DATA

    _OUTPUT_COLS = AbsQuantSampleSheetv10CreationTests._OUTPUT_COLS

    _SAMPLE_CONTEXT = MetagenomicSampleSheetv101CreationTests._SAMPLE_CONTEXT

    def test_AbsQuantSampleSheetv11_instantiate_from_path(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_AbsQuantSampleSheetv11_make_sample_sheet(self):
        self._help_test_make_sample_sheet(self.sheet_class)

    def test_AbsQuantSampleSheetv11_load_sample_sheet(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_AbsQuantSampleSheetv11_roundtrip(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)


class TellseqMetagSampleSheetv10CreationTests(SampleSheetLoadMakeAndLoadTests):
    sheet_class = TellseqMetagSampleSheetv10
    sample_sheet_name = "tellseq_metag_dummy_sample_sheet_2.csv"

    _INPUT_COLS = [
        SS_SAMPLE_ID_KEY, 'Sample', 'Row', 'Col', 'Blank',
        'Well', 'Project Plate', 'barcode_id', 'Project Name'
    ]

    _INPUT_DATA = [
        ['sample_1', 'sample.1', '1', '1', 'sample_plate_1', 'False',
         'A1', 'C501', 'MyProject_99999'],
        ['sample_2', 'sample.2', '2', '1', 'sample_plate_1', 'False',
         'A2', 'C509', 'MyProject_99999'],
        ['sample_3', 'sample.3', '3', '1', 'sample_plate_1', 'False',
         'A3', 'C520', 'MyProject_99999'],
    ]

    _OUTPUT_COLS = [
        'Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
        'barcode_id',
        'Sample_Project', 'Well_description'
    ]

    _SAMPLE_CONTEXT = MetagenomicSampleSheetv101CreationTests._SAMPLE_CONTEXT

    def test_TellseqMetagSampleSheetv10_instantiate_from_path(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_TellseqMetagSampleSheetv10_make_sample_sheet(self):
        self._help_test_make_sample_sheet(self.sheet_class)

    def test_TellseqMetagSampleSheetv10_load_sample_sheet(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_TellseqMetagSampleSheetv10_roundtrip(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)


class TellseqAbsquantMetagSampleSheetv10CreationTests(
        SampleSheetLoadMakeAndLoadTests):
    sheet_class = TellseqAbsquantMetagSampleSheetv10
    sample_sheet_name = "tellseq_absquant_dummy_sample_sheet_2.csv"

    _INPUT_COLS = [
        SS_SAMPLE_ID_KEY, 'Sample', 'Row', 'Col', 'Blank',
        'Well', 'Project Plate', 'barcode_id', 'Project Name',
        'mass_syndna_input_ng', 'extracted_gdna_concentration_ng_ul',
        'vol_extracted_elution_ul', 'syndna_pool_number'
    ]

    _INPUT_DATA = [
        ['sample_1', 'sample.1', '1', '1', 'sample_plate_1', 'False',
         'A1', 'C501', 'MyProject_99999',
         '0.2', '1.0', '1.1', '1'],
        ['sample_2', 'sample.2', '2', '1', 'sample_plate_1', 'False',
         'A2', 'C509', 'MyProject_99999',
         '0.22', '1.0', '1.1', '1'],
        ['sample_3', 'sample.3', '3', '1', 'sample_plate_1', 'False',
         'A3', 'C520', 'MyProject_99999',
         '0.25', '1.0', '1.1', '1'],
    ]

    _OUTPUT_COLS = [
        'Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
        'barcode_id',
        'Sample_Project', 'Well_description',
        'mass_syndna_input_ng', 'extracted_gdna_concentration_ng_ul',
        'vol_extracted_elution_ul', 'syndna_pool_number'
    ]

    _SAMPLE_CONTEXT = MetagenomicSampleSheetv101CreationTests._SAMPLE_CONTEXT

    def test_TellseqAbsquantMetagSampleSheetv10_instantiate_from_path(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_TellseqAbsquantMetagSampleSheetv10_make_sample_sheet(self):
        self._help_test_make_sample_sheet(self.sheet_class)

    def test_TellseqAbsquantMetagSampleSheetv10_load_sample_sheet(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_TellseqAbsquantMetagSampleSheetv10_roundtrip(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)


class MetatranscriptomicSampleSheetv10CreationTests(
        SampleSheetLoadMakeAndLoadTests):
    sheet_class = MetatranscriptomicSampleSheetv10
    sample_sheet_name = "good_standard_metatv10.csv"

    _INPUT_COLS = [
        SS_SAMPLE_ID_KEY, 'Sample', 'Row', 'Col', 'Blank',
        'Well', 'Project Plate', 'i7 name', 'i7 sequence',
        'i5 name', 'i5 sequence', 'Project Name',
        'total_rna_concentration_ng_ul', 'vol_extracted_elution_ul']

    _INPUT_DATA = [
        ['sample_1', 'sample.1', '1', '1', 'False',
         'A1', 'sample_plate_1', 'iTru7_107_07', 'CCGACTAT',
         'iTru5_01_A', 'ACCGACAA', 'MyProject_99999',
         '94.454', '70'],
        ['sample_2', 'sample.2', '2', '1', 'False',
         'A2', 'sample_plate_1', 'iTru7_107_07', 'CCGACTAC',
         'iTru5_01_A', 'ACCGACAT', 'MyProject_99999',
         '5.3', '70'],
        ['sample_3', 'sample.3', '3', '1', 'False',
         'A3', 'sample_plate_1', 'iTru7_107_07', 'CCGACTAG',
         'iTru5_01_A', 'ACCGACAG', 'MyProject_99999',
         '0.5', '70'],
    ]

    _OUTPUT_COLS = [
        'Sample_ID', 'Sample_Name', 'Sample_Plate', 'well_id_384',
        'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
        'Sample_Project', 'total_rna_concentration_ng_ul',
        'vol_extracted_elution_ul', 'Well_description']

    @staticmethod
    def _make_metadata(a_self):
        return MetagenomicSampleSheetv90CreationTests._make_metadata(a_self)

    def test_MetatranscriptomicSampleSheetv10_instantiate_from_path(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_MetatranscriptomicSampleSheetv10_make_sample_sheet(self):
        self._help_test_make_sample_sheet(self.sheet_class)

    def test_MetatranscriptomicSampleSheetv10_load_sample_sheet(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_MetatranscriptomicSampleSheetv10_roundtrip(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)


class MetagenomicSampleSheetv102CreationTests(SampleSheetLoadMakeAndLoadTests):
    sheet_class = MetagenomicSampleSheetv102
    sample_sheet_name = "good_standard_metagv102_wo_katharoseq.csv"
    kath_sheet_name = "good_standard_metagv102_w_katharoseq.csv"

    _KATH_COLS = [
        'Kathseq_RackID', TUBECODE_KEY, 'katharo_description',
        'number_of_cells', 'platemap_generation_date', 'project_abbreviation',
        'vol_extracted_elution_ul', 'well_id_96']

    _KATH_INPUT_COLS = \
        SampleSheetLoadMakeAndLoadTests._INPUT_COLS.copy() + _KATH_COLS.copy()

    _KATH_INPUT_DATA = [
        ['sample_1', 'sample.1', '1', '1', 'False',
         'A1', 'sample_plate_1', 'iTru7_107_07', 'CCGACTAT',
         'iTru5_01_A', 'ACCGACAA', 'MyProject_99999',
         '', '', '', '', '', '', '', ''],
        ['sample_2', 'sample.2', '2', '1', 'False',
         'A2', 'sample_plate_1', 'iTru7_107_07', 'CCGACTAC',
         'iTru5_01_A', 'ACCGACAT', 'MyProject_99999',
         '', '', '', '', '', '', '', ''],
        ['sample_3', 'sample.3', '3', '1', 'False',
         'A3', 'sample_plate_1', 'iTru7_107_07', 'CCGACTAG',
         'iTru5_01_A', 'ACCGACAG', 'MyProject_99999',
         '', '', '', '', '', '', '', ''],
        # added katharoseq control here. apparently katharoseq-specific
        # columns have to exist but don't have to be filled even for the
        # katharoseq control?  Seems like an oversight.
        ['katharo001', 'katharo0001', '4', '1', 'False',
         'A3', 'sample_plate_1', 'iTru7_107_07', 'CCGACTCT',
         'iTru5_01_A', 'ACCGACCG', 'MyProject_99999',
         '', '', '', '', '', '', '', '']
    ]

    _KATH_OUTPUT_COLS = \
        SampleSheetLoadMakeAndLoadTests._OUTPUT_COLS.copy() + _KATH_COLS.copy()

    _SAMPLE_CONTEXT = MetagenomicSampleSheetv101CreationTests._SAMPLE_CONTEXT

    _MISSING_COLS_ERR_LINES = [
        "Sample sheet instantiation failed: ",
        "The Kathseq_RackID column in the Data section is missing\n",
        "The TubeCode column in the Data section is missing\n",
        "The katharo_description column in the Data section is missing\n",
        "The number_of_cells column in the Data section is missing\n",
        "The platemap_generation_date column in the Data section is missing\n",
        "The project_abbreviation column in the Data section is missing\n",
        "The vol_extracted_elution_ul column in the Data section is missing\n",
        "The well_id_96 column in the Data section is missing"]

    def setUp(self):
        self.data_dir = join(dirname(__file__), 'data')
        self.good_wo_katharoseq_sheet_fp = join(
            self.data_dir, self.sample_sheet_name)

        self.good_w_katharoseq_sheet_fp = join(
            self.data_dir, self.kath_sheet_name)

        self.bad_missing_katharoseq_col_sheet_fp = join(
            self.data_dir, 'test_katharoseq_sheet3.csv')

    def test_MetagenomicSampleSheetv102_instantiate_from_path_wo_kath(self):
        self._help_test_instantiate_sample_sheet_from_path(self.sheet_class)

    def test_MetagenomicSampleSheetv102_make_sample_sheet_wo_kath(self):
        self._help_test_make_sample_sheet(self.sheet_class)

    def test_MetagenomicSampleSheetv102_load_sample_sheet_wo_kath(self):
        self._help_test_load_sample_sheet(self.sheet_class)

    def test_MetagenomicSampleSheetv102_roundtrip_wo_kath(self):
        self._help_test_roundtrip_sample_sheet(self.sheet_class)

    def test_MetagenomicSampleSheetv102_instantiate_from_path_w_kath(self):
        self._help_test_instantiate_sample_sheet_from_path(
            self.sheet_class, self.good_w_katharoseq_sheet_fp,
            self._KATH_OUTPUT_COLS)

    def test_MetagenomicSampleSheetv102_make_sample_sheet_w_kath(self):
        self._help_test_make_sample_sheet(
            self.sheet_class, self._KATH_INPUT_COLS, self._KATH_INPUT_DATA,
            self._KATH_OUTPUT_COLS)

    def test_MetagenomicSampleSheetv102_load_sample_sheet_w_kath(self):
        self._help_test_load_sample_sheet(
            self.sheet_class, self.good_w_katharoseq_sheet_fp,
            self._KATH_OUTPUT_COLS)

    def test_MetagenomicSampleSheetv102_roundtrip_w_kath(self):
        self._help_test_roundtrip_sample_sheet(
            self.sheet_class, self.good_w_katharoseq_sheet_fp)

    def test_katharoseq_enabled_sheet_load_no_state_change(self):
        # confirm that class-wide state is not permanently changed by loading
        # a karathoseq-enabled file. Loading sheet1 (no katharoseq columns)
        # after sheet2 (containing katharoseq columns) should yield a sheet1
        # samplesheet with only the shorter set of columns.
        sheet2 = load_sample_sheet(self.good_w_katharoseq_sheet_fp)
        self.assertEqual(type(sheet2), MetagenomicSampleSheetv102)

        sheet1 = load_sample_sheet(self.good_wo_katharoseq_sheet_fp,
                                   defer_validate=False)
        self.assertEqual(type(sheet1), MetagenomicSampleSheetv102)
        obs = sheet1._get_expected_data_columns()
        self.assertEqual(obs, tuple(self._OUTPUT_COLS))
        self.assertTrue(sheet1.validate_and_scrub_sample_sheet())

    def test_katharoseq_enabled_sheet_err_load_missing_col(self):
        err = ("Sample sheet instantiation failed: The number_of_cells column"
               " in the Data section is missing")
        with self.assertRaisesRegex(ValueError, err):
            _ = MetagenomicSampleSheetv102(
                self.bad_missing_katharoseq_col_sheet_fp, defer_validate=False)

    def test_katharoseq_enabled_sheet_creation(self):
        data = self._INPUT_DATA.copy()
        data[-1][0] = data[-1][1] = "katharo0001"
        table = pd.DataFrame(columns=self._INPUT_COLS, data=data)
        err = "".join(self._MISSING_COLS_ERR_LINES)
        with self.assertRaisesRegex(ValueError, err):
            make_sample_sheet(self._make_metadata(self),
                              table, 'iSeq', [1])

    def test_katharoseq_make_sample_sheet_implicit_not_strict(self):
        data = self._KATH_INPUT_DATA.copy()
        data = data[:-1]  # drop the only kath sample
        table = pd.DataFrame(columns=self._KATH_INPUT_COLS, data=data)
        sheet = make_sample_sheet(self._make_metadata(self),
                                  table, 'iSeq', [1])

        # confirm that we get a sample-sheet w/out katharoseq-control-related
        # columns when we allow implicit table remapping and have no kath
        # samples in the data
        self.assertIsNotNone(sheet)
        self.assertIsInstance(sheet, MetagenomicSampleSheetv102)
        self.assertFalse(sheet.contains_katharoseq_samples())
        obs_columns = set(sheet.samples[0].to_json().keys())
        exp_columns = set(self._OUTPUT_COLS) | {'Lane'}
        self.assertEqual(exp_columns, obs_columns)

    def test_katharoseq_make_sample_sheet_kath_cols_wo_kath_sample_ok(self):
        # input data with katharoseq columns but no katharoseq sample
        # will have kath columns dropped when the sample sheet is created.

        # first, delete the last item in the self._KATH_INPUT_DATA list
        # so that the katharoseq control is not present in the data.
        data = self._KATH_INPUT_DATA.copy()
        data = data[:-1]

        table = pd.DataFrame(columns=self._KATH_INPUT_COLS, data=data)
        metadata = self._make_metadata(self)

        # sheet will be created but extended columns will not be present
        # and no error is raised. Kathseq_RackID is silently dropped.
        sheet = make_sample_sheet(metadata, table, 'iSeq', [1])

        self.assertIsNotNone(sheet)
        self.assertIsInstance(sheet, MetagenomicSampleSheetv102)
        self.assertFalse(sheet.contains_katharoseq_samples())
        obs_columns = set(sheet.samples[0].to_json().keys())
        exp_columns = {'Sample_ID', 'Sample_Name', 'Sample_Plate',
                       'well_id_384', 'I7_Index_ID', 'index', 'I5_Index_ID',
                       'index2', 'Sample_Project', 'Well_description',
                       'Lane'}
        self.assertEqual(obs_columns, exp_columns)

    def test_katharoseq_make_sample_sheet_kath_sample_wo_kath_cols_err(self):
        # input data with katharoseq sample but not all katharoseq columns
        # will throw an error during sample sheet creation.

        table = pd.DataFrame(
            columns=self._KATH_INPUT_COLS, data=self._KATH_INPUT_DATA)

        cols_to_remove = set(self._KATH_COLS) - {'Kathseq_RackID'}
        table.drop(columns=list(cols_to_remove), inplace=True)

        err_lines = self._MISSING_COLS_ERR_LINES.copy()
        err_lines.remove(
            "The Kathseq_RackID column in the Data section is missing\n")
        err = "".join(err_lines)
        with self.assertRaisesRegex(ValueError, err):
            make_sample_sheet(self._make_metadata(self),
                              table, 'iSeq', [1])


class MissingValueTests(BaseTests):
    def setUp(self):
        super().setUp()
        self.sheet_fp = join(self.data_dir,
                             'good_standard_metagv102_wo_katharoseq.csv')

    @staticmethod
    def _write(sheet):
        out = StringIO()
        sheet.write(out)
        return out.getvalue().splitlines()

    def test_write_sample_missing_a_data_column(self):
        # a sample with a column the existing samples lack leaves a gap in
        # their rows, which must be written as an empty cell, not 'nan'
        sheet = load_sample_sheet(self.sheet_fp)
        first = sheet.samples[0]
        new_sample = {k: getattr(first, k) for k in sheet.all_sample_keys}
        new_sample.update({
            SS_SAMPLE_ID_KEY: 'NEW_SAMPLE', 'Sample_Name': 'NEW.SAMPLE',
            'index': 'ACGTACGTAC', 'index2': 'TTGGCCAATT',
            'katharo_description': 'a description'})
        sheet.add_sample(sample_sheet.Sample(new_sample))

        obs = [x for x in self._write(sheet)
               if x.startswith(f'{first.Sample_ID},')]
        self.assertEqual(len(obs), 1)
        self.assertNotIn('nan', obs[0].split(','))
        self.assertTrue(obs[0].endswith(',1,'))

    def test_write_none_in_a_section(self):
        sheet = load_sample_sheet(self.sheet_fp)
        contact = sheet.Contact.copy()
        contact.loc[contact.index[0], 'Email'] = None
        sheet.Contact = contact

        lines = self._write(sheet)
        contact_start = next(i for i, x in enumerate(lines)
                             if x.startswith('[Contact]'))
        obs = lines[contact_start + 2]
        self.assertTrue(obs.startswith(',MyProject_99999,'), obs)

    def test_load_sample_sheet_blank_sample_id(self):
        with open(self.sheet_fp) as f:
            lines = f.read().splitlines()
        data_header = next(i for i, x in enumerate(lines)
                           if x.startswith('[Data]')) + 1
        first_row = lines[data_header + 1].split(',')
        first_row[0] = ''
        lines[data_header + 1] = ','.join(first_row)

        with tempfile.NamedTemporaryFile('w', suffix='.csv') as tmp:
            tmp.write('\n'.join(lines) + '\n')
            tmp.flush()
            with self.assertRaisesRegex(
                    ValueError,
                    'The Sample_ID column in the Data section is empty in '
                    r'sample row\(s\): 1'):
                load_sample_sheet(tmp.name)


if __name__ == '__main__':
    assert not hasattr(sys.stdout, "getvalue")
    unittest.main(module=__name__, buffer=True, exit=False)
