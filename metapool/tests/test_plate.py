import sys
import os
import pandas as pd
from pandas.testing import assert_frame_equal

from unittest import TestCase, main
from metapool.mp_strings import (get_short_name_and_id, parse_project_name,
                                 QIITA_ID_KEY, PROJECT_SHORT_NAME_KEY,
                                 PROJECT_FULL_NAME_KEY)
from metapool.plate import (_well_to_row_and_col, _decompress_well,
                            _plate_position, validate_plate_metadata,
                            _validate_plate, Message, ErrorMessage,
                            WarningMessage, requires_dilution, dilute_gDNA,
                            find_threshold, autopool, _validate_well_id_96,
                            PlateReplication)
from metapool.metapool import (read_plate_map_csv, read_pico_csv,
                               calculate_norm_vol, assign_index,
                               compute_pico_concentration)


class DilutionTests(TestCase):
    def test_dilution_test(self):
        plate_df = read_plate_map_csv('notebooks/test_data/Plate_Maps/Finrisk'
                                      ' 33-36_plate_map.tsv')
        sample_concs = read_pico_csv(open('notebooks/test_data/Quant/MiniPico'
                                          '/FinRisk_33-36_gDNA_quant.tsv',
                                          'r'))
        plate_df = pd.merge(plate_df, sample_concs, on='Well')

        self.assertFalse(requires_dilution(plate_df,
                                           threshold=20,
                                           tolerance=.10))
        self.assertTrue(requires_dilution(plate_df,
                                          threshold=10,
                                          tolerance=.10))

        # these values shouldn't change
        self.assertEqual(plate_df['Sample DNA Concentration'][0], 0.044)
        self.assertEqual(plate_df['Sample DNA Concentration'][2], 3.336)
        self.assertEqual(plate_df['Sample DNA Concentration'][3], 1.529)
        self.assertEqual(plate_df['Sample DNA Concentration'][4], 6.705)

        # these should change
        self.assertEqual(plate_df['Sample DNA Concentration'][1], 25.654)
        self.assertEqual(plate_df['Sample DNA Concentration'][8], 12.605)
        self.assertEqual(plate_df['Sample DNA Concentration'][9], 12.378)
        self.assertEqual(plate_df['Sample DNA Concentration'][10], 10.595)

        plate_df = dilute_gDNA(plate_df, threshold=10)

        # these values shouldn't change
        self.assertEqual(plate_df['Sample DNA Concentration'][0], 0.044)
        self.assertEqual(plate_df['Sample DNA Concentration'][2], 3.336)
        self.assertEqual(plate_df['Sample DNA Concentration'][3], 1.529)
        self.assertEqual(plate_df['Sample DNA Concentration'][4], 6.705)

        # these should change
        self.assertEqual(plate_df['Sample DNA Concentration'][1], 2.5654)
        self.assertEqual(plate_df['Sample DNA Concentration'][8], 1.2605)
        self.assertEqual(plate_df['Sample DNA Concentration'][9], 1.2378)
        self.assertEqual(plate_df['Sample DNA Concentration'][10], 1.0595)

    def test_find_threshold_autopool(self):
        plate_map_fp = ('notebooks/test_data/Plate_Maps/'
                        'Finrisk 33-36_plate_map.tsv')
        plate_df = read_plate_map_csv(plate_map_fp)
        sample_concs_fp = ('notebooks/test_data/Quant/MiniPico/'
                           'FinRisk_33-36_gDNA_quant.tsv')
        sample_concs = read_pico_csv(open(sample_concs_fp, 'r'),
                                     plate_reader='Synergy_HT',
                                     max_conc=60)
        plate_df = pd.merge(plate_df, sample_concs, on='Well')
        total_vol = 3500
        dna_vols = calculate_norm_vol(plate_df['Sample DNA Concentration'],
                                      ng=5, min_vol=25, max_vol=3500,
                                      resolution=2.5)
        plate_df['Normalized DNA volume'] = dna_vols
        plate_df['Normalized water volume'] = total_vol - dna_vols
        plate_df['Library Well'] = plate_df['Well']
        index_combos = pd.read_csv('notebooks/test_output/iTru/'
                                   'temp_iTru_combos.csv')
        indices = assign_index(len(plate_df['Sample']),
                               index_combos,
                               start_idx=0).reset_index()
        plate_df = pd.concat([plate_df, indices], axis=1)
        lib_concs_fp = ('notebooks/test_data/Quant/MiniPico/'
                        '10-13-17_FinRisk_33-36_library_quant.txt')
        col_name = 'MiniPico Library DNA Concentration'
        lib_concs = read_pico_csv(open(lib_concs_fp, 'r'),
                                  plate_reader='Synergy_HT',
                                  conc_col_name=col_name)
        plate_df = pd.merge(plate_df, lib_concs, on='Well')
        result = compute_pico_concentration(plate_df[col_name], size=500)
        plate_df['MiniPico Library Concentration'] = result
        col1 = 'Sample DNA Concentration'
        col2 = 'Normalized DNA volume'
        plate_df['Input DNA'] = plate_df[col1] * plate_df[col2] / 1000
        threshold = find_threshold(plate_df['MiniPico Library Concentration'],
                                   plate_df['Blank'])
        self.assertAlmostEqual(5.8818, threshold, places=4)
        plate_df = autopool(plate_df, method='norm', pool_failures='high',
                            automate=True, floor_conc=threshold, offset=0.01)
        self.assertEqual(plate_df['MiniPico Pooled Volume'].tolist()[0:6],
                         [1000.0, 314.5959409755514, 217.8574480616536,
                          295.68518270050976, 236.59171073317077, 1000.0])


class PlateHelperTests(TestCase):
    def test_well_to_row_and_col(self):
        self.assertEqual((1, 1), _well_to_row_and_col('A1'))
        self.assertEqual((2, 1), _well_to_row_and_col('B1'))
        self.assertEqual((6, 10), _well_to_row_and_col('F10'))

        self.assertEqual((1, 1), _well_to_row_and_col('a1'))
        self.assertEqual((2, 1), _well_to_row_and_col('b1'))
        self.assertEqual((6, 10), _well_to_row_and_col('f10'))

    def test_decompress_well(self):
        self.assertEqual('A1', _decompress_well('A1'))
        self.assertEqual('H6', _decompress_well('O12'))

    def test_plate_position(self):
        self.assertEqual('1', _plate_position('A1'))
        self.assertEqual('1', _plate_position('A3'))
        self.assertEqual('1', _plate_position('A5'))

        self.assertEqual('4', _plate_position('B2'))
        self.assertEqual('4', _plate_position('B4'))
        self.assertEqual('4', _plate_position('B6'))

        self.assertEqual('2', _plate_position('C2'))
        self.assertEqual('2', _plate_position('C4'))
        self.assertEqual('2', _plate_position('C6'))

        self.assertEqual('3', _plate_position('D1'))
        self.assertEqual('3', _plate_position('D3'))
        self.assertEqual('3', _plate_position('D5'))

    def test_well_id_96(self):
        tests = []
        for a in 'ABCDEFGH':
            for b in range(1, 13):
                tests.append((f"{a}{b}", (a, b)))
        for a in 'IJK':
            for b in range(1, 13):
                tests.append((f"{a}{b}", None))
        for a in 'ABC':
            for b in [0, 13, 14]:
                tests.append((f"{a}{b}", None))

        for test, exp in tests:
            obs = _validate_well_id_96(test)
            self.assertEqual(obs, exp)


class MessageTests(TestCase):
    def test_message_construction(self):
        m = Message('In a world ...')

        self.assertIsNone(m._color)
        self.assertEqual(m.message, 'In a world ...')
        self.assertTrue(isinstance(m, Message))

    def test_error_construction(self):
        e = ErrorMessage('In a world ...')

        self.assertEqual(e._color, 'red')
        self.assertEqual(e.message, 'In a world ...')
        self.assertTrue(isinstance(e, ErrorMessage))

    def test_warning_construction(self):
        w = WarningMessage('In a world ...')

        self.assertEqual(w._color, 'yellow')
        self.assertEqual(w.message, 'In a world ...')
        self.assertTrue(isinstance(w, WarningMessage))

    def test_equal(self):
        a = Message('na na na')
        b = Message('na na na')
        c = Message('Batman')

        self.assertEqual(a, b)
        self.assertNotEqual(a, c)
        self.assertNotEqual(b, c)

        a._color = 'blue'
        self.assertNotEqual(a, b)

    def test_equal_warnings_and_errors(self):
        e = ErrorMessage('D:')
        w = WarningMessage(':D')
        self.assertNotEqual(e, w)

        self.assertEqual(e, ErrorMessage('D:'))
        self.assertNotEqual(e, ErrorMessage(':0'))

    def test_str(self):
        self.assertEqual('Message: test', str(Message('test')))
        self.assertEqual('ErrorMessage: test', str(ErrorMessage('test')))
        self.assertEqual('WarningMessage: test', str(WarningMessage('test')))

    def test_echo(self):
        e = ErrorMessage('Catch me if: you can')
        e.echo()

        # testing stdout: https://stackoverflow.com/a/12683001
        output = sys.stdout.getvalue().strip()
        self.assertEqual(output, 'ErrorMessage: Catch me if: you can')


class PlateValidationTests(TestCase):
    def setUp(self):
        self.metadata = [
            {
                'Plate Position': '1',
                'Project Abbreviation': 'ADAPT',
                'Plate map file':
                    './data/2022_summer_Celeste_Adaptation_16_plate_map.tsv',
                'Primer Plate #': '1',
                'Plating': 'SF',
                'Extraction Kit Lot': '166032128',
                'Extraction Robot': 'Carmen_HOWE_KF3',
                'TM1000 8 Tool': '109379Z',
                'Plate elution volume': '70',
                'TM300 8 Tool': 'NA',
                'TM50 8 Tool': 'NA',
                'Primer Date': '2021-08-17',
                'MasterMix Lot': '978215',
                'Water Lot': 'RNBJ0628',
                'Processing Robot': 'Echo550',
                'Sample Plate': 'THDMI_UK_Plate_2',
                'Project Name': 'THDMI UK',
                'Original Name': '',    # leave empty
                'TM10 8 Tool': '865HS8',
                'run_date': '2023-03-02',
                'instrument_model': 'Illumina MiSeq',
                'center_project_name': 'Rob ABTX',
                'experiment_design_description': ('### sequencing of anti'
                                                  'biotic time series')
            },
            {
                'Plate Position': '2',
                'Primer Plate #': '2',
                'Plate elution volume': '70',
                'Plate map file':
                    './data/2022_summer_Celeste_Adaptation_17_plate_map.tsv',
                'Project Abbreviation': 'ADAPT',
                'Plating': 'AS',
                'Extraction Kit Lot': '166032128',
                'Extraction Robot': 'Carmen_HOWE_KF4',
                'TM1000 8 Tool': '109379Z',
                'TM300 8 Tool': 'NA',
                'TM50 8 Tool': 'NA',
                'Primer Date': '2021-08-17',
                'MasterMix Lot': '978215',
                'Water Lot': 'RNBJ0628',
                'Processing Robot': 'Echo550',
                'Sample Plate': 'THDMI_UK_Plate_3',
                'Project Name': 'THDMI UK',
                'Original Name': '',
                'TM10 8 Tool': '865HS8',
                'run_date': '2023-03-02',
                'instrument_model': 'Illumina MiSeq',
                'center_project_name': 'Rob ABTX',
                'experiment_design_description': ('### sequencing of anti'
                                                  'biotic time series')
            },
            {
                'Plate Position': '3',
                'Primer Plate #': '3',
                'Plate elution volume': '70',
                'Plate map file':
                    './data/2022_summer_Celeste_Adaptation_18_plate_map.tsv',
                'Project Abbreviation': 'ADAPT',
                'Plating': 'MB_SF',
                'Extraction Kit Lot': '166032128',
                'Extraction Robot': 'Carmen_HOWE_KF3',
                'TM1000 8 Tool': '109379Z',
                'TM300 8 Tool': 'NA',
                'TM50 8 Tool': 'NA',
                'Primer Date': '2021-08-17',
                'MasterMix Lot': '978215',
                'Water Lot': 'RNBJ0628',
                'Processing Robot': 'Echo550',
                'Sample Plate': 'THDMI_UK_Plate_4',
                'Project Name': 'THDMI UK',
                'Original Name': '',
                'TM10 8 Tool': '865HS8',
                'run_date': '2023-03-02',
                'instrument_model': 'Illumina MiSeq',
                'center_project_name': 'Rob ABTX',
                'experiment_design_description': ('### sequencing of anti'
                                                  'biotic time series')
            },
            {
                'Plate Position': '4',
                'Primer Plate #': '4',
                'Plate elution volume': '70',
                'Plate map file':
                    './data/2022_summer_Celeste_Adaptation_21_plate_map.tsv',
                'Project Abbreviation': 'ADAPT',
                'Plating': 'AS',
                'Extraction Kit Lot': '166032128',
                'Extraction Robot': 'Carmen_HOWE_KF4',
                'TM1000 8 Tool': '109379Z',
                'TM300 8 Tool': 'NA',
                'TM50 8 Tool': 'NA',
                'Primer Date': '2021-08-17',
                'MasterMix Lot': '978215',
                'Water Lot': 'RNBJ0628',
                'Processing Robot': 'Echo550',
                'Sample Plate': 'THDMI_US_Plate_6',
                'Project Name': 'THDMI US',
                'Original Name': '',
                'TM10 8 Tool': '865HS8',
                'run_date': '2023-03-02',
                'instrument_model': 'Illumina MiSeq',
                'center_project_name': 'Rob ABTX',
                'experiment_design_description': ('### sequencing of anti'
                                                  'biotic time series')
            },
        ]

        self.context = {
            'primers': [],
            'names': [],
            'positions': []
        }

    def test_validate_plate_metadata_full(self):
        expected = pd.DataFrame(self.metadata)

        pd.testing.assert_frame_equal(validate_plate_metadata(self.metadata),
                                      expected)

    def test_validate_plate_metadata_full_extra_column(self):
        for item in self.metadata:
            item['NewColumn'] = 'NewValue'

        expected = pd.DataFrame(self.metadata)

        pd.testing.assert_frame_equal(validate_plate_metadata(self.metadata),
                                      expected)

    def test_validate_plate_metadata_returns_None(self):
        # add a repeated plate position
        self.metadata[1]['Plate Position'] = '1'
        self.assertIsNone(validate_plate_metadata(self.metadata))

        # testing stdout: https://stackoverflow.com/a/12683001
        output = sys.stdout.getvalue().strip()
        self.assertEqual(output, 'Messages for Plate 1 \n'
                         'ErrorMessage: The plate position "1" is repeated')

    def test_validate_plate_extra_columns(self):
        plate = self.metadata[0]
        plate['New Value'] = 'a deer'
        messages, context = _validate_plate(plate, self.context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         WarningMessage('The following columns are not '
                                        'recognized and may be misspelled '
                                        'column names: New Value'))

        expected = {'primers': ['1'], 'names': ['THDMI_UK_Plate_2'],
                    'positions': ['1']}
        self.assertEqual(context, expected)

    def test_validate_plate_missing_columns(self):
        plate = self.metadata[0]
        del plate['Plating']
        messages, context = _validate_plate(plate, self.context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         ErrorMessage('The following columns are missing: '
                                      'Plating'))

        expected = {'primers': ['1'], 'names': ['THDMI_UK_Plate_2'],
                    'positions': ['1']}
        self.assertEqual(context, expected)

    def test_validate_plate_repeated_primers(self):
        plate = self.metadata[0]
        context = {'primers': ['1'], 'names': ['THDMI_UK_Plate_3'],
                   'positions': ['2']}

        messages, context = _validate_plate(plate, context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         ErrorMessage('The Primer Plate "1" is repeated'))

        expected = {
            'primers': ['1', '1'],
            'names': ['THDMI_UK_Plate_3', 'THDMI_UK_Plate_2'],
            'positions': ['2', '1']
        }
        self.assertEqual(context, expected)

    def test_validate_plate_invalid_primer_plate(self):
        plate = self.metadata[0]
        plate['Primer Plate #'] = '11'
        context = {'primers': ['1'], 'names': ['THDMI_UK_Plate_3'],
                   'positions': ['2']}

        messages, context = _validate_plate(plate, context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         ErrorMessage('The Primer Plate # "11" is not between '
                                      '1-10'))

        expected = {
            'primers': ['1', '11'],
            'names': ['THDMI_UK_Plate_3', 'THDMI_UK_Plate_2'],
            'positions': ['2', '1']
        }
        self.assertEqual(context, expected)

    def test_validate_plate_rare_primer_plate(self):
        plate = self.metadata[0]
        plate['Primer Plate #'] = '9'
        context = {'primers': ['1'], 'names': ['THDMI_UK_Plate_3'],
                   'positions': ['2']}

        messages, context = _validate_plate(plate, context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         WarningMessage('Primer Plate # "9" is unusual, please'
                                        ' verify this value is correct'))

        expected = {
            'primers': ['1', '9'],
            'names': ['THDMI_UK_Plate_3', 'THDMI_UK_Plate_2'],
            'positions': ['2', '1']
        }
        self.assertEqual(context, expected)

    def test_validate_plate_repeated_names(self):
        plate = self.metadata[0]
        context = {'primers': ['2'], 'names': ['THDMI_UK_Plate_2'],
                   'positions': ['2']}

        messages, context = _validate_plate(plate, context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         ErrorMessage('The plate name "THDMI_UK_Plate_2" is '
                                      'repeated'))

        expected = {
            'primers': ['2', '1'],
            'names': ['THDMI_UK_Plate_2', 'THDMI_UK_Plate_2'],
            'positions': ['2', '1']
        }
        self.assertEqual(context, expected)

    def test_validate_plate_bad_position(self):
        plate = self.metadata[0]
        plate['Plate Position'] = '100'
        context = {'primers': ['2'], 'names': ['THDMI_UK_Plate_3'],
                   'positions': ['1']}

        messages, context = _validate_plate(plate, context)
        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0], ErrorMessage("Only the values '1', '2', "
                                                   "'3' and '4' are allowed in"
                                                   " the 'Plate Position' "
                                                   "field, you entered: "
                                                   "100"))

        expected = {
            'primers': ['2', '1'],
            'names': ['THDMI_UK_Plate_3', 'THDMI_UK_Plate_2'],
            'positions': ['1', '100']
        }
        self.assertEqual(context, expected)

    def test_validate_plate_repeated_position(self):
        plate = self.metadata[0]
        context = {'primers': ['2'], 'names': ['THDMI_UK_Plate_3'],
                   'positions': ['1']}

        messages, context = _validate_plate(plate, context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         ErrorMessage('The plate position "1" is repeated'))

        expected = {
            'primers': ['2', '1'],
            'names': ['THDMI_UK_Plate_3', 'THDMI_UK_Plate_2'],
            'positions': ['1', '1']
        }
        self.assertEqual(context, expected)

    def test_validate_no_empty_metadata(self):
        messages, context = _validate_plate({}, self.context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         ErrorMessage("Can't leave the first plate empty"))

        expected = {'primers': [], 'names': [],
                    'positions': []}
        self.assertEqual(context, expected)

    def test_validate_trailing_empty(self):
        context = {
            'primers': ['3', '2'],
            'names': ['THDMI_UK_Plate_3', 'THDMI_UK_Plate_2'],
            'positions': ['3', '2']
        }

        messages, context = _validate_plate({}, context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         WarningMessage("This plate has no metadata"))

        expected = {
            'primers': ['3', '2'],
            'names': ['THDMI_UK_Plate_3', 'THDMI_UK_Plate_2'],
            'positions': ['3', '2']
        }
        self.assertEqual(context, expected)

    def test_correct_date_format(self):
        plate = self.metadata[0]
        plate['Primer Date'] = '2000/01/01'
        messages, context = _validate_plate(plate, self.context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         ErrorMessage('Date format is invalid should be '
                                      'YYYY-mm-dd'))

        expected = {'primers': ['1'], 'names': ['THDMI_UK_Plate_2'],
                    'positions': ['1']}
        self.assertEqual(context, expected)

    def test_date_is_in_the_past(self):
        plate = self.metadata[0]
        plate['Primer Date'] = '2100-01-01'
        messages, context = _validate_plate(plate, self.context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         WarningMessage('The Primer Date is in the future'))

        expected = {'primers': ['1'], 'names': ['THDMI_UK_Plate_2'],
                    'positions': ['1']}
        self.assertEqual(context, expected)

    def test_non_ascii(self):
        plate = self.metadata[0]
        plate['Plating'] = 'é'
        messages, context = _validate_plate(plate, self.context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         ErrorMessage("The value for 'Plating' has non-ASCII "
                                      'characters'))

        expected = {'primers': ['1'], 'names': ['THDMI_UK_Plate_2'],
                    'positions': ['1']}
        self.assertEqual(context, expected)

    def test_sample_plate_has_no_spaces(self):
        plate = self.metadata[0]
        plate['Sample Plate'] = plate['Sample Plate'].replace('_', ' ')
        messages, context = _validate_plate(plate, self.context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         WarningMessage("Spaces are not recommended in the "
                                        "Sample Plate field"))

        expected = {'primers': ['1'], 'names': ['THDMI UK Plate 2'],
                    'positions': ['1']}
        self.assertEqual(context, expected)


class PlateReplicationTests(TestCase):
    def setUp(self):
        data_dir = os.path.dirname(__file__)
        input_plate_fp = os.path.join(data_dir, 'data/input_plate.tsv')
        self.input_df = pd.read_csv(input_plate_fp, sep='\t', dtype=str)

    def test_overwrite_source_quad(self):
        # replicate a valid source to empty sources 2 and 4 plus overwriting
        # source 3.

        pr = PlateReplication('Library Well')

        obs = pr.make_replicates(self.input_df, {1: [2, 3, 4]}, overwrite=True)

        exp = pd.read_csv('metapool/tests/data/file1.tsv',
                          sep='\t', dtype=str)

        obs = obs.set_index('Sample')
        exp = exp.set_index('Sample')

        assert_frame_equal(obs, exp, check_like=True)

        # replicate a valid source to empty sources 2 and 4 plus overwriting
        # source 3 with overwrites not allowed. Should return an Error.

        with self.assertRaisesRegex(ValueError, 'Quadrant 3 is a source '
                                                'quadrant'):
            pr.make_replicates(self.input_df, {1: [2, 3, 4]}, overwrite=False)

    def test_another_source(self):
        # confirm that a source other than quad 1 can replicate to empty
        # sources.
        pr = PlateReplication('Library Well')

        obs = pr.make_replicates(self.input_df, {3: [2, 4]}, overwrite=True)

        exp = pd.read_csv('metapool/tests/data/file2.tsv',
                          sep='\t', dtype=str)

        obs = obs.set_index('Sample')
        exp = exp.set_index('Sample')

        assert_frame_equal(obs, exp, check_like=True)

    def test_two_replications(self):
        # confirm that two sources can be replicated successfully.
        pr = PlateReplication('Library Well')

        obs = pr.make_replicates(self.input_df, {1: [2], 3: [4]},
                                 overwrite=True)

        exp = pd.read_csv('metapool/tests/data/file3.tsv',
                          sep='\t', dtype=str)

        obs = obs.set_index('Sample')
        exp = exp.set_index('Sample')

        assert_frame_equal(obs, exp, check_like=True)

    def test_parameter_conversion(self):
        # confirm conversion to lists works as intended.
        pr = PlateReplication('Library Well')

        obs = pr.make_replicates(self.input_df, {1: 2, 3: 4}, overwrite=True)

        exp = pd.read_csv('metapool/tests/data/file4.tsv',
                          sep='\t', dtype=str)

        obs = obs.set_index('Sample')
        exp = exp.set_index('Sample')

        assert_frame_equal(obs, exp, check_like=True)

    def test_replicate_empty_quad(self):
        # confirm replicating an empty quad to an empty quad raises an Error.

        pr = PlateReplication('Library Well')

        with self.assertRaisesRegex(ValueError, 'Quadrant 2 is not a source '
                                                'quadrant'):
            pr.make_replicates(self.input_df, {2: [4]}, overwrite=False)

        # confirm replicating an empty quad to a source quad raises an Error.
        with self.assertRaisesRegex(ValueError, 'Quadrant 2 is not a source '
                                                'quadrant'):
            pr.make_replicates(self.input_df, {2: [3]}, overwrite=False)

        # confirm two sources replicating to the same quad raises an Error.
        with self.assertRaisesRegex(ValueError, 'Quadrant 2 is already '
                                                'occupied with replicate '
                                                'samples'):
            pr.make_replicates(self.input_df, {1: [2], 3: [2]},
                               overwrite=False)

        # confirm source replicating to empty quad then empty quad attempting
        # to replicate and overwrite original raises an Error.
        # A source that's been replicated once shouldn't be allowed to be
        # overwritten.
        with self.assertRaisesRegex(ValueError, 'Quadrant 2 is not a source '
                                                'quadrant'):
            pr.make_replicates(self.input_df, {1: [2], 2: [1]},
                               overwrite=False)

    def test_get_short_name_and_id(self):
        # normal, non-pathological cases
        obs = get_short_name_and_id('project_1')
        self.assertEqual(2, len(obs))
        self.assertEqual(obs[0], 'project')
        self.assertEqual(obs[1], '1')

        obs = get_short_name_and_id('project_00333333')
        self.assertEqual(2, len(obs))
        self.assertEqual(obs[0], 'project')
        self.assertEqual(obs[1], '00333333')

        # pathological cases
        obs = get_short_name_and_id('project')
        self.assertEqual(2, len(obs))
        self.assertEqual(obs[0], 'project')
        self.assertIsNone(obs[1])

        obs = get_short_name_and_id('project_')
        self.assertEqual(2, len(obs))
        self.assertEqual(obs[0], 'project_')
        self.assertIsNone(obs[1])

    def test_parse_project_name_err_none(self):
        expected_err = "project_name cannot be None or empty string"
        with self.assertRaisesRegex(ValueError, expected_err):
            parse_project_name(None)

        with self.assertRaisesRegex(ValueError, expected_err):
            parse_project_name("")

    def test_parse_project_name_err_qiita_id(self):
        with self.assertRaisesRegex(
                ValueError, "'project' does not contain a Qiita-ID."):
            parse_project_name("project")

        with self.assertRaisesRegex(
                ValueError, "'project_blue' does not contain a Qiita-ID."):
            parse_project_name("project_blue")

    def test_parse_project_name(self):
        exp = {
            QIITA_ID_KEY: '1',
            PROJECT_SHORT_NAME_KEY: "project_green",
            PROJECT_FULL_NAME_KEY: "project_green_1"
        }
        obs = parse_project_name("project_green_1")
        self.assertDictEqual(exp, obs)


if __name__ == '__main__':
    assert not hasattr(sys.stdout, "getvalue")
    main(module=__name__, buffer=True, exit=False)
