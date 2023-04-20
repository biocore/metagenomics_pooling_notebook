import sys
import pandas as pd

from unittest import TestCase, main
from metapool.plate import (_well_to_row_and_col, _decompress_well,
                            _plate_position, validate_plate_metadata,
                            _validate_plate, Message, ErrorMessage,
                            WarningMessage, requires_dilution, dilute_gDNA,
                            find_threshold, autopool, _validate_well_id_96)
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
        # a range of test inputs, including all valid inputs and some invalid
        # ones.
        test_input = ['A0', 'B0', 'C0', 'D0', 'E0', 'F0', 'G0', 'H0', 'I0',
                      'A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1', 'I1',
                      'A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2', 'I2',
                      'A3', 'B3', 'C3', 'D3', 'E3', 'F3', 'G3', 'H3', 'I3',
                      'A4', 'B4', 'C4', 'D4', 'E4', 'F4', 'G4', 'H4', 'I4',
                      'A5', 'B5', 'C5', 'D5', 'E5', 'F5', 'G5', 'H5', 'I5',
                      'A6', 'B6', 'C6', 'D6', 'E6', 'F6', 'G6', 'H6', 'I6',
                      'A7', 'B7', 'C7', 'D7', 'E7', 'F7', 'G7', 'H7', 'I7',
                      'A8', 'B8', 'C8', 'D8', 'E8', 'F8', 'G8', 'H8', 'I8',
                      'A9', 'B9', 'C9', 'D9', 'E9', 'F9', 'G9', 'H9', 'I9',
                      'A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10',
                      'I10', 'A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11',
                      'H11', 'I11', 'A12', 'B12', 'C12', 'D12', 'E12', 'F12',
                      'G12', 'H12', 'I12', 'A13', 'B13', 'C13', 'D13', 'E13',
                      'F13', 'G13', 'H13', 'I13']

        results = []

        for well in test_input:
            results.append(_validate_well_id_96(well))

        # invalid inputs will cause _validate_well_id_96() to return None.
        obs = [x for x in results if x != None]

        exp = [('A', 1), ('B', 1), ('C', 1), ('D', 1), ('E', 1), ('F', 1),
               ('G', 1), ('H', 1), ('A', 2), ('B', 2), ('C', 2), ('D', 2),
               ('E', 2), ('F', 2), ('G', 2), ('H', 2), ('A', 3), ('B', 3),
               ('C', 3), ('D', 3), ('E', 3), ('F', 3), ('G', 3), ('H', 3),
               ('A', 4), ('B', 4), ('C', 4), ('D', 4), ('E', 4), ('F', 4),
               ('G', 4), ('H', 4), ('A', 5), ('B', 5), ('C', 5), ('D', 5),
               ('E', 5), ('F', 5), ('G', 5), ('H', 5), ('A', 6), ('B', 6),
               ('C', 6), ('D', 6), ('E', 6), ('F', 6), ('G', 6), ('H', 6),
               ('A', 7), ('B', 7), ('C', 7), ('D', 7), ('E', 7), ('F', 7),
               ('G', 7), ('H', 7), ('A', 8), ('B', 8), ('C', 8), ('D', 8),
               ('E', 8), ('F', 8), ('G', 8), ('H', 8), ('A', 9), ('B', 9),
               ('C', 9), ('D', 9), ('E', 9), ('F', 9), ('G', 9), ('H', 9),
               ('A', 10), ('B', 10), ('C', 10), ('D', 10), ('E', 10),
               ('F', 10), ('G', 10), ('H', 10), ('A', 11), ('B', 11),
               ('C', 11), ('D', 11), ('E', 11), ('F', 11), ('G', 11),
               ('H', 11), ('A', 12), ('B', 12), ('C', 12), ('D', 12),
               ('E', 12), ('F', 12), ('G', 12), ('H', 12)]

        self.assertEqual(set(obs), set(exp))


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
                'Primer Plate #': '1',
                'Plating': 'SF',
                'Extraction Kit Lot': '166032128',
                'Extraction Robot': 'Carmen_HOWE_KF3',
                'TM1000 8 Tool': '109379Z',
                'TM300 8 Tool': 'NA',
                'TM50 8 Tool': 'NA',
                'Primer Date': '2021-08-17',
                'MasterMix Lot': '978215',
                'Water Lot': 'RNBJ0628',
                'Processing Robot': 'Echo550',
                'Sample Plate': 'THDMI_UK_Plate_2',
                'Project_Name': 'THDMI UK',
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
                'Project_Name': 'THDMI UK',
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
                'Project_Name': 'THDMI UK',
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
                'Project_Name': 'THDMI US',
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


if __name__ == '__main__':
    assert not hasattr(sys.stdout, "getvalue")
    main(module=__name__, buffer=True, exit=False)
