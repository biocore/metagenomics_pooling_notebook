import pandas as pd

from unittest import TestCase, main
from metapool.plate import (_well_to_row_and_col, _decompress_well,
                            _plate_position, validate_plate_metadata,
                            _validate_plate, ErrorMessage, WarningMessage)


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
                'Primer Date': '2021-08-17',
                'MasterMix Lot': '978215',
                'Water Lot': 'RNBJ0628',
                'Processing Robot': 'Echo550',
                'Sample Plate': 'THDMI_UK_Plate_2',
                'Project_Name': 'THDMI UK',
                'Original Name': ''
            },
            {
                'Plate Position': '2',
                'Primer Plate #': '2',
                'Plating': 'AS',
                'Extraction Kit Lot': '166032128',
                'Extraction Robot': 'Carmen_HOWE_KF4',
                'TM1000 8 Tool': '109379Z',
                'Primer Date': '2021-08-17',
                'MasterMix Lot': '978215',
                'Water Lot': 'RNBJ0628',
                'Processing Robot': 'Echo550',
                'Sample Plate': 'THDMI_UK_Plate_3',
                'Project_Name': 'THDMI UK',
                'Original Name': ''
            },
            {
                'Plate Position': '3',
                'Primer Plate #': '3',
                'Plating': 'MB_SF',
                'Extraction Kit Lot': '166032128',
                'Extraction Robot': 'Carmen_HOWE_KF3',
                'TM1000 8 Tool': '109379Z',
                'Primer Date': '2021-08-17',
                'MasterMix Lot': '978215',
                'Water Lot': 'RNBJ0628',
                'Processing Robot': 'Echo550',
                'Sample Plate': 'THDMI_UK_Plate_4',
                'Project_Name': 'THDMI UK',
                'Original Name': ''
            },
            {
                'Plate Position': '4',
                'Primer Plate #': '4',
                'Plating': 'AS',
                'Extraction Kit Lot': '166032128',
                'Extraction Robot': 'Carmen_HOWE_KF4',
                'TM1000 8 Tool': '109379Z',
                'Primer Date': '2021-08-17',
                'MasterMix Lot': '978215',
                'Water Lot': 'RNBJ0628',
                'Processing Robot': 'Echo550',
                'Sample Plate': 'THDMI_US_Plate_6',
                'Project_Name': 'THDMI US',
                'Original Name': ''
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

    def test_validate_plate_extra_columns(self):
        plate = self.metadata[0]
        plate['New Value'] = 'a deer'
        messages, context = _validate_plate(plate, self.context)

        self.assertTrue(len(messages) == 1)
        self.assertEqual(messages[0],
                         ErrorMessage('The following columns are not needed: '
                                      'New Value'))

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
                         ErrorMessage('The primer plate "1" is repeated'))

        expected = {
            'primers': ['1', '1'],
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


if __name__ == '__main__':
    main()
