from unittest import TestCase, main
from metapool.plate import (_well_to_row_and_col, _decompress_well,
                            _plate_position, validate_plate_metadata,
                            _validate_plate)


class PlateTests(TestCase):
    def setUp(self):
        pass

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


if __name__ == '__main__':
    main()
