import unittest

from metapool.metadata import get_platemap_data
class Test1(unittest.TestCase):
    def get_zip(self):
        self.assertEqual(get_platemap_data('PL'), 'CA')