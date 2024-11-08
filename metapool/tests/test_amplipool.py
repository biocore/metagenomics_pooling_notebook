import pandas as pd

from unittest import TestCase, main
from metapool.amplipool import assign_emp_index, _load_emp_indices


class AmplipoolTests(TestCase):
    def setUp(self):
        self.plate_metadata = pd.DataFrame([
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
            }
        ])

        columns = ['Sample', 'Row', 'Col', 'Blank', 'Project Plate',
                   'Project Name', 'Compressed Plate Name', 'Well']
        data = [
            ['X00180471', 'A', 1, False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'A1'],
            ['X00180199', 'C', 1, False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'C1'],
            ['X00179789', 'E', 1, False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'E1'],
            ['X00180201', 'G', 1, False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'G1'],
            ['X00180464', 'I', 1, False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'I1'],
            ['X00179796', 'K', 1, False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'K1']]

        self.df = (pd.DataFrame(columns=columns, data=data).astype(
            dtype={'Col': 'str'}))
        self.seqtype = '16S'

    def test_assign_emp_index_position_one(self):
        obs = assign_emp_index(self.df, self.plate_metadata, self.seqtype)

        data = [
            ['X00180471', 'A', '1', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'A1', '1', '1', 'SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_2', 'THDMI UK', '', '1', 'A1',
             '515rcbc0', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'AGCCTTCGTCGC',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             ('AATGATACGGCGACCACCGAGATCTACACGCTAG'
              'CCTTCGTCGCTATGGTAATTGTGTGYCAGCMGCCGCGGTAA')],
            ['X00180199', 'C', '1', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'C1', '1', '1', 'SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_2', 'THDMI UK', '', '1', 'B1',
             '515rcbc12', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'CGTATAAATGCG',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             ('AATGATACGGCGACCACCGAGATCTACACGCTCG'
              'TATAAATGCGTATGGTAATTGTGTGYCAGCMGCCGCGGTAA')],
            ['X00179789', 'E', '1', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'E1', '1', '1', 'SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_2', 'THDMI UK', '', '1', 'C1',
             '515rcbc24', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'TGACTAATGGCC',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             ('AATGATACGGCGACCACCGAGATCTACACGCTTGAC'
              'TAATGGCCTATGGTAATTGTGTGYCAGCMGCCGCGGTAA')],
            ['X00180201', 'G', '1', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'G1', '1', '1', 'SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_2', 'THDMI UK', '', '1', 'D1',
             '515rcbc36', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'GTGGAGTCTCAT',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             ('AATGATACGGCGACCACCGAGATCTACACGC'
              'TGTGGAGTCTCATTATGGTAATTGTGTGYCAGCMGCCGCGGTAA')],
            ['X00180464', 'I', '1', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'I1', '1', '1', 'SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_2', 'THDMI UK', '', '1', 'E1',
             '515rcbc48', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'TGATGTGCTAAG',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             ('AATGATACGGCGACCACCGAGATCTACACGCTTGATGTG'
              'CTAAGTATGGTAATTGTGTGYCAGCMGCCGCGGTAA')],
            ['X00179796', 'K', '1', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'K1', '1', '1', 'SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_2', 'THDMI UK', '', '1', 'F1',
             '515rcbc60', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'TGTGCACGCCAT',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             ('AATGATACGGCGACCACCGAG'
              'ATCTACACGCTTGTGCACGCCATTATGGTAATTGTGTGYCAGCMGCCGCGGTAA')]
        ]

        exp = pd.DataFrame(
            columns=['Sample', 'Row', 'Col', 'Blank', 'Project Plate',
                     'Project Name', 'Compressed Plate Name', 'Well',
                     'Plate Position', 'Primer Plate #', 'Plating',
                     'Extraction Kit Lot', 'Extraction Robot', 'TM1000 8 Tool',
                     'Primer Date', 'MasterMix Lot', 'Water Lot',
                     'Processing Robot', 'Sample Plate', 'Project_Name',
                     'Original Name', 'Plate', 'EMP Primer Plate Well', 'Name',
                     "Illumina 5prime Adapter", 'Golay Barcode',
                     'Forward Primer Pad', 'Forward Primer Linker',
                     '515FB Forward Primer (Parada)', 'Primer For PCR'],
            data=data
        )

        pd.testing.assert_frame_equal(obs, exp)

    def test_assign_emp_index_multiple_positions(self):
        # change some of the well ids and their primer plates to spot check
        # that correct barcodes are retrieved from the EMP indices file

        # position 1 gets primer plate 5
        self.plate_metadata.loc[0, 'Primer Plate #'] = '5'
        self.df.loc[0, 'Row'] = 'A'
        self.df.loc[0, 'Col'] = '1'
        self.df.loc[0, 'Well'] = 'A1'

        # position 2 gets primer plate 6
        self.plate_metadata.loc[1, 'Primer Plate #'] = '6'
        self.df.loc[1, 'Row'] = 'A'
        self.df.loc[1, 'Col'] = '4'
        self.df.loc[1, 'Well'] = 'A4'

        # position 3 gets primer plate 9
        self.plate_metadata.loc[2, 'Primer Plate #'] = '9'
        self.df.loc[2, 'Row'] = 'D'
        self.df.loc[2, 'Col'] = '7'
        self.df.loc[2, 'Well'] = 'D7'

        self.df.loc[3, 'Row'] = 'F'
        self.df.loc[3, 'Col'] = '9'
        self.df.loc[3, 'Well'] = 'F9'

        # position 4 gets primer plate 10
        self.plate_metadata.loc[3, 'Primer Plate #'] = '10'
        self.df.loc[4, 'Row'] = 'B'
        self.df.loc[4, 'Col'] = '6'
        self.df.loc[4, 'Well'] = 'B6'

        self.df.loc[5, 'Row'] = 'F'
        self.df.loc[5, 'Col'] = '10'
        self.df.loc[5, 'Well'] = 'F10'

        obs1 = assign_emp_index(self.df, self.plate_metadata, '16S')
        obs2 = assign_emp_index(self.df, self.plate_metadata, '18S')
        obs3 = assign_emp_index(self.df, self.plate_metadata, 'ITS')

        data1 = [
            ['X00180471', 'A', '1', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'A1', '1', '5', 'SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_2', 'THDMI UK', '', '5', 'A1',
             '515rcbc384', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'ATGTTAGGGAAT',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             ('AATGATACGGCGACCACCGAGATCTACACGCTATGTTAGG'
              'GAATTATGGTAATTGTGTGYCAGCMGCCGCGGTAA')],
            ['X00180199', 'A', '4', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'A4', '2', '6', 'AS', '166032128',
             'Carmen_HOWE_KF4', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_3', 'THDMI UK', '', '6', 'A2',
             '515rcbc481', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'CTACCGATTGCG',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             ('AATGATACGGCGACCACCGAGATCTACACGCTCTACCGAT'
              'TGCGTATGGTAATTGTGTGYCAGCMGCCGCGGTAA')],
            ['X00179789', 'D', '7', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'D7', '3', '9', 'MB_SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_4', 'THDMI UK', '', '9', 'B4',
             '515rcbc783', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'GGTTACGGTTAC',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             ('AATGATACGGCGACCACCGAGATCTACACGCTGGTTACGG'
              'TTACTATGGTAATTGTGTGYCAGCMGCCGCGGTAA')],
            ['X00180201', 'F', '9', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'F9', '3', '9', 'MB_SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_4', 'THDMI UK', '', '9', 'C5',
             '515rcbc796', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'GATCTGCGATCC',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             ('AATGATACGGCGACCACCGAGATCTACACGCTGATCTGCG'
              'ATCCTATGGTAATTGTGTGYCAGCMGCCGCGGTAA')],
            ['X00180464', 'B', '6', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'B6', '4', '10', 'AS', '166032128',
             'Carmen_HOWE_KF4', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_US_Plate_6', 'THDMI US', '', '10', 'A3',
             '515rcbc866', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'TTGACGACATCG',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             ('AATGATACGGCGACCACCGAGATCTACACGCTTTGACGAC'
              'ATCGTATGGTAATTGTGTGYCAGCMGCCGCGGTAA')],
            ['X00179796', 'F', '10', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'F10', '4', '10', 'AS', '166032128',
             'Carmen_HOWE_KF4', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_US_Plate_6', 'THDMI US', '', '10', 'C5',
             '515rcbc892', 'AATGATACGGCGACCACCGAGATCTACACGCT', 'TAGAGGCGTAGG',
             'TATGGTAATT', 'GT', 'GTGYCAGCMGCCGCGGTAA',
             ('AATGATACGGCGACCACCGAGATCTACACGCTTAGAGGCG'
              'TAGGTATGGTAATTGTGTGYCAGCMGCCGCGGTAA')]
        ]
        data2 = [
            ['X00180471', 'A', '1', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'A1', '1', '5', 'SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_2', 'THDMI UK', '', '5', 'A1',
             'EukBr_Hiseq_0411', 'CAAGCAGAAGACGGCATACGAGAT', 'AACAAACTGCCA',
             'AGTCAGTCAG', 'CA', 'TGATCCTTCTGCAGGTTCACCTAC',
             ('CAAGCAGAAGACGGCATACGAGATAACAAACTGCCAAGTCAG'
              'TCAGCATGATCCTTCTGCAGGTTCACCTAC')],
            ['X00180199', 'A', '4', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'A4', '2', '6', 'AS', '166032128',
             'Carmen_HOWE_KF4', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_3', 'THDMI UK', '', '6', 'A2',
             'EukBr_Hiseq_0508', 'CAAGCAGAAGACGGCATACGAGAT', 'GGAGAGATCACG',
             'AGTCAGTCAG', 'CA', 'TGATCCTTCTGCAGGTTCACCTAC',
             ('CAAGCAGAAGACGGCATACGAGATGGAGAGATCACGAGTCAG'
              'TCAGCATGATCCTTCTGCAGGTTCACCTAC')],
            ['X00179789', 'D', '7', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'D7', '3', '9', 'MB_SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_4', 'THDMI UK', '', '9', 'B4',
             'EukBr_Hiseq_0810', 'CAAGCAGAAGACGGCATACGAGAT', 'CATTTGACGACG',
             'AGTCAGTCAG', 'CA', 'TGATCCTTCTGCAGGTTCACCTAC',
             ('CAAGCAGAAGACGGCATACGAGATCATTTGACGACGAGTCAG'
              'TCAGCATGATCCTTCTGCAGGTTCACCTAC')],
            ['X00180201', 'F', '9', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'F9', '3', '9', 'MB_SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_4', 'THDMI UK', '', '9', 'C5',
             'EukBr_Hiseq_0823', 'CAAGCAGAAGACGGCATACGAGAT', 'CACGTTTATTCC',
             'AGTCAGTCAG', 'CA', 'TGATCCTTCTGCAGGTTCACCTAC',
             ('CAAGCAGAAGACGGCATACGAGATCACGTTTATTCCAGTCAG'
              'TCAGCATGATCCTTCTGCAGGTTCACCTAC')],
            ['X00180464', 'B', '6', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'B6', '4', '10', 'AS', '166032128',
             'Carmen_HOWE_KF4', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_US_Plate_6', 'THDMI US', '', '10', 'A3',
             'EukBr_Hiseq_0893', 'CAAGCAGAAGACGGCATACGAGAT', 'GTATGGAGCTAT',
             'AGTCAGTCAG', 'CA', 'TGATCCTTCTGCAGGTTCACCTAC',
             ('CAAGCAGAAGACGGCATACGAGATGTATGGAGCTATAGTCAG'
              'TCAGCATGATCCTTCTGCAGGTTCACCTAC')],
            ['X00179796', 'F', '10', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'F10', '4', '10', 'AS', '166032128',
             'Carmen_HOWE_KF4', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_US_Plate_6', 'THDMI US', '', '10', 'C5',
             'EukBr_Hiseq_0919', 'CAAGCAGAAGACGGCATACGAGAT', 'CAGCCTGCAAAT',
             'AGTCAGTCAG', 'CA', 'TGATCCTTCTGCAGGTTCACCTAC',
             ('CAAGCAGAAGACGGCATACGAGATCAGCCTGCAAATAGTCAG'
              'TCAGCATGATCCTTCTGCAGGTTCACCTAC')]
        ]
        data3 = [
            ['X00180471', 'A', '1', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'A1', '1', '5', 'SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_2', 'THDMI UK', '', '5', 'A1',
             'kabir_ITS2rcbc384', 'CAAGCAGAAGACGGCATACGAGAT', 'GGAATTATCGGT',
             '', 'CG', 'GCTGCGTTCTTCATCGATGC',
             ('CAAGCAGAAGACGGCATACGAGATGGAATTAT'
              'CGGTCGGCTGCGTTCTTCATCGATGC')],
            ['X00180199', 'A', '4', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'A4', '2', '6', 'AS', '166032128',
             'Carmen_HOWE_KF4', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_3', 'THDMI UK', '', '6', 'A2',
             'kabir_ITS2rcbc481', 'CAAGCAGAAGACGGCATACGAGAT', 'TGACGTAGAACT',
             '', 'CG', 'GCTGCGTTCTTCATCGATGC',
             ('CAAGCAGAAGACGGCATACGAGATTGACGTAG'
              'AACTCGGCTGCGTTCTTCATCGATGC')],
            ['X00179789', 'D', '7', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'D7', '3', '9', 'MB_SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_4', 'THDMI UK', '', '9', 'B4',
             'kabir_ITS2rcbc783', 'CAAGCAGAAGACGGCATACGAGAT', 'ATCGATCCACAG',
             '', 'CG', 'GCTGCGTTCTTCATCGATGC',
             ('CAAGCAGAAGACGGCATACGAGATATCGATCC'
              'ACAGCGGCTGCGTTCTTCATCGATGC')],
            ['X00180201', 'F', '9', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'F9', '3', '9', 'MB_SF', '166032128',
             'Carmen_HOWE_KF3', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_UK_Plate_4', 'THDMI UK', '', '9', 'C5',
             'kabir_ITS2rcbc796', 'CAAGCAGAAGACGGCATACGAGAT', 'TCCAGGGCTATA',
             '', 'CG', 'GCTGCGTTCTTCATCGATGC',
             ('CAAGCAGAAGACGGCATACGAGATTCCAGGGC'
              'TATACGGCTGCGTTCTTCATCGATGC')],
            ['X00180464', 'B', '6', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'B6', '4', '10', 'AS', '166032128',
             'Carmen_HOWE_KF4', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_US_Plate_6', 'THDMI US', '', '10', 'A3',
             'kabir_ITS2rcbc866', 'CAAGCAGAAGACGGCATACGAGAT', 'ACATGTCACGTG',
             '', 'CG', 'GCTGCGTTCTTCATCGATGC',
             ('CAAGCAGAAGACGGCATACGAGATACATGTCA'
              'CGTGCGGCTGCGTTCTTCATCGATGC')],
            ['X00179796', 'F', '10', False, 'THDMI_10317_PUK2', 'THDMI_10317',
             'THDMI_10317_UK2-US6', 'F10', '4', '10', 'AS', '166032128',
             'Carmen_HOWE_KF4', '109379Z', '2021-08-17', '978215', 'RNBJ0628',
             'Echo550', 'THDMI_US_Plate_6', 'THDMI US', '', '10', 'C5',
             'kabir_ITS2rcbc892', 'CAAGCAGAAGACGGCATACGAGAT', 'ACGTGAGGAACG',
             '', 'CG', 'GCTGCGTTCTTCATCGATGC',
             ('CAAGCAGAAGACGGCATACGAGATACGTGAGG'
              'AACGCGGCTGCGTTCTTCATCGATGC')]
        ]

        exp1 = pd.DataFrame(
            columns=['Sample', 'Row', 'Col', 'Blank', 'Project Plate',
                     'Project Name', 'Compressed Plate Name', 'Well',
                     'Plate Position', 'Primer Plate #', 'Plating',
                     'Extraction Kit Lot', 'Extraction Robot',
                     'TM1000 8 Tool',
                     'Primer Date', 'MasterMix Lot', 'Water Lot',
                     'Processing Robot', 'Sample Plate', 'Project_Name',
                     'Original Name', 'Plate',
                     'EMP Primer Plate Well', 'Name',
                     "Illumina 5prime Adapter", 'Golay Barcode',
                     'Forward Primer Pad', 'Forward Primer Linker',
                     '515FB Forward Primer (Parada)', 'Primer For PCR'],
            data=data1
        )
        exp2 = pd.DataFrame(
            columns=['Sample', 'Row', 'Col', 'Blank', 'Project Plate',
                     'Project Name', 'Compressed Plate Name', 'Well',
                     'Plate Position', 'Primer Plate #', 'Plating',
                     'Extraction Kit Lot', 'Extraction Robot',
                     'TM1000 8 Tool',
                     'Primer Date', 'MasterMix Lot', 'Water Lot',
                     'Processing Robot', 'Sample Plate', 'Project_Name',
                     'Original Name', 'Plate',
                     'EMP Primer Plate Well', 'Name',
                     'Reverse complement of 3prime Illumina Adapter',
                     'Golay Barcode',
                     'Reverse Primer Pad', 'Reverse Primer Linker',
                     'Reverse primer (EukBr)', 'Primer For PCR'],
            data=data2
        )
        exp3 = pd.DataFrame(
            columns=['Sample', 'Row', 'Col', 'Blank', 'Project Plate',
                     'Project Name', 'Compressed Plate Name', 'Well',
                     'Plate Position', 'Primer Plate #', 'Plating',
                     'Extraction Kit Lot', 'Extraction Robot',
                     'TM1000 8 Tool',
                     'Primer Date', 'MasterMix Lot', 'Water Lot',
                     'Processing Robot', 'Sample Plate', 'Project_Name',
                     'Original Name', 'Plate',
                     'EMP Primer Plate Well', 'Name',
                     'Reverse complement of 3prime Illumina Adapter',
                     'Golay Barcode',
                     'Reverse Primer Pad', 'Reverse Primer Linker',
                     'ITS2 Reverse Primer', 'Primer For PCR'],
            data=data3
        )

        pd.testing.assert_frame_equal(obs1, exp1)
        pd.testing.assert_frame_equal(obs2, exp2)
        pd.testing.assert_frame_equal(obs3, exp3)

    def test_load_emp_indices(self):
        obs1 = _load_emp_indices('16S')
        obs2 = _load_emp_indices('18S')
        obs3 = _load_emp_indices('ITS')

        # no NaN values of any kind
        pd.testing.assert_frame_equal(obs1, obs1.dropna(how='any'))
        pd.testing.assert_frame_equal(obs2, obs2.dropna(how='any'))
        pd.testing.assert_frame_equal(obs3, obs3.dropna(how='any'))

        self.assertTrue(len(obs1), 961)
        self.assertTrue(len(obs2), 961)
        self.assertTrue(len(obs3), 961)

        expected1 = [
            "Plate", "Well", "Name", "Illumina 5prime Adapter",
            "Golay Barcode",
            "Forward Primer Pad", "Forward Primer Linker",
            "515FB Forward Primer (Parada)", "Primer For PCR"
        ]
        expected2 = [
            "Plate", "Well", "Name",
            "Reverse complement of 3prime Illumina Adapter",
            "Golay Barcode",
            "Reverse Primer Pad", "Reverse Primer Linker",
            "Reverse primer (EukBr)", "Primer For PCR"
        ]
        expected3 = [
            "Plate", "Well", "Name",
            "Reverse complement of 3prime Illumina Adapter",
            "Golay Barcode",
            "Reverse Primer Pad", "Reverse Primer Linker",
            "ITS2 Reverse Primer", "Primer For PCR"
        ]

        self.assertTrue(obs1.columns.tolist(), expected1)
        self.assertTrue(obs2.columns.tolist(), expected2)
        self.assertTrue(obs3.columns.tolist(), expected3)


if __name__ == '__main__':
    main()
