from unittest import TestCase, main

import os
import tempfile
import warnings
import pandas as pd
import numpy as np
import numpy.testing as npt
from io import StringIO

from metapool.metapool import (read_plate_map_csv, read_pico_csv,
            calculate_norm_vol,
            format_dna_norm_picklist, assign_index, format_index_picklist,
            compute_qpcr_concentration, compute_shotgun_pooling_values_eqvol,
            compute_shotgun_pooling_values_qpcr,
            compute_shotgun_pooling_values_qpcr_minvol, estimate_pool_conc_vol,
            format_pooling_echo_pick_list, plot_plate_vals, make_2D_array,
            combine_dfs, parse_dna_conc_csv, add_dna_conc,
            compute_pico_concentration, ss_temp, format_sheet_comments,
            format_sample_sheet, bcl_scrub_name, rc, sequencer_i5_index,
            format_sample_data, reformat_interleaved_to_columns,
            validate_sample_sheet, parse_sample_sheet, write_sample_sheet)


class Tests(TestCase):

    def setUp(self):
        self.maxDiff = None
        self.cp_vals = np.array([[10.14, 7.89, 7.9, 15.48],
                                 [7.86, 8.07, 8.16, 9.64],
                                 [12.29, 7.64, 7.32, 13.74]])

        self.dna_vals = np.array([[10.14, 7.89, 7.9, 15.48],
                                  [7.86, 8.07, 8.16, 9.64],
                                  [12.29, 7.64, 7.32, 13.74]])

        self.qpcr_conc = \
            np.array([[98.14626462, 487.8121413, 484.3480866, 2.183406934],
                      [498.3536649, 429.0839787, 402.4270321, 140.1601735],
                      [21.20533391, 582.9456031, 732.2655041, 7.545145988]])

        self.pico_conc = \
            np.array([[38.4090909, 29.8863636, 29.9242424, 58.6363636],
                      [29.7727273, 30.5681818, 30.9090909, 36.5151515],
                      [46.5530303, 28.9393939, 27.7272727, 52.0454545]])

        data_dir = os.path.join(os.path.dirname(__file__), 'data')

        self.good_ss = os.path.join(data_dir, 'good-sample-sheet.csv')
        self.no_comments_ss = os.path.join(data_dir, 'no-comments-'
                                           'sample-sheet.csv')
        self.no_project_ss = os.path.join(data_dir,
                                          'no-project-name-sample-sheet.csv')

        # "valid" upfront but will have repeated values after scrubbing
        self.ok_ss = os.path.join(data_dir, 'ok-sample-sheet.csv')

        self.scrubbable_ss = os.path.join(data_dir,
                                          'scrubbable-sample-sheet.csv')

        self.bad_project_name_ss = os.path.join(data_dir, 'bad-project'
                                                '-name-sample-sheet.csv')

    # def test_compute_shotgun_normalization_values(self):
    #     input_vol = 3.5
    #     input_dna = 10
    #     plate_layout = []
    #     for i in range(4):
    #         row = []
    #         for j in range(4):
    #             row.append({'dna_concentration': 10,
    #                         'sample_id': "S%s.%s" % (i, j)})
    #         plate_layout.append(row)

    #     obs_sample, obs_water = compute_shotgun_normalization_values(
    #         plate_layout, input_vol, input_dna)

    #     exp_sample = np.zeros((4, 4), dtype=np.float)
    #     exp_water = np.zeros((4, 4), dtype=np.float)
    #     exp_sample.fill(1000)
    #     exp_water.fill(2500)

    #     npt.assert_almost_equal(obs_sample, exp_sample)
    #     npt.assert_almost_equal(obs_water, exp_water)

    #     # Make sure that we don't go above the limit
    #     plate_layout[1][1]['dna_concentration'] = 0.25
    #     obs_sample, obs_water = compute_shotgun_normalization_values(
    #         plate_layout, input_vol, input_dna)

    #     exp_sample[1][1] = 3500
    #     exp_water[1][1] = 0

    #     npt.assert_almost_equal(obs_sample, exp_sample)
    #     npt.assert_almost_equal(obs_water, exp_water)
    def test_read_plate_map_csv(self):
        plate_map_csv = \
        'Sample\tRow\tCol\tBlank\n' + \
        'sam1\tA\t1\tFalse\n' + \
        'sam2\tA\t2\tFalse\n' + \
        'blank1\tB\t1\tTrue\n' + \
        'sam3\tB\t2\tFalse\n' 

        plate_map_f = StringIO(plate_map_csv)

        exp_plate_df = pd.DataFrame({'Sample': ['sam1','sam2','blank1','sam3'],
                                     'Row': ['A','A','B','B'],
                                     'Col': [1,2,1,2],
                                     'Well': ['A1','A2','B1','B2'],
                                     'Blank': [False, False, True, False]})

        obs_plate_df = read_plate_map_csv(plate_map_f)

        pd.testing.assert_frame_equal(obs_plate_df, exp_plate_df, check_like=True)

    def test_read_pico_csv(self):
        # Test a normal sheet
        pico_csv = '''Results

        Well ID\tWell\t[Blanked-RFU]\t[Concentration]
        SPL1\tA1\t5243.000\t3.432
        SPL2\tA2\t4949.000\t3.239
        SPL3\tB1\t15302.000\t10.016
        SPL4\tB2\t4039.000\t2.644

        Curve2 Fitting Results

        Curve Name\tCurve Formula\tA\tB\tR2\tFit F Prob
        Curve2\tY=A*X+B\t1.53E+003\t0\t0.995\t?????
        '''
        exp_pico_df = pd.DataFrame({'Well': ['A1','A2','B1','B2'],
                                    'Sample DNA Concentration': 
                                     [3.432, 3.239, 10.016, 2.644]})

        pico_csv_f = StringIO(pico_csv)

        obs_pico_df = read_pico_csv(pico_csv_f)

        pd.testing.assert_frame_equal(obs_pico_df, exp_pico_df, check_like=True)

        # Test a sheet that has some ???? zero values
        pico_csv = '''Results

        Well ID\tWell\t[Blanked-RFU]\t[Concentration]
        SPL1\tA1\t5243.000\t3.432
        SPL2\tA2\t4949.000\t3.239
        SPL3\tB1\t15302.000\t10.016
        SPL4\tB2\t\t?????

        Curve2 Fitting Results

        Curve Name\tCurve Formula\tA\tB\tR2\tFit F Prob
        Curve2\tY=A*X+B\t1.53E+003\t0\t0.995\t?????
        '''
        exp_pico_df = pd.DataFrame({'Well': ['A1','A2','B1','B2'],
                                    'Sample DNA Concentration': 
                                     [3.432, 3.239, 10.016, np.nan]})

        pico_csv_f = StringIO(pico_csv)

        obs_pico_df = read_pico_csv(pico_csv_f)

        pd.testing.assert_frame_equal(obs_pico_df, exp_pico_df, check_like=True)

    def test_calculate_norm_vol(self):
        dna_concs = np.array([[2, 7.89],
                              [np.nan, .0]])

        exp_vols = np.array([[2500., 632.5],
                              [3500., 3500.]])

        obs_vols = calculate_norm_vol(dna_concs)

        np.testing.assert_allclose(exp_vols, obs_vols)

    def test_format_dna_norm_picklist(self):

        exp_picklist = \
        'Sample\tSource Plate Name\tSource Plate Type\tSource Well\t' + \
        'Concentration\tTransfer Volume\tDestination Plate Name\tDestination Well\n' + \
        'sam1\tWater\t384PP_AQ_BP2_HT\tA1\t2.0\t1000.0\tNormalizedDNA\tA1\n' + \
        'sam2\tWater\t384PP_AQ_BP2_HT\tA2\t7.89\t2867.5\tNormalizedDNA\tA2\n' + \
        'blank1\tWater\t384PP_AQ_BP2_HT\tB1\tnan\t0.0\tNormalizedDNA\tB1\n' + \
        'sam3\tWater\t384PP_AQ_BP2_HT\tB2\t0.0\t0.0\tNormalizedDNA\tB2\n' + \
        'sam1\tSample\t384PP_AQ_BP2_HT\tA1\t2.0\t2500.0\tNormalizedDNA\tA1\n' + \
        'sam2\tSample\t384PP_AQ_BP2_HT\tA2\t7.89\t632.5\tNormalizedDNA\tA2\n' + \
        'blank1\tSample\t384PP_AQ_BP2_HT\tB1\tnan\t3500.0\tNormalizedDNA\tB1\n' + \
        'sam3\tSample\t384PP_AQ_BP2_HT\tB2\t0.0\t3500.0\tNormalizedDNA\tB2'

        dna_vols = np.array([[2500., 632.5],
                              [3500., 3500.]])

        water_vols = 3500 - dna_vols

        wells = np.array([['A1', 'A2'],
                          ['B1', 'B2']])

        sample_names =  np.array([['sam1', 'sam2'],
                          ['blank1', 'sam3']])

        dna_concs = np.array([[2, 7.89],
                              [np.nan, .0]])

        obs_picklist = format_dna_norm_picklist(dna_vols, water_vols, wells,
                                                sample_names = sample_names,
                                                dna_concs = dna_concs)

        self.assertEqual(exp_picklist, obs_picklist)

        # test if switching dest wells
        exp_picklist = \
        'Sample\tSource Plate Name\tSource Plate Type\tSource Well\t' + \
        'Concentration\tTransfer Volume\tDestination Plate Name\tDestination Well\n' + \
        'sam1\tWater\t384PP_AQ_BP2_HT\tA1\t2.0\t1000.0\tNormalizedDNA\tD1\n' + \
        'sam2\tWater\t384PP_AQ_BP2_HT\tA2\t7.89\t2867.5\tNormalizedDNA\tD2\n' + \
        'blank1\tWater\t384PP_AQ_BP2_HT\tB1\tnan\t0.0\tNormalizedDNA\tE1\n' + \
        'sam3\tWater\t384PP_AQ_BP2_HT\tB2\t0.0\t0.0\tNormalizedDNA\tE2\n' + \
        'sam1\tSample\t384PP_AQ_BP2_HT\tA1\t2.0\t2500.0\tNormalizedDNA\tD1\n' + \
        'sam2\tSample\t384PP_AQ_BP2_HT\tA2\t7.89\t632.5\tNormalizedDNA\tD2\n' + \
        'blank1\tSample\t384PP_AQ_BP2_HT\tB1\tnan\t3500.0\tNormalizedDNA\tE1\n' + \
        'sam3\tSample\t384PP_AQ_BP2_HT\tB2\t0.0\t3500.0\tNormalizedDNA\tE2'

        dna_vols = np.array([[2500., 632.5],
                              [3500., 3500.]])

        water_vols = 3500 - dna_vols

        wells = np.array([['A1', 'A2'],
                          ['B1', 'B2']])
        dest_wells = np.array([['D1', 'D2'],
                               ['E1', 'E2']])
        sample_names =  np.array([['sam1', 'sam2'],
                          ['blank1', 'sam3']])

        dna_concs = np.array([[2, 7.89],
                              [np.nan, .0]])

        obs_picklist = format_dna_norm_picklist(dna_vols, water_vols, wells,
                                                dest_wells = dest_wells,
                                                sample_names = sample_names,
                                                dna_concs = dna_concs)

        self.assertEqual(exp_picklist, obs_picklist)

        # test if switching source plates
        exp_picklist = \
        'Sample\tSource Plate Name\tSource Plate Type\tSource Well\t' + \
        'Concentration\tTransfer Volume\tDestination Plate Name\tDestination Well\n' + \
        'sam1\tWater\t384PP_AQ_BP2_HT\tA1\t2.0\t1000.0\tNormalizedDNA\tA1\n' + \
        'sam2\tWater\t384PP_AQ_BP2_HT\tA2\t7.89\t2867.5\tNormalizedDNA\tA2\n' + \
        'blank1\tWater\t384PP_AQ_BP2_HT\tB1\tnan\t0.0\tNormalizedDNA\tB1\n' + \
        'sam3\tWater\t384PP_AQ_BP2_HT\tB2\t0.0\t0.0\tNormalizedDNA\tB2\n' + \
        'sam1\tSample_Plate1\t384PP_AQ_BP2_HT\tA1\t2.0\t2500.0\tNormalizedDNA\tA1\n' + \
        'sam2\tSample_Plate1\t384PP_AQ_BP2_HT\tA2\t7.89\t632.5\tNormalizedDNA\tA2\n' + \
        'blank1\tSample_Plate2\t384PP_AQ_BP2_HT\tB1\tnan\t3500.0\tNormalizedDNA\tB1\n' + \
        'sam3\tSample_Plate2\t384PP_AQ_BP2_HT\tB2\t0.0\t3500.0\tNormalizedDNA\tB2'

        dna_vols = np.array([[2500., 632.5],
                              [3500., 3500.]])

        water_vols = 3500 - dna_vols

        wells = np.array([['A1', 'A2'],
                          ['B1', 'B2']])

        sample_names =  np.array([['sam1', 'sam2'],
                          ['blank1', 'sam3']])

        sample_plates = np.array([['Sample_Plate1', 'Sample_Plate1'],
                          ['Sample_Plate2', 'Sample_Plate2']])

        dna_concs = np.array([[2, 7.89],
                              [np.nan, .0]])

        obs_picklist = format_dna_norm_picklist(dna_vols, water_vols, wells,
                                                sample_names = sample_names,
                                                sample_plates = sample_plates,
                                                dna_concs = dna_concs)

        self.assertEqual(exp_picklist, obs_picklist)

    def test_format_index_picklist(self):
        exp_picklist = \
            'Sample\tSource Plate Name\tSource Plate Type\tSource Well\tTransfer Volume\tIndex Name\t' + \
            'Index Sequence\tIndex Combo\tDestination Plate Name\tDestination Well\n' + \
            'sam1\tiTru5_plate\t384LDV_AQ_B2_HT\tA1\t250\tiTru5_01_A\tACCGACAA\t0\tIndexPCRPlate\tA1\n' + \
            'sam2\tiTru5_plate\t384LDV_AQ_B2_HT\tB1\t250\tiTru5_01_B\tAGTGGCAA\t1\tIndexPCRPlate\tA2\n' + \
            'blank1\tiTru5_plate\t384LDV_AQ_B2_HT\tC1\t250\tiTru5_01_C\tCACAGACT\t2\tIndexPCRPlate\tB1\n' + \
            'sam3\tiTru5_plate\t384LDV_AQ_B2_HT\tD1\t250\tiTru5_01_D\tCGACACTT\t3\tIndexPCRPlate\tB2\n' + \
            'sam1\tiTru7_plate\t384LDV_AQ_B2_HT\tA1\t250\tiTru7_101_01\tACGTTACC\t0\tIndexPCRPlate\tA1\n' + \
            'sam2\tiTru7_plate\t384LDV_AQ_B2_HT\tA2\t250\tiTru7_101_02\tCTGTGTTG\t1\tIndexPCRPlate\tA2\n' + \
            'blank1\tiTru7_plate\t384LDV_AQ_B2_HT\tA3\t250\tiTru7_101_03\tTGAGGTGT\t2\tIndexPCRPlate\tB1\n' + \
            'sam3\tiTru7_plate\t384LDV_AQ_B2_HT\tA4\t250\tiTru7_101_04\tGATCCATG\t3\tIndexPCRPlate\tB2'

        sample_wells = np.array(['A1', 'A2', 'B1', 'B2'])

        sample_names =  np.array(['sam1', 'sam2', 'blank1', 'sam3'])

        indices = pd.DataFrame({'i5 name': {0: 'iTru5_01_A',
                               1: 'iTru5_01_B',
                               2: 'iTru5_01_C',
                               3: 'iTru5_01_D'},
                 'i5 plate': {0: 'iTru5_plate',
                              1: 'iTru5_plate',
                              2: 'iTru5_plate',
                              3: 'iTru5_plate'},
                 'i5 sequence': {0: 'ACCGACAA', 1: 'AGTGGCAA',
                                 2: 'CACAGACT', 3: 'CGACACTT'},
                 'i5 well': {0: 'A1', 1: 'B1', 2: 'C1', 3: 'D1'},
                 'i7 name': {0: 'iTru7_101_01',
                             1: 'iTru7_101_02',
                             2: 'iTru7_101_03',
                             3: 'iTru7_101_04'},
                 'i7 plate': {0: 'iTru7_plate',
                              1: 'iTru7_plate',
                              2: 'iTru7_plate',
                              3: 'iTru7_plate'},
                 'i7 sequence': {0: 'ACGTTACC', 1: 'CTGTGTTG',
                                 2: 'TGAGGTGT', 3: 'GATCCATG'},
                 'i7 well': {0: 'A1', 1: 'A2', 2: 'A3', 3: 'A4'},
                 'index combo': {0: 0, 1: 1, 2: 2, 3: 3},
                 'index combo seq': {0: 'ACCGACAAACGTTACC',
                                     1: 'AGTGGCAACTGTGTTG',
                                     2: 'CACAGACTTGAGGTGT',
                                     3: 'CGACACTTGATCCATG'}})

        obs_picklist = format_index_picklist(sample_names, sample_wells, indices)

        self.assertEqual(exp_picklist, obs_picklist)

    def test_compute_qpcr_concentration(self):
        obs = compute_qpcr_concentration(self.cp_vals)
        exp = self.qpcr_conc

        npt.assert_allclose(obs, exp)

    def test_compute_shotgun_pooling_values_eqvol(self):
        obs_sample_vols = \
            compute_shotgun_pooling_values_eqvol(self.qpcr_conc,
                                                 total_vol=60.0)

        exp_sample_vols = np.zeros([3, 4]) + 60.0/12*1000

        npt.assert_allclose(obs_sample_vols, exp_sample_vols)

    def test_compute_shotgun_pooling_values_eqvol_intvol(self):
        obs_sample_vols = \
            compute_shotgun_pooling_values_eqvol(self.qpcr_conc,
                                                 total_vol=60)

        exp_sample_vols = np.zeros([3, 4]) + 60.0/12*1000

        npt.assert_allclose(obs_sample_vols, exp_sample_vols)

    def test_compute_shotgun_pooling_values_qpcr(self):
        sample_concs = np.array([[1, 12, 400],
                                 [200, 40, 1]])

        exp_vols = np.array([[0, 50000, 6250],
                             [12500, 50000, 0]])

        obs_vols = compute_shotgun_pooling_values_qpcr(sample_concs)

        npt.assert_allclose(exp_vols, obs_vols)

    def test_compute_shotgun_pooling_values_qpcr_minvol(self):
        sample_concs = np.array([[1, 12, 400],
                                 [200, 40, 1]])

        exp_vols = np.array([[100, 100, 4166.6666666666],
                             [8333.33333333333, 41666.666666666, 100]])

        obs_vols = compute_shotgun_pooling_values_qpcr_minvol(sample_concs)

        npt.assert_allclose(exp_vols, obs_vols)

    def test_estimate_pool_conc_vol(self):
        obs_sample_vols = compute_shotgun_pooling_values_eqvol(
                                        self.qpcr_conc, total_vol=60.0)

        obs_pool_conc, obs_pool_vol = estimate_pool_conc_vol(
                                        obs_sample_vols, self.qpcr_conc)

        exp_pool_conc = 323.873027979
        exp_pool_vol = 60000.0

        npt.assert_almost_equal(obs_pool_conc, exp_pool_conc)
        npt.assert_almost_equal(obs_pool_vol, exp_pool_vol)


    def test_format_pooling_echo_pick_list(self):
        vol_sample = np.array([[10.00, 10.00, 5.00, 5.00, 10.00, 10.00]])

        header = ['Source Plate Name,Source Plate Type,Source Well,'
                'Concentration,Transfer Volume,Destination Plate Name,'
                'Destination Well']

        exp_values = ['1,384LDV_AQ_B2_HT,A1,,10.00,NormalizedDNA,A1',
                      '1,384LDV_AQ_B2_HT,A2,,10.00,NormalizedDNA,A1',
                      '1,384LDV_AQ_B2_HT,A3,,5.00,NormalizedDNA,A1',
                      '1,384LDV_AQ_B2_HT,A4,,5.00,NormalizedDNA,A2',
                      '1,384LDV_AQ_B2_HT,A5,,10.00,NormalizedDNA,A2',
                      '1,384LDV_AQ_B2_HT,A6,,10.00,NormalizedDNA,A2']

        exp_str = '\n'.join(header + exp_values)

        obs_str = format_pooling_echo_pick_list(vol_sample,
                                  max_vol_per_well=26,
                                  dest_plate_shape=[16,24])
        self.maxDiff = None
        self.assertEqual(exp_str, obs_str)


    def test_format_pooling_echo_pick_list(self):
        vol_sample = np.array([[10.00, 10.00, np.nan, 5.00, 10.00, 10.00]])

        header = ['Source Plate Name,Source Plate Type,Source Well,'
                'Concentration,Transfer Volume,Destination Plate Name,'
                'Destination Well']

        exp_values = ['1,384LDV_AQ_B2_HT,A1,,10.00,NormalizedDNA,A1',
                      '1,384LDV_AQ_B2_HT,A2,,10.00,NormalizedDNA,A1',
                      '1,384LDV_AQ_B2_HT,A3,,0.00,NormalizedDNA,A1',
                      '1,384LDV_AQ_B2_HT,A4,,5.00,NormalizedDNA,A1',
                      '1,384LDV_AQ_B2_HT,A5,,10.00,NormalizedDNA,A2',
                      '1,384LDV_AQ_B2_HT,A6,,10.00,NormalizedDNA,A2']

        exp_str = '\n'.join(header + exp_values)

        obs_str = format_pooling_echo_pick_list(vol_sample,
                                  max_vol_per_well=26,
                                  dest_plate_shape=[16,24])
        self.maxDiff = None
        self.assertEqual(exp_str, obs_str)
      


    def test_make_2D_array(self):
        example_qpcr_df = pd.DataFrame({'Cp': [12, 0, 5, np.nan],
                                        'Pos': ['A1','A2','A3','A4']})

        exp_cp_array = np.array([[12.0,0.0,5.0,np.nan]])

        np.testing.assert_allclose(make_2D_array(example_qpcr_df, rows=1, cols=4).astype(float), exp_cp_array)


        example2_qpcr_df = pd.DataFrame({'Cp': [12, 0, 1, np.nan,
                                               12, 0, 5, np.nan],
                                        'Pos': ['A1','A2','A3','A4',
                                                'B1','B2','B3','B4']})
        exp2_cp_array = np.array([[12.0,0.0,1.0,np.nan],
                                  [12.0,0.0,5.0,np.nan]])

        np.testing.assert_allclose(make_2D_array(example2_qpcr_df, rows=2, cols=4).astype(float), exp2_cp_array)

    def combine_dfs(self):
        exp_df_f = '''Sample\tWell\tPlate\tCounter\tPrimer_i5\tSource_Well_i5\tIndex_i5\tPrimer_i7\tSource_Well_i7\tIndex_i7\tDNA_concentration\tTransfer_Volume\tCp
        8_29_13_rk_rh\tA1\tABTX_35\t1841.0\tiTru5_01_G\tG1\tGTTCCATG\tiTru7_110_05\tA23\tCGCTTAAC\t12.751753\t80.0\t20.55
        8_29_13_rk_lh\tC1\tABTX_35\t1842.0\tiTru5_01_H\tH1\tTAGCTGAG\tiTru7_110_06\tB23\tCACCACTA\t17.582063\t57.5\t9.15'''

        test_index_picklist_f = '''\tWell Number\tPlate\tSample Name\tSource Plate Name\tSource Plate Type\tCounter\tPrimer\tSource Well\tIndex\tUnnamed: 9\tUnnamed: 10\tUnnamed: 11\tTransfer volume\tDestination Well\tUnnamed: 14
        0\t1\tABTX_35\t8_29_13_rk_rh\ti5 Source Plate\t384LDV_AQ_B2_HT\t1841.0\tiTru5_01_G\tG1\tGTTCCATG\tiTru7_110_05\tA23\tCGCTTAAC\t250\tA1\tNaN
        1\t2\tABTX_35\t8_29_13_rk_lh\ti5 Source Plate\t384LDV_AQ_B2_HT\t1842.0\tiTru5_01_H\tH1\tTAGCTGAG\tiTru7_110_06\tB23\tCACCACTA\t250\tC1\tNaN
        2\t1\tABTX_35\t8_29_13_rk_rh\ti7 Source Plate\t384LDV_AQ_B2_HT\t1841.0\tiTru7_110_05\tA23\tCGCTTAAC\t\t\t\t250\tA1\tNaN
        3\t2\tABTX_35\t8_29_13_rk_lh\ti7 Source Plate\t384LDV_AQ_B2_HT\t1842.0\tiTru7_110_06\tB23\tCACCACTA\t\t\t\t250\tC1\tNaN'''

        test_dna_picklist_f = '''\tSource Plate Name\tSource Plate Type\tSource Well\tConcentration\tTransfer Volume\tDestination Plate Name\tDestination Well
        0\twater\t384LDV_AQ_B2_HT\tA1\tNaN\t3420.0\tNormalizedDNA\tA1
        1\twater\t384LDV_AQ_B2_HT\tC1\tNaN\t3442.5\tNormalizedDNA\tC1
        5\t1\t384LDV_AQ_B2_HT\tA1\t12.751753\t80.0\tNormalizedDNA\tA1
        6\t1\t384LDV_AQ_B2_HT\tC1\t17.582063\t57.5\tNormalizedDNA\tC1'''

        test_qpcr_f = '''\tInclude\tColor\tPos\tName\tCp\tConcentration\tStandard\tStatus
        0\tTRUE\t255\tA1\tSample 1\t20.55\tNaN\t0\tNaN
        1\tTRUE\t255\tC1\tSample 2\t9.15\tNaN\t0\tNaN'''

        exp_out_f = '''Well\tCp\tDNA Concentration\tDNA Transfer Volume\tSample Name\tPlate\tCounter\tPrimer i7\tSource Well i7\tIndex i7\tPrimer i5\tSource Well i5\tIndex i5
        A1\t20.55\t12.751753\t80.0\t8_29_13_rk_rh\tABTX_35\t1841.0\tiTru7_110_05\tA23\tCGCTTAAC\tiTru5_01_G\tG1\tGTTCCATG
        C1\t9.15\t17.582063\t57.5\t8_29_13_rk_lh\tABTX_35\t1842.0\tiTru7_110_06\tB23\tCACCACTA\tiTru5_01_H\tH1\tTAGCTGAG'''

        test_index_picklist_df = pd.read_csv(StringIO(test_index_picklist_f), header=0, sep='\t')
        test_dna_picklist_df = pd.read_csv(StringIO(test_dna_picklist_f), header=0, sep='\t')
        test_qpcr_df = pd.read_csv(StringIO(test_qpcr_f), header=0, sep='\t')

        exp_df = pd.read_csv(StringIO(exp_out_f), header=0, sep='\t')

        combined_df = combine_dfs(test_qpcr_df, test_dna_picklist_df, test_index_picklist_df)

        pd.testing.assert_frame_equal(combined_df, exp_df, check_like=True)

    def test_add_dna_conc(self):
        test_dna = '''Well\tpico_conc
        A1\t2.5
        C1\t20'''

        test_combined = '''Well\tCp\tDNA Concentration\tDNA Transfer Volume\tSample Name\tPlate\tCounter\tPrimer i7\tSource Well i7\tIndex i7\tPrimer i5\tSource Well i5\tIndex i5
        A1\t20.55\t12.751753\t80.0\t8_29_13_rk_rh\tABTX_35\t1841.0\tiTru7_110_05\tA23\tCGCTTAAC\tiTru5_01_G\tG1\tGTTCCATG
        C1\t9.15\t17.582063\t57.5\t8_29_13_rk_lh\tABTX_35\t1842.0\tiTru7_110_06\tB23\tCACCACTA\tiTru5_01_H\tH1\tTAGCTGAG'''

        test_exp_out = '''Well\tCp\tDNA Concentration\tDNA Transfer Volume\tSample Name\tPlate\tCounter\tPrimer i7\tSource Well i7\tIndex i7\tPrimer i5\tSource Well i5\tIndex i5\tpico_conc
        A1\t20.55\t12.751753\t80.0\t8_29_13_rk_rh\tABTX_35\t1841.0\tiTru7_110_05\tA23\tCGCTTAAC\tiTru5_01_G\tG1\tGTTCCATG\t2.5
        C1\t9.15\t17.582063\t57.5\t8_29_13_rk_lh\tABTX_35\t1842.0\tiTru7_110_06\tB23\tCACCACTA\tiTru5_01_H\tH1\tTAGCTGAG\t20'''

        exp_df = pd.read_csv(StringIO(test_exp_out), header=0, sep='\t')
        test_in_df = pd.read_csv(StringIO(test_combined), header=0, sep='\t')
        test_dna_df = pd.read_csv(StringIO(test_dna), header=0, sep='\t')

        obs_df = add_dna_conc(test_in_df, test_dna_df)

        pd.testing.assert_frame_equal(obs_df, exp_df, check_like=True)

    def test_compute_pico_concentration(self):
        obs = compute_pico_concentration(self.dna_vals)
        exp = self.pico_conc

        npt.assert_allclose(obs, exp)

    def test_format_sheet_comments(self):
        contacts = {'Jeff Dereus': 'jdereus@ucsd.edu',
                    'Gail Ackermann': 'ackermag@ucsd.edu',
                    'Jon Sanders': 'jonsan@gmail.com',
                    'Greg Humphrey': 'ghsmu414@gmail.com'}

        PI = {'Knight': 'robknight@ucsd.edu'}

        other = None

        sep = '\t'

        exp_comment = (
            'PI\tKnight\trobknight@ucsd.edu\n'
            'Contact\tGail Ackermann\tGreg Humphrey'
            '\tJeff Dereus\tJon Sanders\n'
            '\tackermag@ucsd.edu\tghsmu414@gmail.com'
            '\tjdereus@ucsd.edu\tjonsan@gmail.com\n'
            )

        obs_comment = format_sheet_comments(PI, contacts, other, sep)

        self.assertEqual(exp_comment, obs_comment)

    def test_format_sample_sheet(self):

        self.maxDiff = None

        exp_sample_sheet = (
            '[Header]\n'
            'IEMFileVersion\t4\n'
            'Investigator Name\tKnight\n'
            'Experiment Name\t\n'
            'Date\t2017-08-13\n'
            'Workflow\tGenerateFASTQ\n'
            'Application\tFASTQ Only\n'
            'Assay\tMetagenomics\n'
            'Description\t\n'
            'Chemistry\tDefault\n\n'
            '[Reads]\n'
            '150\n'
            '150\n\n'
            '[Settings]\n'
            'ReverseComplement\t0\n\n'
            '[Data]\n'
            'Lane\tSample_ID\tSample_Name\tSample_Plate\tSample_Well'
            '\tI7_Index_ID\tindex\tI5_Index_ID\tindex2\tSample_Project'
            '\tDescription\n'
            '1\tsam1\tsam1\texample\tA1\tiTru7_101_01\tACGTTACC\tiTru5_01_A'
            '\tACCGACAA\texample_proj\t\n'
            '1\tsam2\tsam2\texample\tA2\tiTru7_101_02\tCTGTGTTG\tiTru5_01_B'
            '\tAGTGGCAA\texample_proj\t\n'
            '1\tblank1\tblank1\texample\tB1\tiTru7_101_03\tTGAGGTGT\tiTru5_01_C'
            '\tCACAGACT\texample_proj\t\n'
            '1\tsam3\tsam3\texample\tB2\tiTru7_101_04\tGATCCATG\tiTru5_01_D'
            '\tCGACACTT\texample_proj\t'
            )


        exp_sample_sheet_2 = (
            '# PI\tKnight\trobknight@ucsd.edu\t\t\n'
            '# Contact\tJeff Dereus\tGail Ackermann\tJon Sanders\tGreg Humphrey\n'
            '# \tjdereus@ucsd.edu\tackermag@ucsd.edu\tjonsan@gmail.com\tghsmu414@gmail.com\n'
            '[Header]\n'
            'IEMFileVersion\t4\n'
            'Investigator Name\tKnight\n'
            'Experiment Name\t\n'
            'Date\t2017-08-13\n'
            'Workflow\tGenerateFASTQ\n'
            'Application\tFASTQ Only\n'
            'Assay\tMetagenomics\n'
            'Description\t\n'
            'Chemistry\tDefault\n\n'
            '[Reads]\n'
            '150\n'
            '150\n\n'
            '[Settings]\n'
            'ReverseComplement\t0\n\n'
            '[Data]\n'
            'Lane\tSample_ID\tSample_Name\tSample_Plate\t'
            'Sample_Well\tI7_Index_ID\tindex\tI5_Index_ID\t'
            'index2\tSample_Project\tDescription\n'
            '1\tsam1\tsam1\texample\tA1\tiTru7_101_01\tACGTTACC'
            '\tiTru5_01_A\tACCGACAA\texample_proj\t\n'
            '1\tsam2\tsam2\texample\tA2\tiTru7_101_02\tCTGTGTTG'
            '\tiTru5_01_B\tAGTGGCAA\texample_proj\t\n'
            '1\tblank1\tblank1\texample\tB1\tiTru7_101_03\tTGAGGTGT'
            '\tiTru5_01_C\tCACAGACT\texample_proj\t\n'
            '1\tsam3\tsam3\texample\tB2\tiTru7_101_04\tGATCCATG'
            '\tiTru5_01_D\tCGACACTT\texample_proj\t'
            )

        comment = (
            'PI\tKnight\trobknight@ucsd.edu\t\t\n'
            'Contact\tJeff Dereus\tGail Ackermann\t'
            'Jon Sanders\tGreg Humphrey\n'
            '\tjdereus@ucsd.edu\tackermag@ucsd.edu\t'
            'jonsan@gmail.com\tghsmu414@gmail.com\n'
            )


        data = (
            'Lane\tSample_ID\tSample_Name\tSample_Plate\tSample_Well\t'
            'I7_Index_ID\tindex\tI5_Index_ID\tindex2\tSample_Project\t'
            'Description\n'
            '1\tsam1\tsam1\texample\tA1\tiTru7_101_01\tACGTTACC\t'
            'iTru5_01_A\tACCGACAA\texample_proj\t\n'
            '1\tsam2\tsam2\texample\tA2\tiTru7_101_02\tCTGTGTTG\t'
            'iTru5_01_B\tAGTGGCAA\texample_proj\t\n'
            '1\tblank1\tblank1\texample\tB1\tiTru7_101_03\tTGAGGTGT\t'
            'iTru5_01_C\tCACAGACT\texample_proj\t\n'
            '1\tsam3\tsam3\texample\tB2\tiTru7_101_04\tGATCCATG\t'
            'iTru5_01_D\tCGACACTT\texample_proj\t'
            )

        sample_sheet_dict = {'comments': '',
                  'IEMFileVersion': '4',
                  'Investigator Name': 'Knight',
                  'Experiment Name': '',
                  'Date': '2017-08-13',
                  'Workflow': 'GenerateFASTQ',
                  'Application': 'FASTQ Only',
                  'Assay': 'Metagenomics',
                  'Description': '',
                  'Chemistry': 'Default',
                  'read1': 150,
                  'read2': 150,
                  'ReverseComplement': '0',
                  'data': data}

        obs_sample_sheet = format_sample_sheet(sample_sheet_dict, sep='\t')

        self.assertEqual(exp_sample_sheet, obs_sample_sheet)

        sample_sheet_dict_2 = {'comments': comment,
                  'IEMFileVersion': '4',
                  'Investigator Name': 'Knight',
                  'Experiment Name': '',
                  'Date': '2017-08-13',
                  'Workflow': 'GenerateFASTQ',
                  'Application': 'FASTQ Only',
                  'Assay': 'Metagenomics',
                  'Description': '',
                  'Chemistry': 'Default',
                  'read1': 150,
                  'read2': 150,
                  'ReverseComplement': '0',
                  'data': data}

        obs_sample_sheet_2 = format_sample_sheet(sample_sheet_dict_2, sep='\t')

        self.assertEqual(exp_sample_sheet_2, obs_sample_sheet_2)

    def test_bcl_scrub_name(self):
        self.assertEqual('test_1', bcl_scrub_name('test.1'))
        self.assertEqual('test-1', bcl_scrub_name('test-1'))
        self.assertEqual('test_1', bcl_scrub_name('test_1'))

    def test_rc(self):
        self.assertEqual(rc('AGCCT'), 'AGGCT')

    def test_sequencer_i5_index(self):
        indices = ['AGCT','CGGA','TGCC']

        exp_rc = ['AGCT','TCCG','GGCA']

        obs_hiseq4k = sequencer_i5_index('HiSeq4000', indices)
        obs_hiseq25k = sequencer_i5_index('HiSeq2500', indices)
        obs_nextseq = sequencer_i5_index('NextSeq', indices)

        self.assertListEqual(obs_hiseq4k, exp_rc)
        self.assertListEqual(obs_hiseq25k, indices)
        self.assertListEqual(obs_nextseq, exp_rc)

        with self.assertRaises(ValueError):
            sequencer_i5_index('foo', indices)

    def test_format_sample_data(self):
        # test that single lane works
        exp_data = (
            'Lane,Sample_ID,Sample_Name,Sample_Plate'
            ',Sample_Well,I7_Index_ID,index,I5_Index_ID'
            ',index2,Sample_Project,Description\n'
            '1,sam1,sam1,example,A1,iTru7_101_01,ACGTTACC,'
            'iTru5_01_A,ACCGACAA,example_proj,\n'
            '1,sam2,sam2,example,A2,iTru7_101_02,CTGTGTTG,'
            'iTru5_01_B,AGTGGCAA,example_proj,\n'
            '1,blank1,blank1,example,B1,iTru7_101_03,TGAGGTGT,'
            'iTru5_01_C,CACAGACT,example_proj,\n'
            '1,sam3,sam3,example,B2,iTru7_101_04,GATCCATG,'
            'iTru5_01_D,CGACACTT,example_proj,'
            )

        wells = ['A1', 'A2', 'B1', 'B2']
        sample_ids =  ['sam1', 'sam2', 'blank1', 'sam3']
        i5_name = ['iTru5_01_A', 'iTru5_01_B', 'iTru5_01_C', 'iTru5_01_D']
        i5_seq = ['ACCGACAA', 'AGTGGCAA', 'CACAGACT', 'CGACACTT']
        i7_name = ['iTru7_101_01', 'iTru7_101_02',
                   'iTru7_101_03', 'iTru7_101_04']
        i7_seq =['ACGTTACC', 'CTGTGTTG', 'TGAGGTGT', 'GATCCATG']

        obs_data = format_sample_data(sample_ids,i7_name, i7_seq,
                                      i5_name, i5_seq, wells=wells,
                                      sample_plate='example',
                                      sample_proj='example_proj',
                                      lanes=[1])

        self.assertEqual(obs_data, exp_data)

        # test that two lanes works
        exp_data_2 = (
            'Lane,Sample_ID,Sample_Name,Sample_Plate,'
            'Sample_Well,I7_Index_ID,index,I5_Index_ID,'
            'index2,Sample_Project,Description\n'
            '1,sam1,sam1,example,A1,iTru7_101_01,ACGTTACC,'
            'iTru5_01_A,ACCGACAA,example_proj,\n'
            '1,sam2,sam2,example,A2,iTru7_101_02,CTGTGTTG,'
            'iTru5_01_B,AGTGGCAA,example_proj,\n'
            '1,blank1,blank1,example,B1,iTru7_101_03,TGAGGTGT,'
            'iTru5_01_C,CACAGACT,example_proj,\n'
            '1,sam3,sam3,example,B2,iTru7_101_04,GATCCATG,'
            'iTru5_01_D,CGACACTT,example_proj,\n'
            '2,sam1,sam1,example,A1,iTru7_101_01,ACGTTACC,'
            'iTru5_01_A,ACCGACAA,example_proj,\n'
            '2,sam2,sam2,example,A2,iTru7_101_02,CTGTGTTG,'
            'iTru5_01_B,AGTGGCAA,example_proj,\n'
            '2,blank1,blank1,example,B1,iTru7_101_03,TGAGGTGT'
            ',iTru5_01_C,CACAGACT,example_proj,\n'
            '2,sam3,sam3,example,B2,iTru7_101_04,GATCCATG'
            ',iTru5_01_D,CGACACTT,example_proj,'
            )

        obs_data_2 = format_sample_data(sample_ids,i7_name, i7_seq,
                                        i5_name, i5_seq, wells=wells,
                                        sample_plate='example',
                                        sample_proj='example_proj',
                                        lanes=[1,2])

        self.assertEqual(obs_data_2, exp_data_2)

        # test with r/c i5 barcodes
        exp_data = (
            'Lane,Sample_ID,Sample_Name,Sample_Plate'
            ',Sample_Well,I7_Index_ID,index,I5_Index_ID'
            ',index2,Sample_Project,Description\n'
            '1,sam1,sam1,example,A1,iTru7_101_01,ACGTTACC,'
            'iTru5_01_A,ACCGACAA,example_proj,\n'
            '1,sam2,sam2,example,A2,iTru7_101_02,CTGTGTTG,'
            'iTru5_01_B,AGTGGCAA,example_proj,\n'
            '1,blank1,blank1,example,B1,iTru7_101_03,TGAGGTGT,'
            'iTru5_01_C,CACAGACT,example_proj,\n'
            '1,sam3,sam3,example,B2,iTru7_101_04,GATCCATG,'
            'iTru5_01_D,CGACACTT,example_proj,'
            )

        i5_seq = ['ACCGACAA', 'AGTGGCAA', 'CACAGACT', 'CGACACTT']

        obs_data = format_sample_data(sample_ids,i7_name, i7_seq,
                                      i5_name, i5_seq, wells=wells,
                                      sample_plate='example',
                                      sample_proj='example_proj',
                                      lanes=[1])

        self.assertEqual(obs_data, exp_data)

    def test_reformat_interleaved_to_columns(self):
        wells = ['A1','A23','C1','C23',
                 'A2','A24','C2','C24',
                 'B1','B23','D1','D23',
                 'B2','B24','D2','D24']

        exp = ['A1','B6','C1','D6',
               'A7','B12','C7','D12',
               'A13','B18','C13','D18',
               'A19','B24','C19','D24']

        obs = reformat_interleaved_to_columns(wells)
            
        np.testing.assert_array_equal(exp, obs)

    def test_parse_sample_sheet(self):
        sheet = parse_sample_sheet(self.good_ss)

        exp = ('# PI,Knight,robknight@ucsd.edu,,,,,,,,\n'
               '# Contact,Gail Ackermann,Greg Humphrey,Jeff Dereus,'
               'Jon Sanders,,,,,,\n'
               '# ,ackermag@ucsd.edu,ghsmu414@gmail.com,jdereus@ucsd.edu'
               ',jonsan@gmail.com,,,,,,\n')
        self.assertEqual(sheet.comments, exp)

        exp = {'IEMFileVersion': '4',
               'Investigator Name': 'Knight',
               'Experiment Name': 'RKL0042',
               'Date': '2020-02-26',
               'Workflow': 'GenerateFASTQ',
               'Application': 'FASTQ Only',
               'Assay': 'Metagenomics',
               'Description': '',
               'Chemistry': 'Default'}
        self.assertEqual(sheet.Header, exp)
        self.assertEqual(sheet.Reads, [150, 150])
        self.assertEqual(sheet.Settings, {'ReverseComplement': '0'})

        # check that we are parsing the right number of samples
        self.assertEqual(len(sheet.samples), 783)

    def test_parse_sample_sheet_no_comments(self):
        sheet = parse_sample_sheet(self.no_comments_ss)

        self.assertEqual(sheet.comments, '')

        exp = {'IEMFileVersion': '4',
               'Investigator Name': 'Knight',
               'Experiment Name': 'RKL0042',
               'Date': '2020-02-26',
               'Workflow': 'GenerateFASTQ',
               'Application': 'FASTQ Only',
               'Assay': 'Metagenomics',
               'Description': '',
               'Chemistry': 'Default'}
        self.assertEqual(sheet.Header, exp)
        self.assertEqual(sheet.Reads, [150, 150])
        self.assertEqual(sheet.Settings, {'ReverseComplement': '0'})

        # check that we are parsing the right number of samples
        self.assertEqual(len(sheet.samples), 783)

    def test_validate_sample_sheet(self):
        sheet = parse_sample_sheet(self.good_ss)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            validated_sheet = validate_sample_sheet(sheet)

            # check no warnings are raised
            self.assertTrue(len(w) == 0)

    def test_validate_sample_sheet_no_comments(self):
        sheet = parse_sample_sheet(self.no_comments_ss)

        with self.assertRaisesRegex(ValueError, 'This sample sheet does not'):
            validate_sample_sheet(sheet)

    def test_validate_sample_sheet_no_sample_project(self):
        sheet = parse_sample_sheet(self.no_project_ss)

        with self.assertRaisesRegex(ValueError, 'The Sample_project column '):
            validate_sample_sheet(sheet)

    def test_validate_sample_sheet_repeated_sample_ids(self):
        sheet = parse_sample_sheet(self.ok_ss)

        # the sample identifiers are only repeated after scrubbing
        with self.assertRaisesRegex(ValueError, 'After scrubbing samples for '
                                    'bcl2fastq compatibility there are '
                                    'repeated identifiers. The following names'
                                    ' in the Sample_ID column are listed '
                                    'multiple times:\nCDPH-SAL_Salmonella_'
                                    'Typhi_MDL_144: 2'):
            validate_sample_sheet(sheet)

    def test_validate_sample_sheet_scrubbed_names(self):
        sheet = parse_sample_sheet(self.scrubbable_ss)

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
        sheet = parse_sample_sheet(self.bad_project_name_ss)

        with self.assertWarnsRegex(UserWarning, 'The following project names '
                                   'in the Sample_project column are missing a'
                                   ' Qiita study identifier: Feist, Gerwick'):
            validate_sample_sheet(sheet)

    def test_write_sample_sheet(self):
        sheet = parse_sample_sheet(self.good_ss)

        with tempfile.NamedTemporaryFile('w+') as f:
            write_sample_sheet(sheet, f.name)

            f.seek(0)

            contents = f.read()

            # spot check the contents
            self.assertEqual(len(contents.split('\n')), 807)
            print(contents.split('\n')[0])
            print(contents.split('\n')[-1])

            self.assertTrue('# PI,Knight,robknight@ucsd.edu' in contents)
            self.assertTrue('P21_E_coli_ELI352' in contents)

    def test_write_sample_sheet_bad_input(self):
        with self.assertRaisesRegex(ValueError, 'The input sample sheet should'
                                    ' be a SampleSheet instance'):
            write_sample_sheet({}, 'file.txt')


if __name__ == "__main__":
    main()
