from unittest import TestCase

import os
import pandas as pd
import numpy as np
import numpy.testing as npt

from metapool.metapool import *


class Tests(TestCase):
    def setUp(self):
        self.cp_vals = np.array([[10.14, 7.89, 7.9, 15.48],
                                 [7.86, 8.07, 8.16, 9.64],
                                 [12.29, 7.64, 7.32, 13.74]])

        self.qpcr_conc = \
            np.array([[98.14626462, 487.8121413, 484.3480866, 2.183406934],
                      [498.3536649, 429.0839787, 402.4270321, 140.1601735],
                      [21.20533391, 582.9456031, 732.2655041, 7.545145988]])

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

    def test_estimate_pool_conc_vol(self):
        obs_sample_vols = compute_shotgun_pooling_values_eqvol(
                                        self.qpcr_conc, total_vol=60.0)

        obs_pool_conc, obs_pool_vol = estimate_pool_conc_vol(
                                        obs_sample_vols, self.qpcr_conc)

        exp_pool_conc = 323.873027979
        exp_pool_vol = 60000.0

        npt.assert_almost_equal(obs_pool_conc, exp_pool_conc)
        npt.assert_almost_equal(obs_pool_vol, exp_pool_vol)

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


if __name__ == "__main__":
    main()