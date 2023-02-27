import os
import pandas as pd
import unittest
import numpy as np

from diff_dataframes import dataframes_are_equal
from metapool.rescale_counts import convert_read_count_to_cell_count, to_absolute_abundance_read_count, to_absolute_abundance_cell_count, _to_row_percentage, \
    _to_column_percentage


class TestRescaleCounts(unittest.TestCase):

    def test_read_count_to_cell_count(self):
        # Load required data
        table_community = pd.read_csv("data/spike_in/table_community_hits.txt", sep="\t", index_col="OTUID")
        metadata_features = pd.read_csv("data/spike_in/metadata_features.tsv", sep="\t", index_col="OTUID")

        df1 = convert_read_count_to_cell_count(table_community, metadata_features)
        spot_check_gotus = ["G000006785", "G900156305", "G000006725", "G900156885"]

        for i in range(len(spot_check_gotus)):
            for j in range(i+1, len(spot_check_gotus)):
                gotu1 = spot_check_gotus[i]
                gotu2 = spot_check_gotus[j]
                len1 = metadata_features.loc[gotu1, "total_length"]
                len2 = metadata_features.loc[gotu2, "total_length"]
                for sample in table_community.columns:
                    start_ratio = table_community.loc[gotu1, sample] / table_community.loc[gotu2, sample]
                    end_ratio = df1.loc[gotu1, sample] / df1.loc[gotu2, sample]
                    if not np.isnan(start_ratio):
                        self.assertAlmostEqual(end_ratio, start_ratio * len2/len1)

    def test_abundance_read_count(self):
        # Load required data
        table_community = pd.read_csv("data/spike_in/table_community_hits.txt", sep="\t", index_col="OTUID")
        linear_models = pd.read_csv("data/spike_in/results_synDNA_linear_model_variables.txt", sep="\t")
        metadata_features = pd.read_csv("data/spike_in/metadata_features.tsv", sep="\t", index_col="OTUID")

        # Clean sample names (Trims everything after last underscore)
        linear_models["sample_name"] = linear_models["sample_name_pool"].replace("_[^_]+$", "", regex=True)

        df1, failed_cols = to_absolute_abundance_read_count(table_community, linear_models)
        df1_cell_count = convert_read_count_to_cell_count(df1, metadata_features)
        df2, failed_cols = to_absolute_abundance_cell_count(table_community, linear_models, metadata_features)

        self.assertTrue(dataframes_are_equal(df1_cell_count, df2, verbose=False))

    def test_normalizations(self):
        df = np.random.rand(10, 10)
        df = pd.DataFrame(df)
        df_col_pct = _to_column_percentage(df)
        self.assertTrue(((df_col_pct.sum() > 99.99) & (df_col_pct.sum() < 100.01)).all())

        df_row_pct = _to_row_percentage(df)
        self.assertTrue(((df_row_pct.sum(axis=1) > 99.99) & (df_row_pct.sum(axis=1) < 100.01)).all())

    def test_same_as_R(self):
        # Load required data
        table_community = pd.read_csv("data/spike_in/table_community_hits.txt", sep="\t", index_col="OTUID")
        linear_models = pd.read_csv("data/spike_in/results_synDNA_linear_model_variables.txt", sep="\t")
        metadata_features = pd.read_csv("data/spike_in/metadata_features.tsv", sep="\t", index_col="OTUID")

        # Clean sample names (Trims everything after last underscore)
        linear_models["sample_name"] = linear_models["sample_name_pool"].replace("_[^_]+$", "", regex=True)

        df, failed_samples = to_absolute_abundance_cell_count(table_community, linear_models, metadata_features)
        df = df.T  # R version applies absolute abundance transform, cell count transform, and a transpose, Python version you must externally transpose.

        print("WARNING: The following samples don't have applicable linear models:", failed_samples)
        # If you want to check with external tools:
        # df.to_csv("feature_table_cell_counts.txt", sep="\t")
        df_R = pd.read_csv("data/spike_in/feature_table_cell_counts_R.txt", sep="\t", index_col="sample_name")
        self.assertTrue(dataframes_are_equal(df, df_R, verbose=False))

        df_cell_normalized = _to_row_percentage(df)
        # If you want to check with external tools:
        # df_cell_normalized.to_csv("feature_table_cell_counts_normalized.txt", sep="\t")
        df_cell_normalized_R = pd.read_csv("data/spike_in/feature_table_cell_counts_normalized_R.txt", sep="\t", index_col="sample_name")
        self.assertTrue(dataframes_are_equal(df_cell_normalized, df_cell_normalized_R, verbose=False))

        # TODO FIXME HACK:  We would like to normalize cell counts to cell per gram
        #  use pseudocounts (and discuss why some of the files have zero mass!?)