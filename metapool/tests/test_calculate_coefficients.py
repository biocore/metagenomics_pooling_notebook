import unittest
import pandas as pd
from metapool.calculate_coefficients import calculate_coefficients
from diff_dataframes import dataframes_are_equal


class TestRescaleCounts(unittest.TestCase):
    def test_same_as_R(self):
        folder = "metapool/tests/data/spike_in/"
        # Import data
        table_synthetic_hits = pd.read_csv(
            folder + "table_plasmid_sequence_hits.txt",
            sep="\t", index_col="OTUID")
        metadata_pools = pd.read_csv(
            folder + "metadata_samples_plasmid_sequences.txt",
            sep="\t")

        coefall = calculate_coefficients(table_synthetic_hits, metadata_pools,
                                         plot_fit=False)
        coefall = coefall.set_index("sample_name")

        coefall_r = pd.read_csv(
            folder + "results_synDNA_linear_model_variables_R.txt",
            index_col="sample_name", sep="\t")

        self.assertTrue(dataframes_are_equal(coefall, coefall_r, verbose=True))


if __name__ == '__main__':
    unittest.main()
