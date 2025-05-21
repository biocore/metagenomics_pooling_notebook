# test_input_unittest.py
import unittest
import papermill as pm
import pandas as pd
import json
import tempfile
from pathlib import Path

NOTEBOOK = "../tellseq_D_variable_volume_pooling.ipynb"
TEST_DICT = 'test_dict'

PARAM_SETS = [
    ({
        TEST_DICT: {
                "plate_df_set_fp": '../test_output/QC/Tellseq_plate_df_C_set_col19to24.txt',
                "read_counts_fps": ['../test_data/Demux/Tellseq_fastqc_sequence_counts.tsv'],
                'dynamic_range': 5,
                'iseqnormed_picklist_fbase': '../test_output/Pooling/Tellseq_iSeqnormpool'
        }
        },
        "test1"),
    # ({
    #     TEST_DICT: {
    #             "plate_df_set_fp": './test_output/QC/Tellseq_plate_df_C_set_col19to24.txt',
    #             "read_counts_fps": ['./test_data/Demux/Tellseq_fastqc_sequence_counts.tsv'],
    #             'dynamic_range': 5,
    #             'iseqnormed_picklist_fbase': './test_output/Pooling/Tellseq_iSeqnormpool'
    #     },
    #     DATAFRAME_OUTPUT_NAME: "run1_results.csv"
    #     },
    #     "test1"),
]


class TestMyAnalysis(unittest.TestCase):
    def test_my_analysis(self):
        """
        Executes the notebook with two different parameter sets and checks:
          1) the resulting CSV matches expected DataFrame
          2) the JSON of changed variables can be loaded
        """
        for params, expected_fn in PARAM_SETS:
            with self.subTest(test_case=expected_fn):
                # create isolated temp directories
                with tempfile.TemporaryDirectory() as tmp_dir:
                    tmp_path = Path(tmp_dir)
                    tmp_nb_dir = tmp_path / "nb"
                    tmp_csv_dir = tmp_path / "csv"
                    tmp_nb_dir.mkdir()
                    tmp_csv_dir.mkdir()

                    # paths for CSV and dict output
                    dict_output = tmp_nb_dir / f"dict{expected_fn}.txt"

                    run_params = {
                        **params,
                        "test_notebook_output_csv": str(dict_output),
                    }

                    out_nb = tmp_nb_dir / f"{expected_fn}.ipynb"
                    pm.execute_notebook(
                        input_path=NOTEBOOK,
                        output_path=str(out_nb),
                        parameters=run_params,
                        log_output=True,
                    )

                    # Test 1. Check the variables are outputed as expected
                    # ! Currently not actually a test
                    with open(dict_output, 'r') as f:
                        changed_variables = json.load(f)
                        print(changed_variables)


if __name__ == "__main__":
    unittest.main()
