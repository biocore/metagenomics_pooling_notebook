import os

from unittest import TestCase, main
from metapool.metapool import parse_sample_sheet
from metapool.prep import (preparations_for_run, remove_qiita_id,
                           get_run_prefix, is_nonempty_gz_file,
                           sample_sheet_to_dataframe, parse_illumina_run_id)


class Tests(TestCase):
    def setUp(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')

        # "valid" upfront but will have repeated values after scrubbing
        ok_ss = os.path.join(data_dir, 'ok-sample-sheet.csv')
        self.ss = parse_sample_sheet(ok_ss)

    def test_preparations_for_run(self):
        pass

    def test_remove_qiita_id(self):
        obs = remove_qiita_id('project_1')
        self.assertEqual(obs, 'project')

        obs = remove_qiita_id('project_00333333')
        self.assertEqual(obs, 'project')

        obs = remove_qiita_id('project')
        self.assertEqual(obs, 'project')

    def test_get_run_prefix(self):
        pass

    def test_is_non_empty_gz_file(self):
        pass

    def test_sample_sheet_to_dataframe(self):
        pass

    def test_parse_illumina_run_id(self):
        pass


if __name__ == "__main__":
    main()
