#!/usr/bin/env python

import click
import os

from metapool import (preparations_for_run, parse_sample_sheet,
                      sample_sheet_to_dataframe)


@click.command()
@click.argument('run_dir', type=click.Path(exists=True, dir_okay=True,
                                           file_okay=False))
@click.argument('sample_sheet', type=click.Path(exists=True, dir_okay=False,
                                                file_okay=True))
@click.argument('output_dir', type=click.Path(writable=True))
def format_preparation_files(run_dir, sample_sheet, output_dir):
    """Generate the preparation files for the projects in a run

    Preparations are stratified by project and by lane. Only samples with
    non-empty files are included.
    """

    sample_sheet = sample_sheet_to_dataframe(parse_sample_sheet(sample_sheet))

    # returns a map of project_name.lane -> preparation frame
    preps = preparations_for_run(run_dir, sample_sheet)

    os.makedirs(output_dir, exist_ok=True)

    for filename, df in preps.items():
        # qiita requires txt files no tsvs allowed
        filename = os.path.join(output_dir, filename) + '.txt'

        df.to_csv(filename, sep='\t', index=False)


if __name__ == '__main__':
    format_preparation_files()
