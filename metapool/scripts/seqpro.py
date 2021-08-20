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
@click.option('--pipeline', help='Which pipeline generated the data',
              show_default=True, default='fastp-and-minimap2',
              type=click.Choice(['atropos-and-bowtie2', 'fastp-and-minimap2']))
def format_preparation_files(run_dir, sample_sheet, output_dir, pipeline):
    """Generate the preparation files for the projects in a run

    RUN_DIR: should be the directory where the results of running bcl2fastq are
    saved.

    SAMPLE_SHEET: should be a CSV file that includes information for the
    samples and projects in RUN_DIR.

    OUTPUT_DIR: directory where the outputted preparations should be saved to.

    Preparations are stratified by project and by lane. Only samples with
    non-empty files are included.
    """

    sample_sheet = sample_sheet_to_dataframe(parse_sample_sheet(sample_sheet))

    # returns a map of project_name.lane -> preparation frame
    preps = preparations_for_run(run_dir, sample_sheet, pipeline=pipeline)

    os.makedirs(output_dir, exist_ok=True)

    for filename, df in preps.items():
        filename = os.path.join(output_dir, filename) + '.tsv'

        df.to_csv(filename, sep='\t', index=False)


if __name__ == '__main__':
    format_preparation_files()
