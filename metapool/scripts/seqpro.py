#!/usr/bin/env python

import click
import os
import re
from os.path import abspath

from metapool import (preparations_for_run, load_sample_sheet,
                      sample_sheet_to_dataframe, run_counts,
                      get_qiita_id_from_project_name)


@click.command()
@click.argument('run_dir', type=click.Path(exists=True, dir_okay=True,
                                           file_okay=False))
@click.argument('sample_sheet', type=click.Path(exists=True, dir_okay=False,
                                                file_okay=True))
@click.argument('output_dir', type=click.Path(writable=True))
@click.option('--verbose', help='list prep-file output paths, study_ids',
              is_flag=True)
def format_preparation_files(run_dir, sample_sheet, output_dir, verbose):
    """Generate the preparation files for the projects in a run

    RUN_DIR: should be the directory where the results of running bcl2fastq are
    saved.

    SAMPLE_SHEET: should be a CSV file that includes information for the
    samples and projects in RUN_DIR.

    OUTPUT_DIR: directory where the outputted preparations should be saved to.

    Preparations are stratified by project and by lane. Only samples with
    non-empty files are included. If "fastp-and-minimap2" is used, the script
    will collect sequence count stats for each sample and add them as columns
    in the preparation file.
    """
    sample_sheet = load_sample_sheet(sample_sheet)
    df_sheet = sample_sheet_to_dataframe(sample_sheet)

    stats = run_counts(run_dir, sample_sheet)
    stats['sample_name'] = \
        df_sheet.set_index('lane', append=True)['sample_name']

    # sample_sheet_to_dataframe() automatically lowercases the column names
    # before returning df_sheet. Hence, sample_sheet.CARRIED_PREP_COLUMNS also
    # needs to be lowercased for the purposes of tests in
    # preparation_for_run().
    c_prep_columns = [x.lower() for x in sample_sheet.CARRIED_PREP_COLUMNS]
    # returns a map of (run, project_name, lane) -> preparation frame
    preps = preparations_for_run(run_dir,
                                 df_sheet,
                                 sample_sheet.GENERATED_PREP_COLUMNS,
                                 c_prep_columns)

    os.makedirs(output_dir, exist_ok=True)

    for (run, project, lane), df in preps.items():
        fp = os.path.join(output_dir, f'{run}.{project}.{lane}.tsv')

        # stats are indexed by sample name and lane, lane is the first
        # level index. When merging, make sure to select the lane subset
        # that we care about, otherwise we'll end up with repeated rows
        df = df.merge(stats.xs(lane, level=1), how='left', on='sample_name')

        # strip qiita_id from project names in sample_project column
        df['sample_project'] = df['sample_project'].map(
            lambda x: re.sub(r'_\d+$', r'', x))

        # center_project_name is a legacy column that should mirror
        # the values for sample_project.
        df['center_project_name'] = df['sample_project']

        df.to_csv(fp, sep='\t', index=False)

        if verbose:
            # assume qiita_id is extractable and is an integer, given that
            # we have already passed error-checking.
            qiita_id = get_qiita_id_from_project_name(project)
            print("%s\t%s" % (qiita_id, abspath(fp)))


if __name__ == '__main__':
    format_preparation_files()
