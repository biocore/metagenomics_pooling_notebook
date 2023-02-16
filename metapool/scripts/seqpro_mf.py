#!/usr/bin/env python

import click
import os
import re
import pandas as pd

from metapool import preparations_for_run_mapping_file, run_counts


@click.command()
@click.argument('run_dir', type=click.Path(exists=True, dir_okay=True,
                                           file_okay=False))
@click.argument('mapping_file', type=click.Path(exists=True, dir_okay=False,
                                                file_okay=True))
@click.argument('output_dir', type=click.Path(writable=True))
def format_preparation_files(run_dir, mapping_file, output_dir):
    """Generate the preparation files for the projects in a run

    RUN_DIR: should be the directory where the results of running bcl2fastq are
    saved.

    MAPPING_FILE: should be a TSV file that includes information for the
    samples and projects in RUN_DIR.

    OUTPUT_DIR: directory where the outputted preparations should be saved to.

    Preparations are stratified by project and by lane. Only samples with
    non-empty files are included.
    """
    df_mapping_file = pd.read_csv(mapping_file, delimiter='\t')

    # add a faked value for 'lane' to preserve the original logic.
    # lane will always be '1' for amplicon runs.
    df_mapping_file['lane'] = pd.Series(
        [1 for x in range(len(df_mapping_file.index))])

    # add a faked column for 'Sample_ID' to preserve the original logic.
    # count-related code will need to search run_directories based on
    # sample-id, not sample-name.
    # TODO: Note there may be more to change than just '.' to '_'.
    df_mapping_file['Sample_ID'] = df_mapping_file['sample_name'].str.replace('.', '_')

    # add a faked column for 'Sample_Project' to preserve the original logic.
    df_mapping_file['Sample_Project'] = df_mapping_file.loc[:, 'Project_name']

    stats = run_counts(run_dir, df_mapping_file)

    # stats don't include well description which is the primary key to
    # merge with the preps in the loop below
    stats['sample_name'] = \
        df_mapping_file.set_index('lane', append=True)['well_description']

    # returns a map of (run, project_name, lane) -> preparation frame
    preps = preparations_for_run_mapping_file(run_dir, df_mapping_file)

    os.makedirs(output_dir, exist_ok=True)

    for (run, project, lane), df in preps.items():
        filename = os.path.join(output_dir, f'{run}.{project}.{lane}.tsv')

        # stats are indexed by sample name and lane, lane is the first
        # level index. When merging, make sure to select the lane subset
        # that we care about, otherwise we'll end up with repeated rows
        df = df.merge(stats.xs(lane, level=1), how='left',
                      on='sample_name')

        # strip qiita_id from project names in sample_project column
        df['sample_project'] = df['sample_project'].map(
            lambda x: re.sub(r'_\d+$', r'', x))

        # center_project_name is a legacy column that should mirror
        # the values for sample_project.
        df['center_project_name'] = df['sample_project']

        df.to_csv(filename, sep='\t', index=False)


if __name__ == '__main__':
    format_preparation_files()
