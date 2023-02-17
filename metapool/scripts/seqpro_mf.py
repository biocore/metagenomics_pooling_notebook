#!/usr/bin/env python

import click
import os
import pandas as pd

from metapool import preparations_for_run_mapping_file, run_counts
from metapool.metapool import bcl_scrub_name


@click.command()
@click.argument('run_dir', type=click.Path(exists=True, dir_okay=True,
                                           file_okay=False))
@click.argument('mapping_file', type=click.Path(exists=True, dir_okay=False,
                                                file_okay=True))
@click.argument('output_dir', type=click.Path(writable=True))
def format_preparation_files_mf(run_dir, mapping_file, output_dir):
    """Generate the preparation files for the projects in a run

    RUN_DIR: should be the directory where the results of running bcl2fastq are
    saved.

    MAPPING_FILE: should be a TSV file that includes information for the
    samples and projects in RUN_DIR.

    OUTPUT_DIR: directory where the outputted preparations should be saved to.

    Preparations are stratified by project and by lane. Only samples with
    non-empty files are included.
    """
    df_mf = pd.read_csv(mapping_file, delimiter='\t')

    # add a faked value for 'lane' to preserve the original logic.
    # lane will always be '1' for amplicon runs.
    df_mf['lane'] = pd.Series(
        ['1' for x in range(len(df_mf.index))])

    # add a faked column for 'Sample_ID' to preserve the original logic.
    # count-related code will need to search run_directories based on
    # sample-id, not sample-name.
    df_mf['Sample_ID'] = df_mf.apply(
        lambda x: bcl_scrub_name(x['sample_name']), axis=1)

    stats = run_counts(run_dir, df_mf)

    # stats don't include sample_name which is the primary key to merge w/the
    # preps in the loop below. Create a map of sample_names to our generated
    # Sample_IDs. Assume the sample_ids in stats that were ripped from fastq
    # filenames are a match for our Sample_IDs. Map a new column of
    # sample_names to stats using Sample_IDs/stats's index as the key.
    # This will cause sample_names for all samples where a file wasn't found
    # to be 'NaN'.
    stats = stats.join(
        df_mf[["Sample_ID", "sample_name"]].set_index('Sample_ID'))

    # returns a map of (run, project_name, lane) -> preparation frame
    preps = preparations_for_run_mapping_file(run_dir, df_mf)

    os.makedirs(output_dir, exist_ok=True)

    for (run, project, lane), df in preps.items():
        filename = os.path.join(output_dir, f'{run}.{project}.{lane}.tsv')

        stats = stats.xs(lane, level=1)

        # stats are indexed by sample name and lane, lane is the first
        # level index. When merging, make sure to select the lane subset
        # that we care about, otherwise we'll end up with repeated rows
        df = df.merge(stats, how='left', on='sample_name')

        df.to_csv(filename, sep='\t', index=False)


if __name__ == '__main__':
    format_preparation_files_mf()
