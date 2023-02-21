#!/usr/bin/env python

import click
import os
import pandas as pd

from metapool import preparations_for_run_mapping_file


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

    # run_counts cannot be determined for a single multiplexed file.

    # returns a map of (run, project_name, lane) -> preparation frame
    preps = preparations_for_run_mapping_file(run_dir, df_mf)

    os.makedirs(output_dir, exist_ok=True)

    for (run, project, lane), df in preps.items():
        filename = os.path.join(output_dir, f'{run}.{project}.{lane}.tsv')

        df = df.drop(['index', 'index2', 'i5_index_id', 'i7_index_id',
                      'sample_well', 'well_description'], axis=1)

        df.to_csv(filename,
                  sep='\t',
                  index=False,
                  # finalize column order
                  columns=['sample_name',
                           'experiment_design_description',
                           'library_construction_protocol',
                           'platform',
                           'run_center',
                           'run_date',
                           'run_prefix',
                           'sequencing_meth',
                           'center_name',
                           'center_project_name',
                           'instrument_model',
                           'runid',
                           'lane',
                           'sample_project',
                           'sample_plate'])


if __name__ == '__main__':
    format_preparation_files_mf()
