#!/usr/bin/env python

import click
import os
import pandas as pd
from os.path import abspath

from metapool import (preparations_for_run_mapping_file,
                      get_qiita_id_from_project_name)


@click.command()
@click.argument('run_dir', type=click.Path(exists=True, dir_okay=True,
                                           file_okay=False))
@click.argument('mapping_file', type=click.Path(exists=True, dir_okay=False,
                                                file_okay=True))
@click.argument('output_dir', type=click.Path(writable=True))
@click.option('--verbose', help='list prep-file output paths, study_ids',
              is_flag=True)
def format_preparation_files_mf(run_dir, mapping_file, output_dir, verbose):
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
        fp = os.path.join(output_dir, f'{run}.{project}.{lane}.tsv')

        df.to_csv(fp,
                  sep='\t',
                  index=False,
                  # finalize column order

                  columns=['sample_name', 'center_name', 'center_project_name',
                           'experiment_design_description', 'instrument_model',
                           'lane', 'library_construction_protocol',
                           'platform', 'run_center', 'run_date', 'run_prefix',
                           'runid', 'sample_plate', 'sequencing_meth',
                           'barcode', 'linker', 'primer', 'extraction_robot',
                           'extractionkit_lot', 'mastermix_lot', 'orig_name',
                           'pcr_primers', 'plating', 'primer_date',
                           'primer_plate', 'processing_robot', 'project_name',
                           'target_gene', 'target_subfragment',
                           'tm1000_8_tool', 'tm300_8_tool', 'tm50_8_tool',
                           'tm10_8_tool', 'water_lot', 'well_description',
                           'well_id_384', 'well_id_96'])

        if verbose:
            # assume qiita_id is extractable and is an integer, given that
            # we have already passed error-checking.
            qiita_id = get_qiita_id_from_project_name(project)
            print("%s\t%s" % (qiita_id, abspath(fp)))


if __name__ == '__main__':
    format_preparation_files_mf()
