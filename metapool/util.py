from os.path import abspath
from metapool import (preparations_for_run, load_sample_sheet,
                      sample_sheet_to_dataframe, run_counts,
                      remove_qiita_id, preparations_for_run_mapping_file)
import re
import os
import pandas as pd


def generate_prep_from_sample_sheet(run_dir, sample_sheet, output_dir):
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

    results = []

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

        project_name = remove_qiita_id(project)
        # assume qiita_id is extractable and is an integer, given that
        # we have already passed error-checking.
        qiita_id = project.replace(project_name + '_', '')
        results.append((qiita_id, abspath(fp)))

    return results


def generate_prep_from_pre_prep(run_dir, mapping_file, output_dir):
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

    results = []

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

        project_name = remove_qiita_id(project)
        # assume qiita_id is extractable and is an integer, given that
        # we have already passed error-checking.
        qiita_id = project.replace(project_name + '_', '')
        results.append((qiita_id, abspath(fp)))

    return results
