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
    preps = preparations_for_run(os.walk(run_dir), sample_sheet)

    os.makedirs(output_dir, exist_ok=True)

    for filename, df in preps.items():
        filename = os.path.join(output_dir, filename)

        df.to_tsv(filename, sep='\t')


    '/sequencing/ucsd_2/complete_runs/200417_A00953_0098_BHC3MYDSXY/'
    ['Lisa_Mar2020', 'Feist_11661', 'Burk_Gold_II_GMB_Sent_Plates_11666',
            'Reports', 'Hiutung_Isolates', 'Qi_WIHS_GMB_13223', 'Stats',
            'Lara_Zuniga']
    ['Undetermined_S0_L002_R1_001.fastq.gz',
            'Undetermined_S0_L002_R2_001.fastq.gz',
            '1319993.barnacle.ucsd.edu.txt',
            'Undetermined_S0_L003_R1_001.fastq.gz',
            'Undetermined_S0_L003_R2_001.fastq.gz',
            'Undetermined_S0_L003_I1_001.fastq.gz',
            'Undetermined_S0_L003_I2_001.fastq.gz',
            'Undetermined_S0_L002_I1_001.fastq.gz',
            '1320045.barnacle.ucsd.edu.txt',
            '2020-03-27_KNIGHT_RKL0051_Feist_11661_52-55_samplesheet.csv',
            '1320151.barnacle.ucsd.edu.txt',
            'Undetermined_S0_L002_I2_001.fastq.gz',
            '1320157.barnacle.ucsd.edu.txt',
            'Undetermined_S0_L004_I1_001.fastq.gz',
            'Undetermined_S0_L004_I2_001.fastq.gz',
            '1319899.barnacle.ucsd.edu.txt', '1319900.barnacle.ucsd.edu.txt',
            '20200417_Knight_RKL0052_GMB_gDNA_1-2_Lisa_Hiutung_Zuniga_samplesheet.csv',
            'Undetermined_S0_L004_R1_001.fastq.gz',
            '1341457.barnacle.ucsd.edu.txt',
            '20200417_Knight_RKL0052_GMB_1-8_samplesheet.csv',
            'Undetermined_S0_L004_R2_001.fastq.gz',
            '1319793.barnacle.ucsd.edu.txt']


if __name__ == '__main__':
    format_preparation_files()
