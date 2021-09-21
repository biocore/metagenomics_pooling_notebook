import re
import os
import glob
import pandas as pd
import json
import warnings

# first group is the sample name from the sample sheet, second group is the
# cell number, third group is the lane number, and fourth group is the
# forward/reverse/index id.
#
# Here's a few examples for this regular expression: tinyurl.com/filenamepatt
SAMPLE_PATTERN = re.compile(r'(.*)_(S\d{1,4})_(L\d{1,3})_([RI][12]).*')

# first group is the number of sequences written to fwd and reverse files
# Here's a few examples for this regular expression: tinyurl.com/samtoolspatt
SAMTOOLS_PATTERN = re.compile(r'\[.*\] processed (\d+) reads',
                              flags=re.MULTILINE)


def _extract_name_and_lane(filename):
    search = re.match(SAMPLE_PATTERN, filename)
    if search is None:
        raise ValueError(f'Unrecognized filename pattern {filename}')

    name, _, lane, _ = search.groups()

    # remove the leading L and any leading zeroes
    lane = lane[1:].lstrip('0')
    return name, lane


def _parse_fastp_counts(path):
    with open(path) as fp:
        stats = json.load(fp)
        return int(stats['summary']['after_filtering']['total_reads'])


def _parse_samtools_counts(path):
    with open(path, 'r') as f:
        matches = re.match(SAMTOOLS_PATTERN, f.read())

        if matches is None:
            raise ValueError('The samtools output for %s is malformed')

        # divided by 2 because samtools outputs the number of records found
        # in the forward and reverse files
        return int(matches.groups()[0]) / 2.0


def _parsefier(run_dir, sample_sheet, subdir, suffix, name, funk):
    """High order helper to search through a run directory

    Parameters
    ----------
    run_dir: str
        Illumina's run directory.
    sample_sheet: metapool.KLSampleSheet
        Sample sheet for the samples to get counts for.
    subdir: str
        Name of the directory in the project folder.
    suffix: str
        Suffix of the log files.
    name: str
        Column name for the output counts.
    funk: callable
        Function to parse log files.

    Returns
    -------
    pd.DataFrame
        Table with sample, lane and counts.
    """
    out = []

    for project, lane in {(s.Sample_Project, s.Lane) for s in sample_sheet}:
        lane = lane.zfill(3)

        # log files are named after the sequence files themselves, specifically
        # the forward sequences, so we just make sure we match the right lane,
        # the forward file and the suffix. The suffix is something like ".log"
        # or "*.json" if you expect to see other characters before the
        # extension
        logs = glob.glob(os.path.join(run_dir, project, subdir,
                                      f'*_L{lane}_R1_001' + suffix))

        for log in logs:
            out.append([*_extract_name_and_lane(os.path.basename(log)),
                        project, log])

    out = pd.DataFrame(columns=['Sample_ID', 'Lane', 'Sample_Project', 'path'],
                       data=out)

    found = set(out['Sample_ID'])
    expected = {s.Sample_ID for s in sample_sheet}

    if found > expected:
        out = out[out['Sample_ID'].isin(expected)]

        # when running bcl2fastq/bclconvert multiple times you can run into
        # situations where the cell number is the only thing that changes. For
        # those situations, make sure you flag this as a possible error
        if len(out) > len(expected):
            raise ValueError('Multiple matches found for the same samples in'
                             ' the same lane, only one match is expected')

    elif expected > found:
        warnings.warn(f'No {name} log found for these samples: %s' %
                      ', '.join(expected - found))

    out[name] = out.path.apply(funk)
    # TODO: do we need to get rid of sample_project?
    out.drop(columns=['path', 'Sample_Project'], inplace=True)
    out.set_index(['Sample_ID', 'Lane'], inplace=True, verify_integrity=True)
    return out


def bcl2fastq_counts(run_dir, sample_sheet):
    path = os.path.join(os.path.abspath(run_dir), 'Stats/Stats.json')
    with open(path) as fp:
        contents = json.load(fp)

    out = []
    for lane in contents['ConversionResults']:
        table = pd.DataFrame(lane['DemuxResults'])
        table['Lane'] = str(lane['LaneNumber'])
        out.append(table)

    out = pd.concat(out)
    out.rename(columns={'SampleId': 'Sample_ID', 'NumberReads': 'bcl'},
               inplace=True)
    out = out[['Sample_ID', 'Lane', 'bcl']]
    out.set_index(['Sample_ID', 'Lane'], inplace=True, verify_integrity=True)
    return out


def fastp_counts(run_dir, sample_sheet):
    return _parsefier(run_dir, sample_sheet, 'json', '.json', 'fastp',
                      _parse_fastp_counts)


def minimap2_counts(run_dir, sample_sheet):
    return _parsefier(run_dir, sample_sheet, 'samtools', '.log',
                      'minimap2', _parse_samtools_counts)


def run_counts(run_dir, sample_sheet):
    out = bcl2fastq_counts(run_dir, sample_sheet).join([
            fastp_counts(run_dir, sample_sheet),
            minimap2_counts(run_dir, sample_sheet),

        ])

    out.fillna(value='NA', inplace=True)

    return out
