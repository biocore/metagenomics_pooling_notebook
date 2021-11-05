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

        # check all the required keys are present, otherwise the file could be
        # malformed and we would only see a weird KeyError exception
        if ('summary' not in stats or
            'after_filtering' not in stats['summary'] or
           'total_reads' not in stats['summary']['after_filtering']):
            raise ValueError(f'The fastp log for {path} is malformed')

        return int(stats['summary']['after_filtering']['total_reads'])


def _parse_samtools_counts(path):
    with open(path, 'r') as f:
        matches = re.match(SAMTOOLS_PATTERN, f.read())

        if matches is None:
            raise ValueError(f'The samtools log for {path} is malformed')

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

    # ignore the things not present in the sheet
    out = out[out['Sample_ID'].isin(expected)]

    dups = out.duplicated(subset=['Sample_ID', 'Lane'])
    if dups.any():
        pairs = [f"{r['Sample_ID']} in lane {r['Lane']}"
                 for _, r in out[dups].iterrows()]

        # when running bcl2fastq/bclconvert multiple times you can run into
        # situations where the cell number is the only thing that changes. For
        # those situations, make sure you flag this as a possible error
        raise ValueError('Multiple matches found for the same samples in'
                         ' the same lane, only one match is expected: %s' %
                         ', '.join(pairs))

    if expected > found:
        warnings.warn(f'No {name} log found for these samples: %s' %
                      ', '.join(expected - found))

    out[name] = out.path.apply(funk)

    out.drop(columns=['path', 'Sample_Project'], inplace=True)
    out.set_index(['Sample_ID', 'Lane'], inplace=True, verify_integrity=True)
    return out


def _safe_get(_document, _key):
    """Prevent generic KeyError exceptions"""
    if _key not in _document:
        raise KeyError(f'bcl stats file is missing {_key} attribute')
    else:
        return _document[_key]


def bcl2fastq_counts(run_dir, sample_sheet):
    bcl2fastq_path = os.path.join(os.path.abspath(run_dir),
                                  'Stats/Stats.json')
    bclconvert_path = os.path.join(os.path.abspath(run_dir),
                                   'Reports/Demultiplex_Stats.csv')

    if os.path.exists(bcl2fastq_path):
        if os.path.exists(bclconvert_path):
            raise IOError(f"both '{bcl2fastq_path}' and '{bclconvert_path}'"
                          " exist")
        else:
            return _bcl2fastq_counts(bcl2fastq_path)
    elif os.path.exists(bclconvert_path):
        return _bclconvert_counts(bclconvert_path)
    else:
        raise IOError(f"Cannot find Stats.json '{bcl2fastq_path}' or "
                      f"Demultiplex_Stats.csv '{bclconvert_path}' for this"
                      " run")


def _bcl2fastq_counts(path):
    with open(path) as fp:
        contents = json.load(fp)

    out = []
    for lane in _safe_get(contents, 'ConversionResults'):
        table = pd.DataFrame(_safe_get(lane, 'DemuxResults'))
        table['Lane'] = str(_safe_get(lane, 'LaneNumber'))

        out.append(table)

    out = pd.concat(out)
    out.rename(columns={'SampleId': 'Sample_ID', 'NumberReads': 'raw_reads'},
               inplace=True)
    out = out[['Sample_ID', 'Lane', 'raw_reads']]
    out.set_index(['Sample_ID', 'Lane'], inplace=True, verify_integrity=True)
    return out


def _bclconvert_counts(path):
    # read the csv in from file
    df = pd.read_csv(path)
    # subselect only the columns we're concerned with
    df = df[["SampleID", "Lane", "# Reads"]]
    # filter out rows that reference an 'Undetermined' fastq.gz file
    # and create our own copy to return to the user
    df = df.loc[df['SampleID'] != 'Undetermined'].copy()
    # rename columns to standard values for metapool
    df.rename(columns={'SampleID': 'Sample_ID', '# Reads': 'raw_reads'},
              inplace=True)
    # create indexes on these columns
    df['Lane'] = df['Lane'].astype(str)
    df.set_index(['Sample_ID', 'Lane'], inplace=True, verify_integrity=True)

    return df


def fastp_counts(run_dir, sample_sheet):
    return _parsefier(run_dir, sample_sheet, 'json', '.json',
                      'quality_filtered_reads',
                      _parse_fastp_counts)


def minimap2_counts(run_dir, sample_sheet):
    return _parsefier(run_dir, sample_sheet, 'samtools', '.log',
                      'non_host_reads', _parse_samtools_counts)


def run_counts(run_dir, sample_sheet):
    out = bcl2fastq_counts(run_dir, sample_sheet).join([
            fastp_counts(run_dir, sample_sheet),
            minimap2_counts(run_dir, sample_sheet),

        ])

    # convenience columns to assess sample quality
    out['fraction_passing_quality_filter'] = (out['quality_filtered_reads'] /
                                              out['raw_reads'])
    out['fraction_non_human'] = (out['non_host_reads'] /
                                 out['quality_filtered_reads'])

    out.fillna(value='NA', inplace=True)

    return out
