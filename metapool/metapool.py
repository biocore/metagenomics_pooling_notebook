import json
import re
import numpy as np
import pandas as pd
import string
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
from random import choices
from configparser import ConfigParser
from qiita_client import QiitaClient
from .prep import remove_qiita_id


REVCOMP_SEQUENCERS = ['HiSeq4000', 'MiniSeq', 'NextSeq', 'HiSeq3000',
                      'iSeq', 'NovaSeq']
OTHER_SEQUENCERS = ['HiSeq2500', 'HiSeq1500', 'MiSeq']


def extract_stats_metadata(stats_json_fp, lane_numbers):
    """
    Extract metadata from Stats.json file.
    :param stats_json_fp: A file-path to Stats.json file.
    :param lane_numbers: A list of ints (lane-numbers) to extract metadata for.
    :return: A legacy 3-tuple containing:
                A dict() containing metadata,
                A DataFrame() containing conversion-results,
                A DataFrame() containing unknown barcodes found.
    """
    with open(stats_json_fp, 'r') as f:
        stats = json.load(f)

    # gather general metadata
    flowcell = stats['Flowcell']
    run_number = stats['RunNumber']
    run_id = stats['RunId']
    conversion_results = stats['ConversionResults']
    unknown_barcodes = stats['UnknownBarcodes']

    # lane_numbers must be a list with one or more valid lane_numbers.
    valid_lane_numbers = [x['LaneNumber'] for x in conversion_results]
    if not set(valid_lane_numbers).issuperset(set(lane_numbers)):
        raise ValueError(f"Valid lane numbers are: {valid_lane_numbers}")

    # extract conversion_results from a list of nested dicts to a list of
    # dataframes containing only the information required.
    filtered_results = []
    for conversion_result in conversion_results:
        if conversion_result['LaneNumber'] in lane_numbers:
            # process the metadata for this lane.
            lane_number = conversion_result['LaneNumber']
            rows = []
            for demux_result in conversion_result['DemuxResults']:
                # these prefixes are meant to simplify access for deeply
                # nested elements.
                index_metrics = demux_result['IndexMetrics'][0]
                read_metrics_0 = demux_result['ReadMetrics'][0]
                read_metrics_1 = demux_result['ReadMetrics'][1]

                # extract required general metadata.
                number_reads = demux_result['NumberReads']
                sample_id = demux_result['SampleId']
                sample_name = demux_result['SampleName']
                primary_yield = demux_result['Yield']

                # extract required index-metrics.
                index_sequence = index_metrics['IndexSequence']
                mismatch_0 = index_metrics['MismatchCounts']['0']
                mismatch_1 = index_metrics['MismatchCounts']['1']

                # extract required read-metrics.
                yield_r1 = read_metrics_0['Yield']
                yield_q30r1 = read_metrics_0['YieldQ30']
                yield_r2 = read_metrics_1['Yield']
                yield_q30r2 = read_metrics_1['YieldQ30']

                # store data as a tuple until processing is complete.
                # order the tuple in legacy order, for backwards-compatibility.
                rows.append((index_sequence, sample_name, sample_id,
                             lane_number, mismatch_0, mismatch_1, number_reads,
                             yield_r1, yield_q30r1, yield_r2, yield_q30r2,
                             primary_yield))

            # create dataframe at one time using all rows, instead of using
            # concat() or append(). columns listed in legacy order,
            # for backwards-compatibility.
            filtered_results.append(pd.DataFrame(rows,
                                                 columns=['IndexSequence',
                                                          'SampleName',
                                                          'SampleId',
                                                          'Lane',
                                                          'Mismatch0',
                                                          'Mismatch1',
                                                          'NumberReads',
                                                          'YieldR1',
                                                          'YieldQ30R1',
                                                          'YieldR2',
                                                          'YieldQ30R2',
                                                          'Yield']))
    filtered_unknowns = []
    for unknown_barcode in unknown_barcodes:
        if unknown_barcode['Lane'] in lane_numbers:
            # convert dict of (k,v) pairs into list form.
            lst = [(k, unknown_barcode['Lane'], v) for k, v in
                   unknown_barcode['Barcodes'].items()]
            df = pd.DataFrame(lst, columns=['IndexSequence', 'Lane', 'Value'])
            filtered_unknowns.append(df)

    # convert the list of dataframes into a single dataframe spanning all
    # selected lanes.
    filtered_unknowns = pd.concat(filtered_unknowns)
    filtered_results = pd.concat(filtered_results)

    # set index to legacy index + lane-number.
    filtered_results = filtered_results.set_index(['IndexSequence',
                                                   'SampleName',
                                                   'SampleId'])

    filtered_unknowns.set_index(["IndexSequence"], inplace=True)

    return {'Flowcell': flowcell,
            'RunNumber': run_number,
            'RunId': run_id}, filtered_results, filtered_unknowns


def sum_lanes(multi_lane_df, lanes_to_sum):
    """
    Sum the values of a DataFrame across a given set of lanes.
    :param multi_lane_df: A DataFrame containing a 'Lane' column.
    :param lanes_to_sum: A list of lane numbers to sum across.
    :return: A DataFrame containing sums across the given set of lanes.
    """
    if 'Lane' not in multi_lane_df:
        raise ValueError("DataFrame must contain 'Lane' column.")

    if not set(lanes_to_sum).issubset(set(multi_lane_df.Lane.unique())):
        raise ValueError("One or more specified lanes does not occur in "
                         "DataFrame.")

    total = None

    for lane_id in lanes_to_sum:
        # create a dataframe w/only the rows from lane_id
        lane = multi_lane_df[multi_lane_df['Lane'] == lane_id]
        # drop the 'Lane' column
        lane = lane.drop(columns=['Lane'])

        if total is None:
            total = lane
        else:
            # add the values for lane to the running total.

            # at times, a lane will contain an index not present in another
            # lane. Find and add these indexes as empty rows so that they
            # won't disappear from the results.
            rt = set(total.index.values)
            cl = set(lane.index.values)

            for x in cl - rt:
                # for each key in cl not found in the running total, create
                # a new empty row in the running total so that this row of
                # data does not disappear from the final result.

                # one less column assuming that 'Lane' is present.
                total.loc[x] = [0] * (len(multi_lane_df.columns) - 1)

            # now sum the values in current_lane to the running-total.
            total = total.add(lane, fill_value=0)

    return total


def read_plate_map_csv(f, sep='\t', qiita_oauth2_conf_fp=None):
    """
    reads tab-delimited plate map into a Pandas dataframe and if
    qiita_oauth2_conf_fp is passed the code will also check for Sample
    overlap between the plate_map and the Qiita study.

    Parameters
    ----------
    f: fp or open filehandle
        plate map file
    qiita_oauth2_conf_fp: str, filepath
        the path to the oauth2 configuration file to connect to Qiita

    Returns
    -------
    plate_df: pandas DataFrame object
        DataFrame relating sample name, well location, and blank status

    Raises
    ------
    UserWarning
        If there are wells with no sample names associated with them.
        If qiita_oauth2_conf_fp is None, no Qiita sample name validation.
    ValueError
        If plate_map doesn't have a 'Project Name' column.
        If there are repeated sample names.
        If the sample(s) in the plate_map are not a subset of Qiita samples.
    """

    qiita_validate = False
    if qiita_oauth2_conf_fp is not None:
        parser = ConfigParser()
        with open(qiita_oauth2_conf_fp, "r") as qfp:
            parser.read_file(qfp)
        url = parser.get("qiita-oauth2", "URL")
        client_id = parser.get("qiita-oauth2", "CLIENT_ID")
        client_secret = parser.get("qiita-oauth2", "CLIENT_SECRET")
        server_cert = parser.get("qiita-oauth2", "SERVER_CERT")
        qclient = QiitaClient(url, client_id, client_secret, server_cert)
        qiita_validate = True
    else:
        warnings.warn('No qiita_oauth2_conf_fp set so not checking plate_map '
                      'and Qiita study overlap')

    plate_df = pd.read_csv(f, sep=sep)
    if 'Project Name' not in plate_df.columns:
        raise ValueError('Missing `Project Name` column.')
    plate_df['Well'] = plate_df['Row'] + plate_df['Col'].map(str)

    null_samples = plate_df.Sample.isnull()
    if null_samples.any():
        warnings.warn(('This plate map contains %d empty wells, these will be '
                      'ignored') % null_samples.sum())

        # slice to the non-null samples and reset the index so samples are
        # still indexed with a continuous list of integers
        plate_df = plate_df[~null_samples]
        plate_df.reset_index(inplace=True, drop=True)

    duplicated_samples = plate_df.Sample[plate_df.Sample.duplicated()]
    if len(duplicated_samples):
        raise ValueError('The following sample names are duplicated %s' %
                         ', '.join(sorted(duplicated_samples)))

    if qiita_validate:
        errors = []
        for project, _df in plate_df.groupby(['Project Name']):
            project_name = remove_qiita_id(project)
            qiita_id = project.replace(f'{project_name}_', '')
            qurl = f'/api/v1/study/{qiita_id}/samples'

            plate_map_samples = {
                s for s in _df['Sample'] if not s.startswith('BLANK')}
            qsamples = {
                s.replace(f'{qiita_id}.', '') for s in qclient.get(qurl)}
            sample_name_diff = plate_map_samples - set(qsamples)
            if sample_name_diff:
                # before we report as an error, check tube_id
                error_tube_id = 'No tube_id column in Qiita.'
                if 'tube_id' in qclient.get(f'{qurl}/info')['categories']:
                    tids = qclient.get(f'{qurl}/categories=tube_id')['samples']
                    tids = {tid[0] for _, tid in tids.items()}
                    tube_id_diff = plate_map_samples - tids
                    if not tube_id_diff:
                        continue
                    len_tube_id_overlap = len(tube_id_diff)

                    tids = list(tids)
                    if len(tids) > 5:
                        # select five samples at random to display to the user.
                        tids_example = ', '.join(choices(tids, k=5))
                    else:
                        tids_example = ', '.join(tids)

                    error_tube_id = (
                        f'tube_id in Qiita but {len_tube_id_overlap} missing '
                        f'samples. Some samples from tube_id: {tids_example}.')
                len_overlap = len(sample_name_diff)

                qsamples = list(qsamples)
                if len(qsamples) > 5:
                    # select five samples at random to display to the user.
                    samples_example = ', '.join(choices(qsamples, k=5))
                else:
                    samples_example = ', '.join(qsamples)

                missing = ', '.join(sorted(sample_name_diff)[:4])
                errors.append(
                    f'{project} has {len_overlap} missing samples (i.e. '
                    f'{missing}). Some samples from Qiita: {samples_example}. '
                    f'{error_tube_id}')
        if errors:
            raise ValueError('\n'.join(errors))

    return plate_df


# method to read minipico output
def read_pico_csv(f, sep='\t', plate_reader='Synergy_HT',
                  conc_col_name='Sample DNA Concentration'):
    """
    reads tab-delimited pico quant

    Parameters
    ----------
    f: fp or open filehandle
        pico quant file
    sep: str
        sep char used in quant file
    plate_reader: str
        plate reader used to generate quant file ['Synergy_HT',
        'SpectraMax_i3x']
    conc_col_name: str
        name to use for concentration column output

    Returns
    -------
    pico_df: pandas DataFrame object
        DataFrame relating well location and DNA concentration
    """
    if plate_reader == 'Synergy_HT':
        encoding, skipfooter = None, 5
    elif plate_reader == 'SpectraMax_i3x':
        encoding, skipfooter = 'utf-16', 15
    else:
        raise ValueError("Invalid plate reader %s" % plate_reader)
    if not hasattr(f, 'read'):
        f = open(f, encoding=encoding)

    pico_df = pd.read_csv(f, sep=sep, skiprows=2,
                          skipfooter=skipfooter, engine='python')

    # synergy's concentration column is "Concentration", spectramax's is
    # [Concentration]. Rename will ignore any labels not in the dataframe so
    # only one of the two label updates should happen
    pico_df.rename(columns={'Concentration': conc_col_name,
                            '[Concentration]': conc_col_name,
                            'Wells': 'Well'}, inplace=True)

    pico_df = pico_df[['Well', conc_col_name]].copy()

    # coerce oddball concentrations to np.nan
    pico_df[conc_col_name] = \
        pd.to_numeric(pico_df[conc_col_name], errors='coerce')
    if plate_reader == 'SpectraMax_i3x':
        # limit concentration range (0 - 60 )
        pico_df[conc_col_name] = np.clip(pico_df[conc_col_name], 0, 60)

    return pico_df


def calculate_norm_vol(dna_concs, ng=5, min_vol=2.5, max_vol=3500,
                       resolution=2.5):
    """
    Calculates nanoliters of each sample to add to achieve a normalized pool

    Parameters
    ----------
    dna_concs : numpy array of float
        The concentrations calculated via PicoGreen (ng/uL)
    ng : float
        The amount of DNA to pool (ng)
    max_vol : float
        The maximum volume to pool (nL)
    min_vol : float
        The minimum volume to pool (nL)

    Returns
    -------
    sample_vols : numpy array of float
        The volumes to pool (nL)
    """
    sample_vols = ng / np.nan_to_num(dna_concs) * 1000

    sample_vols = np.clip(sample_vols, min_vol, max_vol)

    sample_vols = np.round(sample_vols / resolution) * resolution

    return sample_vols


def format_dna_norm_picklist(dna_vols, water_vols, wells, dest_wells=None,
                             dna_concs=None, sample_names=None,
                             sample_plates=None, water_plate_name='Water',
                             dna_plate_type='384PP_AQ_BP2_HT',
                             water_plate_type='384PP_AQ_BP2_HT',
                             dest_plate_name='NormalizedDNA'):
    """
    Writes Echo-format pick list to achieve a normalized input DNA pool

    Parameters
    ----------
    dna_vols:  numpy array of float
        The volumes of dna to add
    water_vols:  numpy array of float
        The volumes of water to add
    wells: numpy array of str
        The well codes in the same orientation as the DNA concentrations
    dest_wells: numpy array of str
        The well codes, in the same orientation as `wells`,
        in which to place each sample if reformatting
    dna_concs:  numpy array of float
        The concentrations calculated via PicoGreen (ng/uL)
    sample_names: numpy array of str
        The sample names in the same orientation as the DNA concentrations
    sample_plates: numpy array of str
        The sample plates in the same orientation as the DNA concentrations

    Returns
    -------
    picklist : str
        The Echo formatted pick list
    """

    # check that arrays are the right size
    if dna_vols.shape != wells.shape != water_vols.shape:
        raise ValueError(('dna_vols %r has a size different from wells %r or '
                          'water_vols %r') %
                         (dna_vols.shape, wells.shape, water_vols.shape))

    # if destination wells not specified, use source wells
    if dest_wells is None:
        dest_wells = wells

    if sample_names is None:
        sample_names = np.empty(dna_vols.shape) * np.nan
    if sample_plates is None:
        sample_plates = 'Sample'
    if isinstance(sample_plates, str):
        sample_plates = np.full_like(dna_vols, sample_plates, dtype=object)
    if dna_plate_type is None:
        dna_plate_type = '384PP_AQ_BP2_HT'
    if isinstance(dna_plate_type, str):
        dna_plate_type = np.full_like(dna_vols, dna_plate_type, dtype=object)
    if dna_concs is None:
        dna_concs = np.empty(dna_vols.shape) * np.nan
    if (dna_concs.shape != sample_names.shape != dna_vols.shape
       != sample_plates.shape != dna_plate_type.shape):
        raise ValueError(('dna_vols %r has a size different from dna_concs %r'
                          ' or sample_names %r') %
                         (dna_vols.shape, dna_concs.shape, sample_names.shape))

    picklist = ''

    # header
    picklist += ('Sample\tSource Plate Name\tSource Plate Type\tSource Well\t'
                 'Concentration\tTransfer Volume\tDestination Plate Name\t'
                 'Destination Well')

    # water additions
    for index, sample in np.ndenumerate(sample_names):
        picklist += '\n' + '\t'.join([str(sample), water_plate_name,
                                      water_plate_type, str(wells[index]),
                                      str(dna_concs[index]),
                                      str(water_vols[index]),
                                      dest_plate_name, str(dest_wells[index])])
    # DNA additions
    for index, sample in np.ndenumerate(sample_names):
        picklist += '\n' + '\t'.join([str(sample), str(sample_plates[index]),
                                      str(dna_plate_type[index]),
                                      str(wells[index]), str(dna_concs[index]),
                                      str(dna_vols[index]),
                                      dest_plate_name, str(dest_wells[index])])

    return picklist


def assign_index(samples, index_df, start_idx=0):
    """
    Writes Echo-format pick list to achieve a normalized input DNA pool

    Parameters
    ----------
    samples:  int
        The number of samples for which to get indices
    index_df:  pandas DataFrame
        The dataframe of complete index combinations and information
    start_idx: int
        The starting index combo to use

    Returns
    -------
    indices : pandasDataFrame
        The index information for the chosen indices
    """

    indices = index_df.iloc[start_idx:(start_idx + samples)]

    return indices


def format_index_picklist(sample_names, sample_wells, indices,
                          i5_vol=250, i7_vol=250,
                          i5_plate_type='384LDV_AQ_B2_HT',
                          i7_plate_type='384LDV_AQ_B2_HT',
                          dest_plate_name='IndexPCRPlate'):
    """
    Writes Echo-format pick list to achieve a normalized input DNA pool

    Parameters
    ----------
    sample_names:  array-like of str
        The sample names matching index order of indices
    sample_wells:  array-like of str
        The wells matching sample name order
    indices: pandas DataFrame
        The dataframe with index info matching sample_names

    Returns
    -------
    picklist : str
        The Echo formatted pick list
    """

    # check that arrays are the right size
    if len(sample_names) != len(sample_wells) != len(indices):
        raise ValueError(('sample_names (%s) has a size different from '
                          'sample_wells (%s) or index list (%s)') %
                         (len(sample_names), len(sample_wells), len(indices)))

    picklist = ''

    # header
    picklist += ('Sample\tSource Plate Name\tSource Plate Type\tSource Well\t'
                 'Transfer Volume\tIndex Name\tIndex Sequence\tIndex Combo\t'
                 'Destination Plate Name\tDestination Well')

    # i5 additions
    for i, (sample, well) in enumerate(zip(sample_names, sample_wells)):
        picklist += '\n' + '\t'.join([str(sample), indices.iloc[i]['i5 plate'],
                                      i5_plate_type,
                                      indices.iloc[i]['i5 well'], str(
                                          i5_vol), indices.iloc[i]['i5 name'],
                                      indices.iloc[i]['i5 sequence'], str(
                                          indices.iloc[i]['index combo']),
                                      dest_plate_name, well])
    # i7 additions
    for i, (sample, well) in enumerate(zip(sample_names, sample_wells)):
        picklist += '\n' + '\t'.join([str(sample), indices.iloc[i]['i7 plate'],
                                      i7_plate_type,
                                      indices.iloc[i]['i7 well'], str(
                                          i7_vol), indices.iloc[i]['i7 name'],
                                      indices.iloc[i]['i7 sequence'], str(
                                          indices.iloc[i]['index combo']),
                                      dest_plate_name, well])

    return picklist


def compute_qpcr_concentration(cp_vals, m=-3.231, b=12.059, dil_factor=25000):
    """Computes molar concentration of libraries from qPCR Cp values.

    Returns a 2D array of calculated concentrations, in nanomolar units

    Parameters
    ----------
    cp_vals : numpy array of float
        The Cp values parsed from the plate reader
    m : float
        The slope of the qPCR standard curve
    b : float
        The intercept of the qPCR standard curve
    dil_factor: float or int
        The dilution factor of the samples going into the qPCR

    Returns
    -------10
    np.array of floats
        A 2D array of floats
    """
    qpcr_concentration = np.power(10, ((cp_vals - b) / m)) * dil_factor / 1000

    return qpcr_concentration


def compute_shotgun_pooling_values_eqvol(sample_concs, total_vol=60.0):
    """Computes molar concentration of libraries from qPCR Cp values.

    Returns a 2D array of calculated concentrations, in nanomolar units

    Parameters
    ----------
    sample_concs : numpy array of float
        The concentrations calculated via qPCR (nM)
    total_vol : float
        The total volume to pool (uL)

    Returns
    -------
    np.array of floats
        A 2D array of floats
    """
    per_sample_vol = (total_vol / sample_concs.size) * 1000.0

    sample_vols = np.zeros(sample_concs.shape) + per_sample_vol

    return sample_vols


def compute_shotgun_pooling_values_qpcr(sample_concs, sample_fracs=None,
                                        min_conc=10, floor_conc=50,
                                        total_nmol=.01):
    """Computes pooling volumes for samples based on qPCR estimates of
    nM concentrations (`sample_concs`).

    Reads in qPCR values in nM from output of `compute_qpcr_concentration`.
    Samples must be above a minimum concentration threshold (`min_conc`,
    default 10 nM) to be included. Samples above this threshold but below a
    given floor concentration (`floor_conc`, default 50 nM) will be pooled as
    if they were at the floor concentration, to avoid overdiluting the pool.

    Samples can be assigned a target molar fraction in the pool by passing a
    np.array (`sample_fracs`, same shape as `sample_concs`) with fractional
    values per sample. By default, will aim for equal molar pooling.

    Finally, total pooling size is determined by a target nanomolar quantity
    (`total_nmol`, default .01). For a perfect 384 sample library, in which you
    had all samples at a concentration of exactly 400 nM and wanted a total
    volume of 60 uL, this would be 0.024 nmol.

    Parameters
    ----------
    sample_concs: 2D array of float
        nM calculated by compute_qpcr_concentration
    sample_fracs: 2D of float
        fractional value for each sample (default 1/N)
    min_conc: float
        minimum nM concentration to be included in pool
    floor_conc: float
        minimum value for pooling for samples above min_conc
        corresponds to a maximum vol in pool
    total_nmol : float
        total number of nM to have in pool

    Returns
    -------
    sample_vols: np.array of floats
        the volumes in nL per each sample pooled
    """

    if sample_fracs is None:
        sample_fracs = np.ones(sample_concs.shape) / sample_concs.size

    # get samples above threshold
    sample_fracs_pass = sample_fracs.copy()
    sample_fracs_pass[sample_concs <= min_conc] = 0

    # renormalize to exclude lost samples
    sample_fracs_pass *= 1/sample_fracs_pass.sum()

    # floor concentration value
    sample_concs_floor = sample_concs.copy()
    sample_concs_floor[sample_concs < floor_conc] = floor_conc

    # calculate volumetric fractions including floor val
    sample_vols = (total_nmol * sample_fracs_pass) / sample_concs_floor

    # convert L to nL
    sample_vols *= 10**9

    return sample_vols


def compute_shotgun_pooling_values_qpcr_minvol(sample_concs, sample_fracs=None,
                                               floor_vol=100, floor_conc=40,
                                               total_nmol=.01):
    """Computes pooling volumes for samples based on qPCR estimates of
    nM concentrations (`sample_concs`), taking a minimum volume of samples
    below a threshold.

    Reads in qPCR values in nM from output of `compute_qpcr_concentration`.
    Samples below a minimum concentration (`floor_conc`, default 40 nM)
    will be included, but at a decreased volume (`floor_vol`, default 100 nL)
    to avoid overdiluting the pool.

    Samples can be assigned a target molar fraction in the pool by passing a
    np.array (`sample_fracs`, same shape as `sample_concs`) with fractional
    values per sample. By default, will aim for equal molar pooling.

    Finally, total pooling size is determined by a target nanomolar quantity
    (`total_nmol`, default .01). For a perfect 384 sample library, in which you
    had all samples at a concentration of exactly 400 nM and wanted a total
    volume of 60 uL, this would be 0.024 nmol.

    For a Novaseq, we expect to need 150 uL at 4 nM, or about 0.0006 nmol.
    Taking into account sample loss on the pippin prep (1/2) and molar loss
    due to exclusion of primer dimers (1/2), figure we need 4 times that or
    0.0024.

    Parameters
    ----------
    sample_concs: 2D array of float
        nM calculated by compute_qpcr_concentration
    sample_fracs: 2D of float
        fractional value for each sample (default 1/N)
    floor_vol: float
        volume (nL) at which samples below floor_conc will be pooled
    floor_conc: float
        minimum value (nM) for pooling at real estimated value (default 40)
    total_nmol : float
        total number of nM to have in pool

    Returns
    -------
    sample_vols: np.array of floats
        the volumes in nL per each sample pooled
    """

    if sample_fracs is None:
        sample_fracs = np.ones(sample_concs.shape) / sample_concs.size

    # calculate volumetric fractions including floor val
    sample_vols = (total_nmol * sample_fracs) / sample_concs

    # convert L to nL
    sample_vols *= 10**9

    # drop volumes for samples below floor concentration to floor_vol
    sample_vols[sample_concs < floor_conc] = floor_vol

    return sample_vols


def estimate_pool_conc_vol(sample_vols, sample_concs):
    """Estimates the actual molarity and volume of a pool.

    Parameters
    ----------
    sample_concs : numpy array of float
        The concentrations calculated via qPCR (nM)
    sample_vols : numpy array of float
        The calculated pooling volumes (nL)

    Returns
    -------
    pool_conc : float
        The estimated actual concentration of the pool, in nM
    total_vol : float
        The total volume of the pool, in nL
    """
    # scalar to adjust nL to L for molarity calculations
    nl_scalar = 10**-9

    # calc total pool pmols
    total_pmols = np.multiply(sample_concs, sample_vols) * nl_scalar

    # calc total pool vol in nanoliters
    total_vol = sample_vols.sum()

    # pool pM is total pmols divided by total liters
    # (total vol in nL * 1 L / 10^9 nL)
    pool_conc = total_pmols.sum() / (total_vol * nl_scalar)

    return pool_conc, total_vol


def format_pooling_echo_pick_list(vol_sample,
                                  max_vol_per_well=60000,
                                  dest_plate_shape=[16, 24]):
    """Format the contents of an echo pooling pick list

    Parameters
    ----------
    vol_sample : 2d numpy array of floats
        The per well sample volume, in nL
    max_vol_per_well : 2d numpy array of floats
        Maximum destination well volume, in nL
    """
    contents = ['Source Plate Name,Source Plate Type,Source Well,'
                'Concentration,Transfer Volume,Destination Plate Name,'
                'Destination Well']
    # Write the sample transfer volumes
    rows, cols = vol_sample.shape

    # replace NaN values with 0s to leave a trail of unpooled wells
    pool_vols = np.nan_to_num(vol_sample)

    running_tot = 0
    d = 1
    for i in range(rows):
        for j in range(cols):
            well_name = "%s%d" % (chr(ord('A') + i), j+1)
            # Machine will round, so just give it enough info to do the
            # correct rounding.
            val = "%.2f" % pool_vols[i][j]

            # test to see if we will exceed total vol per well
            if running_tot + pool_vols[i][j] > max_vol_per_well:
                d += 1
                running_tot = pool_vols[i][j]
            else:
                running_tot += pool_vols[i][j]

            dest = "%s%d" % (chr(ord('A') +
                                 int(np.floor(d/dest_plate_shape[0]))),
                             (d % dest_plate_shape[1]))

            contents.append(
                ",".join(['1', '384LDV_AQ_B2_HT', well_name, "",
                          val, 'NormalizedDNA', dest]))

    return "\n".join(contents)


def plot_plate_vals(dataset, color_map='YlGnBu', annot_str=None,
                    annot_fmt='.5s'):
    """
    Plots values in a plate format. Returns a heatmap in the shape of the
    plate, with bar graphs aligned to the rows and columns showing the mean and
    spread of each row and column, and a histogram showing the distribution of
    values.

    Optionally can plot an array of names or other annotations on top of the
    heatmap.

    Parameters
    ----------
    dataset: 2D array of numeric
        data to plot
    color_map: str
        matplotlib color map name for heatmap
    annot_str: 2D array of str
        values to write over heatmap values to annotate wells
    annot_fmt: str
        string formatting values for annotations. Defaults to first 5 char per
        well.

    Returns
    -------
    """
    plt.figure(figsize=(20, 20))

    with sns.axes_style("white"):
        ax1 = plt.subplot2grid((40, 20), (20, 0), colspan=18, rowspan=18)
        ax1.xaxis.tick_top()
        if annot_str is None:
            sns.heatmap(dataset,
                        ax=ax1,
                        xticklabels=[x + 1 for x in range(dataset.shape[1])],
                        yticklabels=list(string.ascii_uppercase)[
                            0:dataset.shape[0]],
                        # square = True,
                        annot=True,
                        fmt='.0f',
                        cmap=color_map,
                        cbar=False)
        else:
            sns.heatmap(dataset,
                        ax=ax1,
                        xticklabels=[x + 1 for x in range(dataset.shape[1])],
                        yticklabels=list(string.ascii_uppercase)[
                            0:dataset.shape[0]],
                        # square = True,
                        annot=annot_str,
                        fmt=annot_fmt,
                        cmap=color_map,
                        cbar=False)

    with sns.axes_style("white"):
        ax2 = plt.subplot2grid((40, 20), (38, 0), colspan=18, rowspan=2)
        ax3 = plt.subplot2grid((40, 20), (20, 18), colspan=2, rowspan=18)
        sns.despine()
        sns.barplot(data=dataset, orient='v', ax=ax2, color='grey')
        sns.barplot(data=dataset.transpose(), orient='h', ax=ax3,
                    color='grey')
        ax2.set(xticklabels=[], yticklabels=[])
        ax3.set(xticklabels=[], yticklabels=[])

    with sns.axes_style():
        ax4 = plt.subplot2grid((40, 20), (0, 0), colspan=18, rowspan=18)
        sns.distplot(dataset.flatten()[~np.isnan(dataset.flatten())], ax=ax4,
                     bins=20)

    return


def make_2D_array(qpcr, data_col='Cp', well_col='Pos', rows=16, cols=24):
    """
    Pulls a column of data out of a dataframe and puts into array format
    based on well IDs in another column

    Parameters
    ----------
    qpcr: Pandas DataFrame
        dataframe from which to pull values
    data_col: str
        name of column with data
    well_col: str
        name of column with well IDs, in 'A1,B12' format
    rows: int
        number of rows in array to return
    cols: int
        number of cols in array to return

    Returns
    -------
    """
    # initialize empty Cp array
    cp_array = np.empty((rows, cols), dtype=object)

    # fill Cp array with the post-cleaned values from the right half of the
    # plate
    for record in qpcr.iterrows():
        row = ord(str.upper(record[1][well_col][0])) - ord('A')
        col = int(record[1][well_col][1:]) - 1
        cp_array[row, col] = record[1][data_col]

    return cp_array


def combine_dfs(qpcr_df, dna_picklist, index_picklist):
    """
    Combines information from the three dataframes into a single frame suitable
    for plotting

    Parameters
    ----------
    qpcr_df: Pandas DataFrame
        df from qpcr data import. Expects cols ['Pos','Cp']
    dna_picklist: Pandas DataFrame
        df from DNA picklist import. Expects cols
        ['Destination Well', 'Concentration', 'Transfer Volume']
    index_picklist: Pandas DataFrame
        df from index addition picklist import. Expects cols
        ['Destination Well','Plate','Sample Name',
         'Counter','Primer','Source Well','Index']

    Returns
    -------
    combined_df: Pandas DataFrame
        new DataFrame with the relevant columns
    """
    combined_df = pd.DataFrame({'Well': qpcr_df['Pos'],
                                'Cp': qpcr_df['Cp']})

    combined_df.set_index('Well', inplace=True)

    b = dna_picklist.loc[dna_picklist['Source Plate Name'] != 'water',
                         ].set_index('Destination Well')
    c = index_picklist.loc[index_picklist['Source Plate Name'] ==
                           'i7 Source Plate', ].set_index('Destination Well')
    d = index_picklist.loc[index_picklist['Source Plate Name'] ==
                           'i5 Source Plate', ].set_index('Destination Well')

    # Add DNA conc columns
    combined_df['DNA Concentration'] = b['Concentration']
    combined_df['DNA Transfer Volume'] = b['Transfer Volume']

    # Add Index columns
    combined_df['Sample Name'] = c['Sample Name']
    combined_df['Plate'] = c['Plate']
    combined_df['Counter'] = d['Counter']
    combined_df['Source Well i7'] = c['Source Well']
    combined_df['Index i7'] = c['Index']
    combined_df['Primer i7'] = c['Primer']
    combined_df['Source Well i5'] = d['Source Well']
    combined_df['Index i5'] = d['Index']
    combined_df['Primer i5'] = d['Primer']

    combined_df.reset_index(inplace=True)

    return combined_df


def parse_dna_conc_csv(fp):
    dna_df = pd.read_excel(fp, skiprows=4, parse_cols=[1, 2, 3, 4, 5])

    dna_df = dna_df.loc[list(range(384)), ]

    dna_df['pico_conc'] = pd.to_numeric(
        dna_df['[Concentration]'], errors='Coerce')
    return dna_df


def add_dna_conc(combined_df, dna_df):
    new_df = combined_df.set_index('Well')

    new_df['pico_conc'] = dna_df.set_index('Well')['pico_conc']

    new_df.reset_index(inplace=True)

    return new_df


def compute_pico_concentration(dna_vals, size=400):
    """Computes molar concentration of libraries from library DNA concentration
    values.

    Returns a 2D array of calculated concentrations, in nanomolar units

    Parameters
    ----------
    dna_vals : numpy array of float
        The DNA concentration in ng/uL
    size : int
        The average library molecule size in bp

    Returns
    -------
    np.array of floats
        A 2D array of floats
    """
    lib_concentration = (dna_vals / (660 * float(size))) * 10**6

    return lib_concentration


def bcl_scrub_name(name):
    """Modifies a sample name to be BCL2fastq compatible

    Parameters
    ----------
    name : str
        the sample name

    Returns
    -------
    str
        the sample name, formatted for bcl2fastq
    """

    return re.sub(r'[^0-9a-zA-Z\-\_]+', '_', name)


def rc(seq):
    """
    from http://stackoverflow.com/a/25189185/7146785
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    rev_seq = "".join(complement.get(base, base) for base in reversed(seq))

    return rev_seq


def sequencer_i5_index(sequencer, indices):

    if sequencer in REVCOMP_SEQUENCERS:
        print('%s: i5 barcodes are output as reverse compliments' % sequencer)
        return [rc(x) for x in indices]
    elif sequencer in OTHER_SEQUENCERS:
        print('%s: i5 barcodes are output in standard direction' % sequencer)
        return indices
    else:
        raise ValueError(('Your indicated sequencer [%s] is not recognized.\n'
                          'Recognized sequencers are: \n %s') %
                         (sequencer,
                          ', '.join(REVCOMP_SEQUENCERS + OTHER_SEQUENCERS)))


def reformat_interleaved_to_columns(wells):
    """
    converts condensed 96-to-384 plates in this format:

    plate1 | plate2
    ---------------
    plate3 | plate4

    to this format:

    plate1 | plate2 | plate3 | plate4

    where samples for each of the constituent plates are packed into contiguous
    columns of the 384 well plate.

    This is useful when doing a 192 well plate in order to save Mosquito tips /
    time

    Parameters
    ----------
    wells: array-like of str
        the sample source wells

    Returns
    -------
    new_wells: np array of str
        then new well locations in matching array positions
    """

    wells = np.array(wells)
    new_wells = np.empty(np.shape(wells), dtype='object')

    for i, owell in np.ndenumerate(wells):
        row = ord(str(owell[0]).upper()) - 65
        col = int(owell[1:]) - 1

        # ROWS
        # roffset = ROW % 2
        # row = ROW - roffset + floor(COL / 12)

        roffset = row % 2
        nrow = int(row - roffset + np.floor(col / 12))

        # COLS
        # coffset = COL % 2 + (ROW % 2) * 2
        # col = coffset * 6 + (col / 2) % 6

        coffset = col % 2 + (row % 2) * 2
        ncol = int(coffset * 6 + (col / 2) % 6)

        nwell = '%s%s' % (chr(nrow + 65), ncol + 1)

        new_wells[i] = nwell

    return new_wells


def merge_read_counts(plate_df, counts_df, reads_column_name='Filtered Reads'):
    """ Merges reads counts from FASTQC report or per_sample_FASTQ summary
    :param plate_df: A DataFrame containing the growing plate dataframe.
    :param counts_df: A DataFrame containing the counts.
    :param reads_column_name: A string column header for
    the merged read counts column.
    :return: A DataFrame containing the read counts
    from the input file in a column
    with the reads_column_name.
    """

    # Logic for multiple input_file format support
    if 'Category' in counts_df.columns:
        sample_column = 'Category'
        file_type = 'FastQC'
    elif 'filename' in counts_df.columns:
        sample_column = 'filename'
        file_type = 'per_sample_FASTQ'
    else:
        raise Exception("Unsupported input file type.")

    # Parse table to find sample names, and sum forward and rev reads.
    for i in counts_df.index:
        filename = counts_df.loc[i, sample_column]
        match = re.match(r'^(.*)_S\d+_L00\d', filename)
        if not match:
            raise LookupError(f'id not found in {filename}')
        sample_id = match.group(1)
        counts_df.loc[i, 'Sample'] = sample_id
    counts_df = counts_df.groupby('Sample').sum(numeric_only=True)

    # Logic for multiple input_file format support
    if file_type == 'FastQC':
        counts_df[reads_column_name] = counts_df['Unique Reads'] + \
                                       counts_df['Duplicate Reads']
    elif file_type == 'per_sample_FASTQ':
        counts_df.rename(columns={'reads': reads_column_name},
                         inplace=True)

    # Merge reads with plate_df
    to_merge = counts_df[[reads_column_name]]
    plate_df_w_reads = plate_df.merge(to_merge, left_on='Sample',
                                      right_on='Sample', how='left')

    return(plate_df_w_reads)


def read_survival(reads, label='Remaining', rmin=0, rmax=10 ** 6, rsteps=100):
    """
    Calculates a sample_retention curve from a read distribution
    :param reads: A DataFrame Series of a distribution of reads.
    :param label: A string label for the column that counts samples retained.
    :param rmin: An integer that counts the minimum number of reads
    for the sample_retention curve
    :param rmax: An integer that counts the maximum number of reads
    for the sample_retention curve
    :param rsteps: An integer that counts the steps between
    rmin and rmax to compute the sample_retention.
    :return: A DataFrame containing the sample_retention data
    for a survival curve.
    with the reads_column_name.
    """
    rstep = int((rmax - rmin) / rsteps)
    steps = list(range(rmin, rmax, rstep))
    remaining = np.zeros(rsteps)
    for i, t in enumerate(steps):
        remaining[i] = np.greater_equal(reads, t).sum()

    remaining_df = pd.DataFrame(remaining,
                                index=steps, columns=[label])
    return remaining_df


def linear_transform(input_values, output_min=100, output_max=1000):
    """
    Calculates a linear transformation between a range of
    input_values and output_min and output_max parameters
    :param input_values: A DataFrame Series of a distribution of values.
    :param output_min: An integer that specifies the minimum value
    for the transformed output
    :param output_max: An integer that specifies the maximum value
    for the transformed output
    :return: An array of linearly transformed values.
    """

    input_range = input_values.max() - input_values.min()
    input_min = input_values.min()

    diff1 = input_values - input_min
    diff2 = output_max - output_min
    transformed_values = (((diff1))*(diff2)/input_range) + output_min
    return(transformed_values)


def calculate_iseqnorm_pooling_volumes(plate_df,
                                       normalization_column='Filtered Reads',
                                       dynamic_range=100, lower_bound=1,
                                       min_pool_vol=100, max_pool_vol=1000):
    """
    Calculates optimal iseq normalization pooling volumes
    from a given set of parameters.
    :param plate_df: A DataFrame with ----need to add description----
    :param normalization_column: A string indicating the column name for a
    distribution of read counts used for normalization
    :param dynamic_range: The ratio between the read counts of the
    lower and high end of a clipped distribution of
    Loading Factors for normalization.
    Used to shape normalized distribution
    :param lower_bound: An integer that specifies the lower end of the clipped
    distribution of Loading Factors for normalization.
    Used to shape normalized distribution
    :param min_pool_vol: An integer that specifies the
    minimum pooling volume (nL) per sample for acoustic droplet ejection.
    :param min_pool_vol: An integer that specifies the
    maximum pooling volume (nL) per sample for acoustic droplet ejection.
    :return: A DataFrame with iSeqnorm pooling volumes.
    """
    try:
        norm = plate_df[normalization_column]
        prop = plate_df['proportion']
        plate_df['proportion'] = norm / (norm.sum())
        plate_df['LoadingFactor'] = max(prop) / prop
        plate_df.loc[plate_df['LoadingFactor'].isnull(),
                     'LoadingFactor'] = plate_df['LoadingFactor'].max()
    except ValueError:
        raise ValueError(f'Function expects dataframe with merged, \
                         per_sample read counts [\'{normalization_column}\']')

    # Calculate proportion of samples not normalized (unnorm_proportion)
    n = plate_df.shape[0]
    upper_bound = lower_bound * dynamic_range
    lf = plate_df['LoadingFactor']
    unnorm_proportion = sum((lf < lower_bound) | (lf > upper_bound)) / n

    # Transforming LoadingFactors into Pooling Volumes
    # with linear transformation
    plate_df['LoadingFactor'] = np.clip(plate_df['LoadingFactor'],
                                        lower_bound, upper_bound)
    norm_name = 'iSeq normpool volume'
    plate_df[norm_name] = linear_transform(plate_df['LoadingFactor'],
                                           output_min=min_pool_vol,
                                           output_max=max_pool_vol)

    # Underpooling BLANKS
    nblanks = plate_df.loc[plate_df['Blank']].shape[0]
    if (nblanks == 0):
        warnings.warn("There are no BLANKS in this plate")
    else:
        plate_df.loc[plate_df['Blank'], 'iSeq normpool volume'] = \
            plate_df.loc[plate_df['Blank'], 'iSeq normpool volume'] / 5

    # Plotting
    f, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(5, 8), sharex=True)

    ax1.set_title(f'{unnorm_proportion*100:.2f}% of sample set not normalized')
    sns.scatterplot(x=plate_df[normalization_column],
                    y=plate_df['iSeq normpool volume'],
                    alpha=0.5, ax=ax1)
    plt.xscale('log')
    sns.histplot(plate_df[normalization_column], ax=ax2)
    plt.tight_layout()

    return(plate_df)


def estimate_read_depth(plate_df,
                        estimated_total_output=4e9,
                        on_target_column='Filtered Reads',
                        off_target_column='Raw Reads'):
    """
    Builds a figure to estimate read depth per sample when
    given a DataFrame with iSeq normalization pooling volumes.
    :param plate_df: A DataFrame containing the growing plate dataframe.
    :param estimated_total_output: An integer (reads) that
    estimates the total output of the projected sequencing run.
    :param on_target_column: A string formated column name specifying which
    column was used for normalization, and as the numerator in the proportions.
    :param off_target_column: A string formated column name specifying
    which column to use as denominator in the proportions, typically Raw Reads.
    :return: A figure with proportions of reads on and off target and an
    estimate of average reads on target per sample.
    """

    plate_df['projected_reads'] = \
        plate_df[off_target_column] * plate_df['LoadingFactor']

    plate_df['projected_proportion'] = \
        plate_df['projected_reads'] / (plate_df['projected_reads'].sum())

    plate_df['projected_HO_reads'] = \
        plate_df['projected_proportion'] * estimated_total_output

    plate_df.sort_values(by='projected_HO_reads',
                         ascending=False, inplace=True)

    plate_df['on_target_proportion'] = \
        plate_df[on_target_column] / plate_df[off_target_column]

    plate_df['projected_off_target_reads'] = \
        plate_df['projected_HO_reads'] * (1 - plate_df['on_target_proportion'])

    plate_df['projected_on_target_reads'] = \
        plate_df['projected_HO_reads'] * plate_df['on_target_proportion']

    # PLOT
    plt.subplots(figsize=(11, 6))
    plot_df = plate_df.loc[~plate_df['projected_HO_reads'].isnull()]
    plot_df = plot_df.reset_index()
    plt.bar(range(plot_df.shape[0]),
            plot_df['projected_on_target_reads'], width=1, color='r')
    plt.bar(range(plot_df.shape[0]),
            plot_df['projected_off_target_reads'],
            bottom=plot_df['projected_on_target_reads'],
            width=1, align='center', color='gray')

    plt.title('Average reads on target per sample: ' +
              "{:,}".format(int(plate_df['projected_on_target_reads'].mean())))
    plt.show()

    return(plate_df)
