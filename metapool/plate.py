import click
from math import ceil
from datetime import datetime
import numpy as np
import pandas as pd
import warnings
from scipy.stats import zscore
from sklearn.linear_model import LogisticRegression


EXPECTED_COLUMNS = {
    'Plate Position', 'Primer Plate #', 'Plating', 'Extraction Kit Lot',
    'Extraction Robot', 'TM1000 8 Tool', 'Primer Date', 'MasterMix Lot',
    'Water Lot', 'Processing Robot', 'Sample Plate', 'Project_Name',
    'Original Name'}


class Message(object):
    _color = None

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return '%s: %s' % (self.__class__.__name__, self.message)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.message == other.message and self._color == other._color:
                return True

        return False

    def echo(self):
        prefix, suffix = str(self).split(': ', maxsplit=1)
        click.echo(click.style(prefix + ': ', fg=self._color) + suffix)


class ErrorMessage(Message):
    _color = 'red'


class WarningMessage(Message):
    _color = 'yellow'


def validate_plate_metadata(metadata):
    """Validates a list of plating metadata and outputs error/warning messages

    Parameters
    ----------
    metadata: list of dict
        A list of dictionaries where each dictionary represents plating
        metadata.

    Returns
    -------
    pd.DataFrame or None
        None if there's errors in the metadata.
        pd.DataFrame with the dictionaries as rows and the keys as columns.
    """
    messages = []

    if len(metadata) > 4:
        ErrorMessage("There are more than 4 plates, cannot continue "
                     "validation").echo()
        return

    # used to keep track of attributes to validate across plates
    context = {
        'primers': [],
        'names': [],
        'positions': []
    }

    no_errors = True
    for i, plate in enumerate(metadata):
        messages, context = _validate_plate(plate, context)

        if messages:
            print('Messages for Plate %s ' % plate['Plate Position'])
            for message in messages:
                message.echo()

                # all it takes is one error message in one plate
                if isinstance(message, ErrorMessage):
                    no_errors = False

    # prevent completion of other steps until these issues are fixed
    if no_errors:
        metadata = pd.DataFrame(metadata)
    else:
        metadata = None

    return metadata


def _validate_plate(plate_metadata, context):
    messages = []

    # 1. The first plate cannot be empty, but others can
    if plate_metadata == {}:
        if len(context['names']) == 0:
            messages.append(ErrorMessage("Can't leave the first plate empty"))
        else:
            messages.append(WarningMessage("This plate has no metadata"))

        # terminate the execution here, the rest of the code assumes
        # there's metadata to validate
        return messages, context

    observed = set(plate_metadata.keys())

    # 2. All columns are exactly present, no more no less
    extra = observed - EXPECTED_COLUMNS
    if extra:
        messages.append(
            ErrorMessage('The following columns are not needed: %s'
                         % ', '.join(extra)))
    missing = EXPECTED_COLUMNS - observed
    if missing:
        messages.append(ErrorMessage('The following columns are missing: %s' %
                                     ', '.join(missing)))

    # 3. Plate positions can only be from 1-4
    if plate_metadata['Plate Position'] not in {'1', '2', '3', '4'}:
        messages.append(
            ErrorMessage("Only the values '1', '2', '3' and '4' are allowed in"
                         " the 'Plate Position' field, you entered: "
                         "%s" % plate_metadata['Plate Position']))

    # 4. Are primer plates repeated?
    if plate_metadata['Primer Plate #'] in context['primers']:
        messages.append(ErrorMessage('The Primer Plate "%s" is repeated' %
                                     plate_metadata['Primer Plate #']))
    context['primers'].append(plate_metadata['Primer Plate #'])

    # 5. Are primer plates within the allowed values?
    # The EMP has 10 primer plates but in practice we only use 1-8. So show a
    # warning if the # is 9 or 10 but show an error if the plate is not between
    # 1 and 10.
    if plate_metadata['Primer Plate #'] not in {str(i) for i in range(1, 11)}:
        messages.append(ErrorMessage(('The Primer Plate # "%s" is not between '
                                     '1-10') %
                                     plate_metadata['Primer Plate #']))
    elif plate_metadata['Primer Plate #'] in {'9', '10'}:
        messages.append(WarningMessage(('Primer Plate # "%s" is unusual, '
                                       'please verify this value is correct') %
                                       plate_metadata['Primer Plate #']))

    # 6. Are plate names repeated?
    if plate_metadata['Sample Plate'] in context['names']:
        messages.append(ErrorMessage('The plate name "%s" is repeated' %
                                     plate_metadata['Sample Plate']))
    context['names'].append(plate_metadata['Sample Plate'])

    # 7. Spaces in plate names are not recommended
    if ' ' in plate_metadata['Sample Plate']:
        messages.append(WarningMessage("Spaces are not recommended in the "
                                       "Sample Plate field"))

    # 8. Are positions repeated?
    if plate_metadata['Plate Position'] in context['positions']:
        messages.append(ErrorMessage('The plate position "%s" is repeated' %
                                     plate_metadata['Plate Position']))
    context['positions'].append(plate_metadata['Plate Position'])

    # 9. Check the primer date is not in the future and that it can
    try:
        old = datetime.strptime(plate_metadata['Primer Date'],
                                '%Y-%m-%d').date()
    except ValueError:
        messages.append(ErrorMessage('Date format is invalid should be '
                                     'YYYY-mm-dd'))
    else:
        if old > datetime.now().date():
            messages.append(WarningMessage("The Primer Date is in the future"))

    # 10. Check there's only ASCII characters
    for key, value in plate_metadata.items():
        for c in value:
            if ord(c) > 127:
                messages.append(ErrorMessage(("The value for '%s' has "
                                              "non-ASCII characters") % key))
                break

    return messages, context


def _well_to_row_and_col(well):
    return ord(well[0].upper()) - 64, int(well[1:])


def _decompress_well(well):
    """Returns a 96 well plate ID from a compressed 384 well plate ID"""

    # convert the row name into a number starting at 1 for A, 2 for B, ...
    row, col = _well_to_row_and_col(well)

    if row != 1:
        row = ceil(row/2)
    if col != 1:
        col = ceil(col/2)

    return chr(64 + row) + str(col)


def _plate_position(well):
    """Returns a 96 well plate position from a compressed 394 well plate ID"""
    row, col = _well_to_row_and_col(well)

    if row % 2 == 0 and col % 2 == 0:
        return '4'
    elif row % 2 == 0 and col % 2 == 1:
        return '3'
    elif row % 2 == 1 and col % 2 == 0:
        return '2'
    elif row % 2 == 1 and col % 2 == 1:
        return '1'


def dilute_gDNA(plate_df, threshold=15):
    """
        Reads a plate_df and returns plate_df w/1:10 diluted gDNA samples.
        :param plate_df: processing-plate Pandas dataframe
        :param threshold: Upper gDNA concentration threshold.
        :return: Pandas DataFrame
               If dilution was needed, gDNA concentrations will be adjusted.
               'Diluted' column will be added and contents set to 'True'.
    """
    if 'Diluted' in plate_df.columns:
        # function has already been applied to data
        warnings.warn('Dilution operation was already performed')
        return plate_df

    # return a copy of the input data. Do not overwrite the input data by
    # default.
    df = plate_df.copy()
    df['Diluted'] = True

    df['Project Plate'] = df['Project Plate'].astype(str) + '_diluted'
    df['Compressed Plate Name'] = df['Compressed Plate Name'] + '_dilution'
    df['Sample DNA Concentration'] = df['Sample DNA Concentration'] / 10

    # Picking diluted samples for normalization
    diluted = (df.loc[df['Sample DNA Concentration'] > (threshold / 10)])
    undiluted = plate_df.loc[~plate_df['Sample'].isin(diluted['Sample'])]

    return pd.concat([diluted, undiluted]).sort_index()


def requires_dilution(plate_df, threshold=15, tolerance=0.05):
    """
        Reads a plate_df and determines whether or not it needs diluting.
        :param plate_df: processing-plate Pandas dataframe
        :param threshold: Upper gDNA concentration threshold.
        :param tolerance: Percentage (0.0-1.0) allowed outside threshold.
        :return: True if plate requires dilution. False otherwise.
        """
    # generate a Series of Booleans, based on whether each concentration falls
    # within the threshold or not.
    within_threshold = plate_df['Sample DNA Concentration'] > threshold
    pct_of_total = plate_df.loc[within_threshold].shape[0] / plate_df.shape[0]

    return pct_of_total > tolerance


def find_threshold(concentrations, labels):
    """
    Returns value that best separates concentration data into labels.
    :param concentrations: Pandas series w/concentrations.
    :param labels: Pandas series w/categorical labels.
    :return: Threshold value.
    """
    # Training model
    log_reg = LogisticRegression()
    log_reg.fit(np.array(concentrations).reshape(-1, 1), np.array(labels))

    # Making label predictions
    predictions = pd.DataFrame(
        [log_reg.predict(np.array(concentrations).reshape(-1, 1)),
         np.array(labels)]).transpose().rename(
        columns={0: 'predicted', 1: 'ground_truth'})

    # Appending concentration values to predictions dataframe and sorting.
    predictions['concentration'] = concentrations
    predictions = predictions.sort_values('concentration').reset_index()

    for i in predictions.index:
        # Finds the first sample predicted as non-blank.
        if not predictions.loc[i, 'predicted']:
            # Returns the concentation of the first sample predicted as
            # non-blank. This concentration is optimal boundary between
            # classes passed as labels.
            return predictions.loc[i, 'concentration']


def autopool(plate_df, method='norm', pool_failures='low', automate=True,
             offset=0.01, min_pool=100, total_vol=100, floor_conc=10,
             min_conc=0, total_nmol=0.0040):
    """
    reads a plate_df and calculates pooling volumes based on parameters.
    :param plate_df: processing-plate Pandas dataframe.
    :param method: 'norm' (normalized pooling) or 'evp' (equal volume pooling).
    :param pool_failures: Pool failed samples/blanks w/either 'low' or 'high'
                          volumes. Placing them close to the left or right
                          edge of the pooling volume distribution,
                          respectively.
    :param automate: estimate optimal parameters if True.
                     Relies on manually supplied params if Fase.
    :param offset:
    :param min_pool:
    :param total_vol:
    :param floor_conc: lowest conc. for which sample will be accurately pooled.
    :param min_conc: min conc. for a sample to be considered for pooling.
                     Set to 0 to pool all samples regardless.
    :param total_nmol:
    :return: DataFrame w/calculated pooling volumes (in nL).
    """
    if method not in ['evp', 'norm']:
        raise ErrorMessage('method must be either "evp" or "norm".')

    if pool_failures not in ['low', 'high']:
        raise ErrorMessage('pool_failures must be either "low" or "high".')

    sample_concs = plate_df['MiniPico Library Concentration']

    if method == 'evp':
        return _autopool_evp(plate_df, total_vol, sample_concs)

    return _autopool_norm(pool_failures, total_nmol, min_conc, sample_concs,
                          floor_conc, offset, automate, plate_df, min_pool)


def _autopool_evp(plate_df, total_vol, sample_concs):
    per_sample_vol = (total_vol / sample_concs.size) * 1000.0
    sample_vols = np.zeros(sample_concs.shape) + per_sample_vol
    plate_df['MiniPico Pooled Volume'] = sample_vols

    return plate_df


def _autopool_norm(pool_failures, total_nmol, min_conc, sample_concs,
                   floor_conc, offset, automate, plate_df, min_pool):
    sample_fracs = np.ones(sample_concs.shape) / sample_concs.size
    sample_fracs[sample_concs <= min_conc] = 0
    sample_fracs *= 1 / sample_fracs.sum()

    if pool_failures == 'low':
        use_this = sample_concs
        sort_ascending = True
    else:
        floored_concs = sample_concs.copy()
        floored_concs[sample_concs < floor_conc] = floor_conc
        use_this = floored_concs
        sort_ascending = False

    sample_vols = (total_nmol * sample_fracs) / use_this
    sample_vols *= 10 ** 9
    sample_vols[sample_concs < floor_conc] = \
        sample_vols.sort_values(ascending=sort_ascending, ignore_index=True)[
            int(sample_vols.size * offset)]

    if automate is True:
        sample_vols_zscored = np.clip(zscore(sample_vols, nan_policy='omit'),
                                      -2, 2)
        plate_df['zscore'] = sample_vols_zscored
        zrange = plate_df['zscore'].max() - plate_df['zscore'].min()
        zmin = plate_df['zscore'].min()
        sample_vols_scaled = (((sample_vols_zscored - (zmin)) * (
                    1000 - min_pool)) / zrange) + min_pool
        plate_df['MiniPico Pooled Volume'] = sample_vols_scaled
    else:
        if pool_failures == 'low':
            sample_vols[sample_concs < floor_conc] = \
                sample_vols.sort_values(ignore_index=True)[
                    int(sample_vols.size * offset)]
            my_func = sample_vols.min()
        else:
            my_func = sample_vols.max()

        plate_df['MiniPico Pooled Volume'] = sample_vols
        plate_df.loc[plate_df['MiniPico Pooled Volume'].isnull(),
                     'MiniPico Pooled Volume'] = my_func

    return plate_df
