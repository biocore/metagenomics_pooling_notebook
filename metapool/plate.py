import click
from math import ceil
from datetime import datetime
import numpy as np
import pandas as pd
import re
import warnings
from scipy.stats import zscore
from sklearn.linear_model import LogisticRegression
from collections import OrderedDict
from string import ascii_uppercase
from metapool.mp_strings import EXPT_DESIGN_DESC_KEY, PM_PROJECT_NAME_KEY, \
    PM_PROJECT_PLATE_KEY, PM_PROJECT_ABBREV_KEY, PM_SAMPLE_KEY, \
    PM_WELL_KEY, PM_COMPRESSED_PLATE_NAME_KEY, SAMPLE_DNA_CONC_KEY, \
    ORIG_NAME_KEY, CONTAINS_REPLICATES_KEY, \
    PM_WELL_ID_96_KEY, PM_WELL_ID_384_KEY
from metapool.util import convert_to_bool

EXPECTED_COLUMNS = {
    'Plate Position', 'Plate map file', 'Plate elution volume',
    'Primer Plate #', 'Plating', 'Extraction Kit Lot',
    'Extraction Robot', 'TM1000 8 Tool', 'Primer Date', 'MasterMix Lot',
    'Water Lot', 'Processing Robot', 'Sample Plate', PM_PROJECT_NAME_KEY,
    PM_PROJECT_ABBREV_KEY, 'Original Name', 'TM300 8 Tool', 'TM50 8 Tool',
    'TM10 8 Tool', 'run_date', 'instrument_model', 'center_project_name',
    EXPT_DESIGN_DESC_KEY}


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

    # 2. All expected columns are present. additional columns are now okay.
    extra = observed - EXPECTED_COLUMNS
    if extra:
        messages.append(
            WarningMessage('The following columns are not recognized and may '
                           'be misspelled column names: %s'
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


def _parse_and_validate_well_384(well):
    """Validate input 384-well plate well ID and parse into (row, column).

    Parameters
    ----------
    well : str
        Well identifier like 'A1', 'P24', 'A01', etc.

    Returns
    -------
    tuple
        (row_number, column_number) where row A=1, B=2, etc.

    Raises
    ------
    ValueError
        If well format is invalid or outside 384-well plate bounds.
    """

    VALID_384_WELL_ROWS = set('ABCDEFGHIJKLMNOP')
    MAX_384_WELL_COL = 24

    # Check type and basic format using regex
    # Format: single letter A-Pa-p followed by 1-2 digits
    if not isinstance(well, str) or not re.match(r'^[A-Pa-p]\d{1,2}$', well):
        raise ValueError(
            "Well must be a letter A-P followed by a number 1-24, "
            "e.g., A1 or A01")

    row_letter = well[0].upper()
    col_str = well[1:]

    # Validate row is A-P
    if row_letter not in VALID_384_WELL_ROWS:
        raise ValueError(
            f"Well row '{row_letter}' is invalid. "
            f"Row must be A-P for a 384-well plate.")

    # Parse and validate column
    col = int(col_str)
    if col < 1 or col > MAX_384_WELL_COL:
        raise ValueError(
            f"Well column {col} is invalid. "
            f"Column must be 1-24 for a 384-well plate.")

    row = ord(row_letter) - 64
    return (row, col)


def sort_by_interleaved_plates(df, well_column):
    """Sort DataFrame rows by interleaved 96-well plate order.

    Sorts wells in the order of four interleaved 96-well plates within a
    384-well plate. The order is: quadrant 1, then quadrant 2, then quadrant 3,
    then quadrant 4. Within each quadrant, wells are sorted by row first,
    then by column.

    Quadrant mapping (based on 384-well position):
    - Quadrant 1: Odd rows (A, C, E, ...) + Odd columns (1, 3, 5, ...)
    - Quadrant 2: Odd rows + Even columns (2, 4, 6, ...)
    - Quadrant 3: Even rows (B, D, F, ...) + Odd columns
    - Quadrant 4: Even rows + Even columns

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing well identifiers.
    well_column : str
        Name of the column containing 384-well IDs.

    Returns
    -------
    pd.DataFrame
        DataFrame sorted by quadrant, then row, then column.

    Raises
    ------
    ValueError
        If well_column doesn't exist or contains invalid well IDs.
    """
    if well_column not in df.columns:
        raise ValueError(
            f"Column '{well_column}' not found in DataFrame.")

    def _get_sort_key(well):
        row, col = _parse_and_validate_well_384(well)
        # not necessary to validate plate position output since we know the
        # well is a valid one if we got here
        quadrant = int(_plate_position(well))
        return (quadrant, row, col)

    sort_keys = df[well_column].apply(_get_sort_key)
    sorted_indices = sort_keys.sort_values().index
    return df.loc[sorted_indices].reset_index(drop=True)


def _decompress_well(well):
    """Returns a 96 well plate ID from a compressed 384 well plate ID"""

    # convert the row name into a number starting at 1 for A, 2 for B, ...
    row, col = _well_to_row_and_col(well)

    if row != 1:
        row = ceil(row/2)
    if col != 1:
        col = ceil(col/2)

    return chr(64 + row) + str(col)


def _validate_well_id_96(well):
    VALID_96_WELL_COLUMNS = {i for i in range(1, 13)}
    VALID_96_WELL_ROWS = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'}

    if well in [None, '']:
        return None
    try:
        row = well[0]
        col = well[1:]
    except IndexError:
        return None

    try:
        col = int(col)
    except ValueError:
        return None

    if col not in VALID_96_WELL_COLUMNS:
        return None

    if row.upper() not in VALID_96_WELL_ROWS:
        return None

    return (row, col)


def _plate_position(well):
    """Returns a 96 well plate position from a compressed 384 well plate ID"""
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

    temp_key = 'TEMP_CONC'
    a_df = plate_df.copy()
    diluted_mask = a_df[SAMPLE_DNA_CONC_KEY] > threshold
    a_df[temp_key] = a_df[SAMPLE_DNA_CONC_KEY] / 10
    result = record_gdna_dilution(a_df, diluted_mask, temp_key)
    result.drop(columns=[temp_key], inplace=True)
    return result


def record_gdna_dilution(plate_df, diluted_mask, diluted_conc_key):
    """
        Reads a plate_df and sets info for samples that need dilution.
        :param plate_df: processing-plate Pandas dataframe
        :param diluted_mask: Boolean mask for samples that need dilution.
        :param diluted_conc_key: Key for the diluted concentration column.
        :return: Pandas DataFrame
               Note: 'Diluted' column will be added and contents set to 'True'.
    """
    if 'Diluted' in plate_df.columns:
        # function has already been applied to data
        warnings.warn('Dilution operation was already performed')
        return plate_df

    # return a copy of the input data. Do not overwrite input data by default.
    df = plate_df.copy()
    df['Diluted'] = True
    df[PM_PROJECT_PLATE_KEY] = \
        df[PM_PROJECT_PLATE_KEY].map(str) + '_diluted'
    df[PM_COMPRESSED_PLATE_NAME_KEY] = \
        df[PM_COMPRESSED_PLATE_NAME_KEY] + '_dilution'
    df[SAMPLE_DNA_CONC_KEY] = df[diluted_conc_key]

    # Make output dataframe from the diluted and undiluted samples
    diluted = (df.loc[diluted_mask])
    undiluted = plate_df.loc[~diluted_mask]

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
    within_threshold = plate_df[SAMPLE_DNA_CONC_KEY] > threshold
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
             offset=0.01, min_pool=100, total_vol=190, floor_conc=10,
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
    :param total_vol: (default = 190µL)
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


def merge_plate_dfs(plate_df_one, plate_df_two, wells_col,
                    plate_df_one_name, plate_df_two_name):
    """
    Merges the dilution plate dataframe with the plate dataframe

    Parameters:
    plate_df_one : pandas DataFrame
        A DataFrame of plate info, containing at least the wells_col.  May or
        may not contain rows for all wells in the plate.
    plate_df_two : pandas DataFrame
        A DataFrame of plate info, containing at least the wells_col.  May or
        may not contain rows for all wells in the plate.
    wells_col : str
        The column indicating the well positions, for merging the DataFrames
    plate_df_one_name : str
        The name of the first DataFrame, used only in warning message
    plate_df_two_name : str
        The name of the second DataFrame, used only in warning message

    Returns:
    a_df : pandas DataFrame
        A new dataframe that combines the contents of the two plate dataframes.
    """

    # warn if the wells in the two dataframes do not match;
    # this isn't NECESSARILY an error, as in the case where the lab is
    # prepping only half a plate, but it is worth noting
    wells_one = set(plate_df_one[wells_col])
    wells_two = set(plate_df_two[wells_col])
    overlap_wells = wells_one.intersection(wells_two)
    if wells_one != wells_two:
        warn_msg = (
            f"The two DataFrames do not have the same set of wells. "
            f"DataFrame {plate_df_one_name} has {len(wells_one)} unique wells "
            f"and DataFrame {plate_df_two_name} has {len(wells_two)} unique "
            f"wells, with {len(overlap_wells)} shared unique wells. Only data "
            f"for wells that are shared in both plates will be returned."
        )

        warnings.warn(warn_msg, UserWarning)
    # end if merge col contents sets not the same

    return pd.merge(plate_df_one, plate_df_two, on=wells_col)


class _WellMapper96to384:
    """Mixin specifying interface for 96-well to 384-well mapping"""

    def get_384_well_location(self, well_96_id, well_96_unit_id):
        """Translate 96-well id + 96-well-plate unit id to 384-well id

        Parameters
        ----------
        well_96_id : str
            A 96-well plate ID such as "A1", "H12", etc.
        well_96_unit_id : str
            The identifier for the 96-well unit that the well_id_96 belongs
             to.  Could be e.g. a plate name, or a quadrant number.

        Returns
        -------
        str
            The 384-well location corresponding to the 96-well location on
            the specified plate.

        Raises
        -------
        NotImplementedError
            If the method is not implemented in the subclass.
        """
        raise NotImplementedError("Subclasses must implement this method")

    def _raise_missing_well_error(
            self, well_96_id, well_96_unit_id, well_96_unit_key):
        """Raise an error indicating that the 96-well location was not found

        Parameters
        ----------
        well_96_id : str
            A 96-well plate ID such as "A1", "H12", etc.
        well_96_unit_id : str
            The identifier for the 96-well unit that the well_id_96 belongs
             to.  Could be e.g. a plate name, or a quadrant number.
        """
        raise ValueError(
            f"{PM_WELL_ID_96_KEY} '{well_96_id}' not found in "
            f"{well_96_unit_key} '{well_96_unit_id}'")


class PlateReplication(_WellMapper96to384):
    STATUS_EMPTY = 'empty'
    STATUS_SOURCE = 'source'
    STATUS_DESTINATION = 'destination'

    row_letters = list(ascii_uppercase[:16])

    # aka ['blue', 'green', 'red', 'yellow']
    quadrants = ['1', '2', '3', '4']

    @classmethod
    def _get_96_to_384_mapping_for_quadrant(cls, quadrant):
        """
        Returns a dictionary keyed by 96-well plate location, with value being
        that well's 384-well plate location if placed in the specified
        quadrant using the "interleaved" mapping method.
        corresponding 384-well plate locations for the given quadrant.

        Parameters
        ----------
        quadrant: str
            The quadrant to map, one of '1', '2', '3', or '4'.

        Returns
        -------
        OrderedDict
            An OrderedDict where keys are 96-well plate locations (e.g., 'A1')
            and values are the corresponding 384-well plate locations if the
            96-well plate were placed in the specified quadrant using the
            "interleaved" mapping method.
        """
        if quadrant not in cls.quadrants:
            raise ValueError(f"Invalid quadrant '{quadrant}'. Must be one of "
                             f"{cls.quadrants}.")

        d = OrderedDict()
        for curr_row_num in range(1, 9):
            row_96 = cls.row_letters[curr_row_num - 1]
            if quadrant in ['1', '2']:
                row_384 = cls.row_letters[(2 * curr_row_num) - 2]
            else:
                row_384 = cls.row_letters[(2 * curr_row_num) - 1]

            for curr_col_num in range(1, 13):
                col_96 = curr_col_num
                if quadrant in ['1', '3']:
                    col_384 = (2 * curr_col_num) - 1
                else:
                    col_384 = (2 * curr_col_num)

                curr_96_well_id = "%s%s" % (row_96, col_96)
                curr_384_well_id = "%s%s" % (row_384, col_384)
                d[curr_96_well_id] = curr_384_well_id

        return d

    @staticmethod
    def df_contains_replicates(a_df, df_name):
        if a_df is None or CONTAINS_REPLICATES_KEY not in a_df:
            return False

        contains_rep_vals = a_df[CONTAINS_REPLICATES_KEY].tolist()

        contains_rep_bools = []
        illegal_vals = []
        for val in contains_rep_vals:
            bool_val = convert_to_bool(val)
            if bool_val is None:
                illegal_vals.append(val)
            else:
                contains_rep_bools.append(bool_val)
        # next val

        if illegal_vals:
            err_msg = (f"Values for '{CONTAINS_REPLICATES_KEY}' in "
                       f"'{df_name}' must be True or False. The following "
                       f"values are not valid: {', '.join(illegal_vals)}")
            raise ValueError(err_msg)

        # remove duplicates, so we can check if all values are the same
        contains_replicates = list(set(contains_rep_bools))

        # by convention, all projects described together are either going
        # to be True or False. If some entries are True while others are
        # False, we should raise an Error.
        if len(contains_replicates) > 1:
            raise ValueError(f"All projects in {df_name} "
                             f"must either contain replicates or not.")

        return contains_replicates[0]

    def __init__(self, dest_well_column_name):
        self._mapping_dicts_by_quadrant = {}

        for quadrant in PlateReplication.quadrants:
            self._mapping_dicts_by_quadrant[quadrant] = \
                self._get_96_to_384_mapping_for_quadrant(quadrant)

        if dest_well_column_name is None:
            self.dest_well_column_name = 'Library Well'
        else:
            self.dest_well_column_name = dest_well_column_name

        self._reset()

    def _reset(self):
        # Used for (re)initialization so that make_replicate() calls can
        # be made multiple times using the same object.

        # All quadrants begin at one as they represent a potential source
        # quadrant, even if they are empty. If quad 1 is replicated once,
        # counter '1' will become 2. If quad 2 is replicated three times,
        # counter '1' will become 4. If a quadrant becomes a destination then
        # this value will be overwritten. If empty, this value will be unused.
        self._rep_counters = {'1': 1, '2': 1, '3': 1, '4': 1}

        self._status = {'1': PlateReplication.STATUS_EMPTY,
                        '2': PlateReplication.STATUS_EMPTY,
                        '3': PlateReplication.STATUS_EMPTY,
                        '4': PlateReplication.STATUS_EMPTY}

        # this allows us to store output dataframes for each quadrant,
        # overwrite them on demand, and save final concatenation for when
        # we're done.
        self._data = {'1': None, '2': None, '3': None, '4': None}

    def _map_between_quadrants(self, src_quad, dst_quad):
        """
        Returns a mapping between wells in the src_quad and dst_quad
        (assuming "interleaved" mapping).
        For example, if src_quad is '1' and dst_quad is '2', then
        the mapping will include A1: A2, C1: C2, A3: A4, etc.

        Parameters
        ----------
        src_quad: str
            The source quadrant, one of '1', '2', '3', or '4'.
        dst_quad: str
            The destination quadrant, one of '1', '2', '3', or '4'.

        Returns
        -------
        dict
            A dictionary where keys are 384-well locations in the source
            quadrant and values are the corresponding 384-well locations in
            the destination quadrant.
        """
        # since the src quad well ids and dst quad well ids are generated by
        # the same function and saved in the same order, the 384-well locations
        # from any src quadrant will map directly to those of the dst quad.
        src_quad_384_well_ids = \
            self._mapping_dicts_by_quadrant[str(src_quad)].values()
        dst_quad_384_well_ids = \
            self._mapping_dicts_by_quadrant[str(dst_quad)].values()

        return dict(zip(src_quad_384_well_ids, dst_quad_384_well_ids))

    def _get_all_384_locations(self, quadrant):
        return list(self._mapping_dicts_by_quadrant[quadrant].values())

    def _get_quadrant_and_96_well_location(self, well_384_id):
        """
        Translate a 384-well well ID to a quadrant and a 96-well well ID
        :param well_384_id: A 384-well well ID
        :return: A tuple of (quadrant, 96-well well ID)
        """
        # not an optimal search but mean-search-space is 192 and won't grow.
        for quadrant in PlateReplication.quadrants:
            for k in self._mapping_dicts_by_quadrant[quadrant]:
                if self._mapping_dicts_by_quadrant[quadrant][k] == well_384_id:
                    return quadrant, k

    def _get_quadrants(self, wells_384):
        results = []

        for i in range(1, 5):
            all_wells_in_quadrant = set(self._get_all_384_locations(str(i)))
            if set(wells_384) & all_wells_in_quadrant:
                results.append(str(i))

        return results

    def _replicate(self, plate_384, src_quad, dst_quad, overwrite=False):
        if self._status[src_quad] != PlateReplication.STATUS_SOURCE:
            raise ValueError(f'Quadrant {src_quad} is not a source quadrant')

        if self._status[dst_quad] == PlateReplication.STATUS_SOURCE and \
                overwrite is False:
            raise ValueError(f'Quadrant {dst_quad} is a source quadrant')

        if self._status[dst_quad] == PlateReplication.STATUS_DESTINATION and \
                overwrite is False:
            raise ValueError(f'Quadrant {dst_quad} is already occupied with '
                             f'replicate samples')

        # if there are wells in the prospective destination quadrant that are
        # already occupied in plate_384, raise an Error.
        dst_wells = self._get_all_384_locations(str(dst_quad))
        occupied = plate_384.loc[plate_384[PM_WELL_KEY].isin(dst_wells)].copy()

        if occupied.shape[0] > 0 and overwrite is False:
            raise ValueError(f'Quadrant {dst_quad} contains source samples')

        self._rep_counters[src_quad] += 1

        rows = []
        for src_384_well_id, dst_384_well_id in (
                self._map_between_quadrants(src_quad, dst_quad).items()):
            # rows are clipped one at a time rather than as a subset of Well
            # ids because the order of the subset returned is according to
            # their numeric index id, rather than the order passed through
            # .in().
            row = self._make_replicate_row(
                plate_384, src_quad, src_384_well_id, dst_384_well_id)
            if row is not None:
                rows.append(row)

        self._status[dst_quad] = PlateReplication.STATUS_DESTINATION
        # overwrites any existing values for that quadrant if present.
        self._data[dst_quad] = pd.concat(rows, axis=0, ignore_index=True)

    def _populate_source(self, src_quad, plate_384):
        # note that the number of columns and their type differ in the
        # output from the input. This means that even a source-quadrant can't
        # be copied wholesale from the input document into self.data[src_quad].
        # It must be munged as well.
        rows = []

        for curr_src_384_well_id in self._get_all_384_locations(src_quad):
            row = self._make_replicate_row(
                plate_384, src_quad, curr_src_384_well_id)
            if row is not None:
                rows.append(row)

        self._status[src_quad] = PlateReplication.STATUS_SOURCE
        self._data[src_quad] = pd.concat(rows, axis=0, ignore_index=True)

    def _make_replicate_row(self, plate_df_384, src_quad,
                            src_384_well_id, dst_384_well_id=None):
        """
        Given a 384-well plate df and a source well ID, generate a row for the
        destination well ID. If the destination well ID is None, then the
        source well ID will be used as the destination well ID.

        Parameters
        ----------
        plate_df_384: pandas DataFrame
            A 384-well plate dataframe (one row per sample).
        src_quad: str
            The source quadrant, one of '1', '2', '3', or '4'.
        src_384_well_id: str
            The 384-well location of the sample in the source quadrant
            (e.g., '020').
        dst_384_well_id: str or None
            The 384-well location of the sample in the destination quadrant.
            If None, the source well ID will be used as the destination
            well ID.

        Returns
        -------
        pandas DataFrame or None
            A DataFrame containing the row for the destination well ID, or
            None if the source well ID does not match any row in the
            plate_df_384.

        Raises
        -------
        ValueError
            If the source well ID matches more than one row in plate_df_384.
        """

        well_mask = plate_df_384[PM_WELL_KEY] == src_384_well_id
        row = plate_df_384.loc[well_mask].copy()
        # reset the numeric index, otherwise the row will keep the index
        # from the old plate and row.loc[] will not modify the right row.
        row.reset_index(inplace=True, drop=True)

        if row.shape[0] > 1:
            raise ValueError(f'{src_384_well_id} matched more than one row in '
                             'plate_map')

        # assume row.shape[0] is either 0 (no-match) or 1 (exact match)
        if row.shape[0] == 1:
            relevant_well = dst_384_well_id if dst_384_well_id is not None \
                else src_384_well_id

            row.loc[0, self.dest_well_column_name] = relevant_well
            row.loc[0, ORIG_NAME_KEY] = row.loc[0, PM_SAMPLE_KEY]
            # add dst as suffix to sample_name for uniqueness
            row.loc[0, PM_SAMPLE_KEY] = \
                str(row.loc[0, PM_SAMPLE_KEY]) + '.' + relevant_well
            # convert _rep_counters to string so it doesn't become '1.0.'
            row.loc[0, 'replicate'] = str(self._rep_counters[src_quad])
        elif row.shape[0] == 0:
            row = None
        # endif

        return row

    def make_replicates(self, plate_df_384, replicates=None, overwrite=False):
        """
        Given a 384-well plate and replication orders, generate output.
        :param plate_df_384: A plate dataframe (one row per sample) for a
            384-well plate.
        :param replicates: A dict containing a source quadrant and a list of
        destinations e.g.: {1: [2, 3, 4]}
        :param overwrite: Allow overwriting of an occupied quadrant.
        :return: A Pandas DataFrame containing the output of the replication
        operation. The output will contain the same columns as the input
        plate_df_384, with the addition of the following columns:
            - well_column_name: The 384-well location of the sample.
            - CONTAINS_REPLICATES_KEY: all true if the plate contains
              replicates, all false if it does not.
        If replicates were made, the output will also contain
            - ORIG_NAME_KEY: The original sample name before replication.
            - replicate: The replicate number for the sample in the source
              quadrant.
        and the sample name will be suffixed with the 384-well
        location of the sample in the 384-well plate.
        """

        # re-initialize object state for current call.
        self._reset()

        if replicates is None:
            # This generates output equal to legacy no-replication case.
            # Useful for generating output even when no replications are
            # needed.
            result = plate_df_384.copy()
            result[self.dest_well_column_name] = result[PM_WELL_KEY].copy()
            result[CONTAINS_REPLICATES_KEY] = False
            return result

        for key in replicates:
            if not isinstance(replicates[key], list):
                v = replicates[key]
                replicates[key] = [v]

        # replication writes text into these columns. An input column that is
        # entirely empty reads as float64, which pandas 3 refuses text in, so
        # hold any non-text one as str, as a text read of it would be.
        plate_df_384 = plate_df_384.copy()
        for col in (self.dest_well_column_name, ORIG_NAME_KEY, 'replicate',
                    PM_SAMPLE_KEY):
            if col in plate_df_384.columns and \
                    not pd.api.types.is_string_dtype(plate_df_384[col]):
                plate_df_384[col] = plate_df_384[col].astype("str")

        # discover which quads in plate_df contain samples, mark them as
        # source plates, and populate self.data for src_quad.
        for src_quad in self._get_quadrants(
                list(plate_df_384[PM_WELL_KEY].copy())):
            self._populate_source(src_quad, plate_df_384)

        for src_quad in replicates.keys():
            # every source quadrant (src_quad) is going to have a list in
            # replicates w/one or more destination quadrants (dst_quad) to
            # replicate to. We need to process each one in turn.
            for dst_quad in replicates[src_quad]:
                # the values in replicates are likely to be integers, not
                # strings, as they are from the user.
                self._replicate(plate_df_384, str(src_quad), str(dst_quad),
                                overwrite)

        # Take the final output for each quad after all replications and
        # potential overwrites and concatenate them before returning the
        # result to the user. Reset dataframe index so that iTru index merging
        # doesn't fail on duplicate index integer.
        quads = [self._data[quad] for quad in self._data if
                 self._data[quad] is not None]
        result = pd.concat(quads, axis=0, ignore_index=True)
        # contains_replicates is a document-wide boolean (so not clear why it
        # is being stored in this dataframe, but ...). If there are any
        # replicates, this value should be True in ALL rows--even source rows
        result[CONTAINS_REPLICATES_KEY] = True
        result.reset_index(drop=True, inplace=True)

        return result

    def unmake_replicates(self, combined_df):
        """
        Splits 384-well df into list of 96-well dfs, one for each quadrant used

        Parameters
        ----------
        combined_df: pandas DataFrame
            A DataFrame containing 384-well plate data, with a column with the
            name given in the PlateReplication object's dest_well_column_name
            property.

        Returns
        -------
        per_replicate_dataframes: list of pandas DataFrame
            A list of DataFrames, each containing the records for one used
            quadrant of the 384-well plate. If not every quadrant was used,
            the list will be shorter than four.

        """
        quad_col = "quad"
        per_replicate_dataframes = []

        # check that specified de-replication column exists
        if self.dest_well_column_name not in combined_df.columns:
            raise ValueError(
                f"Column '{self.dest_well_column_name}' not found in "
                f"the input dataframe.")

        # convert each sample's 384 well location into a
        # 96-well location + quadrant. Since replication is performed at the
        # 96-well-plate level, this will identify which replicates belong in
        # which dataframe.
        combined_df[quad_col] = combined_df.apply(
            lambda row: self._get_quadrant_and_96_well_location(
                row[self.dest_well_column_name])[0],
            axis=1)

        for quad in sorted(combined_df[quad_col].unique()):
            # for each unique quadrant found, create a new dataframe that's a
            # subset containing only members of that quadrant. Delete the
            # temporary 'quad' column afterward and reset the index to an int
            # value starting at zero; the current-index will revert to a column
            # named 'sample_id'. Return the list of new dataframes.
            per_replicate_dataframes.append(
                combined_df[combined_df[quad_col] == quad].drop(
                    [quad_col], axis=1)
            )

        return per_replicate_dataframes

    def get_384_well_location(self, well_96_id, quadrant):
        """
        Translate a 96-well plate + a quadrant into a 384-well plate cell
        :param well_96_id: A 96-well plate ID
        :param quadrant: A quadrant of a 384-well plate e.g. '1', '2', '3', '4'
        :return: A 384-well plate ID
        """
        quadrant = str(quadrant)
        if quadrant in self._mapping_dicts_by_quadrant:
            if well_96_id in self._mapping_dicts_by_quadrant[quadrant]:
                return self._mapping_dicts_by_quadrant[quadrant][well_96_id]

        self._raise_missing_well_error(
            well_96_id, quadrant, "quadrant")


class PlateRemapper(_WellMapper96to384):
    def __init__(self, a_df):
        """ Initialize a PlateRemapper object

        Parameters
        ----------
        a_df : pd.DataFrame
            A pandas DataFrame containing
            the mapping between 96-well plate name and position and
            corresponding 384-well location. Must contain columns
            PM_PROJECT_PLATE_KEY, PM_WELL_ID_96_KEY (96-well position), and
            PM_WELL_ID_384_KEY (384-well position).  Note this can only handle
            a *single* 384-well plate.
        """

        # if any of the columns PM_PROJECT_PLATE_KEY, PM_WELL_ID_96_KEY, or
        # PM_WELL_ID_384_KEY are NOT in a_df, raise a ValueError saying which
        # columns are missing
        required_cols = \
            [PM_PROJECT_PLATE_KEY, PM_WELL_ID_96_KEY, PM_WELL_ID_384_KEY]
        missing_required_cols = set(required_cols) - set(a_df.columns)
        if len(missing_required_cols) > 0:
            raise ValueError(
                f"Input dataframe lacks required columns: "
                f"{missing_required_cols}")

        self._df = a_df.copy()

    def get_384_well_location(self, well_96_id, plate_name):
        """ Translate 96-well location + a plate name into 384-well location

        Parameters
        ----------
        well_96_id : str
            A 96-well plate ID such as "A1", "H12", etc.
        plate_name : str
            The name of the plate to search in the dataframe, e.g. "Plate_16"

        Returns
        -------
        str
            The 384-well location corresponding to the 96-well location on
            the specified plate according to the remapper's dataframe.
        """

        well_mask = (self._df[PM_PROJECT_PLATE_KEY] == plate_name) & \
                    (self._df[PM_WELL_ID_96_KEY] == well_96_id)

        if not well_mask.any():
            self._raise_missing_well_error(
                well_96_id, plate_name, PM_PROJECT_PLATE_KEY)

        # get a single value from the dataframe
        result = self._df.loc[well_mask, PM_WELL_ID_384_KEY].values[0]
        return result
