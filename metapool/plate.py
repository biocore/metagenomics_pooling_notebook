import click
from math import ceil
from datetime import datetime
import numpy as np
import pandas as pd
import warnings
from scipy.stats import zscore
from sklearn.linear_model import LogisticRegression
from collections import OrderedDict
from string import ascii_uppercase
from metapool.mp_strings import EXPT_DESIGN_DESC_KEY, PM_PROJECT_NAME_KEY, \
    PM_PROJECT_PLATE_KEY, PM_PROJECT_ABBREV_KEY, PM_COMPRESSED_PLATE_NAME_KEY

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
    if 'Diluted' in plate_df.columns:
        # function has already been applied to data
        warnings.warn('Dilution operation was already performed')
        return plate_df

    # return a copy of the input data. Do not overwrite the input data by
    # default.
    plate_df['extracted_gdna_concentration_ng_ul'] = \
        plate_df['Sample DNA Concentration'].copy()
    df = plate_df.copy()
    df['Diluted'] = True

    df[PM_PROJECT_PLATE_KEY] = \
        df[PM_PROJECT_PLATE_KEY].astype(str) + '_diluted'
    df[PM_COMPRESSED_PLATE_NAME_KEY] = \
        df[PM_COMPRESSED_PLATE_NAME_KEY] + '_dilution'
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


class PlateReplication:
    STATUS_EMPTY = 'empty'
    STATUS_SOURCE = 'source'
    STATUS_DESTINATION = 'destination'

    row_letters = list(ascii_uppercase[:16])

    # aka ['blue', 'green', 'red', 'yellow']
    quadrants = ['1', '2', '3', '4']

    def __init__(self, well_column_name):
        self.map_to_384 = {}

        for quadrant in PlateReplication.quadrants:
            self.map_to_384[quadrant] = self._get_quadrant(quadrant)

        if well_column_name is None:
            self.well_column_name = 'Library Well'
        else:
            self.well_column_name = well_column_name

        self._reset()

    def _reset(self):
        # Used for (re)initialization so that make_replicate() calls can
        # be made multiple times using the same object.

        # All quadrants begin at one as they represent a potential source
        # quadrant, even if they are empty. If quad 1 is replicated once,
        # counter '1' will become 2. If quad 2 is replicated three times,
        # counter '1' will become 4. If a quadrant becomes a destination then
        # this value will be overwritten. If empty, this value will be unused.
        self.rep_counters = {'1': 1, '2': 1, '3': 1, '4': 1}

        self.status = {'1': PlateReplication.STATUS_EMPTY,
                       '2': PlateReplication.STATUS_EMPTY,
                       '3': PlateReplication.STATUS_EMPTY,
                       '4': PlateReplication.STATUS_EMPTY}

        # this allows us to store output dataframes for each quadrant,
        # overwrite them on demand, and save final concatenation for when
        # we're done.
        self.data = {'1': None, '2': None, '3': None, '4': None}

    def _get_quadrant(self, quadrant):
        d = OrderedDict()
        for i in range(1, 9):
            row_96 = PlateReplication.row_letters[i - 1]
            if quadrant in ['1', '2']:
                row_384 = PlateReplication.row_letters[(2 * i) - 2]
            else:
                row_384 = PlateReplication.row_letters[(2 * i) - 1]

            for j in range(1, 13):
                col_96 = j
                if quadrant in ['1', '3']:
                    col_384 = (2 * j) - 1
                else:
                    col_384 = (2 * j)

                k = "%s%s" % (row_96, col_96)
                v = "%s%s" % (row_384, col_384)
                d[k] = v

        return d

    def get_384_well_location(self, well_96_id, quadrant):
        '''
        Translate a 96-well plate + a quadrant into a 384-well plate cell
        :param well_96_id: A 96-well plate ID
        :param quadrant: A quadrant of a 384-well plate e.g. '1', '2', '3', '4'
        :return: A 384-well plate ID
        '''
        quadrant = str(quadrant)
        if quadrant in self.map_to_384:
            if well_96_id in self.map_to_384[quadrant]:
                return self.map_to_384[quadrant][well_96_id]

    def get_96_well_location_and_quadrant(self, well_384_id):
        '''
        Translate a 384-well plate ID to a 96-well plate ID and a quadrant
        :param well_384_id: A 384-well plate ID
        :return: A tuple of (quadrant, 96-well plate ID)
        '''
        # not an optimal search but mean-search-space is 192 and won't grow.
        for quadrant in PlateReplication.quadrants:
            for k in self.map_to_384[quadrant]:
                if self.map_to_384[quadrant][k] == well_384_id:
                    return quadrant, k

    def _map_quadrants(self, src_quad, dst_quad):
        # since all values in self.d are generated in the same order by
        # the same function and saved as OrderedDicts, the 384-well locations
        # from any src quadrant will map directly to those of the dst quad.
        src = self.map_to_384[str(src_quad)].values()
        dst = self.map_to_384[str(dst_quad)].values()

        return dict(map(lambda i, j: (i, j), src, dst))

    def _get_all_384_locations(self, quadrant):
        return list(self.map_to_384[quadrant].values())

    def _get_quadrants(self, wells_384):
        results = []

        for i in range(1, 5):
            all_wells_in_quadrant = set(self._get_all_384_locations(str(i)))
            if set(wells_384) & all_wells_in_quadrant:
                results.append(str(i))

        return results

    def check_bounds_384(self, locations):
        '''
        Check if one or more 384-well plate IDs are valid.
        :param locations: A list of one or more locations
        :return: A list of valid locations.
        '''
        if not isinstance(locations, list):
            # if a single value was passed instead of a list, convert it
            # into a list of one.
            locations = [locations]

        results = []

        for location in locations:
            # ensure location is always a string
            location = str(location)

            row = location[0]
            col = location[1:]
            try:
                col = int(col)
            except ValueError:
                results.append(location)

            if row not in PlateReplication.row_letters:
                results.append(location)

            if col < 1 or col > 24:
                results.append(location)

        return results

    def _replicate(self, plate_384, src_quad, dst_quad, overwrite=False):
        if self.status[src_quad] != PlateReplication.STATUS_SOURCE:
            raise ValueError(f'Quadrant {src_quad} is not a source quadrant')

        if self.status[dst_quad] == PlateReplication.STATUS_SOURCE and \
                overwrite is False:
            raise ValueError(f'Quadrant {dst_quad} is a source quadrant')

        if self.status[dst_quad] == PlateReplication.STATUS_DESTINATION and \
                overwrite is False:
            raise ValueError(f'Quadrant {dst_quad} is already occupied with '
                             f'replicate samples')

        # if there are wells in the prospective destination quadrant that are
        # already occupied in plate_384, raise an Error.
        dst_wells = self._get_all_384_locations(str(dst_quad))
        occupied = plate_384.loc[plate_384['Well'].isin(dst_wells)].copy()

        if occupied.shape[0] > 0 and overwrite is False:
            raise ValueError(f'Quadrant {dst_quad} contains source samples')

        self.rep_counters[src_quad] += 1

        rows = []
        for src, dst in self._map_quadrants(src_quad, dst_quad).items():
            # rows are clipped one at a time rather than as a subset of Well
            # ids because the order of the subset returned is according to
            # their numeric index id, rather than the order passed through
            # .in().
            # row will be a df that is just one row long.
            row = plate_384.loc[plate_384['Well'] == src].copy()
            # reset the numeric index, otherwise the row will keep the index
            # from the old plate and row.loc[] will not modify the right row.
            row.reset_index(inplace=True, drop=True)

            if row.shape[0] > 1:
                raise ValueError(f'{src} matched more than one row in '
                                 'plate_map')

            # assume row.shape[0] is either 0 (no-match) or 1 (exact match)
            if row.shape[0] == 1:
                row.loc[0, self.well_column_name] = dst
                row.loc[0, 'original_sample_name'] = row.loc[0, 'Sample']
                # add dst as suffix to sample_name for uniqueness
                row.loc[0, 'Sample'] = str(row.loc[0, 'Sample']) + '.' + dst
                row.loc[0, 'replicate'] = str(self.rep_counters[src_quad])
                row.loc[0, 'contains_replicates'] = 'True'

                rows.append(row)

        self.status[dst_quad] = PlateReplication.STATUS_DESTINATION
        # overwrites any existing values for that quadrant if present.
        self.data[dst_quad] = pd.concat(rows, axis=0, ignore_index=True)

    def _populate_source(self, src_quad, plate_384):
        # note that the number of columns and their type differ in the
        # output from the input. This means that even a source-quadrant can't
        # be copied wholesale from the input document into self.data[src_quad].
        # It must be munged as well.
        rows = []

        for src in self._get_all_384_locations(src_quad):
            row = plate_384.loc[plate_384['Well'] == src].copy()
            # reset the numeric index, otherwise the row will keep the index
            # from the old plate and row.loc[] will not modify the right row.
            row.reset_index(inplace=True, drop=True)

            if row.shape[0] > 1:
                raise ValueError(f'{src} matched more than one row in '
                                 'plate_map')

            # assume row.shape[0] is either 0 (no-match) or 1 (exact match)
            if row.shape[0] == 1:
                row.loc[0, self.well_column_name] = src
                row.loc[0, 'original_sample_name'] = row.loc[0, 'Sample']
                # add dst as suffix to sample_name for uniqueness
                row.loc[0, 'Sample'] = str(row.loc[0, 'Sample']) + '.' + src
                # convert rep_counters to string so it doesn't become '1.0.'
                row.loc[0, 'replicate'] = str(self.rep_counters[src_quad])
                # contains_replicates is a document-wide boolean. It should
                # be True even in the source row columns.
                row.loc[0, 'contains_replicates'] = 'True'
                rows.append(row)

        self.status[src_quad] = PlateReplication.STATUS_SOURCE
        self.data[src_quad] = pd.concat(rows, axis=0, ignore_index=True)

    def make_replicates(self, plate_384, replicates=None, overwrite=False):
        '''
        Given a 384-well plate and replication orders, generate output.
        :param plate_384: A 384-well plate.
        :param replicates: A dict containing a source quadrant and a list of
        destinations e.g.: {1: [2, 3, 4]}
        :param overwrite: Allow overwriting of an occupied quadrant.
        :return:
        '''

        # re-initialize object state for current call.
        self._reset()

        if replicates is None:
            # This generates output equal to legacy no-replication case.
            # Useful for generating output even when no replications are
            # needed.
            result = plate_384.copy()
            result[self.well_column_name] = result['Well'].copy()
            result['contains_replicates'] = 'False'
            return result

        for key in replicates:
            if not isinstance(replicates[key], list):
                v = replicates[key]
                replicates[key] = [v]

        # discover which quads in plate_df contain samples, mark them as
        # source plates, and populate self.data for src_quad.
        for src_quad in self._get_quadrants(list(plate_384['Well'].copy())):
            self._populate_source(src_quad, plate_384)

        for src_quad in replicates.keys():
            # every source quadrant (src_quad) is going to have a list in
            # replicates w/one or more destination quadrants (dst_quad) to
            # replicate to. We need to process each one in turn.
            for dst_quad in replicates[src_quad]:
                # the values in replicates are likely to be integers, not
                # strings, as they are from the user.
                self._replicate(plate_384, str(src_quad), str(dst_quad),
                                overwrite)

        # Take the final output for each quad after all replications and
        # potential overwrites and concatenate them before returning the
        # result to the user. Reset dataframe index so that iTru index merging
        # doesn't fail on duplicate index integer.
        quads = [self.data[quad] for quad in self.data if
                 self.data[quad] is not None]
        result = pd.concat(quads, axis=0, ignore_index=True)
        result.reset_index(drop=True, inplace=True)

        return result
