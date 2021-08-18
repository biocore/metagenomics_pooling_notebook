import click
import pandas as pd

from math import ceil
from datetime import datetime


class Message(object):
    def __init__(self, message):
        self._color = None
        self.message = message
        return self

    def __str__(self):
        return '%s: %s' % (self.__class__.__name__, self.message)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.message == other.message:
                return True
        else:
            return False

    def echo(self):
        prefix, suffix = str(self).split(': ', maxsplit=1)
        click.echo(click.style(prefix + ': ', fg=self._color) + suffix)


class ErrorMessage(Message):
    def __init__(self, message):
        super().__init__(message)
        self._color = 'red'


class WarningMessage(Message):
    def __init__(self, message):
        super().__init__(message)
        self._color = 'yellow'


def validate_plate_metadata(metadata):
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

    expected = {'Plate Position', 'Primer Plate #', 'Plating',
                'Extraction Kit Lot', 'Extraction Robot', 'TM1000 8 Tool',
                'Primer Date', 'MasterMix Lot', 'Water Lot',
                'Processing Robot', 'Sample Plate', 'Project_Name',
                'Original Name'}
    observed = set(plate_metadata.keys())

    # 2. All columns are exactly present, no more no less
    extra = observed - expected
    if extra:
        messages.append(
            ErrorMessage('The following columns are not needed: %s'
                         % ', '.join(extra)))
    missing = expected - observed
    if missing:
        messages.append(ErrorMessage('The following columns are missing: %s' %
                                     ', '.join(missing)))

    # 3. Are primer plates repeated?
    if plate_metadata['Primer Plate #'] in context['primers']:
        messages.append(ErrorMessage('The primer plate "%s" is repeated' %
                                     plate_metadata['Primer Plate #']))
    context['primers'].append(plate_metadata['Primer Plate #'])

    # 4. Are plate names repeated?
    if plate_metadata['Sample Plate'] in context['names']:
        messages.append(ErrorMessage('The plate name "%s" is repeated' %
                                     plate_metadata['Sample Plate']))
    context['names'].append(plate_metadata['Sample Plate'])

    # 5. Are positions repeated?
    if plate_metadata['Plate Position'] in context['positions']:
        messages.append(ErrorMessage('The plate position "%s" is repeated' %
                                     plate_metadata['Plate Position']))
    context['positions'].append(plate_metadata['Plate Position'])

    # 6. Check the primer date is not in the future and that it can
    try:
        old = datetime.strptime(plate_metadata['Primer Date'],
                                '%Y-%m-%d').date()
    except ValueError:
        messages.append(ErrorMessage('Date format is invalid should be '
                                     'YYYY-mm-dd'))
    else:
        if old > datetime.now().date():
            messages.append(WarningMessage("The Primer Date is in the future"))

    # 7. Check there's only ASCII characters
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
    """Returns a 96 well plate ID from a compressed 394 well plate ID"""

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
