import click
import pandas as pd

from math import ceil


class Message(object):
    def __init__(self, message):
        # TODO: Add a color to format output
        self.color = None
        self.message = message
        return self

    def __str__(self):
        return self.message

    def print(self):
        click.echo(click.style(str(self), fg=self.color))


class Error(Message):
    def __init__(self, message):
        super().__init__(message)
        self.color = 'red'

    def __str__(self):
        return 'Error: %s' % self.message


class Warning(Message):
    def __init__(self, message):
        super().__init__(message)
        self.color = 'yellow'

    def __str__(self):
        return 'Warning: %s' % self.message


def validate_plate_metadata(metadata):
    messages = []

    # used to keep track of attributes to validate across plates
    context = {
        'observed primer plates': [],
        'observed plate names': []
    }

    for plate in metadata:
        messages, context = _validate_plate(plate, context)

        if messages:
            print('Plate 1')
            for message in messages:
                print(message)

    metadata = pd.DataFrame(metadata)
    return metadata


def _validate_plate(plate_metadata, context):
    # things to validate:
    #   Error: Are primer plates repeated? -> check with context
    #   Error: Are sample plate names repeated? -> check with context
    #   Error: Are dates in the past?
    #   Error: Are dates formatted correctly?
    #   Error: Are all the fields present?
    #   Warn: If plates are left empty, check that all plates before are filled
    #   Error: All columns are exactly present, no more no less
    #
    #   Error: When there's more than 4 plates total, let users know when you
    #   are validating the fifth plate.
    #   Error: check for only ascii characters
    #   Check for an exact match of keys

    return [], context


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
