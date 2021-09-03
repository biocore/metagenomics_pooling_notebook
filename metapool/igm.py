import os

from datetime import datetime

from openpyxl import load_workbook
from openpyxl.styles import PatternFill

_IGM_ROW_OFFSET = 25


def _property_maker(obj, name, cell, default):
    def fget(self):
        return getattr(self, '_' + name)

    def fset(self, value):
        if value is not None:
            self._sheet[cell] = value
        setattr(self, '_' + name, value)

    # set the default values
    setattr(obj, '_' + name, default)
    if default is not None:
        getattr(obj, '_sheet')[cell] = default

    setattr(type(obj), name, property(fget, fset))


class IGMManifest(object):
    """Interface to write an IGM manifest file

    Attributes
    ----------
    submission_date: default today's date
        When the sample sheet is submitted.
    institute: default Knight Lab
        Institution or laboratory.
    pi_name: default Dr. Knight
        Name of the PI.
    pi_email: default mackkenzie.m.bryant@gmail.com
        Email for the PI
    contact_name: default MacKenzie Bryant
        Contact's name.
    project_number:
        Project's number.
    task_number:
        Tasks' number
    platform: default NovaSeq S4:
        Sequencing platform.
    run_type: default PE150
        Sequencing run type.
    custom_primer:
        Notes about the primers in the run.
    number_of_samples:
        Samples in the run.
    number_of_lanes:
        Samples in the lane.
    """
    def __init__(self):
        self._workbook = _load_igm_template()
        self._sheet = self._workbook['Sample Information']

        _property_maker(self, 'submission_date', 'B2',
                        datetime.strftime(datetime.today(), '%m/%d/%y'))
        _property_maker(self, 'institute', 'B3', 'Knight Lab')
        _property_maker(self, 'pi_name', 'B4', 'Dr. Knight')
        _property_maker(self, 'pi_email', 'B5', 'mackenzie.m.bryant@gmail.com')
        _property_maker(self, 'contact_name', 'B6', 'MacKenzie Bryant')
        _property_maker(self, 'contact_email', 'B7',
                        'mackenzie.m.bryant@gmail.com')

        _property_maker(self, 'project_number', 'B12', None)
        _property_maker(self, 'task_number', 'B13', None)

        _property_maker(self, 'platform', 'B18', 'NovaSeq S4')
        _property_maker(self, 'run_type', 'B19', 'PE150')
        _property_maker(self, 'custom_primer', 'B20',
                        'No-Standard Illumina Primers are fine')

        _property_maker(self, 'number_of_samples', 'B22', None)
        _property_maker(self, 'number_of_lanes', 'B23', '1')

        self._pools = None

        # default value requested by IGM
        self._sheet['D1'] = '150x8x8x150'
        self._sheet['D1'].fill = PatternFill(fill_type='solid',
                                             fgColor='FFFF00')

    @property
    def pools(self):
        return self._pools

    @pools.setter
    def pools(self, value):
        """List of pools to write in the manifest"""
        if value is None:
            for i, pool in enumerate(self._pools):
                row = str(_IGM_ROW_OFFSET + i)

                del self._sheet['A' + row]
                del self._sheet['B' + row]
                del self._sheet['C' + row]
                del self._sheet['D' + row]
        else:
            for i, pool in enumerate(value):
                row = str(_IGM_ROW_OFFSET + i)

                self._sheet['A' + row] = pool
                self._sheet['B' + row] = pool
                self._sheet['C' + row] = 500
                self._sheet['D' + row] = 'KHP'

        self._pools = value

    def write(self, path=None):
        """Write the manifest to disk

        Parameters
        ----------
        path: str, optional
            Path where to write the manifest. If `None`, IGM's default format
            is used (YYYY_MM_DD_PI_Sequencing_Runs_Manifest_2021.xlsx).

        Raises
        ------
        ValueError
            If a required attribute is missing.
        """
        required = ['project_number', 'task_number', 'platform',
                    'number_of_samples', 'number_of_lanes', 'pools']
        for prop in required:
            if getattr(self, prop) is None:
                raise ValueError('%s cannot be empty, you need to set a value'
                                 % prop)

        path = self._default_path() if path is None else path
        print('Saving manifest to %s' % path)

        self._workbook.save(path)

    def _default_path(self):
        path = (datetime.strftime(datetime.today(), '%Y_%m_%d_') +
                # should be the last name
                self.pi_name.split(' ')[-1] + '_' +
                # all the pools with spaces replaced for underscores
                '_'.join([p.replace(' ', '_') for p in self._pools]) +
                '_Manifest_2021.xlsx')
        return path


def _load_igm_template():
    """Helper function to load IGM's manifest template"""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data',
                        'template.xlsx')
    return load_workbook(path)
