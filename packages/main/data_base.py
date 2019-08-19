# -*- coding: utf-8 -*-
"""
Database isolata


author: eablonet Emeryk Ablonet
    phD Student
copyright: IMFT and Toulouse INP
date: #2019

methods
-------
    get_constant(data=None)
    get_ramp(data=None)

"""
from pandas import concat
import math

try:
    from packages.main.read_online_data import get_data as read_data
except ModuleNotFoundError:
    from read_online_data import get_data as read_data


# get methods
# -----------
def get_data(ref, data=None):
    """Get, generic method.

    inputs
    ------
        ref: string or list of string.
            To choose among:
            - constant: constant wall temperature
            - ramp: ramp wall temperature
            - treated: treated experiments (=/ usable)
            - untreadted: untreated experiments
            - usable: usable experiments
            - inside: inside the chamber experiments
            - outside: outside the chamber experiments
            - horizontal: inclination == 0°
            - inclined: inclination > 0°
            - symmetrical: symmetric (can be inclined)
            - unsymmetrical: unsymmetric (can be inclined) included experiments
                without top view
        data: dataframe or None
            Pre-filtered data base. To use the whole database pass None.
            Default is None

    output
    ------
        out : dataframe
        Contains the filtred database

    """
    data = read_data() if data is None else None

    case = {
        'constant':
            (data.ref_therm == 'Tw-7') |
            (data.ref_therm == 'Tw-10') |
            (data.ref_therm == 'Tw-12') |
            (data.ref_therm == 'Tw-15'),
        'ramp':
            (data.ref_therm == 'qw1') |
            (data.ref_therm == 'qw2') |
            (data.ref_therm == 'qw3') |
            (data.ref_therm == 'qw4'),
        'treated':
            (data.completed_2 == 1) |
            (data.completed_2 == 2),
        'untreated':
            (data.completed_2 == 0),
        'usable':
            (data.completed_2 == 1),
        'inside':
            (data.lieu == 1),
        'outside':
            (data.lieu == 1),
        'horizontal':
            (data.alpha == 0),
        'inclined':
            (data.alpha > 0),
        'symmetrical':
            (data.symmetrical == 'Oui'),
        'unsymmetrical':
            (data.symmetrical == 'Non') |
            (data.symmetrical.isnull()),
    }

    return data[case[ref]]


# add methods
# -----------
def add_data(data, ref):
    """Add ref values to database.

    inputs
    ------
        data: dataframe (pandas)
            Pre-filtered database in which is wanted to add elements.
        ref: string, list string
            What to include, among:
            -constant: constant wall temperature
            -ramp: ramp wall temperature
            -treated: experiments treated (caution =/ usable)
            -untreated: untreated experiments
            -usable: usable experiments
            - inside: inside the chamber experiments
            - outside: outside the chamber experiments
            - horizontal: inclination == 0°
            - inclined: inclination > 0°
            - symmetrical: symmetric (can be inclined)
            - unsymmetrical: unsymmetric (can be inclined) included experiments
                without top view

    output
    ------
        out : dataframe (pandas)
            Contains the filtered database with added elements

    """
    new_data = get_data(ref)

    data = concat(
        [data, new_data],
        ignore_index=False,
    ).drop_duplicates().sort_index()

    return data


if __name__ == '__main__':
    print(get_data.__doc__)
    a = get_data('unsymmetrical')
    b = add_data(a, 'constant')
    print(a.symmetrical)
