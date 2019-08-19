"""


author: eablonet
date #2019

"""
import texttable as tt
import colored as clr
import numpy as np
from random import randint

from packages.main import read_online_data as rd

manip_data = rd.get_data()

manip_data.info()

table = tt.Texttable()
table.set_deco(tt.Texttable.HEADER)
table.set_cols_align(["c", "c", "c", "c", "c", "c"])
table.add_rows([["#", "id", "env controled", "Ta", "phi", "alpha"]])

id_ca = []
id_not = []
id_z0 = []
id_r0 = []

c = 1
chamber = {0: " ", 1: "x"}
for _, row in manip_data.iterrows():
    cond = (
        row['type'] in ['a', 'c', 'e']  # freezing
    )
    if ~np.isnan(row['serie']) and cond:  # si l'experience existe
        if ~np.isnan(row['Ta']):
            Ta = row['Ta']
        elif ~np.isnan(row['Ta_int']):
            Ta = row['Ta_int']
        elif ~np.isnan(row['Ta_ext']):
            Ta = row['Ta_ext']
        else:
            Ta = -1

        # get humidity
        if ~np.isnan(row['phi']):
            phi = row['phi']
        elif ~np.isnan(row['phi_int']):
            phi = row['phi_int']
        elif ~np.isnan(row['phi_ext']):
            phi = row['phi_ext']
        else:
            phi = -1

        if row['completed_2'] == 1:  # treated and usable
            id = row['date']+'n'+str(int(row['serie']))
        elif row['completed_2'] == 2:  # treated but unusable
            id = '_'+row['date']+'n'+str(int(row['serie']))
        else:
            id = '*'+row['date']+'n'+str(int(row['serie']))
            id_not.append(id)

        if np.isnan(row['ca_left']):
            if id[0] not in ['*', '_']:
                id_ca.append(id)
            id = '#' + id

        if np.isnan(row['z0']):
            if (id[0] not in ['*', '_', '#'] and id[1] not in ['*', '_']):
                id_z0.append(id)

        if np.isnan(row['r0']):
            if (id[0] not in ['*', '_', '#'] and id[1] not in ['*', '_']):
                id_r0.append(id)

        print(id, type(id), len(id))

        table.add_row([
            c,
            id,
            chamber[row['lieu']],
            Ta,
            phi,
            row['alpha'],
        ])

        c += 1

print(table.draw())
print('===========================')
print('* for untreated experiments')
print('# means thats there is not contact angle measurement')
print('_ have been treated but is unusable')

print('===========================')
print(len(id_r0))
print('Contact angle to measure : ', id_ca[randint(0, len(id_ca)-1)] if len(id_ca) > 0 else id_ca)
print('Geometry(z0) to measure : ', id_z0[randint(0, len(id_z0)-1)] if len(id_z0) > 0 else id_z0)
print('Geometry(r0) to measure : ', id_r0[randint(0, len(id_r0)-1)] if len(id_r0) > 0 else id_r0)
print('Profile to identifiy : ', id_not[randint(0, len(id_not)-1)])
