"""


author: eablonet
date #2019

"""
import texttable as tt
import numpy as np

from packages.main import read_only_data as rd

manip_data = rd.get_data()

manip_data.info()

table = tt.Texttable()
table.set_deco(tt.Texttable.HEADER)
table.set_cols_align(["c", "c", "c", "c", "c", "c"])
table.add_rows([["#", "id", "env controled", "Ta", "phi", "alpha"]])

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
print('* for untreated experiments')
