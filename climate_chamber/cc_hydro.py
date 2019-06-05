# eablonet
# date : 2019

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# define colors :
green = np.array([62, 150, 81])/255

cc_data = pd.read_csv(
    './climate_chamber/data_hydro_cc.csv', delimiter='\t', header=None,
    decimal=","
)
manip_data = pd.read_csv(
    'manips.csv', skiprows=[0, 2], decimal=","
)

x = []
y = []

for _, row in cc_data.iterrows():
    x.append(row.values[0])
    y.append(row.values[1])


x_manip_1 = []
y_manip_1 = []
x_manip_2 = []
y_manip_2 = []
x_manip_3 = []
y_manip_3 = []
x_manip_4 = []
y_manip_4 = []

plt.figure(figsize=[8, 5])

# plot the area available by the chamber
plt.fill(
    x, y,
    color=green, alpha=.5,
    label='Usefull area of climate chamber'
)

treated = {0: 's', 1: 'o'}
inside = {0: 'tab:blue', 1: 'tab:orange'}

for _, row in manip_data.iterrows():
    # get ambiant temperature
    if ~np.isnan(row['serie']):
        if ~np.isnan(row['Ta']):
            Ta = row['Ta']
        elif ~np.isnan(row['Ta_int']):
            Ta = row['Ta_int']
        else:
            Ta = row['Ta_ext']

        # get humidity
        if ~np.isnan(row['phi']):
            phi = row['phi']
        elif ~np.isnan(row['phi_int']):
            phi = row['phi_int']
        else:
            phi = row['phi_ext']

        print(row['date'], 'n', row['serie'], ' : lieu :', row['lieu'])
        plt.scatter(
            phi, Ta,
            color=inside[row['lieu']], facecolor='',
            marker=treated[row['completed_2']],
        )



plt.xlabel(r'$\phi$ (%)')
plt.ylabel(r'$T_a$ (°C)')

plt.title('Climate chamber hygrometry')
plt.ylim([1, 40])
plt.legend(fancybox=True)
plt.grid(True, which='major', axis='both', linestyle='-')
plt.minorticks_on()
plt.grid(True, which='minor', axis='both', linestyle='--')
plt.show()
