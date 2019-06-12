# eablonet
# date : 2019

import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import pandas as pd

from packages.main import read_online_data as rd

rcParams.update({'figure.autolayout': True})  # enable tight layout

# define colors :
green = np.array([62, 150, 81])/255

cc_data = pd.read_csv(
    './climate_chamber/data_hydro_cc.csv', delimiter='\t', header=None,
    decimal=","
)
manip_data = rd.get_data()

x = []
y = []

for _, row in cc_data.iterrows():
    x.append(row.values[0])
    y.append(row.values[1])

fig0 = plt.figure(figsize=[16, 9]) # freezing
ax00 = fig0.add_subplot(2, 2, 1)  # inside / outside
ax01 = fig0.add_subplot(2, 2, 2)  # treated / not treated
ax02 = fig0.add_subplot(2, 2, 3)  # substrate
ax03 = fig0.add_subplot(2, 2, 4)  # drop /

fig1 = plt.figure(figsize=[16, 9]) # unfreezing
ax10 = fig1.add_subplot(2, 2, 1)  # inside / outside
ax11 = fig1.add_subplot(2, 2, 2)  # treated / not treated
ax12 = fig1.add_subplot(2, 2, 3)  # substrate
ax13 = fig1.add_subplot(2, 2, 4)  # drop / cuve / condensation

axes = [
    ax00, ax01, ax02, ax03,
    ax10, ax11, ax12, ax13,
]

# plot the area available by the chamber
for a in axes:
    a.fill(
        x, y,
        color=green, alpha=.2,
        label='Usefull area of climate chamber'
    )


def color_met(n):
    x = []
    a = n // 5
    b = n % 5
    for i in range(5):
        if i < b:
            x.extend(i*4 + np.arange(a+1))
        else:
            x.extend(i*4 + np.arange(a))
    print(x)
    c = plt.cm.tab20c(np.array(x).astype(int))
    print(c)
    return c

c = color_met(3)
inside = {0: c[0], 1: c[1]}
treated = {0: c[0], 1: c[1], 2: c[2]}
drop = {
    'a': c[0], 'b': c[0],
    'c': c[1], 'd': c[1],
    'e': c[2], 'f': c[2],
}
c = color_met(6)
substrate = {
    'Cu10a': c[0],
    'Cu10a +PTFE': c[1],
    'Cu10b': c[2],
    'Cu20a': c[3],
    'Cu20b': c[4],
    'Cuve20a': c[5],
}

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

        if row['type'] in ['a', 'c', 'e']:
            ax00.scatter(
                phi, Ta,
                color=inside[row['lieu']],  s=30
            )
            ax01.scatter(
                phi, Ta,
                color=treated[row['completed_2']], s=30
            )
            ax02.scatter(
                phi, Ta,
                color=substrate[row['ref']], s=30
            )
            ax03.scatter(
                phi, Ta,
                color=drop[row['type']], s=30
            )
        elif row['type'] in ['b', 'd', 'f']:
            ax10.scatter(
                phi, Ta,
                color=inside[row['lieu']], s=30
            )
            ax11.scatter(
                phi, Ta,
                color=treated[row['completed_2']], s=30
            )
            ax12.scatter(
                phi, Ta,
                color=substrate[row['ref']], s=30
            )
            ax13.scatter(
                phi, Ta,
                color=drop[row['type']], s=30
            )
phi = np.arange(5, 80)

for a in axes:
    a.plot(
        phi, 50-.5*phi, '--', color='tab:red',
        label='Limit of manipulation'
    )
    a.set_xlabel(r'$\phi$ (%)')
    a.set_ylabel(r'$T_a$ (°C)')
    a.set_ylim([1, 40])
    a.legend(fancybox=True)
    a.grid(True, which='major', axis='both', linestyle='-')
    a.minorticks_on()
    a.grid(True, which='minor', axis='both', linestyle='--')

fig0.suptitle('Hygrometry for Freezing', fontsize=16)
fig1.suptitle('Hygrometry for Unfreezing', fontsize=16)
ax00.set_title('Inside / Outside')
ax10.set_title('Inside / Outside')
ax01.set_title('Treated / Not Treated')
ax11.set_title('Treated / Not Treated')
ax02.set_title('Substrate')
ax12.set_title('Substrate')
ax03.set_title('Drop / condensation / cuve')
ax13.set_title('Drop / condensation / cuve')
c = color_met(3)
c0 = plt.Line2D(
    (0, 1), (0, 0), color=c[0], linestyle='', marker='.'
)
c1 = plt.Line2D(
    (0, 1), (0, 0), color=c[1], linestyle='', marker='.'
)
c2 = plt.Line2D(
    (0, 1), (0, 0), color=c[2], linestyle='', marker='.'
)
ax00.legend(
    [
        c0, c1,
    ],
    [
        'Inside',
        'Outside'
    ],
    fancybox=True, fontsize=10
)
ax10.legend(
    [
        c0, c1,
    ],
    [
        'Inside',
        'Outside'
    ],
    fancybox=True, fontsize=10
)
ax01.legend(
    [
        c0, c1, c2,
    ],
    [
        'Treated',
        'Not Treated',
        'Exploitable under condition'
    ],
    fancybox=True, fontsize=10
)
ax11.legend(
    [
        c0, c1, c2,
    ],
    [
        'Treated',
        'Not Treated',
        'Exploitable under condition'
    ],
    fancybox=True, fontsize=10
)
ax03.legend(
    [
        c0, c1, c2
    ],
    [
        'Drop',
        'Condensation',
        'Cuve'
    ],
    fancybox=True, fontsize=10
)
ax13.legend(
    [
        c0, c1, c2
    ],
    [
        'Drop',
        'Condensation',
        'Cuve'
    ],
    fancybox=True, fontsize=10
)

c = color_met(6)
c0 = plt.Line2D(
    (0, 1), (0, 0), color=c[0], linestyle='', marker='.'
)
c1 = plt.Line2D(
    (0, 1), (0, 0), color=c[1], linestyle='', marker='.'
)
c2 = plt.Line2D(
    (0, 1), (0, 0), color=c[2], linestyle='', marker='.'
)
c3 = plt.Line2D(
    (0, 1), (0, 0), color=c[3], linestyle='', marker='.'
)
c4 = plt.Line2D(
    (0, 1), (0, 0), color=c[4], linestyle='', marker='.'
)
c5 = plt.Line2D(
    (0, 1), (0, 0), color=c[5], linestyle='', marker='.'
)
ax02.legend(
    [
        c0, c1, c2, c3, c4, c5
    ],
    [
        'Cu10a',
        'Cu10a +PTFE',
        'Cu10b',
        'Cu20a',
        'Cu20b',
        'Cuve20a',
    ],
    fancybox=True, fontsize=10
)
ax12.legend(
    [
        c0, c1, c2, c3, c4, c5
    ],
    [
        'Cu10a',
        'Cu10a +PTFE',
        'Cu10b',
        'Cu20a',
        'Cu20b',
        'Cuve20a',
    ],
    fancybox=True, fontsize=10
)

plt.show()
