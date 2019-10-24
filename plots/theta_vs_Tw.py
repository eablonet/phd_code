"""
author : eablonet
date : 2019
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

from packages.main import Stack as st
from packages.main import Stefan as ste
from packages.main import Material as ma
from packages.main import data_base as db

from packages.design import color as ea_color
rcParams.update({'figure.autolayout': True})  # enable tight layout

# init physical params
# --------------------
# ice
ice = ma.Water('solid')
rho = ice.get_rho()
k = ice.get_k()
cp = ice.get_cp()
Lf = ice.get_Lf()

# water
water = ma.Water('liquid')
rho_l = water.get_rho()
k_l = water.get_k()
cp_l = water.get_cp()

# init figures params
# -------------------
inside = {True: 'tab:gray', False: 'none'}
marker = {'constant': '_', 'ramp': '^', 'other': 'x'}


# load database
# -------------
manip_data = db.get_data('usable')

# init stack manager
# ------------------
s = st.Stack()

# init figures
# ------------
"""Figure 0: Geometric data

name: fig0

axes: all represents theta over one data
----
    nb: 6
    0: versus T_w_set
    1: versus phi
    2: versus T_nuc
    3: versus Ta
"""
fig0 = plt.figure(figsize=[16, 9])
fig0.canvas.set_window_title('t_z0 for Geometric')
axes0 = []
for i in range(4):
    axes0.append(fig0.add_subplot(2, 2, i+1))


"""Figure 1: Thermal data

name: fig1

axes: all represents tz0_exp/t_z0_th over one data
----
    nb: 9
    0: versus Ta_set (only for inside)
    1: versus Ta_cc (only for inside)
    2: versus Ta
    3: versus phi_set (only for inside)
    4: versus phi_cc (only for inside)
    5: versus phi
    6: versus DT_nuc (read on cryostat)
    7: versus DT_set (set point of cryostat)
    8: versus DT_set-DT_nuc

"""
fig1 = plt.figure(figsize=[16, 9])
axes1 = []
for i in range(9):
    axes1.append(fig1.add_subplot(3, 3, i+1))

data = {
    'theta_l': [],  # contact angle left
    'theta_r': [],  # contact angle right
    'Tnuc': [],
    'Tset': [],
    'inside': [],  # bool
    'ref_therm': [],  # constant or ramp
    'Ta_cc': [],  # nan for outside
    'Ta': [],  # could be nan
    'phi_cc': [],  # nan for outside
    'phi': [],  # could be nan
}
for _, row in manip_data.iterrows():
    cond = (
        row['type'] == 'a',
        row['alpha'] < 1,
    )
    if cond:

        # load geometric data
        # -------------------
        theta_l = row['ca_left']
        theta_r = row['ca_right']
        data['theta_l'].append(theta_l)
        data['theta_r'].append(theta_r)

        # load thermal data
        # -----------------
        Tc_nuc = row['Tc_nuc']
        Tc_set = row['Tc_set']
        Ta_cc = row['Ta_int']
        Ta = row['Ta_ext']
        phi_cc = row['phi_int']
        phi = row['phi_ext']

        data['Tnuc'].append(Tc_nuc)
        data['Tset'].append(Tc_set)
        data['Ta_cc'].append(Ta_cc)
        data['Ta'].append(Ta)
        data['phi_cc'].append(phi_cc)
        data['phi'].append(phi)

        # other data
        # ----------
        data['inside'].append(True if row['lieu'] else False)

# add figure design
# -----------------
for k in range(len(data['theta_l'])):
    # fig1 plotting
    # -------------
    axes0[0].plot(
        data['Tnuc'][k],
        data['theta_l'][k],
        color='g',
        marker='o',
    )
    axes0[0].plot(
        data['Tnuc'][k],
        data['theta_r'][k],
        color='r',
        marker='o',
    )

    axes0[1].plot(
        data['phi'][k],
        data['theta_l'][k],
        color='g',
        marker='o',
    )
    axes0[1].plot(
        data['phi'][k],
        data['theta_r'][k],
        color='r',
        marker='o',
    )

    axes0[2].plot(
        data['Tset'][k],
        data['theta_l'][k],
        color='g',
        marker='o',
    )
    axes0[2].plot(
        data['Tset'][k],
        data['theta_r'][k],
        color='r',
        marker='o',
    )

    axes0[3].plot(
        data['Ta'][k],
        data['theta_l'][k],
        color='g',
        marker='o',
    )
    axes0[3].plot(
        data['Ta'][k],
        data['theta_r'][k],
        color='r',
        marker='o',
    )


# labler fig1
# -----------
axes0[0].set_xlabel(r'$T_{nuc}$ (°C)')
axes0[1].set_xlabel(r'$\phi$ (%)')
axes0[2].set_xlabel(r'$T_{set}$ (°C)')
axes0[3].set_xlabel(r'$T_a$ (°C)')

[a.grid(True) for a in axes0+axes1]
[a.set_ylabel(r'$\theta$ (°)') for a in axes0]

plt.figure(figsize=(8, 4.5))
theta = [
    (data['theta_l'][k] + data['theta_r'][k])/2
    for k in range(len(data['theta_l']))
]
plt.plot(
    data['Tset'],
    theta,
    ls='none', color='tab:blue',
    marker='o', label='Angle à gauche', mfc='none', mew=2
)
# plt.plot(
#     data['Tset'],
#     data['theta_r'],
#     ls='none', color='r',
#     marker='o', label='Angle à droite'
# )
plt.xlabel(r'$T_{w,setpoint}$ (°C)', fontsize=14)
plt.ylabel(r'$\theta$ (°)', fontsize=14)
plt.grid(True)
# plt.legend(fancybox=True, shadow=True)
plt.xlim((-22, -5))
# show
plt.show()
