"""
Plot inside the chamber cinetic adim by t12.

author : eablonet
date : 2019
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd

from packages.main import Stack as st
from packages.main import Stefan as ste
from packages.main import Material as ma
from packages.main import read_online_data as rd

rcParams.update({'figure.autolayout': True})  # enable tight layout


folder = '/Users/eablonet/Documents/0_phd/4_reports/reu_26_03_2019/images/'

z_lookfor = .5

# ice
ice = ma.Ice()
rho = ice.get_rho()
k = ice.get_k()
cp = ice.get_cp()
Lf = 333500

# water (check values)
water = ma.Water()
rho_l = water.get_rho()
k_l = water.get_k()
cp_l = water.get_cp()

# solver
manip_data = rd.get_data()
manip_data.info()
s = st.Stack()

# horizontal
fig = plt.figure(figsize=[8, 4.5])
ax0 = fig.add_subplot(1, 1, 1)

# inclined
fig = plt.figure(figsize=[8, 4.5])
ax1 = fig.add_subplot(1, 1, 1)

# none symmetrical
fig = plt.figure(figsize=[8, 4.5])
ax2 = fig.add_subplot(1, 1, 1)

# all
fig = plt.figure(figsize=[8, 4.5])
ax3 = fig.add_subplot(1, 1, 1)

colors = {
    5: 'tab:blue',
    15: 'tab:orange',
    20: 'tab:green'
}  # symbols horizontal, inclined, random
symbols = {40: 'o', 50: '^', 60: 'v', 70: 's', 'nan': 'd'}  # if inclined 10:b, 20:r, 30:m, 45:g
nivel = {5: 1, 15: 2, 20: 3}
c = {5: 0, 15: 0, 20: 0}

Tamb = []
phi_amb = []

for _, row in manip_data.iterrows():
    cond = (
        row['completed_2'] == 1 and
        row['lieu'] == 1 and
        row['type'] == 'e'
    )
    if cond:
        print('Found :', row['date'], 'n', row['serie'])
        # load data
        px_mm = row['px_mm']
        z0 = row['z0']/px_mm
        r0 = row['r0']/px_mm
        theta = 1/2*(row['ca_left'] + row['ca_right'])
        T_nuc = row['Tc_nuc']
        Ta = row['Ta_int']
        fps = row['fps_real']
        t_nuc = row['t_nuc_calc']
        t_end = row['t_end_calc']
        t_ref = row['t_ref_calc']
        dt_ini = t_nuc - t_ref

        print('Ta consigne :', row['Ta'])
        print(r'$\phi$ consigne :', row['phi'])

        print('Ta interne :', row['Ta_int'])
        print(r'$\phi$ interne :', row['phi_int'])

        print('Ta externe :', row['Ta_ext'])
        print(r'$\phi$ externe :', row['phi_ext'])

        if row['Ta'] not in Tamb:
            Tamb.append(row['Ta'])
        if row['phi'] not in phi_amb:
            phi_amb.append(row['phi'])

        if np.isnan(row['phi']):
            phi = 'nan'
        else:
            phi = int(row['phi'])

        print(phi)

        c[row['Ta']] += 1

        # read front information

        s.read_by_path(row['date'], int(row['serie']))

        s.read_image(2)  # to load an image else none might be load...to correct
        # al, zf = s.read_manual_front()
        zf = s.read_data()
        # s.view_all_profil(zf, 15)
        t_space, zf, zf_mean, zf_std = s.get_dynamic_front(zf)
        zf = np.nanmax(zf_mean) / px_mm

        time = np.arange(len(zf_mean))/fps

        tf = t_end - t_nuc
        time_scale = time - dt_ini/fps
        t12 = time_scale[np.argmin(abs(zf_mean / px_mm / z0 - z_lookfor))]
        time_scale /= t12

        if row['Ta'] == 5:
            # tf = t_end - t_nuc
            ax0.plot(
                time_scale,
                zf_mean / px_mm / z0,
                ls='none', c=colors[row['Ta']],
                marker=symbols[phi],
                mec=colors[row['Ta']], mfc='none', ms=6,
                zorder=nivel[row['Ta']]
            )
        elif row['Ta'] == 15:
            ax1.plot(
                time_scale,
                zf_mean / px_mm / z0,
                ls='none', c=colors[row['Ta']],
                marker=symbols[phi],
                mec=colors[row['Ta']], mfc='none', ms=6,
                zorder=nivel[row['Ta']]
            )
        elif row['Ta'] == 20:
            ax2.plot(
                time_scale,
                zf_mean / px_mm / z0,
                ls='none', c=colors[row['Ta']],
                marker=symbols[phi],
                mec=colors[row['Ta']], mfc='none', ms=6,
                zorder=nivel[row['Ta']]
            )
        ax3.plot(
            time_scale,
            zf_mean / px_mm / z0,
            ls='none', c=colors[row['Ta']],
            marker=symbols[phi],
            mec=colors[row['Ta']], mfc='none', ms=6,
            zorder=nivel[row['Ta']]
        )

print('All temperature', Tamb)
print('All phi', phi_amb)

# theorical data
mon = ste.Stefan()
mon.import_material('ice')
mon.import_geometry(z0=2e-3, Nz=100, T_down=-10, T_up=0, T_ini=0)
mon.define_stop_condition('zf', Nt=100)

# zs_comp = np.sqrt(2 * 2.2 * (0-T_nuc) * time / (916.2*333500))

dip = ste.Diphasique()
dip.import_geometry(z0=2e-3, Nz=100, T_down=-10, T_up=0, T_ini=0, T_m=0)
dip.import_material('ice', 'water')
dip.define_stop_condition('zf', Nt=100)

zs_mon = mon.front_position(adim=True)
zs_dip = dip.front_position(adim=True)

t12 = mon.time[np.argmin(abs(zs_mon - z_lookfor))]
ax0.plot(
    mon.time/t12, zs_mon,
    '-k', linewidth=2,
    label='Stefan Model',
    zorder=5
)
ax1.plot(
    mon.time/t12, zs_mon,
    '-k', linewidth=2,
    label='Stefan Model',
    zorder=5
)
ax2.plot(
    mon.time/t12, zs_mon,
    '-k', linewidth=2,
    label='Stefan Model',
    zorder=5
)
ax3.plot(
    mon.time/t12, zs_mon,
    '-k', linewidth=2,
    label='Stefan Model',
    zorder=5
)
t12 = dip.time[np.argmin(abs(zs_dip - z_lookfor))]
ax0.plot(
    dip.time/t12, zs_dip,
    '--k', linewidth=2,
    label='Stefan Model with Dilatation',
    zorder=5
)
ax1.plot(
    dip.time/t12, zs_dip,
    '--k', linewidth=2,
    label='Stefan Model with Dilatation',
    zorder=5
)
ax2.plot(
    dip.time/t12, zs_dip,
    '--k', linewidth=2,
    label='Stefan Model with Dilatation',
    zorder=5
)
ax3.plot(
    dip.time/t12, zs_dip,
    '--k', linewidth=2,
    label='Stefan Model with Dilatation',
    zorder=5
)

# legends
horizontal = plt.Line2D(
    (0, 1), (0, 0), color=colors[5], linestyle='', marker='.'
)
inclined = plt.Line2D(
    (0, 1), (0, 0), color=colors[15], linestyle='', marker='.'
)
random = plt.Line2D(
    (0, 1), (0, 0), color=colors[20], linestyle='', marker='.'
)
deg0 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[40], c='k', markerfacecolor=None
)
deg10 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[50], c='k', markerfacecolor=None
)
deg20 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[60], c='k', markerfacecolor=None
)
deg30 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[70], c='k', markerfacecolor=None
)
deg45 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols['nan'], c='k', markerfacecolor=None
)
stefanline = plt.Line2D(
    (0, 1), (0, 0), ls='-', c='k'
)
stefandilatationline = plt.Line2D(
    (0, 1), (0, 0), ls='--', c='k'
)

ax0.legend(
    [
        horizontal, deg0, deg10, deg20, deg30, deg45,
        stefanline, stefandilatationline,
    ],
    [
        'Ta = 5°C',
        r'$\phi$ = 40%', r'$\phi$ = 50%', r'$\phi$ = 60%',
        r'$\phi$ = 70%', r'$\phi$ = n.a.',
        'Stefan Model',
        'Stefan Model with dilatation',
    ],
    fancybox=True, fontsize=10, ncol=2,
)
ax1.legend(
    [
        inclined, deg0, deg10, deg20, deg30, deg45,
        stefanline, stefandilatationline,
    ],
    [
        'Ta = 15°C',
        r'$\phi$ = 40%', r'$\phi$ = 50%', r'$\phi$ = 60%',
        r'$\phi$ = 70%', r'$\phi$ = n.a.',
        'Stefan Model',
        'Stefan Model with dilatation',
    ],
    fancybox=True, fontsize=10, ncol=2,
)
ax2.legend(
    [
        random, deg0, deg10, deg20, deg30, deg45,
        stefanline, stefandilatationline,
    ],
    [
        'Ta = 20°C',
        r'$\phi$ = 40%', r'$\phi$ = 50%', r'$\phi$ = 60%',
        r'$\phi$ = 70%', r'$\phi$ = n.a.',
        'Stefan Model',
        'Stefan Model with dilatation',
    ],
    fancybox=True, fontsize=10, ncol=2,
)
ax3.legend(
    [
        horizontal, inclined, random, deg0, deg10, deg20, deg30, deg45,
        stefanline, stefandilatationline,
    ],
    [
        'Ta = 5°C',
        'Ta = 15°C',
        'Ta = 20°C',
        r'$\phi$ = 40%', r'$\phi$ = 50%', r'$\phi$ = 60%',
        r'$\phi$ = 70%', r'$\phi$ = n.a.',
        'Stefan Model',
        'Stefan Model with dilatation',
    ],
    fancybox=True, fontsize=10, ncol=2,
)

# labels
ax0.set_xlabel(r'$t/t_{1/2}$', fontsize=18)
ax0.set_ylabel(r'$z_f/z_0$', fontsize=18)

ax1.set_xlabel(r'$t/t_{1/2}$', fontsize=18)
ax1.set_ylabel(r'$z_f/z_0$', fontsize=18)

ax2.set_xlabel(r'$t/t_{1/2}$', fontsize=18)
ax2.set_ylabel(r'$z_f/z_0$', fontsize=18)

ax3.set_xlabel(r'$t/t_{1/2}$', fontsize=18)
ax3.set_ylabel(r'$z_f/z_0$', fontsize=18)

# grids
ax0.grid(True)
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)

# x/y limits
ax0.set_xlim([-.5, 5])
ax0.set_ylim([-.1, 1.5])

ax1.set_xlim([-.5, 5])
ax1.set_ylim([-.1, 1.5])

ax2.set_xlim([-.5, 5])
ax2.set_ylim([-.1, 1.5])

ax3.set_xlim([-.5, 5])
ax3.set_ylim([-.1, 1.5])

# show
plt.show()
