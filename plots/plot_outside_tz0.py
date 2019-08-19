"""
author : eablonet
date : 2019
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
# from scipy.special import erf
# import scipy.optimize as op

from packages.main import Stack as st
from packages.main import Stefan as ste
from packages.main import Material as ma
from packages.main import read_online_data as rd

rcParams.update({'figure.autolayout': True})  # enable tight layout

folder = '/Users/eablonet/Documents/0_phd/4_reports/reu_26_03_2019/images/'

z_lookfor = 1
fs = [8, 8]  # generic figsize
x_max = 1.5

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
fig = plt.figure(figsize=fs)
ax0 = fig.add_subplot(1, 1, 1)

# inclined
fig = plt.figure(figsize=fs)
ax1 = fig.add_subplot(1, 1, 1)

# none symmetrical
fig = plt.figure(figsize=fs)
ax2 = fig.add_subplot(1, 1, 1)

# all
fig = plt.figure(figsize=fs)
ax3 = fig.add_subplot(1, 1, 1)

colors = {
    'horizontal': 'tab:blue',
    'inclined': 'tab:orange',
    'random': 'tab:green'
}  # symbols horizontal, inclined, random
symbols = {0: 'o', 10: '^', 20: 'v', 30: 's', 45: 'd'}  # if inclined 10:b, 20:r, 30:m, 45:g
nivel = {'horizontal': 1, 'inclined': 2, 'random': 3}
c = {'horizontal': 0, 'inclined': 0, 'random': 0}

for _, row in manip_data.iterrows():
    cond = (
        row['completed_2'] == 1 and
        row['lieu'] == 0
        # row['completed_2'] == 1 and  # front detected
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

        print(row['symmetrical'])

        if row['alpha'] == 0 and row['symmetrical'] == 'Oui':
            type = 'horizontal'
        elif row['alpha'] > 0 and row['symmetrical'] == 'Oui':
            type = 'inclined'
        elif row['symmetrical'] == 'Non':
            type = 'random'

        c[type] += 1

        s.read_by_date(row['date'], int(row['serie']))

        zf = s.read_data()
        t_space, zf, zf_mean, zf_std = s.get_dynamic_front(zf)
        zf = np.max(zf_mean) / px_mm

        time = np.arange(len(zf_mean))/fps

        tf = t_end - t_nuc
        time_scale = time - dt_ini/fps
        t12 = time_scale[np.argmin(abs(zf_mean / px_mm / z0 - z_lookfor))]
        print(t12)
        time_scale /= t12

        if type == 'horizontal':
            # tf = t_end - t_nuc
            ax0.plot(
                time_scale,
                zf_mean / px_mm / z0,
                ls='none', c=colors[type],
                marker=symbols[int(row['alpha'])],
                mec=colors[type], mfc='none', ms=6,
                zorder=nivel[type]
            )
        elif type == 'inclined':
            ax1.plot(
                time_scale,
                zf_mean / px_mm / z0,
                ls='none', c=colors[type],
                marker=symbols[int(row['alpha'])],
                mec=colors[type], mfc='none', ms=6,
                zorder=nivel[type]
            )
        elif type == 'random':
            ax2.plot(
                time_scale,
                zf_mean / px_mm / z0,
                ls='none', c=colors[type],
                marker=symbols[int(row['alpha'])],
                mec=colors[type], mfc='none', ms=6,
                zorder=nivel[type]
            )
        ax3.plot(
            time_scale,
            zf_mean / px_mm / z0,
            ls='none', c=colors[type],
            marker=symbols[int(row['alpha'])],
            mec=colors[type], mfc='none', ms=6,
            zorder=nivel[type]
        )

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
    (0, 1), (0, 0), color=colors['horizontal'], linestyle='', marker='.'
)
inclined = plt.Line2D(
    (0, 1), (0, 0), color=colors['inclined'], linestyle='', marker='.'
)
random = plt.Line2D(
    (0, 1), (0, 0), color=colors['random'], linestyle='', marker='.'
)
deg0 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[0], c='k', markerfacecolor=None
)
deg10 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[10], c='k', markerfacecolor=None
)
deg20 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[20], c='k', markerfacecolor=None
)
deg30 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[30], c='k', markerfacecolor=None
)
deg45 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[45], c='k', markerfacecolor=None
)
stefanline = plt.Line2D(
    (0, 1), (0, 0), ls='-', c='k'
)
stefandilatationline = plt.Line2D(
    (0, 1), (0, 0), ls='--', c='k'
)

ax0.legend(
    [
        horizontal, deg0,
        stefanline, stefandilatationline,
    ],
    [
        'Horizontal',
        r'$\alpha$ = 0°',
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
        'Inclined',
        r'$\alpha$ = 0°', r'$\alpha$ = 10°', r'$\alpha$ = 20°',
        r'$\alpha$ = 30°', r'$\alpha$ = 45°',
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
        'Unsymmetrical',
        r'$\alpha$ = 0°', r'$\alpha$ = 10°', r'$\alpha$ = 20°',
        r'$\alpha$ = 30°', r'$\alpha$ = 45°',
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
        'Horizontal',
        'Inclined',
        'Unsymmetrical',
        r'$\alpha$ = 0°', r'$\alpha$ = 10°', r'$\alpha$ = 20°',
        r'$\alpha$ = 30°', r'$\alpha$ = 45°',
        'Stefan Model',
        'Stefan Model with dilatation',
    ],
    fancybox=True, fontsize=10, ncol=2,
)

# labels
ax0.set_xlabel(r'$t/t_{z_0}$', fontsize=18)
ax0.set_ylabel(r'$z_f/z_0$', fontsize=18)

ax1.set_xlabel(r'$t/t_{z_0}$', fontsize=18)
ax1.set_ylabel(r'$z_f/z_0$', fontsize=18)

ax2.set_xlabel(r'$t/t_{z_0}$', fontsize=18)
ax2.set_ylabel(r'$z_f/z_0$', fontsize=18)

ax3.set_xlabel(r'$t/t_{z_0}$', fontsize=18)
ax3.set_ylabel(r'$z_f/z_0$', fontsize=18)

# grids
ax0.grid(True)
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)

# x/y limits
ax0.set_xlim([-.1, x_max])
ax0.set_ylim([-.1, 1.5])

ax1.set_xlim([-.1, x_max])
ax1.set_ylim([-.1, 1.5])

ax2.set_xlim([-.1, x_max])
ax2.set_ylim([-.1, 1.5])

ax3.set_xlim([-.1, x_max])
ax3.set_ylim([-.1, 1.5])



# show
plt.show()
