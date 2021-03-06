"""
author : eablonet
date : 2019
"""
import numpy as np
import matplotlib.pyplot as plt
# from scipy.special import erf
# import scipy.optimize as op
import pandas as pd
import Stack as st
import Stefan as ste
import Material as ma

folder = '/Users/eablonet/Documents/0_phd/4_reports/reu_26_03_2019/images/'

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
manip_data = pd.read_csv('manips.csv', skiprows=[0, 2], decimal=",")
manip_data.info()
s = st.Stack()

fig = plt.figure(figsize=[8, 4.5])
ax0 = fig.add_subplot(1, 1, 1)

fig = plt.figure(figsize=[8, 4.5])
ax1 = fig.add_subplot(1, 1, 1)

fig = plt.figure(figsize=[8, 4.5])
ax2 = fig.add_subplot(1, 1, 1)

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

        # read front information

        s.read_by_path(row['date'], int(row['serie']))

        s.read_image(2)  # to load an image else none might be load...to correct
        # al, zf = s.read_manual_front()
        zf = s.read_data()
        # s.view_all_profil(zf, 15)
        t_space, zf, zf_mean, zf_std = s.get_dynamic_front(zf)
        zf = np.max(zf_mean) / px_mm

        time = np.arange(len(zf_mean))/fps

        tf = t_end - t_nuc
        time_scale = time - dt_ini/fps
        ax0.plot(
            time_scale / tf*fps,
            zf_mean / px_mm / z0,
            ls='none', c=colors[type],
            marker=symbols[int(row['alpha'])],
            mec=colors[type], mfc='none', ms=6,
            zorder=nivel[type]
        )

        if type == 'horizontal':
            tf = t_end - t_nuc
            ax1.plot(
                (time - dt_ini/fps) / tf*fps,
                zf_mean / px_mm / z0,
                ls='none', c=colors[type],
                marker=symbols[int(row['alpha'])],
                mec=colors[type], mfc='none', ms=6,
                zorder=nivel[type]
            )
            ax2.plot(
                (time - dt_ini/fps) / tf*fps,
                zf_mean / px_mm / z0,
                ls='none', c=colors[type],
                marker=symbols[int(row['alpha'])],
                mec=colors[type], mfc='none', ms=6,
                zorder=nivel[type]
            )
        if type == 'inclined':
            ax2.plot(
                (time - dt_ini/fps) / tf*fps,
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

dip = ste.Diphasique()
dip.import_geometry(z0=2e-3, Nz=100, T_down=-10, T_up=0, T_ini=0, T_m=0)
dip.import_material('ice', 'water')
dip.define_stop_condition('zf', Nt=100)

ptd = ste.PseudoTwoDim()
ptd.import_geometry(z0=2e-3, Nz=100, T_down=-10, T_up=0, T_ini=0, r0=4e-3)
ptd.import_material('ice')
ptd.define_stop_condition('tf', dt=.01)


zs_mon = mon.front_position(adim=True)
zs_dip = dip.front_position(adim=True)
z_ptd = ptd.front_position(adim=True)  # to generate the time

ax0.plot(
    mon.time/max(mon.time), zs_mon,
    '-k', linewidth=2,
    label='Stefan Model',
    zorder=5
)
ax0.plot(
    dip.time/max(dip.time), zs_dip,
    '--k', linewidth=2,
    label='Stefan Model with Dilatation',
    zorder=5
)

ax0.plot(
    np.array(ptd.time)/max(ptd.time), z_ptd,
    '--k', linewidth=2,
    label='Stefan Model pseudo 2D-axi',
    zorder=5
)

# ax0.legend(fancybox=True, fontsize=18)
ax0.set_xlabel(r'$t/t_f$', fontsize=18)
ax0.set_ylabel(r'$x_f/x_0$', fontsize=18)
ax0.grid(True)

horizontal = plt.Line2D((0,1),(0,0), color=colors['horizontal'], linestyle='', marker='.')
inclined = plt.Line2D((0,1),(0,0), color=colors['inclined'], linestyle='', marker='.')
random = plt.Line2D((0,1),(0,0), color=colors['random'], linestyle='', marker='.')
deg0 = plt.Line2D((0,1),(0,0), ls='', marker=symbols[0], c='k', markerfacecolor=None)
deg10 = plt.Line2D((0,1),(0,0), ls='', marker=symbols[10], c='k', markerfacecolor=None)
deg20 = plt.Line2D((0,1),(0,0), ls='', marker=symbols[20], c='k', markerfacecolor=None)
deg30 = plt.Line2D((0,1),(0,0), ls='', marker=symbols[30], c='k', markerfacecolor=None)
deg45 = plt.Line2D((0,1),(0,0), ls='', marker=symbols[45], c='k', markerfacecolor=None)
stefanline = plt.Line2D((0,1),(0,0), ls='-', c='k')
stefandilatationline = plt.Line2D((0,1),(0,0), ls='--', c='k')

print(c)

ax0.legend(
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
ax1.legend(
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
ax2.legend(
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

ax1.set_xlabel(r'$t/t_f$', fontsize=18)
ax1.set_ylabel(r'$x_f/x_0$', fontsize=18)
ax1.grid(True)

ax2.set_xlabel(r'$t/t_f$', fontsize=18)
ax2.set_ylabel(r'$x_f/x_0$', fontsize=18)
ax2.grid(True)

ax0.set_xlim([-.1, 1.1])
ax0.set_ylim([-.1, 1.5])
ax1.set_xlim([-.1, 1.1])
ax1.set_ylim([-.1, 1.5])
ax2.set_xlim([-.1, 1.1])
ax2.set_ylim([-.1, 1.5])

# ax0.tight_layout()
# ax1.tight_layout()
# ax2.tight_layout()

plt.show()
