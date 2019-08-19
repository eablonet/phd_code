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
from packages.main import read_online_data as rd

rcParams.update({'figure.autolayout': True})  # enable tight layout

folder = '/Users/eablonet/Documents/0_phd/4_reports/reu_26_03_2019/images/'

z_lookfor12 = .5
z_lookforz0 = 1

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

# time scale exp vs theorical
fig0 = plt.figure(figsize=[16, 9])
ax00 = fig0.add_subplot(2, 2, 1)  # tz0_1/2 vs tz0_th
ax01 = fig0.add_subplot(2, 2, 2)  # t12_exp vs t12_th
ax02 = fig0.add_subplot(2, 2, 3)  # tf_exp vs tf_th
ax03 = fig0.add_subplot(2, 2, 4)  # labels

# tz0_exp / tz0_th geometric parameters
fig1 = plt.figure(figsize=[16, 9])
ax10 = fig1.add_subplot(2, 2, 1)  # tz0_exp / tz0_th vs a = z0/r0
ax11 = fig1.add_subplot(2, 2, 2)  # tz0_exp / tz0_th vs ca
ax12 = fig1.add_subplot(2, 2, 3)  # tz0_exp / tz0_th vs volume
ax13 = fig1.add_subplot(2, 2, 4)  # tz0_exp / tz0_th vs Ste

# zf/z0 for geometry
fig2 = plt.figure(figsize=[16, 9])
ax20 = fig2.add_subplot(2, 2, 1)  # zf/z0 vs a = z0/r0
ax21 = fig2.add_subplot(2, 2, 2)  # zf/z0 vs ca
ax22 = fig2.add_subplot(2, 2, 3)  # zf/z0 vs volume

fig3 = plt.figure(figsize=[16, 9])
ax30 = fig3.add_subplot(1, 1, 1)

axes = [
    ax00, ax01, ax02,
    ax10, ax11, ax12, ax13,
    ax20, ax21, ax22,
    ax30,
]

colors = {
    'horizontal': 'tab:blue',
    'inclined': 'tab:orange',
    'random': 'tab:green'
}  # symbols horizontal, inclined, random
symbols = {0: 'o', 10: '^', 20: 'v', 30: 's', 45: 'd'}  # if inclined 10:b, 20:r, 30:m, 45:g
nivel = {'horizontal': 1, 'inclined': 2, 'random': 3}
facecolor = {1: 'tab:blue', 0: 'none'}
c = {'horizontal': 0, 'inclined': 0, 'random': 0}

compteur_sec = 0

for _, row in manip_data.iterrows():
    cond = (
        row['completed_2'] == 1 and
        row['type'] == 'a'
        # row['lieu'] == 0
        # row['completed_2'] == 1 and  # front detected
    )
    if cond:
        # load data
        px_mm = row['px_mm']
        z0 = row['z0']/px_mm
        r0 = row['r0']/px_mm
        theta = 1/2*(row['ca_left'] + row['ca_right'])
        Tw = row['Tc_nuc'] if ~np.isnan(row['Tc_nuc']) else row['Tc_set']
        # Ta = row['Ta_int']
        fps = row['fps_real']
        t_nuc = row['t_nuc_calc']
        t_end = row['t_end_calc']
        t_ref = row['t_ref_calc']
        dt_ini = t_nuc - t_ref

        if row['alpha'] == 0 and row['symmetrical'] == 'Oui':
            type = 'horizontal'
        elif row['alpha'] > 0 and row['symmetrical'] == 'Oui':
            type = 'inclined'
        elif row['symmetrical'] == 'Non':
            type = 'random'
        elif np.isnan(row['symmetrical']) and row['lieu']:
            type = 'horizontal'  # cas des goutte à l'intérieur de la chambre, pas de visu de haut

        c[type] += 1

        # read front information
        s.read_by_date(row['date'], int(row['serie']))
        zf = s.read_data()
        t_space, zf, zf_mean, zf_std = s.get_dynamic_front(zf)
        try:
            zf = zf_mean[int(t_end-t_ref)] / px_mm
        except IndexError:
            zf = zf_mean[-1] / px_mm

        time = np.arange(len(zf_mean))/fps
        time_scale = time - dt_ini/fps  # on enlève le temps avant tnuc

        t12_exp = time_scale[np.argmin(abs(zf_mean / px_mm / z0 - .5))]
        tz0_exp = time_scale[np.argmin(abs(zf_mean / px_mm / z0 - 1))]
        tf_exp = (t_end - t_nuc)/fps

        # theorical data
        mon = ste.Stefan()
        mon.import_material('ice')
        mon.import_geometry(
            z0=z0, Nz=100, T_down=Tw, T_up=0, T_ini=0
        )
        mon.define_stop_condition('zf', Nt=100)
        zs_mon = mon.front_position(adim=True)

        t12_th = mon.time[np.argmin(abs(zs_mon - .5))]
        tz0_th = mon.time[np.argmin(abs(zs_mon - 1))]
        tf_th = max(mon.time)

        Ste = mon.calc_St()

        ax00.plot(
            t12_th, t12_exp,
            ls='none', c=colors[type],
            marker=symbols[int(row['alpha'])],
            mec=colors[type], ms=6,
            markerfacecolor=facecolor[row['lieu']],
            zorder=nivel[type]
        )
        ax01.plot(
            tz0_th, tz0_exp,
            ls='none', c=colors[type],
            marker=symbols[int(row['alpha'])],
            mec=colors[type], ms=6,
            markerfacecolor=facecolor[row['lieu']],
            zorder=nivel[type]
        )
        ax02.plot(
            tf_th, tf_exp,
            ls='none', c=colors[type],
            marker=symbols[int(row['alpha'])],
            mec=colors[type], ms=6,
            markerfacecolor=facecolor[row['lieu']],
            zorder=nivel[type]
        )
        ax10.plot(
            z0/r0, tz0_exp/tz0_th,
            ls='none', c=colors[type],
            marker=symbols[int(row['alpha'])],
            mec=colors[type], ms=6,
            markerfacecolor=facecolor[row['lieu']],
            zorder=nivel[type]
        )
        ax11.plot(
            theta, tz0_exp/tz0_th,
            ls='none', c=colors[type],
            marker=symbols[int(row['alpha'])],
            mec=colors[type], ms=6,
            markerfacecolor=facecolor[row['lieu']],
            zorder=nivel[type]
        )
        ax12.plot(
            row['v_theo_mm'], tz0_exp/tz0_th,
            ls='none', c=colors[type],
            marker=symbols[int(row['alpha'])],
            mec=colors[type], ms=6,
            markerfacecolor=facecolor[row['lieu']],
            zorder=nivel[type]
        )
        ax13.plot(
            Ste, tz0_exp/tz0_th,
            ls='none', c=colors[type],
            marker=symbols[int(row['alpha'])],
            mec=colors[type], ms=6,
            markerfacecolor=facecolor[row['lieu']],
            zorder=nivel[type]
        )

# dip = ste.Diphasique()
# dip.import_geometry(z0=2e-3, Nz=100, T_down=-10, T_up=0, T_ini=0, T_m=0)
# dip.import_material('ice', 'water')
# dip.define_stop_condition('zf', Nt=100)
#
# zs_dip = dip.front_position(adim=True)
#
# t12 = dip.time[np.argmin(abs(zs_dip - z_lookfor))]
# ax0.plot(
#     dip.time/t12, zs_dip,
#     '--k', linewidth=2,
#     label='Stefan Model with Dilatation',
#     zorder=5
# )
# ax1.plot(
#     dip.time/t12, zs_dip,
#     '--k', linewidth=2,
#     label='Stefan Model with Dilatation',
#     zorder=5
# )
# ax2.plot(
#     dip.time/t12, zs_dip,
#     '--k', linewidth=2,
#     label='Stefan Model with Dilatation',
#     zorder=5
# )
# ax3.plot(
#     dip.time/t12, zs_dip,
#     '--k', linewidth=2,
#     label='Stefan Model with Dilatation',
#     zorder=5
# )

# legends
horizontal = plt.Line2D(
    (0, 1), (0, 0), color=colors['horizontal'], linestyle='-',
)
inclined = plt.Line2D(
    (0, 1), (0, 0), color=colors['inclined'], linestyle='-',
)
random = plt.Line2D(
    (0, 1), (0, 0), color=colors['random'], linestyle='-',
)
outside = plt.Line2D(
    (0, 1), (0, 0), linestyle='',
    marker='o', markerfacecolor='none', markeredgecolor=colors['horizontal'],
)
inside = plt.Line2D(
    (0, 1), (0, 0), linestyle='',
    marker='o', c=colors['horizontal'],
)
deg0 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[0], c='k', markerfacecolor='none'
)
deg10 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[10], c='k', markerfacecolor='none'
)
deg20 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[20], c='k', markerfacecolor='none'
)
deg30 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[30], c='k', markerfacecolor='none'
)
deg45 = plt.Line2D(
    (0, 1), (0, 0), ls='', marker=symbols[45], c='k', markerfacecolor='none'
)
stefanline = plt.Line2D(
    (0, 1), (0, 0), ls='-', c='k'
)
stefandilatationline = plt.Line2D(
    (0, 1), (0, 0), ls='--', c='k'
)

# ax0.legend(
#     [
#         horizontal, deg0,
#         stefanline, stefandilatationline,
#     ],
#     [
#         'Horizontal',
#         r'$\alpha$ = 0°',
#         'Stefan Model',
#         'Stefan Model with dilatation',
#     ],
#     fancybox=True, fontsize=10, ncol=2,
# )
# ax1.legend(
#     [
#         inclined, deg0, deg10, deg20, deg30, deg45,
#         stefanline, stefandilatationline,
#     ],
#     [
#         'Inclined',
#         r'$\alpha$ = 0°', r'$\alpha$ = 10°', r'$\alpha$ = 20°',
#         r'$\alpha$ = 30°', r'$\alpha$ = 45°',
#         'Stefan Model',
#         'Stefan Model with dilatation',
#     ],
#     fancybox=True, fontsize=10, ncol=2,
# )
# ax2.legend(
#     [
#         random, deg0, deg10, deg20, deg30, deg45,
#         stefanline, stefandilatationline,
#     ],
#     [
#         'Unsymmetrical',
#         r'$\alpha$ = 0°', r'$\alpha$ = 10°', r'$\alpha$ = 20°',
#         r'$\alpha$ = 30°', r'$\alpha$ = 45°',
#         'Stefan Model',
#         'Stefan Model with dilatation',
#     ],
#     fancybox=True, fontsize=10, ncol=2,
# )
# ax3.legend(
#     [
#         horizontal, inclined, random, deg0, deg10, deg20, deg30, deg45,
#         stefanline, stefandilatationline,
#     ],
#     [
#         'Horizontal',
#         'Inclined',
#         'Unsymmetrical',
#         r'$\alpha$ = 0°', r'$\alpha$ = 10°', r'$\alpha$ = 20°',
#         r'$\alpha$ = 30°', r'$\alpha$ = 45°',
#         'Stefan Model',
#         'Stefan Model with dilatation',
#     ],
#     fancybox=True, fontsize=10, ncol=2,
# )

# # labels
# ax0.set_xlabel(r'$t/t_{z_0}$', fontsize=18)
# ax0.set_ylabel(r'$z_f/z_0$', fontsize=18)
#
# ax1.set_xlabel(r'$t/t_{z_0}$', fontsize=18)
# ax1.set_ylabel(r'$z_f/z_0$', fontsize=18)
#
# ax2.set_xlabel(r'$t/t_{z_0}$', fontsize=18)
# ax2.set_ylabel(r'$z_f/z_0$', fontsize=18)
#
# ax3.set_xlabel(r'$t/t_{z_0}$', fontsize=18)
# ax3.set_ylabel(r'$z_f/z_0$', fontsize=18)

ax00.plot([1, 10], [1, 10], '--k')
ax01.plot([5, 30], [5, 30], '--k')
ax02.plot([5, 30], [5, 30], '--k')

ax10.plot([0, 1], [1, 1], '--k')
ax11.plot([20, 90], [1, 1], '--k')
ax12.plot([10, 30], [1, 1], '--k')
ax12.plot([.01, .09], [1, 1], '--k')

[a.grid(True) for a in axes]

ax03.legend(
    [
        horizontal, inclined, random, deg0, deg10, deg20, deg30, deg45,
        inside, outside,
        stefanline,
    ],
    [
        'Horizontal',
        'Inclined',
        'Unsymmetrical',
        r'$\alpha$ = 0°', r'$\alpha$ = 10°', r'$\alpha$ = 20°',
        r'$\alpha$ = 30°', r'$\alpha$ = 45°',
        'Open : inside',
        'Close : outside',
        'Stefan Model',
    ],
    fancybox=True, fontsize=10, ncol=2,
)

ax00.set_xlabel(r'$t_{{1/2}_{th}}$')
ax01.set_xlabel(r'$t_{{z0}_{th}}$')
ax02.set_xlabel(r'$t_{{f}_{th}}$')
ax10.set_xlabel(r'$a = z_0/r_0$')
ax11.set_xlabel(r'$\theta (°)$')
ax12.set_xlabel(r'$V (mm^3)$')
ax13.set_xlabel(r'$S_{te}$')

ax00.set_ylabel(r'$t_{{1/2}_{exp}}$')
ax01.set_ylabel(r'$t_{{z0}_{exp}}$')
ax02.set_ylabel(r'$t_{{f}_{exp}}$')
ax10.set_ylabel(r'$\frac{t_{{z0}_{exp}} }{ t_{{z0}_{ths}}}$')
ax11.set_ylabel(r'$\frac{t_{{z0}_{exp}} }{ t_{{z0}_{ths}}}$')
ax12.set_ylabel(r'$\frac{t_{{z0}_{exp}} }{ t_{{z0}_{ths}}}$')
ax12.set_ylabel(r'$\frac{t_{{z0}_{exp}} }{ t_{{z0}_{ths}}}$')

# show
print(c)
print(compteur_sec)
plt.show()
