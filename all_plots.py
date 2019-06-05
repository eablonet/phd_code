"""
author : eablonet
date : 2019
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import scipy.optimize as op
import pandas as pd
import StackEA as st
import psat_tsat as psts
import Stefan as ste
import Material as ma
"""
Physical Parameters
"""
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

"""
StefanSolver
"""


def calc_tf(
    z0, dT, dTl, phi=0,
    x0=1
):
    """
    Calculate final time of solidication with Stefan theory.

    tf : monophasic 1D
    tf2 : diphasic 1D with dilatation
    """
    dT = np.array(dT)
    dTl = np.array(dTl)

    def f(x, St):
        return x*erf(x)*np.exp(x**2) - St/np.sqrt(np.pi)

    def f2(x, Sts, Stl, rho, Dl):
        return (
            x -
            Sts/np.sqrt(np.pi)*np.exp(-x**2)/erf(x) -
            Stl/np.sqrt(Dl*np.pi)*np.exp(-x**2*(1+rho)**2/(4*Dl)) /
            erf(x*(1+rho)/(2*np.sqrt(Dl)))
        )

    if dT.size > 1:
        dT, dTl = np.meshgrid(dT, dTl)

        St = cp*dT/Lf
        Stl = k_l/k * cp*dTl/Lf
        Dl = k_l * rho * cp / (k * rho_l * cp_l)
        # if dTl.size > 1:
            # T_int, _ = psts.table(Ta, phi)
        delta = np.empty_like(dT)
        delta2 = np.empty_like(dT)

        for row, val0 in enumerate(St):
            for col, val1 in enumerate(val0):
                delta[row, col] = op.fsolve(
                    f, x0, args=(val1)
                )
                delta2[row, col] = op.fsolve(
                    f2, x0,
                    args=(val1, Stl[row, col], rho/rho_l, Dl)
                )
    else:
        St = cp*dT/Lf
        Stl = k_l/k * cp*dTl/Lf
        Dl = k_l * rho * cp / (k * rho_l * cp_l)

        delta = op.fsolve(f, x0, args=(St))
        delta2 = op.fsolve(f2, x0, args=(St, Stl, rho/rho_l, Dl))

    tf = rho*cp*z0**2/(4*k*delta**2)
    tf2 = rho_l*cp*z0**2/(4*k*delta2**2)

    return tf, tf2, delta, delta2, St, Stl


def fit(
    x, a
):
    """Fit data."""
    return 4*a**2*k/(rho*cp)*(x)


manip_data = pd.read_csv('manips.csv', skiprows=[0, 2], decimal=",")
# import datas
manip_data.info()
# print info of table

s = st.StackEA()
# Create stack player

data = {
    'St': [],
    't_tot': [], 'tf': [], 'tf2': [], 't_nuc': [],
    'T_nuc': [],
    'delta': [], 'delta2': [], 'deltafit': [],
    'inside': [], 'cst': [],
    'zf': [],
    'phi': [], 'phi_ref': [],
}

# fig zs vs tf
fig0 = plt.figure(figsize=[8, 4.5])
ax0 = plt.axes()

# fig zs vs tf2
fig01 = plt.figure(figsize=[8, 4.5])
ax01 = plt.axes()

# fig zs vs t_tot
fig02 = plt.figure(figsize=[8, 4.5])
ax02 = plt.axes()

for _, row in manip_data.iterrows():
    cond = (
        row['completed_2'] == 1 and  # front detected
        row['alpha'] < 1 and   # horizontal
        row['type'] == 'a'  # drop solidification
    )
    if cond:
        ref = row['date'] + 'n' + str(int(row['serie']))

        # load recording data
        inside = bool(row['lieu'])
        fps = row['fps_real']
        px_mm = row['px_mm']

        # load thermals data
        T_nuc = row['Tc_nuc']
        Ta = row['Ta_int']
        T_set = row['Tc_set']
        T_sta = row['Tc_sta']
        phi = row['phi_int']
        phi_ref = row['phi']
        if np.isnan(T_nuc):
            T_nuc = T_set
        if np.isnan(Ta):
            Ta = row['Ta_ext']
        if np.isnan(phi):
            phi = row['phi_ext']

        # load times data
        t_nuc = row['t_nuc_calc']
        t_end = row['t_end_calc']
        t_ref = row['t_ref_calc']+1
        dt_ini = t_nuc - t_ref
        t_tot = t_end - t_nuc
        t_nuc_real = row['t_nuc_mes']

        # load experimental shape
        s.read_by_path(row['date'], int(row['serie']))
        s.read_image(100)
        al, zf = s.read_manual_front()
        t_space, zf, zf_mean, zf_std = s.get_dynamic_front(zf)
        zf = np.max(zf_mean) / px_mm

        # load geometricals data
        z0 = row['z0']/px_mm
        r0 = row['r0']/px_mm
        volume_data = row['v_theo_mm']
        volume = np.pi/6 * z0 * (z0**2 + 3*r0**2)/(1e-3)**3
        theta = 1/2*(row['ca_left'] + row['ca_right'])
        if np.isnan(theta):
            theta = 60.
        print(theta)

        # tf, tf2, delta, delta2, St, Stl = calc_tf(z0, -T_nuc, -Ta)
        mon = ste.Stefan()
        mon.import_material('ice')
        mon.import_geometry(z0=z0, Nz=100, T_down=T_nuc, T_up=0, T_ini=20)
        mon.define_stop_condition('z0', Nt=len(zf_mean))

        dip = ste.Diphasique()
        dip.import_geometry(z0=z0, Nz=100, T_down=T_nuc, T_up=Ta, T_ini=Ta, T_m=0)
        dip.import_material('ice', 'water')
        dip.define_stop_condition('z0')

        ptd = ste.PseudoTwoDim()
        ptd.import_geometry(z0=z0, Nz=100, T_down=T_nuc, T_up=0, T_ini=0, r0=r0, theta=theta)
        ptd.import_material('ice')
        ptd.define_stop_condition('z0', dt=.001)
        ptd.front_position()

        delta = mon.solve_delta()
        tf = mon.calc_tf(adim=False)
        tf12 = mon.calc_t12(adim=False)
        delta2 = dip.solve_delta()
        tf2 = dip.calc_tf(adim=False)/4
        St = dip.calc_St()
        Stl = dip.calc_St_liquid()

        t = np.linspace(0, t_tot/fps, 100)

        # real Stefan
        zs_ste_real = mon.front_position(adim=False)
        zs_ste2_real = dip.front_position(adim=False)
        zs_steptd_real = ptd.front_position(adim=False)

        # adimensional Stefan
        zs_ste_adim = mon.front_position(adim=True)
        zs_ste2_adim = dip.front_position(adim=True)
        zs_steptd_adim = ptd.front_position(adim=True)

        # fit of the experiment
        popt, pcov = op.curve_fit(
            fit,
            t_space[int(dt_ini):int(.85*t_tot)]/fps,
            (zf_mean[int(dt_ini):int(.85*t_tot)]/px_mm)**2
        )
        zs_fit_real = 2*popt[0]*np.sqrt(k*(t)/(rho*cp))
        zs_fit_adim = 2*popt[0]*np.sqrt(k*(t)/(rho*cp))/z0

        # store all information
        t0 = z0**2 * rho*cp/k
        data['St'].append(St)
        data['t_tot'].append(t_tot/fps/t0)
        data['tf'].append(tf/t0)
        data['tf2'].append(tf2/t0)
        data['delta'].append(delta)
        data['delta2'].append(delta2)
        data['deltafit'].append(popt[0])
        data['inside'].append(inside)
        data['cst'].append(1 if abs(T_set - T_sta) < 1e-5 else 0)
        data['zf'].append(zf)
        data['t_nuc'].append(t_nuc_real)
        data['T_nuc'].append(T_nuc)
        data['phi'].append(phi)
        data['phi_ref'].append(phi_ref)

        # create text information
        text = (
                "Lieu : {:}" + ", " +
                r"$T_{{nuc}}$ : " + "{:.1f}°C" + ", " +
                r"$T_a$ : " + "{:.1f}°C" + ", " +
                r"$T_{{sta}}$ : " + "{:.1f}°C" + ", " +
                r"$T_{{set}}$ : " + "{:.1f}°C" + ", " +
                r"$\phi$ : " + "{:.1f}%" + ",\n" +
                r"$t_{{f}}$ : " + "{:.1f}s" + ", " +
                r"$t_{{f2}}$ : " + "{:.1f}s" + ", " +
                r"$\delta$ : " + "{:.1e}" + ", " +
                r"$\delta_2$ : " + "{:.1e}" + ", " +
                r"$S_t$ : " + "{:.1e}" + ", " +
                r"$S_{{t_l}}$ : " + "{:.1e}" + ",\n" +
                r"V : " + "{:2.1e}" + r"$mm^3$" + ", " +
                r"$z_0$ : " + "{:.1f}mm" + ", " +
                r"$z_{{0}}/r_0$ : " + "{:.2f}"
            ).format(
                'Int' if inside else 'Ext',
                T_nuc,
                Ta,
                T_sta,
                T_set,
                phi,
                tf,
                tf2[0],
                delta,
                delta2[0],
                St,
                -Stl,
                volume,
                z0*1e3,
                z0/r0
            )

        img_txt = (
                folder +
                "i{:}_" +
                "{:}_" +
                'St{:.1e}_'
                'id{:}'
            ).format(
                '0' if inside else '1',
                'cst' if abs(T_set - T_sta) < 1e-5 else 'var',
                St,
                ref
            )

        # ax : plot in real coordinates
        fig = plt.figure(figsize=(8, 4.5))
        fig.canvas.set_window_title(ref)
        ax = plt.axes()

        _, caplines, _ = ax.errorbar(
            t_space/fps - dt_ini/fps,
            zf_mean/px_mm * 1e3,
            yerr=(zf_std+5)/px_mm * 1e3,
            ls='none', c='k',
            marker='o', mec='k', mfc='none', ms=4, markevery=1,
            ecolor='k', barsabove=True, elinewidth=.8,
            lolims=True, uplims=True,
            label='Experimental'
        )
        ax.plot(
            t_space, zs_ste_real * 1e3,
            ls='--', c='b',
            label='Stefan 1D transient-monophasic',
        )
        ax.plot(
            t, zs_ste2_real * 1e3,
            ls='--', c='c',
            label='Stefan 1D transient-diphasic with dilatation',
        )
        ax.plot(
            t, zs_fit_real * 1e3,
            ls='--', color='r',
            label='Best fit'
        )
        ax.plot(
            ptd.time, zs_steptd_real * 1e3,
            ls='--', color='y',
            label='Pseudo-2D'
        )

        caplines[0].set_marker('_')
        caplines[1].set_marker('_')
        ax.grid(True)
        ax.set_xlabel('t (s)')
        ax.set_ylabel('z (mm)')

        ax.legend(fancybox=True, shadow=True)

        fig.suptitle(text)

        if True:
            fig.savefig(
                img_txt+'.png', dpi=300, facecolor='w', edgecolor='w',
                format='png',
            )

        # ax1 : plot in adimensional coordinates
        fig1 = plt.figure(figsize=(8, 4.5))
        fig1.canvas.set_window_title(ref)
        ax1 = plt.axes()
        idx = np.argmin(abs(zf_mean/px_mm/z0 - .5))
        t_adim = t_space[idx] - dt_ini
        data['t_tot'].append(t_tot/fps/t_adim)
        _, caplines2, _ = ax1.errorbar(
            t_space/t_adim - dt_ini/t_adim, zf_mean/px_mm/z0,
            yerr=(zf_std+5)/px_mm/z0,
            ls='none', c='k',
            marker='o', mec='k', mfc='none', ms=4, markevery=1,
            ecolor='k', barsabove=True, elinewidth=.8,
            lolims=True, uplims=True,
            label='Experimental'
        )
        ax1.plot(
            [0, 1, 1],
            [.5, .5, 0],
            '-g'
        )
        if inside:
            color = 'orange'
        else:
            if abs(T_set - T_sta):
                color = 'red'
            else:
                color = 'm'
        _, caplines0, _ = ax0.errorbar(
            t_space/t_adim - dt_ini/t_adim, zf_mean/px_mm/z0,
            yerr=(zf_std)/px_mm/z0,
            ls='none', mec=color,
            marker='o', mfc='none', ms=4, markevery=1,
            ecolor=color, barsabove=True, elinewidth=.8,
            lolims=True, uplims=True,
        )
        _, caplines01, _ = ax01.errorbar(
            t_space/t_adim - dt_ini/t_adim, zf_mean/px_mm/z0,
            yerr=(zf_std)/px_mm/z0,
            ls='none', mec=color,
            ecolor=color, marker='o', mfc='none', ms=4, markevery=1,
            barsabove=True, elinewidth=.8,
            lolims=True, uplims=True,
        )
        _, caplines02, _ = ax02.errorbar(
            t_space/t_adim - dt_ini/t_adim, zf_mean/px_mm/z0,
            yerr=(zf_std)/px_mm/z0,
            ls='none', mec=color,
            ecolor=color, marker='o', mfc='none', ms=4, markevery=1,
            barsabove=True, elinewidth=.8,
            lolims=True, uplims=True,
        )

        ax1.plot(
            t_space/t_adim, zs_ste_adim,
            ls='--', c='b',
            label='Stefan 1D transient-monophasic',
        )
        ax1.plot(
            t/t_adim*fps, zs_ste2_adim,
            ls='--', c='c',
            label='Stefan 1D transient-diphasic with dilatation',
        )
        ax1.plot(
            t/t_adim*fps, zs_fit_adim,
            ls='--', color='r',
            label='Best fit'
        )
        ax1.plot(
            ptd.time/t_adim, zs_steptd_adim,
            '--y',
            label='Pseudo-2D correction'
        )

        caplines2[0].set_marker('_')
        caplines2[1].set_marker('_')
        caplines0[0].set_marker('_')
        caplines0[1].set_marker('_')
        caplines01[0].set_marker('_')
        caplines01[1].set_marker('_')
        caplines02[0].set_marker('_')
        caplines02[1].set_marker('_')
        ax1.grid(True)
        ax1.set_xlabel(r'$t/t_{1/2}$')
        ax1.set_ylabel('z/z0')

        ax1.legend(fancybox=True, shadow=True)

        fig1.suptitle(text)

        if True:
            fig1.savefig(
                img_txt + '_adim.png', dpi=300, facecolor='w', edgecolor='w',
                format='png',
            )

        plt.close(fig)
        plt.close(fig1)


idx = sorted(range(len(data['St'])), key=data['St'].__getitem__)

# plot t_tot vs St
fig = plt.figure(figsize=[8, 4.5])
ax = plt.axes()

"""
sort data:
0 - outside - all phi - 'sk'
1 - inside - phi = 40 - 'sb'
2 - inside - phi = 50 - '^b'
3 - inside - phi = 60 - 'vb'
4 - inside - phi = 70 - 'ob'
"""
ls = np.empty_like(idx)
so = {
    'St0': [], 'St1': [], 'St2': [], 'St3': [], 'St4': [], 'St5': [], 'St6': [],
    't0': [], 't1': [], 't2': [], 't3': [], 't4': [], 't5': [], 't6': [],
}
for i in idx:
    if data['inside'][i]:
        if data['phi_ref'][i] == 40:
            so['St1'].append(data['St'][i])
            so['t1'].append(data['t_tot'][i])
        elif data['phi_ref'][i] == 50:
            so['St2'].append(data['St'][i])
            so['t2'].append(data['t_tot'][i])
        elif data['phi_ref'][i] == 60:
            so['St3'].append(data['St'][i])
            so['t3'].append(data['t_tot'][i])
        elif data['phi_ref'][i] == 70:
            so['St4'].append(data['St'][i])
            so['t4'].append(data['t_tot'][i])
    elif not data['inside'][i]:
        if data['phi'][i] >= 35 and data['phi'][i] < 45:
            so['St0'].append(data['St'][i])
            so['t0'].append(data['t_tot'][i])
        elif data['phi'][i] >= 45 and data['phi'][i] < 55:
            so['St5'].append(data['St'][i])
            so['t5'].append(data['t_tot'][i])
        elif data['phi'][i] >= 55 and data['phi'][i] < 65:
            so['St6'].append(data['St'][i])
            so['t6'].append(data['t_tot'][i])

print(so)
ls = ['sk', 'sb', '^b', 'vb', 'ob', '^k', 'vk']
lab = [
    r'Outside - $\phi \in$ $[35:45]%$',
    r'Inside - $\phi = 40%$',
    r'Inside - $\phi = 50%$',
    r'Inside - $\phi = 60%$',
    r'Inside - $\phi = 70%$',
    '_nolegend_',  # r'Outside - $\phi \in$ $[45:55]%$',
    '_nolegend_',  # r'Outside - $\phi \in$ $[55:65]%$',
]
for i in range(len(ls)):
    ax.plot(
        so['St'+str(i)], so['t'+str(i)],
        ls[i], mfc='none',
        label=lab[i],
    )

# theoritical values without phi
Tw = np.arange(-7, -16, -.5)
Ta = np.arange(0, 31, 5)

z0 = 2e-3

# tf, tf2, delta, delta2, St, Stl = calc_tf(z0, -Tw, -Ta)
tf2 = np.zeros([len(Tw), len(Ta)])
St = np.zeros([len(Tw), len(Ta)])
Stl = np.zeros([len(Tw), len(Ta)])

for i in range(len(Tw)):
    for j in range(len(Ta)):
        dip.import_geometry(
            z0=2e-3, Nz=100,
            T_down=Tw[i], T_up=Ta[j], T_ini=Ta[j], T_m=0
        )
        dip.import_material('ice', 'water')
        dip.define_stop_condition('z0')

        tf2[i, j] = dip.calc_tf(adim=True)
        St[i, j] = dip.calc_St()
        Stl[i, j] = dip.calc_St_liquid()

Tw, Ta = np.meshgrid(Ta, Tw)

CS = ax.contour(
    St, tf2, -Stl,
    6,
    colors='k',
    linestyles='dashed', linewidths=.8,
)
CS.collections[0].set_label(r'Iso-$St_l$')
manual_locations = [
    (St[2, 0], tf2[2, 0]),
    (St[2, 0], tf2[2, 1]),
    (St[2, 0], tf2[2, 2]),
    (St[2, 0], tf2[2, 3]),
    (St[2, 0], tf2[2, 4]),
    (St[2, 0], tf2[2, 5]),
]
ax.clabel(CS, fontsize=8, inline=True, manual=manual_locations, fmt='%.1e')

ax.plot(St[:, 0], tf2[:, 0], 'b', label='Limit value of Stefan solution with dilatation')
# ax.plot(St[0], tf[0]/(rho*cp*z0**2/k), 'r', label='Limit value of Stefan solution without dilatation')
ax.set_xlabel(r'$S_t$')
ax.set_ylabel(r'$t_f/t_0$')
ax.grid(True)
ax.legend(fancybox=True, shadow=True)

# plot delta vs St
plt.figure(figsize=[8, 4.5])
plt.plot(data['St'], data['delta'], '--sb', mfc='none')
plt.plot(data['St'], data['delta2'], '--sc', mfc='none')
plt.plot(data['St'], data['deltafit'], 'or', mfc='none')
plt.grid(True)
plt.xlabel('St')
plt.ylabel(r'$\delta$')

# plot T_nuc vs St
plt.figure(figsize=[8, 4.5])
plt.plot(data['St'], data['T_nuc'], 'sk', mfc='none')
plt.grid(True)
plt.xlabel(r'$S_t$')
plt.ylabel(r'$T_{nuc}$')

# plot t_nuc vs St
t_nuc = []
for txt in data['t_nuc']:
    if type(txt) is str:
        temp = txt.split("'")
        t_nuc.append(int(temp[0])*60 + int(temp[1]))
    else:
        t_nuc.append(np.nan)
plt.figure(figsize=[8, 4.5])
plt.plot(data['St'], t_nuc, '*b', mfc='none')
plt.grid(True)
plt.xlabel(r'$S_t$')
plt.ylabel(r'$t_{nuc}$')

ax0.plot(
    np.linspace(0, 1, 100), np.sqrt(np.linspace(0, 1, 100)), '--k',
    label='1D monophasic Stefan without dilatation'
)
ax0.grid(True)
ax0.set_xlabel(r't/$t_{f}$')
ax0.set_ylabel(r'z/$z_0$')
ax0.legend(fancybox=True)

ax02.plot(
    np.linspace(0, 1, 100), np.sqrt(np.linspace(0, 1, 100)), '--k',
    label='1D monophasic Stefan without dilatation'
)
ax02.grid(True)
ax02.set_xlabel(r't/$t_{t_{tot}}$')
ax02.set_ylabel(r'z/$z_0$')
ax02.legend(fancybox=True)

plt.show()
