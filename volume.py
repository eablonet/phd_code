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
"""
Physical Parameters
"""
date = '29-03-2018'
serie = 9

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
# manip_data.info()
s = st.Stack()

for _, row in manip_data.iterrows():
    cond = (
        row['date'] == date and
        row['serie'] == serie
        # row['completed_2'] == 1 and  # front detected
    )
    if cond:
        print('Found')
        # load data
        px_mm = row['px_mm']
        z0 = row['z0']/px_mm
        print(z0)
        r0 = row['r0']/px_mm
        rc = int(row['rc'])
        theta = 1/2*(row['ca_left'] + row['ca_right'])
        T_nuc = row['Tc_nuc']
        Ta = row['Ta_int']
        fps = row['fps_real']
        t_nuc = row['t_nuc_calc']
        t_end = row['t_end_calc']
        t_ref = row['t_ref_calc']+1
        dt_ini = t_nuc - t_ref

        # read front information
        s.read_by_path(row['date'], int(row['serie']))
        s.read_image(2)  # to load an image else none might be load...to correct
        # al, zf = s.read_manual_front()
        zf = s.read_data()
        zf_id = zf[int(t_ref-1):, :]
        zf_loc = np.zeros(len(zf[:, 0]))

        for i in range(1, len(zf_id[:, 0])):
            if zf_id[i, rc] != 0:
                zf_loc[i+int(t_ref-1)] = zf_id[i, rc]
            else:
                zf_loc[i+int(t_ref-1)] = np.nan
        zf_loc[zf_loc == np.inf] = np.nan

        t_space, zf, zf_mean, zf_std = s.get_dynamic_front(zf)
        zf = np.max(zf_mean) / px_mm

        time = np.arange(len(zf_mean))/fps

        volume, Cl, Cr = s.read_sa_contour()
        # rc_l = min(Cl[t_nuc, :])
        # rc_r = max(Cr[t_nuc, :])

        volume /= px_mm**3

        volume_ice = np.zeros_like(volume)

        for t in range(s.n_image_tot):
            for i in range(s.current_image.size[0]):
                if i > zf_loc[t] and Cl[t, i] > 0:
                    volume_ice[t] += (rc - Cl[t, i])**2
                if i > zf_loc[t] and Cr[t, i] > 0:
                    volume_ice[t] += (Cr[t, i] - rc)**2

        volume_ice *= np.pi/2
        volume_ice /= px_mm**3

        # theorical data
        mon = ste.Stefan()
        mon.import_material('ice')
        mon.import_geometry(z0=z0, Nz=100, T_down=T_nuc, T_up=0, T_ini=0)
        mon.define_stop_condition('zf', Nt=len(zf_mean))

        zs_comp = np.sqrt(2 * 2.2 * (0-T_nuc) * time / (916.2*333500))

        dip = ste.Diphasique()
        dip.import_geometry(z0=z0, Nz=100, T_down=T_nuc, T_up=Ta, T_ini=0, T_m=0)
        dip.import_material('ice', 'water')
        dip.define_stop_condition('zf', Nt=len(zf_mean))

        # ptd = ste.PseudoTwoDim()
        # ptd.import_geometry(z0=z0, Nz=100, T_down=T_nuc, T_up=0, T_ini=0, r0=r0)
        # ptd.import_material('ice')
        # ptd.define_stop_condition('zf', dt=.1)
        # ptd.front_position()

        # print('Spherical error :theta_mes - theta_th = ', abs(theta-ptd.theta))

        zs_ste_real = mon.front_position(adim=False)
        zs_ste2_real = dip.front_position(adim=False)
        # zs_steptd_real = ptd.front_position(adim=False)

        # ax : plot in real coordinates
        fig = plt.figure(figsize=(8, 4.5))
        fig.canvas.set_window_title(date+str(serie))
        ax = plt.axes()

        _, caplines, _ = ax.errorbar(
            time - dt_ini/fps,
            zf_mean/px_mm * 1e3,
            yerr=(zf_std+5)/px_mm * 1e3,
            ls='none', c='k',
            marker='o', mec='k', mfc='none', ms=4, markevery=1,
            ecolor='k', barsabove=True, elinewidth=.8,
            lolims=True, uplims=True,
            label='Experimental'
        )
        ax.plot(
            mon.time, zs_ste_real * 1e3,
            ls='--', c='b',
            label='Stefan model',
        )
        """
        ax.plot(
            dip.time, zs_ste2_real * 1e3,
            ls='--', c='c',
            label='Stefan 1D transient-diphasic with dilatation',
        )
        """

        """
        ax.plot(
            ptd.time, zs_steptd_real * 1e3,
            ls='--', color='y',
            label='Pseudo-2D'
        )
        """

plt.legend(fancybox=True)
plt.xlabel('t (s)')
plt.ylabel(r'$x_f$ (mm)')
plt.grid(True)

V0 = volume[volume > 0][0]
print(V0)

t = np.arange(s.n_image_tot)/fps-t_nuc

plt.figure(figsize=[8, 6])
plt.plot(
    np.arange(s.n_image_tot)/fps-t_nuc, volume,
    '+b', ms=8,
    label='Contour'
)
plt.plot(
    t[:int(t_end)-3], V0 + volume_ice[:int(t_end)-3]*(1 - 916/998),
    '+r', ms=8,
    label='Front'
)
plt.grid(True)
plt.legend(fancybox=True,  fontsize=18)
plt.ylim([.98*V0, 1.02*max(volume)])
plt.xlim([.98*t_ref - t_nuc, 1.02*t_end - t_nuc])
plt.xlabel('time (s)', fontsize=18)
plt.ylabel(r'Volume ($m^3$)',  fontsize=18)
plt.show()
