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

axes: all represents tz0_exp/t_z0_th over one data
----
    nb: 6
    0: versus z0
    1: versus r0
    2: versus volume
    3: versus a = z0/r0
    4: versus \theta (contact angle)
    5: versus \alpha (inclination)

"""
fig0 = plt.figure(figsize=[16, 9])
fig0.canvas.set_window_title('t_z0 for Geometric')
axes0 = []
for i in range(6):
    axes0.append(fig0.add_subplot(2, 3, i+1))


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
    't_z0_exp': [],
    't_z0_th': [],
    'z0': [],  # mm
    'r0': [],  # mm
    'theta': [],  # contact angle
    'alpha': [],  # inclination
    'Tnuc': [],
    'Tset': [],
    'inside': [],  # bool
    'ref_therm': [],  # constant or ramp
    'Ta_set': [],  # nan for outside
    'Ta_cc': [],  # nan for outside
    'Ta': [],  # could be nan
    'phi_set': [],  # nan for outside
    'phi_cc': [],  # nan for outside
    'phi': [],  # could be nan
    'Ste': []
}
for _, row in manip_data.iterrows():
    cond = (
        row['type'] == 'a'
    )
    if cond:

        # load geometric data
        # -------------------
        px_mm = row['px_mm']

        z0 = row['z0']/px_mm
        r0 = row['r0']/px_mm
        theta = 1/2*(row['ca_left'] + row['ca_right'])
        alpha = row['alpha']

        data['z0'].append(z0)
        data['r0'].append(r0)
        data['theta'].append(theta)
        data['alpha'].append(alpha)

        # load thermal data
        # -----------------
        if row['ref_therm'] in ['Tw-7', 'Tw-10', 'Tw-12', 'Tw-15']:
            data['ref_therm'].append('constant')
        elif row['ref_therm'] in ['qw1', 'qw2', 'qw3', 'qw4']:
            data['ref_therm'].append('ramp')
        else:
            data['ref_therm'].append('other')
        Tc_nuc = row['Tc_nuc']
        Tc_set = row['Tc_set']
        Ta_set = row['Ta']
        Ta_cc = row['Ta_int']
        Ta = row['Ta_ext']
        phi_set = row['phi']
        phi_cc = row['phi_int']
        phi = row['phi_ext']

        data['Tnuc'].append(Tc_nuc)
        data['Tset'].append(Tc_set)
        data['Ta_set'].append(Ta_set)
        data['Ta_cc'].append(Ta_cc)
        data['Ta'].append(Ta)
        data['phi_set'].append(phi_set)
        data['phi_cc'].append(phi_cc)
        data['phi'].append(phi)

        # load temporal data
        # ------------------
        fps = row['fps_real']
        t_nuc = row['t_nuc_calc']
        t_end = row['t_end_calc']
        t_ref = row['t_ref_calc']
        dt_ini = t_nuc - t_ref

        # other data
        # ----------
        data['inside'].append(True if row['lieu'] else False)

        # calculate t_z0_exp
        # ------------------
        s.read_by_date(row['date'], int(row['serie']))
        zf = s.read_data()
        _, _, zf_mean, _ = s.get_dynamic_front(zf)

        time = np.arange(len(zf_mean))/fps
        time_scale = time - dt_ini/fps  # on enlève le temps avant tnuc

        t_z0_exp = time_scale[np.argmin(abs(zf_mean / px_mm / z0 - 1))]
        data['t_z0_exp'].append(t_z0_exp)

        # calculate t_z0_th
        # -----------------
        mon = ste.Stefan()
        mon.import_material('ice')
        mon.import_geometry(
            z0=z0, Nz=100, T_down=Tc_nuc if ~np.isnan(Tc_nuc) else Tc_set,
            T_up=0, T_ini=0
        )
        mon.define_stop_condition('zf', Nt=100)
        zs_mon = mon.front_position(adim=True)

        t_z0_th = mon.time[np.argmin(abs(zs_mon - 1))]

        data['Ste'].append(mon.calc_St())
        data['t_z0_th'].append(t_z0_th)

# add figure design
# -----------------
r_ranges = np.arange(1.5, 5, .2)*1e-3
c = ea_color.color_grad(len(r_ranges), 'Reds')
each_r_range = []
for r in data['r0']:
    j = 0
    while j < len(r_ranges) and r > r_ranges[j]:
        j += 1
    each_r_range.append(j)


def make_selection(data, crit='phi', min=None, max=None):
    c = 0
    for i in range(len(data[crit])):
        if (
            (max is not None and data[crit][i-c] > max) or
            (min is not None and data[crit][i-c] < min)
        ):
            for a in data:
                del(data[a][i-c])
            c += 1
    return data


def make_range(data, crit='Ta', value=[10, 15, 20, 25, 30], ecart=2.5):
    """Create range for value."""
    cat = np.zeros(len(data[crit]))
    for k in range(len(data[crit])):
        j = 0
        while j < len(value) and (
            data[crit][k] < value[j]-ecart or data[crit][k] > value[j]+ecart
        ):
            j += 1
        if j == len(value):
            cat[k] = 0
        else:
            cat[k] = value[j]
    return cat


# data = make_selection(data, 'phi', max=45)
# data = make_selection(data, 'phi_set', max=45)
# data = make_selection(data, 'r0', min=2.7e-3, max=3.2e-3)
# data = make_selection(data, 'alpha', max=5)

T_range = make_range(data, crit='Ta', value=[10, 15, 20, 25, 30], ecart=2.5)
color = ea_color.color_met(6)
T_color = {
    0: color[0], 10: color[1], 15: color[2],
    20: color[3], 25: color[4], 30: color[5]
}

phi_range = make_range(data, crit='phi', value=[30, 40, 50, 60, 70], ecart=5)
marker_unfill = {30: '+', 40: '2', 50: '1', 60: '3', 70: '4', 0: 'x'}
marker_fill = {30: 'P', 40: '^', 50: 'v', 60: '<', 70: '>', 0: 'X'}

mkr = []
for r in range(len(phi_range)):
    if data['ref_therm'][r] == 'ramp':
        mkr.append(marker_fill[phi_range[r]])
    elif data['ref_therm'][r] == 'constant':
        mkr.append(marker_unfill[phi_range[r]])
    else:
        mkr.append('s')

print(mkr)

# mkr = []
# for i, a in enumerate(data['inside']):
#     if a:
#         mkr.append('o')
#     else:
#         mkr.append(marker[data['ref_therm'][i]])

for k in range(len(data['t_z0_exp'])):
    # fig1 plotting
    # -------------
    axes0[0].plot(
        data['z0'][k]*1e3,
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    axes0[1].plot(
        data['r0'][k]*1e3,
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    axes0[2].plot(
        np.pi/6*data['z0'][k]*1e3 * (
            3*np.square(data['r0'][k]*1e3) +
            np.square(data['z0'][k]*1e3)
        ),
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    axes0[3].plot(
        data['z0'][k]/data['r0'][k],
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    axes0[4].plot(
        data['theta'][k],
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    axes0[5].plot(
        data['alpha'][k],
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    # fig2 plotting
    # -------------
    axes1[0].plot(
        data['Ta_set'][k],
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
        ms=6,
    )

    axes1[1].plot(
        data['Ta_cc'][k],
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    axes1[2].plot(
        data['Ta'][k],
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    axes1[3].plot(
        data['phi_set'][k],
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    axes1[4].plot(
        data['phi_cc'][k],
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    axes1[5].plot(
        data['phi'][k],
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    axes1[6].plot(
        -data['Tnuc'][k],
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    axes1[7].plot(
        -data['Tset'][k],
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

    axes1[8].plot(
        data['Tnuc'][k] - data['Tset'][k],
        data['t_z0_exp'][k]/data['t_z0_th'][k],
        color=T_color[T_range[k]],
        marker=mkr[k],
        mfc=T_color[T_range[k]] if data['inside'][k] else 'none',
    )

# labler fig1
# -----------
axes0[0].set_xlabel(r'$z_0$ (mm)')
axes0[1].set_xlabel(r'$r_0$ (mm)')
axes0[2].set_xlabel(r'volume ($mm^3$)')
axes0[3].set_xlabel(r'a = $z_0/r_0$ ($mm^3$)')
axes0[4].set_xlabel(r'$\theta$ (°)')
axes0[5].set_xlabel(r' $\alpha$ (°)')

# label fig2
# ----------
axes1[0].set_xlabel(r'$T_{a_{setpoint}}$ (°C)')
axes1[1].set_xlabel(r'$T_{a_{cc}}$ (°C)')
axes1[2].set_xlabel(r'$T_{a}$ (°C)')
axes1[3].set_xlabel(r'$\phi_{setpoint}$ (%)')
axes1[4].set_xlabel(r'$\phi_{cc}$ (%)')
axes1[5].set_xlabel(r'$\phi$ (%)')
axes1[6].set_xlabel(r'$\Delta T_{c_{nuc}}$ (°C)')
axes1[7].set_xlabel(r'$\Delta T_{c_{set}}$ (°C)')
axes1[8].set_xlabel(r'$\Delta T_{c_{set}} - \Delta T_{c_{nuc}}$ (°C)')

[a.grid(True) for a in axes0+axes1]
[a.set_ylabel(r'$t_{z0_{exp}} / t_{z0_{th}}$') for a in axes0+axes1]
[a.plot(a.get_xlim(), [1.5, 1], '-k') for a in axes0]

# show
plt.show()
