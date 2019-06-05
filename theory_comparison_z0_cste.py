# -*- coding: utf-8 -*-
"""
Create a rectangle on an image.

Authors
-------
eablonet

"""

import Stefan as ste
from matplotlib import pyplot as plt
from numpy import arange, argmin

# paramaters
# ----------
z0 = 2e-3
Tsub = -15
Tm = 0
Tini = 0
Tinf = 20
z0 = 2e-3
r0 = arange(1e-3, 3.1e-3, 5e-4)
idx = argmin(abs(z0-r0))
print('Demi-sphere : (z0,r0) = ({:.1e},{:.1e})'.format(
    z0, r0[idx]
))

classic = ste.Stefan()
classic.import_material('ice')
classic.import_geometry(z0=z0, Nz=100, T_down=-15, T_up=0, T_ini=20)
classic.define_stop_condition('zf', Nt=1000)

dip = ste.Diphasique()
dip.import_material('ice', 'water')
dip.import_geometry(z0=z0, Nz=100, T_down=-15, T_up=20, T_ini=20, T_m=0)
dip.define_stop_condition('zf')

ptd = ste.PseudoTwoDim()
ptd.import_material('ice')


plt.figure(figsize=[8, 4.5])
plt.plot(classic.time, classic.front_position(), '-k', label='Quasi-statique')
plt.plot(dip.time, dip.front_position(), '-r', label='Full Stefan')
plt.grid(True)
plt.xlabel('t (s)')
plt.ylabel(r'$z_s$ (m)')
plt.grid(True)

theta = []
for i in range(len(r0)):
    ptd.import_geometry(z0=z0, Nz=100, T_down=-15, T_up=0, T_ini=20, r0=r0[i])
    ptd.define_stop_condition('z0', dt=.001)
    z = ptd.front_position()  # to generate the time
    theta.append(ptd.theta)
    plt.plot(ptd.time, z, '--g')
    plt.text(
        ptd.time[-1], z[-1], str(r'$r_0 =$ {:.1e}'.format(r0[i])),
        fontsize=8, color='k', horizontalalignment='center',
    )
plt.legend(fancybox=True)

plt.figure()
plt.plot(r0, theta, '--k')
plt.xlabel(r'$r_0$ (m)')
plt.ylabel(r'$\theta$ (°)')
plt.grid(True)

# plt.figure()
# plt.plot(st.z, st.temp_field(True)[:, 30])
# plt.plot(dip.z, dip.temp_field(True)[:, 30])
# plt.plot(
#     [dip.front_position(adim=True)[30], dip.front_position(adim=True)[30]],
#     [0, 1],
#     '--b'
# )
#
# plt.figure()
# plt.plot(st.time, st.temp_field(True)[30, :])
# plt.plot(dip.time, dip.temp_field(True)[30, :])
plt.show()
