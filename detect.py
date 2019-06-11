"""
# authors: eablonet
# date : 2019
"""
# global packages
import matplotlib.pyplot as plt
import numpy as np

# local packages
import packages.main.Stack as st

s = st.Stack()

s.read_by_path('20-11-2018', 5)
# s.current_image.show_image()

# s.remove_treatments()
# s.crop()
# s.clahe()

# s.tracker()
# s.contour_tracker2()

# zf, pts = s.read_data('contour')
# t, zf, zf_loc, _ = s.get_dynamic_front(zf, 300)
#
# plt.plot(range(0, len(zf_loc)), zf_loc)
# plt.show()
#
#
zf = s.read_data('front')
t, zf_loc, zf_mean, zf_std = s.get_dynamic_front(zf)
Lf = 333550
rs = 916
ks = 2.16
a = 7/100
Tnuc = -7
Tc = -12.2
t = np.arange(0, len(zf_mean))/1.43
t_frame = np.arange(0, len(zf_mean))

zs = np.sqrt(
    2*ks*(0-Tnuc)*t / (
        (1-a)*rs*Lf
    )
)
zs_o = np.sqrt(
    2*ks*(0-Tnuc)*t / (
        rs*Lf
    )
)
zs_p = np.sqrt(
    2*ks*(0-Tc)*t / (
        rs*Lf
    )
)

plt.figure(figsize=[8, 4.5])

plt.plot((t_frame - 3)/1.43, zf_mean/53250, 'oc')
plt.plot(t, zs, '--g', label='Corrected Stefan')
plt.plot(t, zs_o, '--m', label=r'$\Delta T = T_m - T_{nuc}$')
plt.plot(t, zs_p, '--r', label=r'$\Delta T = T_m - T_{c}$')

plt.legend(fancybox=True)
plt.grid(True)
plt.show()


# s.current_image.show_image()
# al, zf = s.read_manual_front()
# s.view_all_profil()
# s.view_dynamic_front(zf)
# s.view_propagation_front(zf)

#
# t0 = float(input('Enter a t0 for fit : '))
#
# s.view_dynamic_front(zf, fit_win=[t0, 10])
