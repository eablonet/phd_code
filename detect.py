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

s.read_by_path('29-03-2018', 9)
# s.select_folder()
s.read_image(1)
# s.current_image.show_image()

# s.remove_treatments()
# s.crop()
# s.clahe()

# s.tracker('test')
s.contour_tracker2()

# zf, pts = s.read_data('contour')
# t, zf, zf_loc, _ = s.get_dynamic_front(zf, 300)
#
# plt.plot(range(0, len(zf_loc)), zf_loc)
# plt.show()
#
#
# zf, pts = s.read_data('front')
# zf_loc = s.get_dynamic_front_by_loc(zf, 300)
# Lf = 333550
# rs = 916
# ks = 2.16
# a = 7/100
# Tnuc = -7
# Tc = -12.2
# t = np.arange(0, len(zf_loc))/.83
# t_frame = np.arange(0, len(zf_loc))
# zs = np.sqrt(
#     2*ks*(0-Tnuc)*t / (
#         (1-a)*rs*Lf
#     )
# )
# zs_o = np.sqrt(
#     2*ks*(0-Tnuc)*t / (
#         rs*Lf
#     )
# )
# zs_p = np.sqrt(
#     2*ks*(0-Tc)*t / (
#         rs*Lf
#     )
# )
#
# plt.figure(figsize=[8, 4.5])
#
# plt.plot((t_frame - 4)/.83, zf_loc/65833, 'oc')
# plt.plot(t, zs, '--g', label='Corrected Stefan')
# plt.plot(t, zs_o, '--m', label=r'$\Delta T = T_m - T_{nuc}$')
# plt.plot(t, zs_p, '--r', label=r'$\Delta T = T_m - T_{c}$')
#
# plt.legend(fancybox=True)
# plt.grid(True)
# plt.show()
# s.current_image.show_image()
# al, zf = s.read_manual_front()
# s.view_all_profil()
# s.view_dynamic_front(zf)
# s.view_propagation_front(zf)
#
# t0 = float(input('Enter a t0 for fit : '))
#
# s.view_dynamic_front(zf, fit_win=[t0, 10])
