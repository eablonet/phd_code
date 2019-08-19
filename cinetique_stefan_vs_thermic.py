import numpy as np
from scipy.special import erf

from matplotlib import pyplot as plt

t = np.linspace(1e-2, 60, 1000)
delta = np.sqrt(.04/2)

alpha_wa = 1.317119e-07
alpha_ice = 1.181976e-06
Ta = 1
Tm = 0


Vs = delta*np.sqrt(alpha_ice/t)

"""VT = (Ta - Tm) / np.sqrt(np.pi) * np.sqrt(alpha_ice/alpha_wa) * \
    delta * np.exp(-delta**2 * alpha_ice/alpha_wa) /\
    erf(delta) * 1/t
"""
VT = np.sqrt(alpha_wa / (2*delta*np.sqrt(alpha_ice*t)))
plt.figure(figsize=[8, 4.5])
plt.plot(
    t, Vs, '--', c='tab:orange', label='Solidification velocity'
)

plt.plot(
    t, VT, '--', c='tab:blue', label='Thermal cinetic'
)

plt.legend(fancybox=True, shadow=True)

plt.show()
